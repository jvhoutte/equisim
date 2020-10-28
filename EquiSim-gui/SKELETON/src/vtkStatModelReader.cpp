/*=========================================================================

  Program: Bender

  Copyright (c) Kitware Inc.

  Licensed under the Apache License, Version 2.0 (the "License");
  you may not use this file except in compliance with the License.
  You may obtain a copy of the License at

      http://www.apache.org/licenses/LICENSE-2.0.txt

  Unless required by applicable law or agreed to in writing, software
  distributed under the License is distributed on an "AS IS" BASIS,
  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
  See the License for the specific language governing permissions and
  limitations under the License.

=========================================================================*/

#include "SKELETON/vtkStatModelReader.h"


#include "SKELETON/vtkSkeleton.h"


#include "vtkAbstractTransform.h"
#include "vtkCollection.h"
#include "vtkExecutive.h"
#include "vtkInformation.h"
#include "vtkInformationVector.h"
#include "vtkNew.h"
#include "vtkObjectFactory.h"
#include "vtkPolyData.h"
#include "vtkTransform.h"

#include <algorithm>
#include <cctype>
#include <fstream>
#include <sstream>
#include <stdexcept>

namespace {
//----------------------------------------------------------------------------
const std::string GetKeyword(const std::string& line)
{
  std::string keyword = "";
  bool wordHasStarted = false;
  for (std::string::const_iterator it = line.begin(); it < line.end(); ++it)
    {
    // stop at the first space after the keyword
    if (keyword != "" && isspace(*it))
      {
      return keyword;
      }

    if (*it == '{')
      {
      return "{";
      }
    else if(*it == '}')
      {
      return "}";
      }
    else if(*it == '#')
      {
      return wordHasStarted ? keyword : "#";
      }

    if (!wordHasStarted)
      {
      wordHasStarted = isupper(*it);
      }

    if (wordHasStarted)
      {
      keyword += *it;
      }
    }
  return keyword;
}

//----------------------------------------------------------------------------
const std::string
MoveToNextKeyword(std::ifstream& file, std::string& line)
{
  std::string keyword = GetKeyword(line);
  while(keyword == "#")
    {
    if (!std::getline(file, line))
      {
      return "";
      }
    keyword = GetKeyword(line);
    }

  return keyword;
}

//----------------------------------------------------------------------------
template<typename T>
void GetValues(const std::string& line, std::vector<T>& values)
{
  std::stringstream ss;
  size_t valueStart = 0;
  while (valueStart != std::string::npos)
    {
    size_t valueStop = line.find(" ", valueStart+1);
    ss << line.substr(valueStart, valueStop - valueStart);

    T value;
    ss >> value;
    values.push_back(value);

    ss.clear();
    valueStart = valueStop;
    }
}

//----------------------------------------------------------------------------
template<typename T>
void GetValues(const std::string& line,
               std::vector<T>& values,
               const std::string keyword)
{
  size_t valueStart = line.rfind(keyword) + keyword.size() + 1;
  GetValues<T>(line.substr(valueStart), values);
}

//----------------------------------------------------------------------------
void GetHead(const std::string& line, double* offset)
{
  std::vector<double> offsetVect;
  GetValues<double>(line, offsetVect, "HEAD");
  std::copy(offsetVect.begin(), offsetVect.end(), offset);
}

//----------------------------------------------------------------------------
void GetTail(const std::string& line, double* offset)
{
  std::vector<double> offsetVect;
  GetValues<double>(line, offsetVect, "TAIL");
  std::copy(offsetVect.begin(), offsetVect.end(), offset);
}

//----------------------------------------------------------------------------
void GetRollAxis(const std::string& line, double* offset)
{
  std::vector<double> offsetVect;
  GetValues<double>(line, offsetVect, "ROLLAXIS");
  std::copy(offsetVect.begin(), offsetVect.end(), offset);
}

//----------------------------------------------------------------------------
void GetRange(const std::string& line, double* offset)
{
  std::vector<double> offsetVect;
  GetValues<double>(line, offsetVect, "RANGE");
  std::copy(offsetVect.begin(), offsetVect.end(), offset);
}

//----------------------------------------------------------------------------
void GetChannels(const std::string& line, std::vector<std::string>& channels)
{
  GetValues<std::string>(line, channels, "CHANNELS");
}

//----------------------------------------------------------------------------
template<typename T> T GetValue(const std::string& line, std::string keyword)
{
  std::vector<T> values;
  GetValues<T>(line, values, keyword);
  assert(values.size() == 1);
  return values[0];
}

//----------------------------------------------------------------------------
std::string
GetBoneName(const std::string& line, std::string keyword = "JOINT")
{
  return GetValue<std::string>(line, keyword);
}
} // end namespace


vtkStandardNewMacro(vtkStatModelReader);

//----------------------------------------------------------------------------
vtkStatModelReader::vtkStatModelReader()
{  
  this->Armature = NULL;
  this->SetNumberOfInputPorts(0);
  
  Shape = vtkMultiBlockDataSet::New();
}

//----------------------------------------------------------------------------
vtkStatModelReader::~vtkStatModelReader()
{
  if (this->Armature)
  {
    this->Armature->Delete();
  }
  if (this->Shape)
  {
    this->Shape->Delete();
  }
    
}

//----------------------------------------------------------------------------
void vtkStatModelReader::SetBVHFileName(const char* filename)
{
  if( CanReadFile(filename, "bvh")==0 )
    return;
  
  vtkSmartPointer<vtkBVHReader> bvhreader = vtkSmartPointer<vtkBVHReader>::New();
  bvhreader->SetFileName(filename);
  bvhreader->Update();
  Armature = bvhreader->GetArmature();
  
  this->Modified();
}

//----------------------------------------------------------------------------
void vtkStatModelReader::SetSeparateFileNames(std::vector<std::string> filenames){
  for(int f=0; f<filenames.size(); f++){
    // Retrieve the modelname and extension from the filename
    std::string filename = filenames[f];
    std::string ext = filename.substr(filename.find_last_of(".") + 1);
    std::string modelname = filename.substr(filename.find_last_of("/\\") + 1);
    modelname = modelname.substr(0, modelname.find("."));
    
    // Read the polydata
    vtkSmartPointer<vtkPolyData> polydata = vtkSmartPointer<vtkPolyData>::New();
    if( ext == "stl"){
      vtkSmartPointer<vtkSTLReader> reader = vtkSmartPointer<vtkSTLReader>::New();
      reader->SetFileName(filename.c_str());
      reader->Update();
      polydata->DeepCopy( reader->GetOutput() );
    }
    else if( ext == "vtk"){
      vtkSmartPointer<vtkPolyDataReader> reader = vtkSmartPointer<vtkPolyDataReader>::New();
      reader->SetFileName(filename.c_str());
      reader->Update();
      polydata->DeepCopy( reader->GetOutput() );
    }
    
    // Save the bone name together with the mesh
    vtkSmartPointer<vtkStringArray> namearray = vtkSmartPointer<vtkStringArray>::New();   
    namearray->SetNumberOfValues(1);
    namearray->SetName("ModelName");
    namearray->SetValue(0, modelname.c_str());
    polydata->GetFieldData()->AddArray(namearray);
    Shape->SetBlock(f, polydata);
  }
  
  statmode = false;
}

// Checks if the multiblockdataset is a statistical model. 
// It therefore checksif each subblock has an array with title "Blockname"
bool vtkStatModelReader::IsStatModel(vtkSmartPointer<vtkMultiBlockDataSet> multiblock){
    for(int blockidx=0; blockidx<multiblock->GetNumberOfBlocks(); blockidx++){
        if( multiblock->GetBlock(blockidx)->IsA("vtkMultiBlockDataSet")){
            vtkSmartPointer<vtkMultiBlockDataSet> block = vtkMultiBlockDataSet::SafeDownCast(multiblock->GetBlock(blockidx));
            for(int subblockidx=0; subblockidx<block->GetNumberOfBlocks(); subblockidx++){
            vtkSmartPointer<vtkPolyData> poly = vtkPolyData::SafeDownCast(block->GetBlock(subblockidx));
            if( poly->GetFieldData()->GetAbstractArray("BlockName")==NULL){return false;}
            }
        }
        else{
            return false;
        }
    }
    return true;
}

vtkSmartPointer<vtkMultiBlockDataSet> vtkStatModelReader::MergeNoBones(vtkSmartPointer<vtkMultiBlockDataSet> multiblock){
    for(int blockidx=0; blockidx<multiblock->GetNumberOfBlocks(); blockidx++){
        if( multiblock->GetBlock(blockidx)->IsA("vtkMultiBlockDataSet")){
            vtkSmartPointer<vtkMultiBlockDataSet> block = vtkMultiBlockDataSet::SafeDownCast(multiblock->GetBlock(blockidx)); 
            vtkSmartPointer<vtkAppendPolyData> appendFilter = vtkSmartPointer<vtkAppendPolyData>::New();
            for(int subblockidx=0; subblockidx<block->GetNumberOfBlocks(); subblockidx++){
                  appendFilter->AddInputData(vtkPolyData::SafeDownCast(block->GetBlock(subblockidx)));
            }
            appendFilter->Update();
            vtkSmartPointer<vtkPolyData> merged = appendFilter->GetOutput();
            
            // Transfer the bonename array
            merged->GetFieldData()->AddArray( vtkPolyData::SafeDownCast(block->GetBlock(0))->GetFieldData()->GetAbstractArray("ModelName") );
            
            multiblock->SetBlock(blockidx, merged);
        }
    }
    return multiblock;    
}



//----------------------------------------------------------------------------
void vtkStatModelReader::SetVTMFileName(const char* filename)
{
  if( CanReadFile(filename, "vtm")==1 ){
  
  vtkSmartPointer<vtkXMLMultiBlockDataReader> reader = vtkSmartPointer<vtkXMLMultiBlockDataReader>::New();
  reader->SetFileName(filename);
  reader->Update();
  
  
  
  
  vtkSmartPointer<vtkMultiBlockDataSet> multiblock = vtkMultiBlockDataSet::SafeDownCast( reader->GetOutput() );
  Ncompounds = multiblock->GetNumberOfBlocks();
  
  if( !IsStatModel(multiblock) ){
        Shape->DeepCopy( MergeNoBones(vtkMultiBlockDataSet::SafeDownCast(reader->GetOutput())) ); 
        statmode = false;
  }
  else{
        aSSM = TensorBlocks(Ncompounds);
        
        for(int blockidx=0; blockidx<Ncompounds; blockidx++){
            vtkSmartPointer<vtkMultiBlockDataSet> block = vtkMultiBlockDataSet::SafeDownCast(multiblock->GetBlock(blockidx));
                
            std::vector<bool> checks(3, false);
            
            for(int subblockidx=0; subblockidx<block->GetNumberOfBlocks(); subblockidx++){
            vtkSmartPointer<vtkPolyData> poly = vtkPolyData::SafeDownCast(block->GetBlock(subblockidx));
            
            if( poly->GetFieldData()->GetAbstractArray("ModelName")!=NULL){
            std::string modelname = vtkStringArray::SafeDownCast(poly->GetFieldData()->GetAbstractArray("ModelName"))->GetValue(0);
            aSSM[blockidx].name = modelname;
            }
            
            std::string blockname="";
            if( poly->GetFieldData()->GetAbstractArray("BlockName")!=NULL){
            blockname = vtkStringArray::SafeDownCast(poly->GetFieldData()->GetAbstractArray("BlockName"))->GetValue(0);
            }
            else{
            throw std::runtime_error("Block in the aSSM does not contain a BlockName."); 
            }
            
            if(blockname=="geom"){
            Shape->SetBlock( blockidx, poly );
            aSSM[blockidx].geometry.average = poly;
            aSSM[blockidx].geometry.eigenvectors.SetArrays( poly );
            checks[0]=true;
            }
            else if(blockname=="cor"){
            aSSM[blockidx].cor.average = poly;
            aSSM[blockidx].cor.eigenvectors.SetArrays( poly );
            checks[1]=true;
            }
            else if(blockname=="axis"){
            aSSM[blockidx].axis.average = poly;
            aSSM[blockidx].axis.eigenvectors.SetArrays( poly );
            aSSM[blockidx].axis.pga.SetArrays( poly );
            checks[2]=true;
            }
            
            } // end for-loop over subblocks
            
            // Check if everything was provided for this compound 
            if(aSSM[blockidx].name=="" || !checks[0]  || !checks[1] || !checks[2]){
            throw std::runtime_error("The provided aSSM did not provide all necessary information."); 
            }
            
        } // end for-loop over blocks 
        
        Nmodes = aSSM[0].geometry.eigenvectors.GetNumberOfModes();
        pcweights = std::vector<double>(Nmodes, 0.);
        
        statmode = true;
        
  }
  
  this->Modified();
  
  }
  else{ return; } // in case we cannot read the vtm file 
}

std::vector<int> vtkStatModelReader::GetNonZeroModes(){
    std::vector<int> selection;
    for(int i=0; i<Nmodes; i++){ if(std::abs(pcweights[i])>1e-18){selection.push_back(i);} }
    return selection;
}

vtkSmartPointer<vtkFloatArray> vtkStatModelReader::GetCurrentVectorField(Tensor* tensor){
  
    
  vtkSmartPointer<vtkFloatArray> vectorarray = vtkSmartPointer<vtkFloatArray>::New();
  int Ntuples = (*tensor)[0]->GetNumberOfTuples();
  vectorarray->SetNumberOfComponents(3); // Important to first set Ncomp and next the Ntuples
  vectorarray->SetNumberOfTuples(Ntuples);
  
  vectorarray->FillComponent(0,0.);
  vectorarray->FillComponent(1,0.);
  vectorarray->FillComponent(2,0.);
  
  // get list of non-zero pc weights, such that we dont lose time multiplying arrays with 0
  std::vector<int> selection = GetNonZeroModes();
  
  for(int pointidx=0; pointidx<Ntuples; pointidx++){
    double currentvec[3]={0,0,0};
    //for(int mode=0; mode<tensor->GetNumberOfModes(); mode++){
    for(int mode_id=0; mode_id<selection.size(); mode_id++){
      int mode = selection[mode_id];
      double newvec[3]={0,0,0};
      (*tensor)[mode]->GetTuple(pointidx, newvec);
      vtkMath::MultiplyScalar(newvec, pcweights[mode]);
      vtkMath::Add(currentvec, newvec, currentvec);
    }
    vectorarray->SetTuple(pointidx, currentvec);
  }
  
  return vectorarray;
}

vtkSmartPointer<vtkPoints> vtkStatModelReader::ChangePCWeights(PcModel pcmodel){
  
  vtkSmartPointer<vtkPoints> newpoints = vtkSmartPointer<vtkPoints>::New();
  
  // Get the total displacement vectors for each point
  vtkSmartPointer<vtkFloatArray> vectorfield = GetCurrentVectorField(&pcmodel.eigenvectors);
  
  for(int pointidx=0; pointidx<pcmodel.average->GetNumberOfPoints(); pointidx++){
    double vec[3]; double point[3];
    vectorfield->GetTuple(pointidx, vec);
    pcmodel.average->GetPoint(pointidx, point);
    vtkMath::Add(point, vec, point);
    newpoints->InsertNextPoint(point);
  }
  
  
  /*
   * Using the vtk wrapper filter is almost 3 times slower...
  int Nmodes = pcmodel.eigenvectors.GetNumberOfModes();
  vtkSmartPointer<vtkPolyData> poly = pcmodel.average;
  for(int mode=0; mode<Nmodes; mode++){
    vtkSmartPointer<vtkWarpVector> warp = vtkSmartPointer<vtkWarpVector>::New();
    warp->SetInputData(poly);
    warp->SetInputArrayToProcess(0, 0, 0, vtkDataObject::FIELD_ASSOCIATION_POINTS, pcmodel.eigenvectors[mode]->GetName());
    warp->SetScaleFactor(pcweights[mode]);
    warp->Update();
    poly = vtkPolyData::SafeDownCast(warp->GetOutput());      
  }
  return poly->GetPoints();
  */
  
  return newpoints;
}

vtkSmartPointer<vtkStringArray> vtkStatModelReader::GetBlockNameArray(vtkSmartPointer<vtkMultiBlockDataSet> block){
  return vtkStringArray::SafeDownCast(block->GetFieldData()->GetAbstractArray("ModelName"));
}

std::string vtkStatModelReader::GetBlockName(vtkSmartPointer<vtkMultiBlockDataSet> block){
  return GetBlockNameArray(block)->GetValue(0);
}

int vtkStatModelReader::GetModelIndex(int boneidx){
  if ( bonetomodelmap.find(boneidx) == bonetomodelmap.end() ) { // not added yet
    vtkStringArray* BoneNames = vtkStringArray::SafeDownCast(Armature->GetVertexData()->GetAbstractArray("BoneName"));
    std::string bonename = BoneNames->GetValue(boneidx);
    
    for(int modelidx=0; modelidx<Ncompounds; modelidx++){
      std::string modelname = aSSM[modelidx].name;
      if(modelname == bonename){
	bonetomodelmap[boneidx]=modelidx;
      }
    }
    
    // check again
    if ( bonetomodelmap.find(boneidx) == bonetomodelmap.end() ) {
      return Ncompounds; // means that there is no corresponding shape model
    }
  }
  
  return bonetomodelmap[boneidx];
}

void vtkStatModelReader::SetPcWeights(int mode, double weight){
  if(mode<Nmodes){
    pcweights[mode]=weight;
  }  
  this->Modified();
}



int vtkStatModelReader::GetBoneIndex(int modelidx){
  if ( modeltobonemap.find(modelidx) == modeltobonemap.end() ) { // not added yet
    
    std::string modelname = aSSM[modelidx].name;
    
    if( modelname == "Skin" ){
      modeltobonemap[modelidx] = Armature->GetNumberOfVertices(); // give it a random index which exceeds the bounds of the bone ids in the armature 
    }
    else{
      vtkStringArray* BoneNames = vtkStringArray::SafeDownCast(Armature->GetVertexData()->GetAbstractArray("BoneName"));
      for(int counter=0; counter<BoneNames->GetNumberOfValues(); counter++){
	if(modelname == BoneNames->GetValue(counter)){
	  modeltobonemap[modelidx]=counter;
	}
      }
    }
    
    if ( modeltobonemap.find(modelidx) == modeltobonemap.end() ) {
      throw std::runtime_error("Names of the provided surface files need to match the names of the bones in the vtk file."); 
    }
  }
  
  return modeltobonemap[modelidx];
}

void vtkStatModelReader::ExpMap(double (*v)[3], double (*mean)[3])
{
  double norm = vtkMath::Normalize(*v);
  vtkMath::MultiplyScalar(*mean, std::cos(norm));
  vtkMath::MultiplyScalar(*v    , std::sin(norm));
  vtkMath::Add(*mean, *v, *v);
}



void vtkStatModelReader::UpdateBone(){
  vtkSmartPointer<vtkTreeDFSIterator> childreniterator = vtkSmartPointer<vtkTreeDFSIterator>::New();
  childreniterator->SetMode(0); // discover-mode. Top-down iterator 
  childreniterator->SetTree(Armature);
  childreniterator->SetStartVertex(Armature->GetRoot());
  while(childreniterator->HasNext()){
    int boneidx = childreniterator->Next();
	
    int modelidx = GetModelIndex(boneidx);
    
    
    if(modelidx < Ncompounds){
      // transform everything to world coordinate system 
      vtkSmartPointer<vtkTransform> localtoworldtransform = Armature->GetTransformToParentsLocalFrame(boneidx);
      //localtoworldtransform->Inverse();
      
      // Update geometry
      
      vtkSmartPointer<vtkPoints> geompoints = ChangePCWeights(aSSM[modelidx].geometry);
      vtkPolyData::SafeDownCast(Shape->GetBlock(modelidx))->SetPoints(geompoints);
      
      vtkSmartPointer<vtkTransformPolyDataFilter> polytransform = vtkSmartPointer<vtkTransformPolyDataFilter>::New();
      polytransform->SetTransform(localtoworldtransform);
      polytransform->SetInputData(vtkPolyData::SafeDownCast(Shape->GetBlock(modelidx)));
      polytransform->Update();
      Shape->SetBlock(modelidx, polytransform->GetOutput());
      
      // Update COR
      vtkSmartPointer<vtkPoints> corpoints = ChangePCWeights(aSSM[modelidx].cor);
      
      // Update axis 
      vtkSmartPointer<vtkPoints> axispoints = ChangePCWeights(aSSM[modelidx].axis);
      double pgamean[3]; 
              
      if (!aSSM[modelidx].axis.pga.GetPgaMean(&pgamean)){std::runtime_error("No PGA mean loaded");}
      double axisz[3]; axispoints->GetPoint(0, axisz); ExpMap(&axisz, &pgamean);
      
      vtkSmartPointer<vtkPoints> corpoints_transformed = vtkSmartPointer<vtkPoints>::New();
      localtoworldtransform->TransformPoints(corpoints, corpoints_transformed);
      double head[3] = {0,0,0}; localtoworldtransform->TransformPoint(head, head);
      double tail[3]; corpoints_transformed->GetPoint(0, tail);  
      localtoworldtransform->TransformVector(axisz, axisz);
      
      
      vtkMath::Normalize(axisz);
      
      // Update properties in the skeleton 
      Armature->ChangeRestArmatureWithoutSurfaceChange(boneidx, axisz, head, tail);
      
    }
    
    
  }
    
}


  
//----------------------------------------------------------------------------
int vtkStatModelReader::RequestData(vtkInformation *vtkNotUsed(request),
                              vtkInformationVector **vtkNotUsed(inputVector),
                              vtkInformationVector *outputVector)
{
  if (!this->Armature)
    {
      vtkErrorMacro("Didnt provide a bvh file.");
      return 0;
    }

    
  if(statmode == true){
    UpdateBone();
  }
        
  return 1;
}


//----------------------------------------------------------------------------
int vtkStatModelReader::CanReadFile(const char *filename, const char *extension)
{
  std::ifstream file(filename);
  if (!file.good())
    {
    return 0;
    }

  std::string name = filename;
  if (name.size() < 3 || name.substr(name.size() - 3) != extension)
    {
    return 0;
    }

  return 1;
}


//----------------------------------------------------------------------------
void vtkStatModelReader::PrintSelf(ostream& os, vtkIndent indent)
{
  this->Superclass::PrintSelf(os,indent);

}
