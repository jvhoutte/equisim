#include "SKELETON/DeformRigidSurface.h"
 
#include "vtkObjectFactory.h"
#include "vtkStreamingDemandDrivenPipeline.h"
#include "vtkInformationVector.h"
#include "vtkInformation.h"
#include "vtkDataObject.h"
#include "vtkSmartPointer.h"
#include "vtkDoubleArray.h"
#include "vtkCellData.h"

#include <thread>

vtkStandardNewMacro(vtkRigidSurfaceTransformFilter);

vtkRigidSurfaceTransformFilter::vtkRigidSurfaceTransformFilter()
{
  this->SetNumberOfInputPorts(3);
  //this->SetNumberOfOutputPorts(1);
}
 
vtkRigidSurfaceTransformFilter::~vtkRigidSurfaceTransformFilter()
{
}
 
int vtkRigidSurfaceTransformFilter::FillInputPortInformation( int port, vtkInformation* info )
{
  if ( port == 0 )
    {
    info->Set(vtkDataObject::DATA_TYPE_NAME(), "vtkMultiBlockDataSet" );
    return 1;
    }
  else if(port == 1)
    {
    info->Set(vtkDataObject::DATA_TYPE_NAME(), "vtkSkeleton" );
    return 1;
    }
  else if(port == 2)
    {
    info->Set(vtkDataObject::DATA_TYPE_NAME(), "vtkUnstructuredGrid" );
    info->Set(vtkAlgorithm::INPUT_IS_OPTIONAL(), 1);
    return 1;
    }
  return 0;
}
 
 
int vtkRigidSurfaceTransformFilter::FillOutputPortInformation(int vtkNotUsed(portNumber), vtkInformation *info)
{
    info->Set(vtkDataObject::DATA_TYPE_NAME(), "vtkMultiBlockDataSet");
    return 1;
}
 
int vtkRigidSurfaceTransformFilter::RequestData(vtkInformation *vtkNotUsed(request),
                                          vtkInformationVector **inputVector,
                                          vtkInformationVector *outputVector)
{
    // get the input objects 
    
    this->surfaceRest = vtkMultiBlockDataSet::GetData(inputVector[0],0);
    
    this->armature = vtkTree::GetData(inputVector[1],0);

    // get the output objects 
    vtkInformation *outInfo = outputVector->GetInformationObject(0);
    this->surfacePose = vtkMultiBlockDataSet::SafeDownCast(outInfo->Get(vtkDataObject::DATA_OBJECT()));
    this->surfacePose->DeepCopy(surfaceRest);
    
    // Get the optional input object 
    if(this->GetNumberOfInputConnections(2)>0){ 
        this->tetmodel = vtkUnstructuredGrid::GetData(inputVector[2],0); 
        this->tetmodelPose = vtkSmartPointer<vtkUnstructuredGrid>::New();
        this->tetmodelPose->DeepCopy(tetmodel);
        //PBS_SetReferenceInfo();
    }
    else{
        this->tetmodel = NULL;
        this->tetmodelPose = NULL;
    }
        
    // Check inputs 
    
    assert(surfaceRest != NULL);
    assert(armature != NULL);
    assert(surfacePose != NULL);

    // Apply the skinning deformation 
    this->TransformModels();
    
    // Update cell and point normals
    /*
    for(int b=0; b<surfacePose->GetNumberOfBlocks(); b++){
      
      int ncells = vtkPolyData::SafeDownCast( surfacePose->GetBlock(b) )->GetNumberOfCells();
      
      vtkSmartPointer<vtkPolyDataNormals> normalGenerator = vtkSmartPointer<vtkPolyDataNormals>::New();
      normalGenerator->SetInputData(vtkPolyData::SafeDownCast( surfacePose->GetBlock(b) ));
      normalGenerator->ComputePointNormalsOff();
      normalGenerator->ComputeCellNormalsOn();
      normalGenerator->Update();
      vtkPolyData::SafeDownCast( surfacePose->GetBlock(b) )->GetCellData()->SetNormals( normalGenerator->GetOutput()->GetCellData()->GetNormals() );
      assert(ncells == vtkPolyData::SafeDownCast( surfacePose->GetBlock(b) )->GetNumberOfCells());
      
    }
    */
    
  
    
    //this->surfacePose->GetPointData()->DeepCopy( normalGenerator->GetOutput()->GetPointData() );
    //this->surfacePose = normalGenerator->GetOutput();

    return 1;
}
 
//----------------------------------------------------------------------------
void vtkRigidSurfaceTransformFilter::PrintSelf(ostream& os, vtkIndent indent)
{
  this->Superclass::PrintSelf(os,indent);
}

/*
void vtkRigidSurfaceTransformFilter::SetJacobian(){
  
  vtkStringArray* BoneNames = vtkStringArray::SafeDownCast(armature->GetVertexData()->GetAbstractArray("BoneName"));
  vtkDoubleArray* HeadPose = vtkDoubleArray::SafeDownCast(armature->GetVertexData()->GetArray("HeadPose"));
  vtkDoubleArray* TailPose = vtkDoubleArray::SafeDownCast(armature->GetVertexData()->GetArray("TailPose"));
  
  
  int Nbones = BoneNames->GetNumberOfValues();
  

  
    // Loop over all DOF's
    for(int boneidx = 0; boneidx < Nbones; boneidx++){
      
      
      // Get COR
      double* cor; 
      if(armature->GetParent(boneidx)){
	cor = TailPose->GetTuple(armature->GetParent(boneidx));
      }
      else{
	cor = HeadPose->GetTuple(armature->GetRoot());
      }
            
      double ax[9]; vtkDoubleArray::SafeDownCast(armature->GetVertexData()->GetArray("BoneRotAxes"))->GetTuple(boneidx, ax);
      for(int a=0; a<3; a++){
	  double axis[3] = {ax[3*a], ax[3*a+1], ax[3*a+2]};
	  
	  // Loop over blocks in multiblockdataset
	  for(int blockid = 0; blockid < surfacePose->GetNumberOfBlocks(); blockid++){
	    if(blockid >= boneidx){ // vertices connected to bones higher up in the kinematical tree are not influenced by the lower bones 
		
	      vtkSmartPointer<vtkPolyData> polydatablock = vtkPolyData::SafeDownCast(surfacePose->GetBlock(blockid));
	      
	      // Create array with jacobians for this DOF and for the points of this block
	      vtkSmartPointer<vtkDoubleArray> JacobianArray = vtkSmartPointer<vtkDoubleArray>::New();
	      //JacobianArray->SetNumberOfTuples(polydatablock->GetNumberOfPoints());
	      JacobianArray->SetNumberOfComponents(3);
	      JacobianArray->SetName(("Jac_"+BoneNames->GetValue(boneidx)+"_"+std::to_string(a)).c_str());
	      
	      // Loop over vertices of this block
	      for(int vid = 0; vid < polydatablock->GetNumberOfPoints(); vid++){
		double point[3]; polydatablock->GetPoint(vid, point);
		double res[3]; vtkMath::Subtract(point, cor, res);
		
		double jacvec[3];
		vtkMath::Cross(axis, res, jacvec);
		JacobianArray->InsertNextTuple(jacvec);
	      }
	      
	      // Add array
	      polydatablock->GetPointData()->AddArray(JacobianArray);
	    }
	    else{ // skip this block
	      // do nothing
	    }  
	} // end loop over blockids
      } // end loop over axis ids
    } // end loop over boneidx
}
*/

vtkIdType vtkRigidSurfaceTransformFilter::GetBoneIndex(vtkIdType modelidx){
    
    return modelidx;
    /*
    
  if ( modelbonemap.find(modelidx) == modelbonemap.end() ) { // not added yet
    vtkStringArray* ModelNames = vtkStringArray::SafeDownCast(surfaceRest->GetBlock(modelidx)->GetFieldData()->GetAbstractArray("ModelName"));
    std::cout << "ModelName in skinning " << ModelNames->GetValue(0) << std::endl;
    
    if( ModelNames->GetValue(0) == "Skin" ){
      modelbonemap[modelidx] = armature->GetNumberOfVertices(); // give it a random index which exceeds the bounds of the bone ids in the armature 
    }
    else{
      vtkStringArray* BoneNames = vtkStringArray::SafeDownCast(armature->GetVertexData()->GetAbstractArray("BoneName"));
      for(int counter=0; counter<BoneNames->GetNumberOfValues(); counter++){
        if(ModelNames->GetValue(0) == BoneNames->GetValue(counter)){
        modelbonemap[modelidx]=counter;
        }
      }
    }
    
    if ( modelbonemap.find(modelidx) == modelbonemap.end() ) {
      throw std::runtime_error("Names of the provided surface files need to match the names of the bones in the vtk file."); 
    }
  }
  
  return modelbonemap[modelidx];
  */
}


void vtkRigidSurfaceTransformFilter::SetTransforms(){
  // For articulation-based deformations
  vtkDoubleArray* HeadPose = vtkDoubleArray::SafeDownCast(armature->GetVertexData()->GetArray("HeadPose"));
  vtkDoubleArray* HeadRest = vtkDoubleArray::SafeDownCast(armature->GetVertexData()->GetArray("HeadRest"));
  vtkDoubleArray* RestToPoseRotation = vtkDoubleArray::SafeDownCast(armature->GetVertexData()->GetArray("RestToPoseRotation"));
  vtkDoubleArray* WorldToPoseRotation = vtkDoubleArray::SafeDownCast(armature->GetVertexData()->GetArray("WorldToBoneRotation"));
  
  // For anthropometric deformations
  vtkDoubleArray* ReferenceToPoseTransform = vtkDoubleArray::SafeDownCast(armature->GetVertexData()->GetArray("ReferenceToPoseTransform"));
  
  if (HeadPose == NULL){throw std::runtime_error("Armature has no associated transformations.");}
  if (HeadRest == NULL){throw std::runtime_error("Armature has no associated transformations.");}
  if (RestToPoseRotation == NULL){throw std::runtime_error("Armature has no associated transformations.");}
  if (WorldToPoseRotation == NULL){throw std::runtime_error("Armature has no associated transformations.");}

  for(vtkIdType b = 0; b < armature->GetNumberOfVertices(); b++){
	  // Anthropometric deformation
	  double transf[9];
	  ReferenceToPoseTransform->GetTuple(b, transf);
	  double (*rotation)[3] = (double (*)[3])&transf[0];
	  
	  // Articulation (rest to pose rotation)
	  double quatarray[4];
	  RestToPoseRotation->GetTuple(b, quatarray); vtkQuaterniond quat(quatarray);
	  double resttopose[3][3]; quat.ToMatrix3x3(resttopose);
	  
	  // Combine both transformations
	  vtkMath::Multiply3x3(resttopose, rotation, rotation);
	  
	  // Get translation vector  
	  double headpose[3]; double headrest[3]={0,0,0}; double transl[3]={0,0,0};
	  HeadPose->GetTuple(b, headpose); 
	  HeadRest->GetTuple(b, headrest);
	  double rothead[3];
	  vtkMath::Multiply3x3(rotation, headrest, rothead);
	  vtkMath::Subtract(headpose, rothead, transl);
	  
	  vtkSmartPointer<vtkMatrix4x4> totaltransformmatrix = vtkSmartPointer<vtkMatrix4x4>::New();
	  totaltransformmatrix->Identity();
	  totaltransformmatrix->SetElement(0,0,rotation[0][0]);
	  totaltransformmatrix->SetElement(0,1,rotation[0][1]);
	  totaltransformmatrix->SetElement(0,2,rotation[0][2]);
	  totaltransformmatrix->SetElement(1,0,rotation[1][0]);
	  totaltransformmatrix->SetElement(1,1,rotation[1][1]);
	  totaltransformmatrix->SetElement(1,2,rotation[1][2]);
	  totaltransformmatrix->SetElement(2,0,rotation[2][0]);
	  totaltransformmatrix->SetElement(2,1,rotation[2][1]);
	  totaltransformmatrix->SetElement(2,2,rotation[2][2]);
	  totaltransformmatrix->SetElement(0,3,transl[0]);
	  totaltransformmatrix->SetElement(1,3,transl[1]);
	  totaltransformmatrix->SetElement(2,3,transl[2]);

	  transformmap[b] = vtkSmartPointer<vtkTransform>::New();
	  transformmap[b]->SetMatrix(totaltransformmatrix);
	}
}

//----------------------------------------------------------------------------
vtkSmartPointer<vtkTransform> vtkRigidSurfaceTransformFilter::GetBlockTransform(int blockid){
  return transformmap[GetBoneIndex(blockid)];
}

vtkSmartPointer<vtkPointSet> vtkRigidSurfaceTransformFilter::WeightedTransform(vtkSmartPointer<vtkPointSet> surf){
    int Ncomp = surf->GetPointData()->GetArray("skinning_weights_normalised")->GetNumberOfComponents();
    
    vtkSmartPointer<vtkWeightedTransformFilter> weightedTrans = vtkSmartPointer<vtkWeightedTransformFilter>::New();
    weightedTrans->SetNumberOfTransforms(Ncomp);
    
    for(int c=0; c<Ncomp; c++){
        weightedTrans->SetTransform(GetBlockTransform(c), c);
    }
    
    weightedTrans->SetWeightArray("skinning_weights_normalised");  
    weightedTrans->SetInputData(surf);
    weightedTrans->Update();
    
    return weightedTrans->GetOutput();
}

vtkSmartPointer<vtkPolyData> vtkRigidSurfaceTransformFilter::RigidTransform(vtkSmartPointer<vtkPolyData> surf, int block_id){
    vtkSmartPointer<vtkTransform> transform = GetBlockTransform(block_id);
    vtkSmartPointer<vtkTransformPolyDataFilter> polydatatransformer = vtkSmartPointer<vtkTransformPolyDataFilter>::New();
    polydatatransformer->SetTransform(transform);
    polydatatransformer->SetInputData(surf);
    polydatatransformer->Update();
    return RecalculateNormals(polydatatransformer->GetOutput());
}

vtkSmartPointer<vtkPolyData> vtkRigidSurfaceTransformFilter::RecalculateNormals(vtkSmartPointer<vtkPolyData> surf){
    // Recalculate normals
	vtkSmartPointer<vtkPolyDataNormals> normalGenerator = vtkSmartPointer<vtkPolyDataNormals>::New();
	normalGenerator->SetInputData(surf);
	normalGenerator->ComputePointNormalsOn();
	normalGenerator->ComputeCellNormalsOn();
	normalGenerator->SetFeatureAngle(180.0);
	normalGenerator->SetConsistency( true );
	normalGenerator->SetSplitting( false );
	normalGenerator->Update();

	surf->GetPointData()->RemoveArray("Normals");
	vtkDataArray* normalData = vtkDataArray::SafeDownCast(normalGenerator->GetOutput()->GetPointData()->GetNormals());
	surf->GetPointData()->SetNormals(normalData);
    
    surf->GetCellData()->RemoveArray("Normals");
	vtkDataArray* cellnormalData = vtkDataArray::SafeDownCast(normalGenerator->GetOutput()->GetCellData()->GetNormals());
	surf->GetCellData()->SetNormals(cellnormalData);
    
    return surf;
}

template <typename modelType>
vtkSmartPointer<modelType> vtkRigidSurfaceTransformFilter::ApplyTransform(vtkDataObject* surf){
    vtkSmartPointer<modelType> newsurf = vtkSmartPointer<modelType>::New();
    newsurf->DeepCopy( modelType::SafeDownCast(surf) );
    
    if( newsurf->GetPointData()->GetArray("skinning_weights_normalised") != NULL ){
        newsurf->SetPoints( WeightedTransform(newsurf)->GetPoints() );
    }
    else{throw std::runtime_error("Surface model does not contain skinning weights.");}
	
    return newsurf;
}

// Check if model is a surface model, by verifying if it has a skinning weight array associated to it. 
bool vtkRigidSurfaceTransformFilter::IsBoneModel(vtkDataObject* model){
    vtkSmartPointer<vtkDataArray> array; 
    if(model->IsA("vtkPolyData")){
        array = vtkPolyData::SafeDownCast(model)->GetPointData()->GetArray("skinning_weights_normalised");
    }
    else if(model->IsA("vtkUnstructuredGrid")){
        array = vtkUnstructuredGrid::SafeDownCast(model)->GetPointData()->GetArray("skinning_weights_normalised");
    }
    
    if(array==NULL){return true;}
    else{return false;}    
}


void vtkRigidSurfaceTransformFilter::TransformModels(){
    SetTransforms(); // stores the transformations of ALL bones in the skeleton, irrespective of which surfacemodels are provided. 
    
    // Transform the surface models in the multiblock object "surfaceRest"
    SurfaceTransform();
    
    /*
    vtkSmartPointer<vtkXMLMultiBlockDataWriter> writer = vtkSmartPointer<vtkXMLMultiBlockDataWriter>::New();
    writer->SetInputData(surfacePose);
    writer->SetFileName("surfacePose.vtm");
    writer->Update();
    */
    
    // Transform the volumetric tetrahedral model if available
    VolumeTransform();
    
}

//----------------------------------------------------------------------------
void vtkRigidSurfaceTransformFilter::VolumeTransform(){
    if(this->tetmodel==NULL){return;}
    
    tetmodelPose = ApplyTransform<vtkUnstructuredGrid>(this->tetmodel);
        
    //boneposedistfunc = PBS_GetBoneDistFunc();
    
    // refine deformation 
    //SolveConstraints();
}

int vtkRigidSurfaceTransformFilter::PBS_GetMaxElement(vtkSmartPointer<vtkIntArray> array, int comp){
    int max=0;
    for(int i=0; i<array->GetNumberOfTuples(); i++){
        max = std::max(max,int(array->GetComponent(i,comp)));
    }
    return max;
}

std::vector<vtkSmartPointer<vtkImplicitPolyDataDistance>> vtkRigidSurfaceTransformFilter::PBS_GetBoneDistFunc(){
    int Nbones = this->armature->GetNumberOfVertices(); 
    std::vector<vtkSmartPointer<vtkImplicitPolyDataDistance>> bonedistfunc(Nbones);
    
    for(int b=0; b<Nbones; b++){
        bonedistfunc[b] = vtkSmartPointer<vtkImplicitPolyDataDistance>::New(); 
        bonedistfunc[b]->SetInput(vtkPolyData::SafeDownCast(this->surfacePose->GetBlock(b)));
    }
    
    return bonedistfunc;
}

//----------------------------------------------------------------------------
void vtkRigidSurfaceTransformFilter::PBS_SetReferenceInfo(){
    vtkSmartPointer<vtkIntArray> edgeinfo_array = vtkIntArray::SafeDownCast( this->tetmodel->GetFieldData()->GetArray("EdgeInfo") );
    vtkSmartPointer<vtkIntArray> cellinfo_array = vtkIntArray::SafeDownCast( this->tetmodel->GetCellData()->GetArray("Color") );
    vtkSmartPointer<vtkIntArray> bone_array = vtkIntArray::SafeDownCast( this->tetmodel->GetCellData()->GetArray("Color") );
    vtkSmartPointer<vtkIntArray> boneid_cells = vtkIntArray::SafeDownCast( this->tetmodel->GetCellData()->GetArray("BoneId") );
    vtkSmartPointer<vtkIntArray> boneid_points = vtkIntArray::SafeDownCast( this->tetmodel->GetPointData()->GetArray("BoneId") );
    vtkSmartPointer<vtkIntArray> closestbone_array = vtkIntArray::SafeDownCast( this->tetmodel->GetPointData()->GetArray("ClosestBone") );
    
    if(edgeinfo_array == NULL || cellinfo_array == NULL ){throw std::runtime_error("Volumetric model is missing edge and/or cell info.");}
    
    int NumberOfEdgeColors = PBS_GetMaxElement(edgeinfo_array,2) + 1;
    int NumberOfCellColors = PBS_GetMaxElement(cellinfo_array,0) + 1;
    int NumberOfThreads = std::thread::hardware_concurrency();
    
    edgeinfo = std::vector<std::vector<CellStructList>>( NumberOfEdgeColors );
    cellinfo = std::vector<std::vector<CellStructList>>( NumberOfCellColors );
    bindinfo = std::vector<CellStructList>(NumberOfThreads); 
    
    for(int edgecolor=0; edgecolor<NumberOfEdgeColors; edgecolor++){edgeinfo[edgecolor]=std::vector<CellStructList>(NumberOfThreads);}
    for(int cellcolor=0; cellcolor<NumberOfCellColors; cellcolor++){cellinfo[cellcolor]=std::vector<CellStructList>(NumberOfThreads);}
    
    std::vector<int> LastThreadPerCellColor(NumberOfCellColors,0);
    std::vector<int> LastThreadPerEdgeColor(NumberOfEdgeColors,0);
    
    
    totalVolume = 0; 
    modelvolume = std::vector<double>(tetmodel->GetNumberOfCells());
    
    for(int cell_id=0; cell_id<cellinfo_array->GetNumberOfTuples(); cell_id++){
        
        //if(boneid_cells->GetComponent(cell_id,0)==0){ // ie cell belongs to soft tissue
        
        int c = cellinfo_array->GetComponent(cell_id,0);
        
        vtkSmartPointer<vtkIdList> cellpoints = vtkSmartPointer<vtkIdList>::New();
        tetmodel->GetCellPoints(cell_id, cellpoints);
        
        CellStruct newcell;
        
        newcell.cellid = cell_id;
        
        double p1[3]; tetmodel->GetPoint(cellpoints->GetId(0), p1);
        double p2[3]; tetmodel->GetPoint(cellpoints->GetId(1), p2);
        double p3[3]; tetmodel->GetPoint(cellpoints->GetId(2), p3);
        double p4[3]; tetmodel->GetPoint(cellpoints->GetId(3), p4);
        
        double edge21[3]; vtkMath::Subtract(p2,p1,edge21); vtkMath::Normalize(edge21);
        double edge31[3]; vtkMath::Subtract(p3,p1,edge31); vtkMath::Normalize(edge31);
        double edge41[3]; vtkMath::Subtract(p4,p1,edge41); vtkMath::Normalize(edge41);
        vtkMath::Cross(edge21, edge31, edge31);
        if(vtkMath::Dot(edge31, edge41)>0){
            newcell.pointids->InsertNextId(cellpoints->GetId(0));  
            newcell.pointids->InsertNextId(cellpoints->GetId(1));  
            newcell.pointids->InsertNextId(cellpoints->GetId(2));  
            newcell.pointids->InsertNextId(cellpoints->GetId(3));  
            
            vtkSmartPointer<vtkTetra> tetra = vtkSmartPointer<vtkTetra>::New();
            newcell.vol = tetra->ComputeVolume(p1,p2,p3,p4);
        }
        else{
            newcell.pointids->InsertNextId(cellpoints->GetId(0));  
            newcell.pointids->InsertNextId(cellpoints->GetId(2));  
            newcell.pointids->InsertNextId(cellpoints->GetId(1));  
            newcell.pointids->InsertNextId(cellpoints->GetId(3));  
            
            vtkSmartPointer<vtkTetra> tetra = vtkSmartPointer<vtkTetra>::New();
            newcell.vol = tetra->ComputeVolume(p1,p3,p2,p4);
        }
        
        totalVolume += std::abs(newcell.vol);
        
        cellinfo[c][LastThreadPerCellColor[c]].push_back(newcell);
        
        LastThreadPerCellColor[c] = (LastThreadPerCellColor[c]+1)%NumberOfThreads;
        
        //}
    }
    
    for(int edge_id=0; edge_id<edgeinfo_array->GetNumberOfTuples(); edge_id++){
        int id1 = edgeinfo_array->GetComponent(edge_id,0);
        int id2 = edgeinfo_array->GetComponent(edge_id,1);
        int c  = edgeinfo_array->GetComponent(edge_id,2);
        
        //if(boneid_points->GetComponent(id1,0)==0 || boneid_points->GetComponent(id2,0)==0){ // ie both points belongs to soft tissue
        
        CellStruct newcell;
        
        newcell.cellid = edge_id;
        
        newcell.pointids->InsertNextId(id1);
        newcell.pointids->InsertNextId(id2);  
        
        double p1[3]; tetmodel->GetPoint(id1, p1);
        double p2[3]; tetmodel->GetPoint(id2, p2);
        
        double dist = std::pow(vtkMath::Distance2BetweenPoints(p1,p2),0.5);
        
        newcell.vol = dist; 
        
        edgeinfo[c][LastThreadPerEdgeColor[c]].push_back(newcell);
        
        LastThreadPerEdgeColor[c] = (LastThreadPerEdgeColor[c]+1)%NumberOfThreads;
        //}
    }
    
      
    std::vector<vtkSmartPointer<vtkImplicitPolyDataDistance>> bonedistfunc = PBS_GetBoneDistFunc();
      
    vtkSmartPointer<vtkDoubleArray> distarray = vtkSmartPointer<vtkDoubleArray>::New();
    distarray->SetNumberOfComponents(1);
    distarray->SetNumberOfTuples(tetmodel->GetNumberOfPoints());
    distarray->SetName("DistArray");
    
    for(int p=0; p<tetmodel->GetNumberOfPoints(); p++){
        double point[3]; tetmodel->GetPoint(p, point);
        int boneid = closestbone_array->GetComponent(p,0);
        
        double dist = bonedistfunc[boneid]->FunctionValue(point);
        
        std::cout << "Point " << p << ", boneid " << boneid << ", dist " << dist << std::endl;
        
        CellStruct newcell;
        newcell.cellid = p;
        newcell.vol = dist;
        
        bindinfo[0].push_back(newcell);
        
        distarray->SetComponent(p,0,dist);
    }
    
    this->tetmodel->GetPointData()->AddArray(distarray);
    
    vtkSmartPointer<vtkXMLUnstructuredGridWriter> wr = vtkSmartPointer<vtkXMLUnstructuredGridWriter>::New();
    wr->SetInputData(tetmodel);
    wr->SetFileName("tetmodel_pose_.vtu");
    wr->Update();
    
    
}


//----------------------------------------------------------------------------
void vtkRigidSurfaceTransformFilter::SolveBindConstraint(CellStructList celllist){
    double ks = 0.8; 
    
    vtkSmartPointer<vtkIntArray> pointboneid_array = vtkIntArray::SafeDownCast( tetmodel->GetPointData()->GetArray("BoneId") );
    vtkSmartPointer<vtkIntArray> closestbone_array = vtkIntArray::SafeDownCast( this->tetmodel->GetPointData()->GetArray("ClosestBone") );
    
    bool only_softtissue = true; // modify only soft tissue vertices    
    
    for(int id=0; id< celllist.size(); id++){
        int point_idx = celllist[id].cellid;
        double dist   = celllist[id].vol;
        
        double point[3]; this->tetmodelPose->GetPoint(point_idx, point);
        int boneid = closestbone_array->GetComponent(point_idx,0);
        
        double closestpoint[3];
        double newdist = boneposedistfunc[boneid]->EvaluateFunctionAndGetClosestPoint(point, closestpoint);
        double diff[3]; vtkMath::Subtract(closestpoint, point, diff); vtkMath::Normalize(diff);
        
        //boneposedistfunc[boneid]->FunctionGradient(point, diff); vtkMath::Normalize(diff);
        //vtkMath::MultiplyScalar(diff,-1);
        
        
        
        std::cout << "eval distfunc " << dist << "\t" << newdist << std::endl; 
        
        //if(std::abs((newdist-dist)/dist)>0.01){
            double fac;
            if(newdist*dist<0){ fac = 0.5 * ks*(std::abs(newdist)+std::abs(dist)); } // two points lie on opposite sides of the surface
            else{ fac = 0.5 * ks*(std::abs(newdist)-std::abs(dist)); } // two points lie at the same side of the surface
            vtkMath::MultiplyScalar(diff, fac);
            
            double diffpoint[3]; tetmodelPose->GetPoint(point_idx, diffpoint);
            
            vtkMath::Add(diffpoint, diff, diffpoint);
            
            if((pointboneid_array->GetComponent(point_idx,0)==0 && only_softtissue) || !only_softtissue){ tetmodelPose->GetPoints()->SetPoint(point_idx, diffpoint); }
        //}
    }
}

//----------------------------------------------------------------------------
void vtkRigidSurfaceTransformFilter::SolveEdgeConstraint(CellStructList celllist){
    double ks = 0.8; //0.4
    
    vtkSmartPointer<vtkIntArray> pointboneid_array = vtkIntArray::SafeDownCast( tetmodel->GetPointData()->GetArray("BoneId") );
    
    bool only_softtissue = true; // modify only soft tissue vertices    
    
    for(int id=0; id< celllist.size(); id++){
        int p1_idx = celllist[id].pointids->GetId(0);
        int p2_idx = celllist[id].pointids->GetId(1);
        double dist   = celllist[id].vol;
        
        double p1[3]; this->tetmodelPose->GetPoint(p1_idx, p1);
        double p2[3]; this->tetmodelPose->GetPoint(p2_idx, p2);
        
        double diff[3]; 
        vtkMath::Subtract(p1,p2,diff); 
        double newdist = vtkMath::Normalize(diff);
        
        if(std::abs(newdist-dist)/dist>0.01){
            double fac = 0.5*ks*(newdist-dist);
            vtkMath::MultiplyScalar(diff, fac);
            
            double diff1[3]; tetmodelPose->GetPoint(p1_idx, diff1);
            double diff2[3]; tetmodelPose->GetPoint(p2_idx, diff2);
            
            vtkMath::Subtract(diff1, diff, diff1);      // tetmodelPose->GetPoints()->SetPoint(p1_idx, diff1);
            vtkMath::Add(diff2, diff, diff2);           // tetmodelPose->GetPoints()->SetPoint(p2_idx, diff2);
            
            if((pointboneid_array->GetComponent(p1_idx,0)==0 && only_softtissue) || !only_softtissue){ tetmodelPose->GetPoints()->SetPoint(p1_idx, diff1); }
            if((pointboneid_array->GetComponent(p2_idx,0)==0 && only_softtissue) || !only_softtissue){ tetmodelPose->GetPoints()->SetPoint(p2_idx, diff2); }
        }
    }
}

//----------------------------------------------------------------------------
void vtkRigidSurfaceTransformFilter::SolveVolumeConstraint(CellStructList celllist){
    double kvol = 0.8; //0.4
    
    vtkSmartPointer<vtkIntArray> cellboneid_array = vtkIntArray::SafeDownCast( tetmodel->GetCellData()->GetArray("BoneId") );
    vtkSmartPointer<vtkIntArray> pointboneid_array = vtkIntArray::SafeDownCast( tetmodel->GetPointData()->GetArray("BoneId") );
    
    bool only_softtissue = true; // modify only soft tissue vertices
    
    for(int id=0; id<celllist.size(); id++){
        
        if((cellboneid_array->GetComponent(celllist[id].cellid,0)==0 && only_softtissue) || !only_softtissue){
        
        int p1_idx = celllist[id].pointids->GetId(0);
        int p2_idx = celllist[id].pointids->GetId(1);
        int p3_idx = celllist[id].pointids->GetId(2);
        int p4_idx = celllist[id].pointids->GetId(3);
        double vol   = celllist[id].vol;
                
        double p1[3]; this->tetmodelPose->GetPoint(p1_idx, p1);
        double p2[3]; this->tetmodelPose->GetPoint(p2_idx, p2);
        double p3[3]; this->tetmodelPose->GetPoint(p3_idx, p3);
        double p4[3]; this->tetmodelPose->GetPoint(p4_idx, p4);
        
        vtkSmartPointer<vtkTetra> tetra = vtkSmartPointer<vtkTetra>::New();
        double newvol = tetra->ComputeVolume(p1,p2,p3,p4);
        
        modelvolume[celllist[id].cellid] = newvol;
        
        if(std::abs(std::abs(newvol)-std::abs(vol))/std::abs(vol)>0){
        
        //std::cout <<celllist[id].cellid << " " << newvol << " " << vol << std::endl;
            
                
        double edge21[3]; vtkMath::Subtract(p2,p1,edge21);
        double edge31[3]; vtkMath::Subtract(p3,p1,edge31);
        double edge41[3]; vtkMath::Subtract(p4,p1,edge41);
        
        double norm21 = vtkMath::Normalize(edge21);
        double norm31 = vtkMath::Normalize(edge31);
        double norm41 = vtkMath::Normalize(edge41);
        
        double grad2[3]; vtkMath::Cross(edge31,edge41, grad2); vtkMath::MultiplyScalar(grad2, norm31*norm41/6);
        double grad3[3]; vtkMath::Cross(edge41,edge21, grad3); vtkMath::MultiplyScalar(grad3, norm41*norm21/6);
        double grad4[3]; vtkMath::Cross(edge21,edge31, grad4); vtkMath::MultiplyScalar(grad4, norm21*norm31/6);
        double grad1[3]; vtkMath::Add(grad2, grad3, grad1); 
        vtkMath::Add(grad1, grad4, grad1); vtkMath::MultiplyScalar(grad1, -1);
        double norm = std::pow(vtkMath::Norm(grad1),2) + std::pow(vtkMath::Norm(grad2),2) + std::pow(vtkMath::Norm(grad3),2) + std::pow(vtkMath::Norm(grad4),2);
        double s = -(std::abs(newvol) - std::abs(vol))/norm;
        
        
        //if(newvol<0){s = -(std::abs(newvol) - std::abs(vol))/norm;}
        
        
        
        //std::cout << newvol << " " << vol << std::endl; 
        
        //if(newvol>0){ // stop it from escalating
        
        double diff1[3]; tetmodelPose->GetPoint(p1_idx, diff1);
        double diff2[3]; tetmodelPose->GetPoint(p2_idx, diff2);
        double diff3[3]; tetmodelPose->GetPoint(p3_idx, diff3);
        double diff4[3]; tetmodelPose->GetPoint(p4_idx, diff4);        
        
        vtkMath::MultiplyScalar(grad1, s*kvol);
        vtkMath::MultiplyScalar(grad2, s*kvol);
        vtkMath::MultiplyScalar(grad3, s*kvol);
        vtkMath::MultiplyScalar(grad4, s*kvol);
        
        vtkMath::Add(diff1, grad1, diff1);
        vtkMath::Add(diff2, grad2, diff2);
        vtkMath::Add(diff3, grad3, diff3);
        vtkMath::Add(diff4, grad4, diff4);
        
        if((pointboneid_array->GetComponent(p1_idx,0)==0 && only_softtissue) || !only_softtissue){ tetmodelPose->GetPoints()->SetPoint(p1_idx, diff1); }
        if((pointboneid_array->GetComponent(p2_idx,0)==0 && only_softtissue) || !only_softtissue){ tetmodelPose->GetPoints()->SetPoint(p2_idx, diff2); }
        if((pointboneid_array->GetComponent(p3_idx,0)==0 && only_softtissue) || !only_softtissue){ tetmodelPose->GetPoints()->SetPoint(p3_idx, diff3); }
        if((pointboneid_array->GetComponent(p4_idx,0)==0 && only_softtissue) || !only_softtissue){ tetmodelPose->GetPoints()->SetPoint(p4_idx, diff4); }
        
        //tetmodelPose->GetPoints()->SetPoint(p1_idx, diff1);
        //tetmodelPose->GetPoints()->SetPoint(p2_idx, diff2);
        //tetmodelPose->GetPoints()->SetPoint(p3_idx, diff3);
        //tetmodelPose->GetPoints()->SetPoint(p4_idx, diff4);
        
        //}
        }
        }
    }
}

//----------------------------------------------------------------------------
void vtkRigidSurfaceTransformFilter::SolveConstraints(){
    
    //unsigned maxThreads = std::thread::hardware_concurrency();
    
    int MaxNumberOfIterations = 100;
    int it = 0;
    double newvol = 0;
    
    std::ofstream myfile;
    myfile.open ("pbs_volume.txt");
    
    while(it<MaxNumberOfIterations && std::abs(newvol-totalVolume)/totalVolume > 0.001){ // 0.001
        
        std::cout << "Iteration "<< it << std::endl;
        
        // Solve edge constraint 
        for(int edgecolor=0; edgecolor<edgeinfo.size(); edgecolor++){
            std::vector<std::thread> threads;
            
            for (int threadid = 0; threadid < edgeinfo[edgecolor].size(); threadid++){
                threads.push_back(std::thread(&vtkRigidSurfaceTransformFilter::SolveEdgeConstraint, this, edgeinfo[edgecolor][threadid]));
            }
            
            for (auto& t : threads)
                t.join();
        }
        
        // Solve volume constraint 
        for(int volcolor=0; volcolor<cellinfo.size(); volcolor++){
            std::vector<std::thread> threads;
                        
            for (int threadid = 0; threadid < cellinfo[volcolor].size(); threadid++){
                threads.push_back(std::thread(&vtkRigidSurfaceTransformFilter::SolveVolumeConstraint, this, cellinfo[volcolor][threadid]));
            }
            
            for (auto& t : threads)
                t.join();
        }
        
        // Solve bind constraint 
        std::vector<std::thread> threads;
                    
        for (int threadid = 0; threadid < 1; threadid++){
            threads.push_back(std::thread(&vtkRigidSurfaceTransformFilter::SolveBindConstraint, this, bindinfo[threadid]));
        }
        
        for (auto& t : threads)
            t.join();
        
        
        // Calculate new model volume 
        newvol=0;
        for(int i=0; i<tetmodel->GetNumberOfCells(); i++){newvol+=modelvolume[i];}
        
        std::cout << "Current volume " << newvol << "\t" << totalVolume << "\t" << std::abs(newvol-totalVolume)/totalVolume << std::endl;
        
        myfile << newvol << "\n";
        
        vtkSmartPointer<vtkXMLUnstructuredGridWriter> wr = vtkSmartPointer<vtkXMLUnstructuredGridWriter>::New();
        wr->SetInputData(tetmodelPose);
        wr->SetFileName(("tetmodel_pose_"+std::to_string(it)+".vtu").c_str());
        wr->Update();
        
        it++;
    }
    
    myfile.close();
}
    
//----------------------------------------------------------------------------
void vtkRigidSurfaceTransformFilter::SurfaceTransform(){
    for(vtkIdType j = 0; j < surfacePose->GetNumberOfBlocks(); j++)
    {        
        if(IsBoneModel(surfaceRest->GetBlock(j))==true){              // block belongs to a bone model   
            vtkSmartPointer<vtkPolyData> newsurf = vtkSmartPointer<vtkPolyData>::New();
            newsurf->DeepCopy( vtkPolyData::SafeDownCast(surfaceRest->GetBlock(j)) );
            surfacePose->SetBlock( j, RigidTransform(newsurf, j)); 
        }
        else{    // block belongs to the surface model 
            if(surfaceRest->GetBlock(j)->IsA("vtkPolyData")){
                surfacePose->SetBlock( j, RecalculateNormals(ApplyTransform<vtkPolyData>(surfaceRest->GetBlock(j)))); 
            }
            if(surfaceRest->GetBlock(j)->IsA("vtkUnstructuredGrid")){
                std::cout << "ugrid in multiblock" << std::endl;
                surfacePose->SetBlock( j, ApplyTransform<vtkUnstructuredGrid>(surfaceRest->GetBlock(j))); 
            }
        }
    }
}

//----------------------------------------------------------------------------
// Specify a source object at a specified table location.
void vtkRigidSurfaceTransformFilter::SetArmature(int id, vtkSkeleton *pd)
{
  int numConnections = this->GetNumberOfInputConnections(0);
  
  if (id < 0 || id > numConnections)
  {
    vtkErrorMacro("Bad index " << id << " for source.");
    return;
  }

  vtkTrivialProducer* tp = 0;
  if (pd)
  {
    tp = vtkTrivialProducer::New();
    tp->SetOutput(pd);
  }

  if (id < numConnections)
  {
    if (tp)
    {
      this->SetNthInputConnection(1, id, tp->GetOutputPort());
    }
    else
    {
      this->SetNthInputConnection(1, id, 0);
    }
  }
  else if (id == numConnections && tp)
  {
    this->AddInputConnection(1, tp->GetOutputPort());
  }

  if (tp)
  {
    tp->Delete();
  }
}

