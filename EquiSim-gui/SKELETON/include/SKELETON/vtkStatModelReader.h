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

// .NAME vtkStatModelReader - Reads BVH files.
// .SECTION Description
// Reads a BVH (motion capture) file and outputs the corresponding polydata.
//
// Using an armature, the reader creates the rest position of the armature
// from the HIERARCHY part of the BVH. The reader assumes that there is
// only one root.
//
// Since an armature can only have one pose, the SetFrame property allow
// to choose from the different motion frame. The movement data is gathered
// under The MOTION part of the file. Upon reading, the animation information
// is stored. This provides faster look up when changing the frame.
//
// The polydata is obtained from the armature widget. To learn more about
// its structure, see the vtkArmatureWidget->GetPolyData().
//
// .SECTION See Also
// vtkStatModelReader, vtkArmatureWidget,

#ifndef __vtkStatModelReader_h
#define __vtkStatModelReader_h


#include "SKELETON/vtkSkeleton.h"
#include "SKELETON/vtkBVHReader.h"
#include "SKELETON/Defs.h"

#include "vtkMultiBlockDataSet.h"
#include "vtkQuaternion.h"
#include "DLLDefines.h" 
#include <vector>


// Convenient class to more easily access the different PC modes of an object. 
class Tensor{
  public:
    
    Tensor(){}
    
    void SetArrays(vtkSmartPointer<vtkPolyData> polydata){
      // maximum number of modes possible, the polydata can also contain arrays which to not contain the eigenmodes 
      Nmodes = polydata->GetPointData()->GetNumberOfArrays(); 
      for(unsigned int modeNr = 0; modeNr < Nmodes+1; ++modeNr)
      {
	  //std::ostringstream name;
	  //std::ostringstream num;
	  //num << std::setw( log10( Nmodes+1 ) + 1 ) << std::setfill( '0' ) << modeNr + 1;
	  //name << "Model_Eigenmode_" << num.str();
      //std::string name = name.str()    
          
      std::string name = "Model_Eigenmode_" + std::to_string(modeNr+1);
          
	  vtkSmartPointer<vtkFloatArray> array = vtkFloatArray::SafeDownCast( polydata->GetPointData()->GetArray(name.c_str()) );
            
	  if(array != NULL){
	    if(array->GetNumberOfTuples() != polydata->GetNumberOfPoints()){
            throw std::runtime_error("Difference detected between the number of PC-vectors and number of modelpoints."); 
        }
	    modelEigenmodes.push_back( array );
	  }
      }
      Nmodes = modelEigenmodes.size();
    }
  
    
  
  
    int GetNumberOfModes(){return Nmodes;}
    
  
    vtkSmartPointer<vtkFloatArray> operator[](int idx) { return modelEigenmodes[idx];}
    
    std::vector< vtkSmartPointer<vtkFloatArray> > modelEigenmodes;
  
  private:
    int Nmodes=0;
  
};

class PGAinfo{
  public:
    PGAinfo(){}
    
    void SetArrays(vtkSmartPointer<vtkPolyData> polydata){
        std::string arrayname = "pgamean";
        vtkSmartPointer<vtkDoubleArray> array = vtkDoubleArray::SafeDownCast( polydata->GetFieldData()->GetArray(arrayname.c_str()) );
        if(array != NULL){
            if(array->GetNumberOfTuples() != 1 || array->GetNumberOfComponents() != 3){
                throw std::runtime_error("Size of pga mean array is different than expected"); 
            }   
            pgamean = vtkSmartPointer<vtkPoints>::New();
            double mean[3];
            array->GetTuple(0,mean);
            pgamean->InsertNextPoint(mean);
        }
        
    }
    
    bool GetPgaMean(double (*mean)[3]){
        if(pgamean != NULL){
            pgamean->GetPoint(0,*mean);
            return true;
        }
        else{
            return false;
        }
    }
    
  private:
    vtkSmartPointer<vtkPoints> pgamean;    
};


struct PcModel{
  vtkSmartPointer<vtkPolyData> average;
  Tensor eigenvectors = Tensor(); 
  PGAinfo pga = PGAinfo(); 
};

struct BlockModel{
  PcModel geometry;
  PcModel axis;
  PcModel cor;  
  std::string name = "";
};



typedef std::vector<BlockModel> TensorBlocks;

class vtkStatModelReader : public vtkPolyDataAlgorithm
{
public:

  vtkTypeMacro(vtkStatModelReader, vtkPolyDataAlgorithm);
  void PrintSelf(ostream& os, vtkIndent indent);

  // Description:
  static vtkStatModelReader *New();

  void SetBVHFileName(const char* filename);
  void SetVTMFileName(const char* filename);
  void SetSeparateFileNames(std::vector<std::string> filenames); // For reading a template subject
  
  void SetPcWeights(int mode, double weight);
  
  int GetNumberOfModes(){return Nmodes;}
  
  // Description:
  // Get the armature from which the polydata is obtained.
  vtkGetObjectMacro(Armature, vtkSkeleton);
  vtkGetObjectMacro(Shape, vtkMultiBlockDataSet);
  // Description:
  // A simple, non-exhaustive check to see if a file is a valid file.
  

  // Description:
  // Apply the frame to the given armature. The armature must be the
  // exact same armature than the armature read.
  // This method is meant for exterior application to be able to drive
  // which pose the armature has. Return if the operation succeeded.
  //bool ApplyFrameToArmature(vtkArmatureWidget* armature,  unsigned int frame);

  // Description:
  // Access method to the frame rotation data.
  // No check is performed on the frame nor the boneId for perfomance reasons.
  // vtkQuaterniond GetParentToBoneRotation(unsigned int frame, unsigned int boneId);


  vtkStatModelReader();
  ~vtkStatModelReader();

  
protected:
  
  static int CanReadFile(const char *filename, const char *extension);
  
  vtkSkeleton* Armature;
  vtkMultiBlockDataSet* Shape; // shape tranformed to world coord system, consistent with Armature
  vtkSmartPointer<vtkMultiBlockDataSet> Shape_local; // shape in local coord frame, as read from the vtm file
  
  TensorBlocks aSSM;
  
  virtual int RequestData(vtkInformation*, vtkInformationVector**, vtkInformationVector*);

  vtkSmartPointer<vtkFloatArray> GetCurrentVectorField(Tensor* tensor);
  
  vtkSmartPointer<vtkStringArray> GetBlockNameArray(vtkSmartPointer<vtkMultiBlockDataSet> block);
  std::string GetBlockName(vtkSmartPointer<vtkMultiBlockDataSet> block);
  
  int GetBoneIndex(int modelidx);
  int GetModelIndex(int boneidx);
  
  void UpdateBone();
  vtkSmartPointer<vtkPoints> ChangePCWeights(PcModel pcmodel);
  
  void ExpMap(double (*v)[3], double (*mean)[3]);
  
  std::vector<int> GetNonZeroModes();
  
  vtkSmartPointer<vtkMultiBlockDataSet> MergeNoBones(vtkSmartPointer<vtkMultiBlockDataSet> multiblock);
  bool IsStatModel(vtkSmartPointer<vtkMultiBlockDataSet> multiblock);
  
private:
  vtkStatModelReader(const vtkStatModelReader&);  // Not implemented.
  void operator=(const vtkStatModelReader&);  // Not implemented.

  std::vector<double> pcweights;
  
  std::map<int, int> modeltobonemap;
  std::map<int, int> bonetomodelmap;
  
  int Ncompounds;
  int Nmodes;
  
  bool statmode = false;
};

#endif
