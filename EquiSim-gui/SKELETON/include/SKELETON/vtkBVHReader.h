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

// .NAME vtkBVHReader - Reads BVH files.
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
// vtkBVHReader, vtkArmatureWidget,

#ifndef __vtkBVHReader_h
#define __vtkBVHReader_h

#include <vtkPolyDataAlgorithm.h>

#include "SKELETON/vtkSkeleton.h"

#include "vtkQuaternion.h"
#include "DLLDefines.h" 
#include <vector>

class vtkBVHReader : public vtkPolyDataAlgorithm
{
public:

  vtkTypeMacro(vtkBVHReader, vtkPolyDataAlgorithm);
  void PrintSelf(ostream& os, vtkIndent indent);

  // Description:
  static vtkBVHReader *New();

  // Description:
  // Set the motion capture file's filename to read.
  // Setting a new filename invalids the current armature (if any).
  void SetFileName(const char* filename);
  const char* GetFileName() const;


  // Description:
  // When linking to the first child, the first child of a bone will
  // always start from its parent tail.
  // When this option is off, if the parent has multiple child, the
  // parent's tail position will be given by the average position
  // of its children's head.
  // Default is false.
  //void SetLinkToFirstChild(bool link);
  //vtkGetMacro(LinkToFirstChild, bool);

  // Description:
  // Get the armature from which the polydata is obtained.
  vtkGetObjectMacro(Armature, vtkSkeleton);

  // Description:
  // A simple, non-exhaustive check to see if a file is a valid file.
  static int CanReadFile(const char *filename);

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


protected:
  vtkBVHReader();
  ~vtkBVHReader();

  std::string FileName;

  vtkSkeleton* Armature;
  //bool LinkToFirstChild;

  virtual int RequestData(
    vtkInformation*, vtkInformationVector**, vtkInformationVector*);

  int Parse(std::ifstream& file);
  void ParseMotions(
    std::ifstream& file, std::vector< std::vector<std::string> >& channels);
  void ParseRestArmature(
    std::ifstream& file,
    vtkIdType parentId);

  void LinkBonesToFirstChild();
  void UnlinkBonesFromFirstChild();
  
  void InvalidReader();

private:
  vtkBVHReader(const vtkBVHReader&);  // Not implemented.
  void operator=(const vtkBVHReader&);  // Not implemented.

};

#endif
