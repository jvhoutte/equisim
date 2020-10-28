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

#include "SKELETON/vtkBVHReader.h"


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


vtkStandardNewMacro(vtkBVHReader);

//----------------------------------------------------------------------------
vtkBVHReader::vtkBVHReader()
{
  this->FileName = "";
  this->Armature = NULL;
  this->SetNumberOfInputPorts(0);
}

//----------------------------------------------------------------------------
vtkBVHReader::~vtkBVHReader()
{
  /* Use vtkSmartPointer when using vtkBVHReader
  if (this->Armature)
    {
    this->Armature->Delete();
    }
  */
}

//----------------------------------------------------------------------------
void vtkBVHReader::SetFileName(const char* filename)
{
  if (!filename && this->FileName == ""
    || filename && strcmp(this->FileName.c_str(), filename) == 0)
    {
    return;
    }

  this->FileName = filename ? filename : "";
  this->Modified();
}

//----------------------------------------------------------------------------
const char* vtkBVHReader::GetFileName() const
{
  return this->FileName.c_str();
}

//----------------------------------------------------------------------------
int vtkBVHReader::RequestData(vtkInformation *vtkNotUsed(request),
                              vtkInformationVector **vtkNotUsed(inputVector),
                              vtkInformationVector *outputVector)
{

  // Create new armature if necessary
  if (!this->Armature)
    {
    this->Armature = vtkSkeleton::New();
    if (this->FileName == "")
      {
      vtkErrorMacro("A file name must be specified.");
      return 0;
      }
    }

  // open a BVH file for reading
   
    std::ifstream file(this->FileName.c_str());
    if (!file.good())
        {
        vtkErrorMacro("Cannot open the given file.");
        return 0;
        }

    if (!this->Parse(file))
        {
        vtkErrorMacro("Error when parsing the file.");
        return 0;
        }
    file.close();

  
    
    
  //outputVector->SetNumberOfInformationObjects(1);
 
  //vtkInformation* outInfo0 = outputVector->GetInformationObject(0);
  //outInfo0->Get(vtkDataObject::DATA_OBJECT())->DeepCopy(this->Armature->GetPolyData());

    
  return 1;
}

//----------------------------------------------------------------------------
int vtkBVHReader::Parse(std::ifstream& file)
{
  // Make sure first line is hierarchy
  std::string line = "";
  std::getline(file, line);
  if (MoveToNextKeyword(file, line) != "HIERARCHY")
    {
    std::cerr<<"Invalid BVH file, no hierarchy was specified." << std::endl;
    return 0;
    }

  this->ParseRestArmature(file, -1);
  
  
  
  // Calculate other basic information based on the parsed data
  if( this->Armature->ParseArmature()){
      this->Armature->AddArmatureArrays(); 
      this->Armature->RebuildTransformations();
  }
  else{ throw std::runtime_error("Incompatible graph.");}
  
  
  
  
  return 1;
}

//----------------------------------------------------------------------------
void vtkBVHReader::ParseRestArmature(std::ifstream& file, vtkIdType parentId)
{
    vtkIdType boneid;
  std::string line;
  while (std::getline(file, line))
    {
    std::string keyword = MoveToNextKeyword(file, line);
    
    if (keyword == "}")
      {
      return;
      }
    else if (keyword == "ROLLAXIS")
      {
        double rollaxis[3];
        GetRollAxis(line, rollaxis);
        this->Armature->SetZAxis(rollaxis);
      }
    else if (keyword == "RANGE")
    {
        double range[8];
        GetRange(line, range);
        this->Armature->SetAngleRanges(range);
        
    }
    else if (keyword == "HEAD")
      {
        double head[3];
        GetHead(line, head);
        this->Armature->SetWorldHead(head);
      }
    else if (keyword == "TAIL")
      {
        double tail[3];
        GetTail(line, tail);
        this->Armature->SetWorldTail(tail);
      }
    else if (keyword == "ROOT")
      {
        boneid = this->Armature->AddNodeToGraph(parentId);  
        
        std::string name = GetBoneName(line, "ROOT");
        this->Armature->SetBoneName(name);
      }
    else if (keyword == "JOINT")
        {
        if (parentId == -1) // first bone after root
        {      
            parentId = boneid;
        }

        if(!this->Armature->CheckCompleteness()){throw std::runtime_error("Missing information in the BVH file.");}
        boneid = this->Armature->AddNodeToGraph(parentId);  

        std::string name = GetBoneName(line);
        this->Armature->SetBoneName(name);
    
        ParseRestArmature(file, boneid);

        }
    }
}


//----------------------------------------------------------------------------
int vtkBVHReader::CanReadFile(const char *filename)
{
  std::ifstream file(filename);
  if (!file.good())
    {
    return 0;
    }

  std::string name = filename;
  if (name.size() < 3 || name.substr(name.size() - 3) != "bvh")
    {
    return 0;
    }

  return 1;
}


//----------------------------------------------------------------------------
void vtkBVHReader::PrintSelf(ostream& os, vtkIndent indent)
{
  this->Superclass::PrintSelf(os,indent);

  os << indent << "File Name: " << this->FileName << "\n";
}
