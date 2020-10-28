#ifndef __vtkBlockColorSelector_cpp
#define __vtkBlockColorSelector_cpp

#include "SKELETON/vtkBlockColorSelector.h"

vtkStandardNewMacro(vtkBlockColorSelector);

vtkBlockColorSelector::vtkBlockColorSelector()
{
  this->SetNumberOfInputPorts(1);
  //this->SetNumberOfOutputPorts(1);
}
 
vtkBlockColorSelector::~vtkBlockColorSelector()
{
}
 
int vtkBlockColorSelector::FillInputPortInformation( int port, vtkInformation* info )
{
  if ( port == 0 )
    {
    info->Set(vtkDataObject::DATA_TYPE_NAME(), "vtkMultiBlockDataSet" );
    return 1;
    }
  return 0;
}
 
 
int vtkBlockColorSelector::FillOutputPortInformation(int vtkNotUsed(portNumber), vtkInformation *info)
{
    info->Set(vtkDataObject::DATA_TYPE_NAME(), "vtkMultiBlockDataSet");
    return 1;
}
 
int vtkBlockColorSelector::RequestData(vtkInformation *vtkNotUsed(request),
                                          vtkInformationVector **inputVector,
                                          vtkInformationVector *outputVector)
{
    // get the input object
    this->surfaceRest = vtkMultiBlockDataSet::GetData(inputVector[0],0);
    
    // get the output object
    vtkInformation *outInfo = outputVector->GetInformationObject(0);
    this->surfacePose = vtkMultiBlockDataSet::SafeDownCast(outInfo->Get(vtkDataObject::DATA_OBJECT()));
    this->surfacePose->DeepCopy(surfaceRest);
    
    SetColor();
    
    return 1;
}



void vtkBlockColorSelector::SetColor(){
    
    
    //std::ostringstream ss;
    //ss << std::setw(2) << std::setfill('0') << std::to_string(array_id);
    //std::string str = ss.str();
    
    for(int bidx=0; bidx<surfaceRest->GetNumberOfBlocks(); bidx++){
        vtkSmartPointer<vtkDataArray> magn_arr = vtkPolyData::SafeDownCast(surfaceRest->GetBlock(bidx))->GetPointData()->GetArray(("Norm_Model_Eigenmode_"+std::to_string(array_id)).c_str());
                
        if(magn_arr == NULL || array_id == -1){
            magn_arr = vtkSmartPointer<vtkDoubleArray>::New();
            magn_arr->SetNumberOfComponents(1);
            magn_arr->SetNumberOfTuples(vtkPolyData::SafeDownCast(surfaceRest->GetBlock(bidx))->GetNumberOfPoints());
            magn_arr->FillComponent(0,bidx);
        }
        
        vtkPolyData::SafeDownCast(surfacePose->GetBlock(bidx))->GetPointData()->SetScalars(magn_arr); 
    }
}

void vtkBlockColorSelector::PrintSelf(ostream& os, vtkIndent indent)
{
  this->Superclass::PrintSelf(os,indent);
}

#endif
