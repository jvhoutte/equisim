#ifndef __vtkBlockColorSelector_h
#define __vtkBlockColorSelector_h
 
#include <numeric>
#include "vtkDataObjectAlgorithm.h"
#include "vtkTrivialProducer.h"
#include "vtkPolyDataNormals.h"
#include "vtkMultiBlockDataSet.h"
#include "SKELETON/Defs.h"
#include "DLLDefines.h"  

class vtkBlockColorSelector : public vtkDataObjectAlgorithm 
{
public:
  vtkTypeMacro(vtkBlockColorSelector,vtkDataObjectAlgorithm);
  void PrintSelf(ostream& os, vtkIndent indent);
 
  static vtkBlockColorSelector *New();
  void SetColorId(int c){array_id = c;}

  
protected:
  vtkBlockColorSelector();
  ~vtkBlockColorSelector();
 
  int FillInputPortInformation( int port, vtkInformation* info ) VTK_OVERRIDE;
  int FillOutputPortInformation(int portNumber, vtkInformation *info) VTK_OVERRIDE;
  
  int RequestData(vtkInformation *, vtkInformationVector **, vtkInformationVector *);
  
  
private:
    vtkBlockColorSelector(const vtkBlockColorSelector&); // not implemented
    void operator=(const vtkBlockColorSelector&);  // Not implemented.
    
    void SetColor();
    int array_id = 1;
    
    
    vtkMultiBlockDataSet* surfaceRest;
    vtkMultiBlockDataSet* surfacePose;
};

#endif
