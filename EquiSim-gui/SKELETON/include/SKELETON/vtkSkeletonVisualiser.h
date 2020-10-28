
#ifndef __vtkSkeletonVisualiser_h
#define __vtkSkeletonVisualiser_h
 
#include <numeric>
#include "vtkDataObjectAlgorithm.h"
#include "vtkTrivialProducer.h"
#include "vtkPolyDataNormals.h"
#include "vtkMultiBlockDataSet.h"
#include "DLLDefines.h"  
#include "SKELETON/vtkSkeleton.h" 


// For compatibility with new VTK generic data arrays
#ifdef vtkGenericDataArray_h
#define InsertNextTupleValue InsertNextTypedTuple
#endif

class vtkSkeletonVisualiser : public vtkDataObjectAlgorithm 
{
public:
  vtkTypeMacro(vtkSkeletonVisualiser,vtkDataObjectAlgorithm);
  void PrintSelf(ostream& os, vtkIndent indent);
 
  static vtkSkeletonVisualiser *New();
  void SetArmature(int id, vtkSkeleton *pd);

  
protected:
  vtkSkeletonVisualiser();
  ~vtkSkeletonVisualiser();
 
  int FillInputPortInformation( int port, vtkInformation* info ) VTK_OVERRIDE;
  int FillOutputPortInformation(int portNumber, vtkInformation *info) VTK_OVERRIDE;
  
  int RequestData(vtkInformation *, vtkInformationVector **, vtkInformationVector *);
  
  vtkSmartPointer<vtkPoints> GetSkeletonGeometry();
  vtkSmartPointer<vtkUnsignedCharArray> GetSkeletonColors();
  vtkSmartPointer<vtkCellArray> GetSkeletonTopology();
  
private:
  vtkSkeletonVisualiser(const vtkSkeletonVisualiser&);  // Not implemented.
  void operator=(const vtkSkeletonVisualiser&);  // Not implemented.

  
  vtkTree* armature;
  vtkPolyData* polyskeleton;
  int NUMBER_OF_BONES;
};

#endif
