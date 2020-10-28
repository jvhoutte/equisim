#ifndef __vtkTestPolyDataRigidFilter_h
#define __vtkTestPolyDataRigidFilter_h
 
#include <numeric>
#include "vtkDataObjectAlgorithm.h"
#include "vtkTrivialProducer.h"
#include "vtkPolyDataNormals.h"
#include "vtkMultiBlockDataSet.h"
#include "DLLDefines.h"  
#include "SKELETON/vtkSkeleton.h" 

struct CellStruct{
    CellStruct(){pointids = vtkSmartPointer<vtkIdList>::New();}

    vtkSmartPointer<vtkIdList> pointids;
    vtkIdType cellid;
    double vol;    
};
 
typedef std::vector<CellStruct> CellStructList;


class vtkRigidSurfaceTransformFilter : public vtkDataObjectAlgorithm 
{
public:
  vtkTypeMacro(vtkRigidSurfaceTransformFilter,vtkDataObjectAlgorithm);
  void PrintSelf(ostream& os, vtkIndent indent);
 
  static vtkRigidSurfaceTransformFilter *New();
  void SetArmature(int id, vtkSkeleton *pd);
  
  vtkSmartPointer<vtkTransform> GetBlockTransform(int blockid);
  
  vtkSmartPointer<vtkUnstructuredGrid> GetTetModel(){return tetmodelPose;}
  
protected:
  vtkRigidSurfaceTransformFilter();
  ~vtkRigidSurfaceTransformFilter();
 
  int FillInputPortInformation( int port, vtkInformation* info ) VTK_OVERRIDE;
  int FillOutputPortInformation(int portNumber, vtkInformation *info) VTK_OVERRIDE;
  
  int RequestData(vtkInformation *, vtkInformationVector **, vtkInformationVector *);
  
  void SetTransforms();
  //void SetJacobian();
  
  
private:
  vtkRigidSurfaceTransformFilter(const vtkRigidSurfaceTransformFilter&);  // Not implemented.
  void operator=(const vtkRigidSurfaceTransformFilter&);  // Not implemented.
  void TransformModels();
  void SurfaceTransform();
  void VolumeTransform();
  vtkSmartPointer<vtkPointSet> WeightedTransform(vtkSmartPointer<vtkPointSet> surf);
  vtkSmartPointer<vtkPolyData> RigidTransform(vtkSmartPointer<vtkPolyData> surf, int transform_id);
  template <typename modelType> vtkSmartPointer<modelType> ApplyTransform(vtkDataObject* surf);
  vtkSmartPointer<vtkPolyData> RecalculateNormals(vtkSmartPointer<vtkPolyData> surf);
  bool IsBoneModel(vtkDataObject* model);
  
  
  vtkIdType GetBoneIndex(vtkIdType modelidx);
  vtkMultiBlockDataSet* surfaceRest;
  vtkMultiBlockDataSet* surfacePose;
  vtkUnstructuredGrid* tetmodel;
  vtkSmartPointer<vtkUnstructuredGrid> tetmodelPose;
  vtkTree* armature;
  std::map<int,vtkSmartPointer<vtkTransform>> transformmap;
  std::map<vtkIdType, vtkIdType> modelbonemap;
 
  
  // Position based skinning functions
  int PBS_GetMaxElement(vtkSmartPointer<vtkIntArray> array, int comp);
  void SolveConstraints();
  void SolveEdgeConstraint(CellStructList celllist);
  void SolveVolumeConstraint(CellStructList celllist);
  void SolveBindConstraint(CellStructList celllist);
  void PBS_SetReferenceInfo();
  std::vector<vtkSmartPointer<vtkImplicitPolyDataDistance>> PBS_GetBoneDistFunc();
  
  double totalVolume = 0;
  std::vector<double> modelvolume; 
  
  std::vector<vtkSmartPointer<vtkImplicitPolyDataDistance>> boneposedistfunc;
  
  std::vector<std::vector<CellStructList>> edgeinfo;
  std::vector<std::vector<CellStructList>> cellinfo;
  std::vector<CellStructList> bindinfo;
  
};

#endif
