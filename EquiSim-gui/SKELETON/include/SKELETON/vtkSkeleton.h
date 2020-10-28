/**
 * @class   vtkSkeleton
 * @brief   A rooted tree data structure.
 *
 *
 * vtkTree is a connected directed graph with no cycles. A tree is a type of
 * directed graph, so works with all graph algorithms.
 *
 * vtkTree is a read-only data structure.
 * To construct a tree, create an instance of vtkMutableDirectedGraph.
 * Add vertices and edges with AddVertex() and AddEdge(). You may alternately
 * start by adding a single vertex as the root then call graph->AddChild(parent)
 * which adds a new vertex and connects the parent to the child.
 * The tree MUST have all edges in the proper direction, from parent to child.
 * After building the tree, call tree->CheckedShallowCopy(graph) to copy the
 * structure into a vtkTree. This method will return false if the graph is
 * an invalid tree.
 *
 * vtkTree provides some convenience methods for obtaining the parent and
 * children of a vertex, for finding the root, and determining if a vertex
 * is a leaf (a vertex with no children).
 *
 * @sa
 * vtkDirectedGraph vtkMutableDirectedGraph vtkGraph
*/

#ifndef vtkSkeleton_h
#define vtkSkeleton_h

#include <sstream>
#include <fstream>
#include <iostream>
#include <string>
#include <map>

#include "vtkTree.h"
#include "DLLDefines.h" // For export macro

#include "vtkMutableDirectedGraph.h"
#include "vtkAdjacentVertexIterator.h"
#include "vtkVertexListIterator.h"
#include "SKELETON/Defs.h"



//class vtkQuaternionArray;
class vtkVectorArray;
class vtkScalarArray;

class vtkQuaternionArray: public vtkDoubleArray {
public:
    vtkQuaternionArray() : helper(this,-1)
    {}

    // Nested class
    class Helper{
        public:
            // Initialisator
            Helper(vtkQuaternionArray* quatarray, vtkIdType _index) : m_index(_index), m_self(quatarray) 
            {}
            
            // Getter
            operator vtkQuaterniond () const
            { double array[4]; m_self->GetTuple(m_index, array); vtkQuaterniond quat(array); return quat; }
            
            // Setter
            void operator = (vtkQuaterniond quat)
            { double array[4]; quat.Get(array); m_self->SetTuple(m_index, array); }
            
        private:
            vtkIdType m_index;
            vtkQuaternionArray *m_self;
    } helper;
    
    Helper operator[](vtkIdType id){ return Helper(this, id); }
};




class vtkSkeleton: public vtkTree
{
public:
  static vtkSkeleton *New();
  vtkTypeMacro(vtkSkeleton, vtkTree);
  
  void PrintBoneAngles();
  
  vtkIdType AddNodeToGraph(vtkIdType parentid);
  void SetWorldTail(double tail[3]);
    void SetWorldHead(double head[3]);
    void SetBoneName(std::string name);
    void SetZAxis(double axis[3]);
    void SetAngleRanges(double axis[8]);
  bool CheckCompleteness();

  bool CheckROM(vtkIdType boneid, double (*angles)[3]);
    
    bool ParseArmature();

  
    vtkIdType GetBoneId(std::string bonename);
    std::string GetBoneName(vtkIdType id);
    
    
    void AddBvhArmatureArrays();
    void AddArmatureArrays();
    
    
    void RebuildTransformations();
    
    void ConvertGlobalRotationMatrixToEulerAngles(vtkQuaterniond R, double (*angles)[3]);
    
    void ModelGlobalTranform(vtkSmartPointer<vtkTransform> R);
    void SetModelTranslation(double (*modelorigin)[3]);
    void SetGlobalScale(double scale);
    void GetOrigin(double (*origin)[3]);
    void ScaleArmature(vtkIdType id);
    void SetBoneThickness(std::string bonename, double factor);
    void SetBoneThicknessById(vtkIdType boneid, double factor);
    void SetBoneLength(std::string bonename, double factor);
    void SetBoneLengthById(vtkIdType boneid, double factor);
    void SetBoneLengthRest(std::string bonename, double length); 
    void SetBoneLengthRestById(vtkIdType boneid, double length);
    
    void ScaleOffsetWidth(std::string bonename, double alpha);
    void ScaleOffsetWidth(vtkIdType boneid, double alpha);
    
    // Absolute angle setters
    bool SetBoneAngles(std::string bonename, double (*angles)[3]); // most absolute angle setter -> needs to set all angles of a bone
    bool SetBoneAnglesById(vtkIdType boneid, double (*angles)[3]);
    
    void UpdatePosition(vtkIdType parent);
    void UpdateRestPosition(vtkIdType parent);
     
    template <class ArrayType> ArrayType* GetVertexArray(std::string name);
    vtkSmartPointer<vtkDoubleArray> AddVertexArray(std::string name, size_t ncomp);
    template <class ArrayType> ArrayType* ConnectArrayToTree(vtkSmartPointer<ArrayType> array);
    
    //vtkDoubleArray* ConnectArrayToTree(vtkSmartPointer<vtkDoubleArray> array);
    vtkDoubleArray* AddAndConnectArray(std::string name, size_t ncomp);
    
    vtkSmartPointer<vtkPolyData> GetSimpleSkeleton();

    
    vtkSmartPointer<vtkTransform> GetTransformToParentsLocalFrame(vtkIdType boneid, bool wrtparent = true, bool getinv = false);
    
    // This method resets the skeleton for boneid, like you would do by reading a bvh file for this bone.
    // This method is used to change PC weights associated to the skeleton
    void ChangeRestArmatureWithoutSurfaceChange(vtkIdType boneid, double axis[3], double head[3], double tail[3]);
    
    double ConvertDummyToPhys(vtkIdType boneid, vtkIdType angleid, double angle);
    double ConvertPhysToDummy(vtkIdType boneid, vtkIdType angleid, double angle);
    
    double GetBoneLength(std::string bonename);
    double GetBoneAngle(std::string bonename, vtkIdType angleid);
    
    double GetGlobalScale();
    
    // Swap rest and pose skeleton 
    void SwapBoneRestPose(vtkIdType boneid);
    
    // Export function 
    void ExportSkeletonToTXT(std::string filename);
    void ExportSkeletonToBVH(std::string filename);
    
    
    void GetRotationUnitAxis(vtkIdType boneid, vtkIdType angleid, double (*ax)[3]);
    void GetRotationCenter(vtkIdType boneid, double (*ax)[3]);
    void SetRotationAxes(vtkIdType boneid);
    void InitRotationAxes(vtkIdType boneid);
    
    void ShallowCopy(vtkSkeleton* t);
    
    
    
    vtkGetMacro(ParentToBoneRotation, vtkQuaternionArray*);
    vtkGetMacro(WorldToBoneRestRotation, vtkQuaternionArray*);
    vtkGetMacro(ParentToOffsetRotation, vtkQuaternionArray*);
    vtkGetMacro(WorldToOffsetRotation, vtkQuaternionArray*);
    vtkGetMacro(RestToPoseRotation, vtkQuaternionArray*);
    
    vtkGetMacro(AxisZ, vtkVectorArray*);
    vtkGetMacro(TailPose, vtkVectorArray*);
    vtkGetMacro(HeadPose, vtkVectorArray*);
    vtkGetMacro(TailRest, vtkVectorArray*);
    vtkGetMacro(HeadRest, vtkVectorArray*);
    vtkGetMacro(ReferenceToPoseTransform, vtkVectorArray*);
    vtkGetMacro(AnthroFactor, vtkVectorArray*);
    vtkGetMacro(BoneAngles, vtkVectorArray*);
    vtkGetMacro(BoneAnglesRest, vtkVectorArray*);
    vtkGetMacro(BoneAngleRanges, vtkVectorArray*);
    vtkGetMacro(BoneRotAxes, vtkVectorArray*);
    
    vtkGetMacro(BoneScale_length, vtkScalarArray*);
    vtkGetMacro(BoneScale_thick, vtkScalarArray*);
    vtkGetMacro(BoneLength, vtkScalarArray*);
    vtkGetMacro(BoneLengthRest, vtkScalarArray*);
    vtkGetMacro(BoneOffsetScale_radius, vtkScalarArray*);
    vtkGetMacro(BoneOffsetScale_length, vtkScalarArray*);
    vtkGetMacro(BoneOffset, vtkScalarArray*);
    vtkGetMacro(BoneOffsetRest, vtkScalarArray*);
    
    vtkGetMacro(BoneName, vtkStringArray*);

    vtkSmartPointer<vtkMutableDirectedGraph> graph;
    
    void GetXAxis(vtkIdType boneid, double (*axisX)[3]);
    void GetYAxis(vtkIdType boneid, double (*axisY)[3]);
    
protected:
  vtkSkeleton();
  ~vtkSkeleton() VTK_OVERRIDE;

    
private:
  vtkSkeleton(const vtkSkeleton&) VTK_DELETE_FUNCTION;
  void operator=(const vtkSkeleton&) VTK_DELETE_FUNCTION;
  

    void GetLocalRotationAxis(vtkIdType boneid, vtkIdType axisid, double (*rotaxis)[3]);
    void RotateVectorByQuaternion(vtkQuaterniond* quat, double (*vec)[3], double (*result)[3]);
    void GetRelativeBonePosition(vtkIdType boneid, double (*pos)[3]);
    
    double GetLength(vtkIdType boneid);
    void RebuildWorldToBoneRestRotation(vtkIdType boneid);
    void InitRestToPoseRotation(int boneid);
    vtkQuaterniond RotationFromReferenceAxis(double referenceAxis[3], double axis[3]);
    vtkQuaterniond GetWorldToOffsetRotation(vtkIdType boneid, double* offsetlength);
    void SetParentToBoneOffset(vtkIdType boneid);
    vtkQuaterniond GetParentToBoneRotation(vtkIdType boneid);
    vtkQuaterniond GetParentToOffsetRotation(vtkIdType boneid);
    void SetParentToBoneRotations(vtkIdType boneid);
    void ProjectPointOnPlane(double orig[3], double evec[3], double point[3], double (*proj)[3]);
    void GetEulerAngles(vtkIdType boneid, double (*angles)[3]);
    void SetEulerAngles(vtkIdType boneid);
    void UpdateEulerAngles(vtkIdType boneid);
    
    void RebuildScaleBoneById(vtkIdType boneid);
    
    void SetAnthroExtVector(int boneid);
    
    void UpdateParentToOffsetRotation(vtkIdType boneid);
    void UpdateReferenceToPoseTransform(vtkIdType boneid);
    void UpdateParentToBonePoseTransform(vtkIdType boneid, double angles[3]);
    void UpdateWorldToBoneRotation(vtkIdType boneid);
    void MultiplyVectorByQuaternion(vtkQuaterniond rot, double (*vec)[3], double (*res)[3]);
    void RotateOffset(vtkIdType boneid);
    void RotateRestOffset(vtkIdType boneid);
    void RotateBone(vtkIdType boneid);
    void RotateRestBone(vtkIdType boneid);
    void PositionUpdate(vtkIdType boneid);
    void RestPositionUpdate(vtkIdType boneid);
    
    double GlobalScale;
    
    double Origin[3]={0,0,0};
    double GlobalAxisX[3]={1,0,0};
    double GlobalAxisY[3]={0,1,0};
    double GlobalAxisZ[3]={0,0,1};
  
    vtkQuaternionArray* ParentToBoneRotation;
    vtkQuaternionArray* WorldToBoneRotation;
    vtkQuaternionArray* WorldToBoneRestRotation;
    vtkQuaternionArray* ParentToOffsetRotation;
    vtkQuaternionArray* WorldToOffsetRotation;
    vtkQuaternionArray* RestToPoseRotation;
    
    vtkVectorArray* AxisZ;
    vtkVectorArray* TailPose;
    vtkVectorArray* HeadPose;
    vtkVectorArray* TailRest;
    vtkVectorArray* HeadRest;
    vtkVectorArray* ReferenceToPoseTransform;
    vtkVectorArray* AnthroFactor;
    vtkVectorArray* BoneAngles;
    vtkVectorArray* BoneAnglesRest;
    vtkVectorArray* BoneAngleRanges;
    vtkVectorArray* BoneRotAxes;
    
    vtkScalarArray* BoneScale_length;
    vtkScalarArray* BoneScale_thick;
    vtkScalarArray* BoneLength;
    vtkScalarArray* BoneLengthRest;
    vtkScalarArray* BoneOffsetScale_radius;
    vtkScalarArray* BoneOffsetScale_length;
    vtkScalarArray* BoneOffset;
    vtkScalarArray* BoneOffsetRest;
    
    vtkStringArray* BoneName;
    
    vtkSmartPointer<vtkDoubleArray> AxisZ_bvh;
    vtkSmartPointer<vtkDoubleArray> HeadPose_bvh;
    vtkSmartPointer<vtkDoubleArray> TailPose_bvh;
    vtkSmartPointer<vtkDoubleArray> HeadRest_bvh;
    vtkSmartPointer<vtkDoubleArray> TailRest_bvh;
    vtkSmartPointer<vtkDoubleArray> BoneAngleRanges_bvh;
    vtkSmartPointer<vtkStringArray> BoneName_bvh;
    
    void ExportBones(vtkIdType parentid, std::fstream& file);
    void ExportBonesBVH(vtkIdType parentid, std::fstream& file);
    
    const double PI = 3.14159265359;
    
};

#endif
