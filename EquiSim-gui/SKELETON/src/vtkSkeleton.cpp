#include "vtkInformation.h"
#include "vtkInformationVector.h"
#include "vtkObjectFactory.h"
#include "vtkOutEdgeIterator.h"
#include "vtkSmartPointer.h"
#include "vtkStdString.h"

#include <vector>
#include <SKELETON/vtkSkeleton.h>
#include "SKELETON/Defs.h"



vtkStandardNewMacro(vtkSkeleton);

namespace{
    void CopyVector3(const double* vec, double* copyVec)
    {
        copyVec[0] = vec[0];
        copyVec[1] = vec[1];
        copyVec[2] = vec[2];
    }
    
    bool CompareVector3(const double* v1, const double* v2)
    {
        double diff[3];
        vtkMath::Subtract(v1, v2, diff);
        if (vtkMath::Dot(diff, diff) < 1e-6)
        {
        return true;
        }
        return false;
    }
    
    double DotProd(double (*v1)[3], double (*v2)[3])
    {
        vtkMath::Normalize(*v1);
        vtkMath::Normalize(*v2);
        
        double dot = vtkMath::Dot(*v1, *v2);

        // Expect v1 and v2 to be normalized unit vectors, but sometimes the dot product
        // exceeds 1 by just enough to cause errors when calling acos.
        // If it's greater than 1.0, or less than -1.0, clamp it:
        if (dot > 1.0)
        {dot = 1.0;}
        else if (dot < -1.0)
        {dot = -1.0;}
        return dot;
    } 
}


// construction based on https://social.msdn.microsoft.com/Forums/vstudio/en-US/517a601d-5e1f-4d72-8927-32e80a094da9/overloading-subscript-operator-separating-setting-and-getting-values?forum=vclanguage



class vtkVectorArray: public vtkDoubleArray {
public:
    vtkVectorArray() : helper(this,-1)
    {}

    // Nested class
    class Helper{
        public:
            // Initialisator
            Helper(vtkVectorArray* vecarray, vtkIdType _index) : m_index(_index), m_self(vecarray) 
            {}
            
            // Getter (assumes that vec is allocated well!)
            void operator > (double* vec) const
            { 
                //int Ncomp = m_self->GetNumberOfComponents();
                //double array[3]; 
                m_self->GetTuple(m_index, vec); 
                //CopyVector3(array, vec); 
                
            }
            
            // Setter
            void operator = (double* newvec)
            { m_self->SetTuple(m_index, newvec); }
            
        private:
            vtkIdType m_index;
            vtkVectorArray *m_self;
    } helper;
    
    Helper operator[](vtkIdType id){ return Helper(this, id); }
};

class vtkScalarArray: public vtkDoubleArray {
public:
    vtkScalarArray() : helper(this,-1)
    {}

    // Nested class
    class Helper{
        public:
            // Initialisator
            Helper(vtkScalarArray* scalarray, vtkIdType _index) : m_index(_index), m_self(scalarray) 
            {}
            
            // Getter
            operator double () const
            { double val; m_self->GetTuple(m_index, &val); return val; }
            
            // Setter
            void operator = (double newval)
            { m_self->SetTuple(m_index, &newval); }
            
        private:
            vtkIdType m_index;
            vtkScalarArray *m_self;
    } helper;
    
    Helper operator[](vtkIdType id){ return Helper(this, id); }
};

//----------------------------------------------------------------------------
vtkSkeleton::vtkSkeleton()
{
    // initialisation
    GlobalScale = 1;
    
    AddBvhArmatureArrays();
        
    graph = vtkSmartPointer< vtkMutableDirectedGraph>::New();
}

vtkIdType vtkSkeleton::AddNodeToGraph(vtkIdType parentid){
    vtkIdType newbone;
    
    if (parentid == -1){
        newbone = graph->AddVertex(); // Add root vertex
    }
    else{
        newbone = graph->AddChild(parentid);
    }
    return newbone;
}

bool vtkSkeleton::ParseArmature(){
    // Copies the structure of the graph in the vtkTree. 
    bool ok = this->CheckedShallowCopy(graph);
    
    // connect arrays to the vtktree
    AxisZ = (vtkVectorArray*) ConnectArrayToTree(AxisZ_bvh);
    HeadPose = (vtkVectorArray*) ConnectArrayToTree(HeadPose_bvh);
    HeadRest = (vtkVectorArray*) ConnectArrayToTree(HeadRest_bvh);
    TailPose = (vtkVectorArray*) ConnectArrayToTree(TailPose_bvh);
    TailRest = (vtkVectorArray*) ConnectArrayToTree(TailRest_bvh);

    BoneAngleRanges = (vtkVectorArray*) ConnectArrayToTree(BoneAngleRanges_bvh);
    
    this->GetVertexData()->AddArray(BoneName_bvh);

    return ok;
}
//----------------------------------------------------------------------------
vtkSkeleton::~vtkSkeleton()
{
}

template <class ArrayType>
ArrayType* vtkSkeleton::GetVertexArray(std::string name)
{
    return ArrayType::SafeDownCast(this->GetVertexData()->GetArray(name.c_str()));
}
template vtkDoubleArray* vtkSkeleton::GetVertexArray<vtkDoubleArray>(std::string name);
template vtkStringArray* vtkSkeleton::GetVertexArray<vtkStringArray>(std::string name);

vtkSmartPointer<vtkDoubleArray> vtkSkeleton::AddVertexArray(std::string name, size_t ncomp){
    vtkSmartPointer<vtkDoubleArray> array = vtkSmartPointer<vtkDoubleArray>::New();
    array->SetNumberOfComponents(ncomp);
    array->SetNumberOfTuples(GetNumberOfVertices()); // not needed if you use InsertNextTuple to fill the array. This would be needed if you use SetTuple. 
    array->SetName(name.c_str());
    return array;
}

// This functions adds an array to the vertex data of the tree and returns a pointer this data array.
template <class ArrayType>
ArrayType* vtkSkeleton::ConnectArrayToTree(vtkSmartPointer<ArrayType> array){
    this->GetVertexData()->AddArray(array);
    return GetVertexArray<ArrayType>(array->GetName());
}
template vtkDoubleArray* vtkSkeleton::ConnectArrayToTree<vtkDoubleArray>(vtkSmartPointer<vtkDoubleArray> array);
template vtkStringArray* vtkSkeleton::ConnectArrayToTree<vtkStringArray>(vtkSmartPointer<vtkStringArray> array);

vtkDoubleArray* vtkSkeleton::AddAndConnectArray(std::string name, size_t ncomp){
    ConnectArrayToTree(AddVertexArray(name, ncomp));    
    return GetVertexArray<vtkDoubleArray>(name);
}



void vtkSkeleton::AddBvhArmatureArrays(){
    AxisZ_bvh = AddVertexArray("AxisZ", 3); // coming from bvh file, only used to intialise the world to bone rotation
    TailPose_bvh = AddVertexArray("TailPose", 3); 
    HeadPose_bvh = AddVertexArray("HeadPose", 3); 
    TailRest_bvh = AddVertexArray("TailRest", 3); 
    HeadRest_bvh = AddVertexArray("HeadRest", 3); 
    
    BoneAngleRanges_bvh = AddVertexArray("AngleRange", 6);
    
    BoneName_bvh = vtkSmartPointer<vtkStringArray>::New();
    BoneName_bvh->SetName("BoneName");
}


void vtkSkeleton::AddArmatureArrays()
{
  // Set arrays to the vtktree
    
    // Make the arrays easy accesible by storing the pointers as private data members. 
    ParentToBoneRotation = (vtkQuaternionArray*) AddAndConnectArray("ParentToBoneRotation", 4);
    WorldToBoneRotation = (vtkQuaternionArray*) AddAndConnectArray("WorldToBoneRotation", 4);
    WorldToBoneRestRotation = (vtkQuaternionArray*) AddAndConnectArray("WorldToBoneRestRotation", 4);
    ParentToOffsetRotation = (vtkQuaternionArray*) AddAndConnectArray("ParentToOffsetRotation", 4);
    WorldToOffsetRotation = (vtkQuaternionArray*) AddAndConnectArray("WorldToOffsetRotation", 4);
    RestToPoseRotation = (vtkQuaternionArray*) AddAndConnectArray("RestToPoseRotation", 4);
    
    ReferenceToPoseTransform = (vtkVectorArray*) AddAndConnectArray("ReferenceToPoseTransform", 9);
    AnthroFactor = (vtkVectorArray*) AddAndConnectArray("AnthroFactor", 3);
    BoneAngles = (vtkVectorArray*) AddAndConnectArray("BoneAngles", 3);
    BoneAnglesRest = (vtkVectorArray*) AddAndConnectArray("BoneAnglesRest", 3);
    BoneRotAxes = (vtkVectorArray*) AddAndConnectArray("BoneRotAxes", 9);
    
    BoneScale_length = (vtkScalarArray*) AddAndConnectArray("BoneScale_length", 1);
    BoneScale_thick = (vtkScalarArray*) AddAndConnectArray("BoneScale_thick", 1);
    BoneLength = (vtkScalarArray*) AddAndConnectArray("BoneLength", 1);
    BoneLengthRest = (vtkScalarArray*) AddAndConnectArray("BoneLengthRest", 1);
    BoneOffsetScale_radius = (vtkScalarArray*) AddAndConnectArray("BoneOffsetScale_radius", 1);
    BoneOffsetScale_length = (vtkScalarArray*) AddAndConnectArray("BoneOffsetScale_length", 1);
    BoneOffset = (vtkScalarArray*) AddAndConnectArray("BoneOffset", 1);
    BoneOffsetRest = (vtkScalarArray*) AddAndConnectArray("BoneOffsetRest", 1);
    
}


vtkSmartPointer<vtkPolyData> vtkSkeleton::GetSimpleSkeleton(){
    vtkSmartPointer<vtkPolyData> polydata = vtkSmartPointer<vtkPolyData>::New();
    vtkSmartPointer<vtkPoints> pts = vtkSmartPointer<vtkPoints>::New();
    vtkSmartPointer<vtkCellArray> cells = vtkSmartPointer<vtkCellArray>::New();
    
    for (vtkIdType boneid=0; boneid<this->GetNumberOfVertices(); boneid++){
        double head[3];
        double tail[3];
        (*HeadPose)[boneid] > head;
        (*TailPose)[boneid] > tail;
      
        pts->InsertNextPoint(head);
        pts->InsertNextPoint(tail);
        
        // Join the two points in one cell
        vtkSmartPointer<vtkLine> line1 = vtkSmartPointer<vtkLine>::New();
        line1->GetPointIds()->SetId(0, 2*boneid); 
        line1->GetPointIds()->SetId(1, 2*boneid+1); 
        
        cells->InsertNextCell(line1);
    }
    
    polydata->SetPoints(pts);
    polydata->SetLines(cells);
    
    return polydata;
}


// This function calculates the world->bone and the parent->bone transformations and sets the kinematical angles. Only used to initialise the skeleton
// Based on head, tail and Z-Axis provided by the bvh-file, this function calculates all the other arrays. 
void vtkSkeleton::RebuildTransformations(){

    vtkSmartPointer<vtkVertexListIterator> vertexit = vtkSmartPointer<vtkVertexListIterator>::New();
    this->GetVertices(vertexit);
   
    while(vertexit->HasNext()){
        vtkIdType boneid = vertexit->Next();
        InitRestToPoseRotation(boneid); // RestToPoseRotation
        SetAnthroExtVector(boneid); //AnthroExtVector
        RebuildWorldToBoneRestRotation(boneid); // WorldToBoneRotation, WorldToBoneRestRotation
        SetParentToBoneOffset(boneid); // WorldToOffsetRotation, BoneOffset, BoneLength
        SetParentToBoneRotations(boneid); // ParentToBoneRotation, ParentToOffsetRotation
        SetEulerAngles(boneid); // BoneAngles
	InitRotationAxes(boneid); // Store the rotation axes (meant as output values used for Inverse Kinematics applications)
    }
}

void vtkSkeleton::SetAnthroExtVector(int boneid){
    double vec[3]={1,1,1};
    AnthroFactor->SetTuple(boneid, vec);
}

bool vtkSkeleton::CheckROM(vtkIdType boneid, double (*angles)[3]){
    if(boneid==0){ return true; } // do not check for the root bone
    double ROM[6];
    (*BoneAngleRanges)[boneid] > ROM;
    bool check = true;
    for(vtkIdType it=0; it<3; it++){
        check = check && ((*angles)[it] >= ROM[2*it]) && ((*angles)[it] <= ROM[2*it+1]);
        (*angles)[it] = std::max(std::min((*angles)[it], ROM[2*it+1]), ROM[2*it]);
    }
    return check;
}

void vtkSkeleton::GetOrigin(double (*origin)[3]){
    double orig[3];
    (*HeadPose)[this->GetRoot()] > orig;
    CopyVector3(orig, *origin);
}

void vtkSkeleton::ScaleOffsetWidth(std::string bonename, double alpha){
  vtkIdType boneid = BoneName_bvh->LookupValue(bonename);
  if(boneid == -1){throw std::runtime_error("Bone not found in the model.");}
  ScaleOffsetWidth(boneid, alpha);
}
  
void vtkSkeleton::ScaleOffsetWidth(vtkIdType boneid, double alpha){ 
  double alpha_prev = (*BoneOffsetScale_radius)[boneid];
  double alpha_new = alpha/alpha_prev;
  double scalingmat[3][3] = {{alpha_new,0,0},{0,1,0},{0,0,alpha_new}};
  (*BoneOffsetScale_radius)[boneid] = alpha;
  
  double head_parent[3];
  (*HeadPose)[boneid] > head_parent;
  
  vtkQuaterniond quat = (*WorldToBoneRotation)[boneid];
  double rotation[3][3]; double rotation_inverse[3][3];
  quat.ToMatrix3x3(rotation); quat.Inverse().ToMatrix3x3(rotation_inverse);
  
  vtkMath::Multiply3x3(rotation, scalingmat, rotation);
  vtkMath::Multiply3x3(rotation, rotation_inverse, rotation);
  
  vtkSmartPointer<vtkAdjacentVertexIterator> childreniterator = vtkSmartPointer<vtkAdjacentVertexIterator>::New();
  GetChildren(boneid, childreniterator);
  while(childreniterator->HasNext()){
    vtkIdType childid = childreniterator->Next();
    double head[3];
    (*HeadPose)[childid] > head;
    vtkMath::Subtract(head, head_parent, head);
    vtkMath::Multiply3x3(rotation, head, head);
    vtkMath::Add(head, head_parent, head);
    (*HeadPose)[childid] = head;
    
    UpdateParentToOffsetRotation(childid);
  }
  
  UpdatePosition(boneid); 
  this->Modified();
  
}

void vtkSkeleton::ModelGlobalTranform(vtkSmartPointer<vtkTransform> R){
    // Translation
    double modelorigin[3]; 
    (*HeadRest)[this->GetRoot()] > modelorigin;
    R->TransformPoint(modelorigin, modelorigin);
    
    // Scaling    
    double scale[3];
    R->GetScale(scale);
    if(scale[0] - scale[1] > 1e-9 || scale[2] - scale[1] > 1e-9){throw std::runtime_error("Initial transformation involves non-uniform scaling.");}
    
    // Rotation
    double rotation[4]; double rotangles[3]; 
    R->GetOrientationWXYZ(rotation);
    vtkQuaterniond quat = vtkQuaterniond(cos(rotation[0]/2 *PI/180), rotation[1] * sin(rotation[0]/2 *PI/180), rotation[2] * sin(rotation[0]/2 *PI/180), rotation[3] * sin(rotation[0]/2 *PI/180));
    this->ConvertGlobalRotationMatrixToEulerAngles(quat, &rotangles);
    double angles[3] = {rotangles[0], rotangles[1], rotangles[2]};   
    
    // Apply rotation 
    //vtkQuaterniond quat_orig = (*ParentToBoneRotation)[this->GetRoot()];
    //(*ParentToBoneRotation)[this->GetRoot()] = quat * quat_orig;
    
    //UpdateEulerAngles(0); //updates BoneAngles
    // Updates for the bone and its children
    //UpdatePosition(0); 
    //this->Modified();
    
    // Change the model parameters to these values
    this->SetGlobalScale(scale[0]);
    this->SetBoneAnglesById(this->GetRoot(), &angles);
    this->SetModelTranslation(&modelorigin);
}

// Converts a global rotation matrix to euler angles of the root bone. This makes it easy to rotate the model based on landmarks with a third-party registration method and to interpret its result with model parameters
void vtkSkeleton::ConvertGlobalRotationMatrixToEulerAngles(vtkQuaterniond R, double (*angles)[3]){
  vtkQuaterniond quat_orig = (*ParentToBoneRotation)[this->GetRoot()];
  (*ParentToBoneRotation)[this->GetRoot()] = R * quat_orig;
  GetEulerAngles(this->GetRoot(), angles);
  (*ParentToBoneRotation)[this->GetRoot()] = quat_orig;
}


void vtkSkeleton::SetModelTranslation(double (*modelorigin)[3]){ 
  (*HeadPose)[this->GetRoot()] = *modelorigin;  
  UpdatePosition(this->GetRoot());   
  this->Modified();
}

void vtkSkeleton::SetGlobalScale(double scale){
    GlobalScale = scale;
    ScaleArmature(this->GetRoot());
    UpdatePosition(this->GetRoot()); 
    this->Modified();
}

double vtkSkeleton::GetGlobalScale(){
    return GlobalScale;
}

void vtkSkeleton::ScaleArmature(vtkIdType id){
  
    RebuildScaleBoneById(id);
  
    vtkSmartPointer<vtkAdjacentVertexIterator> childreniterator = vtkSmartPointer<vtkAdjacentVertexIterator>::New();
    GetChildren(id, childreniterator);
    while(childreniterator->HasNext()){
	vtkIdType childid = childreniterator->Next();
	ScaleArmature(childid);
    }
  
}


void vtkSkeleton::RebuildScaleBoneById(vtkIdType boneid){
  
  double fac[3] = {GlobalScale * (*BoneScale_thick)[boneid], GlobalScale * (*BoneScale_length)[boneid], GlobalScale * (*BoneScale_thick)[boneid] };
  (*AnthroFactor)[boneid] = fac;
  (*BoneLength)[boneid] = fac[1] * (*BoneLengthRest)[boneid];
  
  // scale offsets to its children
  vtkSmartPointer<vtkAdjacentVertexIterator> childreniterator = vtkSmartPointer<vtkAdjacentVertexIterator>::New();
  GetChildren(boneid, childreniterator);
  while(childreniterator->HasNext()){
    vtkIdType childid = childreniterator->Next();
    //assert((*BoneOffsetScale_length)[childid] == 1);
    (*BoneOffset)[childid] = fac[1] * (*BoneOffsetRest)[childid] * (*BoneOffsetScale_length)[childid];
  }
  
}

// Change rest bonelength without modifying surface. 
void vtkSkeleton::SetBoneLengthRest(std::string bonename, double factor){
    vtkIdType boneid = BoneName_bvh->LookupValue(bonename);
    if(boneid == -1){throw std::runtime_error("Bone not found in the model.");}
    SetBoneLengthRestById(boneid, factor);
}

// Change rest bonelength without modifying surface. 
void vtkSkeleton::SetBoneLengthRestById(vtkIdType boneid, double length){
  
  // scale offsets to its children
  vtkSmartPointer<vtkAdjacentVertexIterator> childreniterator = vtkSmartPointer<vtkAdjacentVertexIterator>::New();
  GetChildren(boneid, childreniterator);
  while(childreniterator->HasNext()){
    vtkIdType childid = childreniterator->Next();
    (*BoneOffsetRest)[childid] = (*BoneOffsetRest)[childid] * length / (*BoneLengthRest)[boneid];
  }
  
  (*BoneLengthRest)[boneid] = length;
  RebuildScaleBoneById(boneid); // this changes the bonelength, taking prev scaling factors into account
  UpdateRestPosition(boneid); //
  UpdatePosition(boneid);  
  /*
  for(vtkIdType id=0; id<GetNumberOfVertices(); id++){
    double head[3];
    (*HeadPose)[boneid] > head;
    (*HeadRest)[boneid] = head;
  }
  */
  this->Modified(); 
}

// Change the bonelength with a fraction of the original length (absolute setter!)
void vtkSkeleton::SetBoneLength(std::string bonename, double factor){
    vtkIdType boneid = BoneName_bvh->LookupValue(bonename);
    if(boneid == -1){throw std::runtime_error("Bone not found in the model.");}
    SetBoneLengthById(boneid, factor);
}

void vtkSkeleton::SetBoneThickness(std::string bonename, double factor){
    vtkIdType boneid = BoneName_bvh->LookupValue(bonename);
    if(boneid == -1){throw std::runtime_error("Bone not found in the model.");}
    SetBoneThicknessById(boneid, factor);
}

void vtkSkeleton::SetBoneThicknessById(vtkIdType boneid, double factor){
    (*BoneScale_thick)[boneid] = factor;
    RebuildScaleBoneById(boneid);
    UpdateReferenceToPoseTransform(boneid);
    this->Modified();
}

void vtkSkeleton::SetBoneLengthById(vtkIdType boneid, double factor){
    (*BoneScale_length)[boneid] = factor;
    RebuildScaleBoneById(boneid);
    UpdatePosition(boneid);   
    this->Modified();
}

// Set all angles of a bone at once. Needed to define absolute angles, irrespective of previous pose.
bool vtkSkeleton::SetBoneAngles(std::string bonename, double (*angles)[3]){
    // Find the id of the bone
    vtkIdType boneid = BoneName_bvh->LookupValue(bonename);
    if(boneid == -1){throw std::runtime_error("Bone not found in the model.");}
    
    // Check angle is within ROM
    bool check = true;
    CheckROM(boneid, angles);
    //if(check == false){std::cout << bonename << std::endl;throw std::runtime_error("Not in ROM");}
    
    // Updates only for the transformed bone
    UpdateParentToBonePoseTransform(boneid, *angles); // updates ParentToBone
    
    //UpdateEulerAngles(boneid); //updates BoneAngles
    (*BoneAngles)[boneid] = *angles;
    
    // Updates for the bone and its children
    UpdatePosition(boneid);   
    
    this->Modified();
    
    return check;
}

bool vtkSkeleton::SetBoneAnglesById(vtkIdType boneid, double (*angles)[3]){
    // Check angle is within ROM
    // bool check = true; 
    // CheckROM(boneid, angles);
    
    
    // Updates only for the transformed bone
    UpdateParentToBonePoseTransform(boneid, *angles); // updates ParentToBone
    
    
    UpdateEulerAngles(boneid); //updates BoneAngles
    
    double calcangles[3];
    (*BoneAngles)[boneid] > calcangles; 
    
    //std::cout << "Bone angle ("<< (*angles)[0] << " " << (*angles)[1] <<  " " << (*angles)[2] << ") = (" << calcangles[0] << " " << calcangles[1] << " " << calcangles[2] <<")" <<  std::endl;
    
    if(std::abs((*angles)[0]-calcangles[0])>1e-8 || std::abs((*angles)[1]-calcangles[1])>1e-8 || std::abs((*angles)[2]-calcangles[2])>1e-8){
        throw std::runtime_error("Cross-check failed. The requested bone angle does not correspond with the recomputed value.");
    }
    
    (*BoneAngles)[boneid] = *angles;
    
    // Updates for the bone and its children
    UpdatePosition(boneid); 
    
    this->Modified();
    
    return true;
}

// convert dummy variable used in registration program to phys angle
double vtkSkeleton::ConvertDummyToPhys(vtkIdType boneid, vtkIdType angleid, double angle){
    double angles_range[6];
    (*BoneAngleRanges)[boneid] > angles_range;
    double min = angles_range[2*angleid];
    double max = angles_range[2*angleid+1];
    double angle_phys = min + (sin(angle)+1)*(max-min)/2;
    return std::max(std::min(angle_phys, max), min); // just to be sure that number stays within range to avoid nans later. Computation rounding errors possible.  
}

// convert phys angle to dummy variable used in registration program
double vtkSkeleton::ConvertPhysToDummy(vtkIdType boneid, vtkIdType angleid, double angle){
    double angles_range[6];
    (*BoneAngleRanges)[boneid] > angles_range;
    double min = angles_range[2*angleid];
    double max = angles_range[2*angleid+1];
    
    double angle_corrected = std::max(std::min(angle, max), min); // just to be sure that angle lies in between the boundaries. Otherwise asin will give nan.
    
    double angle_dummy = asin(2*(angle_corrected-min)/(max-min)-1);
    return std::max(std::min(angle_dummy, PI/2), -PI/2);
}




void vtkSkeleton::UpdatePosition(vtkIdType boneid){
    
    UpdateWorldToBoneRotation(boneid); // Updates WorldToBone, RestToPose, WorldToOffset
    PositionUpdate(boneid); // Updates head and tail position
    SetRotationAxes(boneid); // Updates the rotation axes (for output)
    
    vtkSmartPointer<vtkAdjacentVertexIterator> childreniterator = vtkSmartPointer<vtkAdjacentVertexIterator>::New();
    GetChildren(boneid, childreniterator);
    while(childreniterator->HasNext()){
        vtkIdType childid = childreniterator->Next();
        UpdatePosition(childid);
    }
}


void vtkSkeleton::UpdateRestPosition(vtkIdType boneid){
    
    RestPositionUpdate(boneid); // Updates head and tail position
    
    vtkSmartPointer<vtkAdjacentVertexIterator> childreniterator = vtkSmartPointer<vtkAdjacentVertexIterator>::New();
    GetChildren(boneid, childreniterator);
    while(childreniterator->HasNext()){
        vtkIdType childid = childreniterator->Next();
        UpdateRestPosition(childid);
    }
}






vtkIdType vtkSkeleton::GetBoneId(std::string bonename){
  return BoneName_bvh->LookupValue(bonename);
}

double vtkSkeleton::GetBoneLength(std::string bonename){
    vtkIdType boneid = BoneName_bvh->LookupValue(bonename);
    if(boneid == -1){throw std::runtime_error("Bone not found in the model.");}
    return (*BoneLength)[boneid];
}


double vtkSkeleton::GetBoneAngle(std::string bonename, vtkIdType angleid){
    vtkIdType boneid = BoneName_bvh->LookupValue(bonename);
    if(boneid == -1){throw std::runtime_error("Bone not found in the model.");}
    double triplet[3];
    (*BoneAngles)[boneid] > triplet;
    return triplet[angleid];
}

// Rotation axes are defined in local coordinate system of parent. 
void vtkSkeleton::GetLocalRotationAxis(vtkIdType boneid, vtkIdType axisid, double (*rotaxis)[3]){

    double relbonepos[3];
    GetRelativeBonePosition(boneid, &relbonepos); // the orientation of the bone in its parent local coordinate system
    
    if(axisid==0){
        vtkMath::Cross(relbonepos, GlobalAxisZ, *rotaxis);
    }
    else if(axisid==1){
        CopyVector3(GlobalAxisZ, *rotaxis);
    }
    else if(axisid==2){ // roll axis
        CopyVector3(relbonepos, *rotaxis);
    }
    
    vtkMath::Normalize(*rotaxis);
}

void vtkSkeleton::RotateVectorByQuaternion(vtkQuaterniond* quat, double (*vec)[3], double (*result)[3]){
    double rotation[3][3];
    quat->ToMatrix3x3(rotation);
    double translation[3];
    vtkMath::Multiply3x3(rotation, *vec, translation);
    CopyVector3(translation, *result);
}

void vtkSkeleton::GetRelativeBonePosition(vtkIdType boneid, double (*pos)[3]){
    double relpos[3];    
    vtkQuaterniond quat = (*ParentToBoneRotation)[boneid];
    RotateVectorByQuaternion(&quat, &GlobalAxisY, &relpos);
    vtkMath::Normalize(relpos);
    CopyVector3(relpos, *pos);
}

// Set-functions, used by the BVH reader

void vtkSkeleton::SetWorldTail(double tail[3])
{
    TailPose_bvh->InsertNextTuple(tail);
    TailRest_bvh->InsertNextTuple(tail);
}

void vtkSkeleton::SetWorldHead(double head[3])
{
    HeadPose_bvh->InsertNextTuple(head);
    HeadRest_bvh->InsertNextTuple(head);
}

void vtkSkeleton::SetBoneName(std::string name){
    BoneName_bvh->InsertNextValue(name);
}

void vtkSkeleton::SetZAxis(double axis[3]) // roll axis comin from bvh file
{
    AxisZ_bvh->InsertNextTuple(axis);
}

void vtkSkeleton::SetAngleRanges(double axis[6])
{
    BoneAngleRanges_bvh->InsertNextTuple(axis);
}

void vtkSkeleton::ChangeRestArmatureWithoutSurfaceChange(vtkIdType boneid, double axis[3], double head[3], double tail[3]){
    
    (*AxisZ)[boneid] = axis;
    (*HeadPose)[boneid] = head;
    (*HeadRest)[boneid] = head;
    (*TailPose)[boneid] = tail;
    (*TailRest)[boneid] = tail;
    
    InitRestToPoseRotation(boneid); // RestToPoseRotation
    SetAnthroExtVector(boneid); //AnthroExtVector
    RebuildWorldToBoneRestRotation(boneid); // WorldToBoneRotation, WorldToBoneRestRotation
    SetParentToBoneOffset(boneid); // WorldToOffsetRotation, BoneOffset, BoneLength
    SetParentToBoneRotations(boneid); // ParentToBoneRotation, ParentToOffsetRotation
    SetEulerAngles(boneid); // BoneAngles
    InitRotationAxes(boneid); // Store the rotation axes (meant as output values used for Inverse Kinematics applications)
  
    UpdatePosition(boneid);
    
        
    if (this->GetNumberOfChildren(boneid) > 0){
        vtkSmartPointer<vtkAdjacentVertexIterator> childreniterator = vtkSmartPointer<vtkAdjacentVertexIterator>::New();
        GetChildren(boneid, childreniterator);
        while(childreniterator->HasNext()){
            vtkIdType childid = childreniterator->Next();
	    
            double headpose[3]; double tailpose[3]; double axispose[3];
            (*HeadPose)[childid] > headpose;
            (*TailPose)[childid] > tailpose;
            (*AxisZ)[childid] > axispose;
	    
            ChangeRestArmatureWithoutSurfaceChange(childid, axispose, headpose, tailpose);
        }
    }
    
    
}


// The arrays head, tail, axesZ are extended based on the information provided in the bvh file
bool vtkSkeleton::CheckCompleteness(){
    int Ntails = TailPose_bvh->GetNumberOfTuples();
    int Nheads = HeadPose_bvh->GetNumberOfTuples();
    int NAxes = AxisZ_bvh->GetNumberOfTuples();
    int NRanges = BoneAngleRanges_bvh->GetNumberOfTuples();
    int Nvert = graph->GetNumberOfVertices();
    return (Ntails == Nheads) && (NAxes == Nvert) && (Nheads == NAxes) && (Nheads == NRanges);
}

void vtkSkeleton::GetXAxis(vtkIdType boneid, double (*axisX)[3]){
    double axisZ[3]; (*AxisZ)[boneid] > axisZ;
    double axisY[3]; GetYAxis(boneid, &axisY);
    vtkMath::Cross(axisY, axisZ, *axisX);
}

void vtkSkeleton::GetYAxis(vtkIdType boneid, double (*axisY)[3]){
    double head[3]; (*HeadPose)[boneid] > head;
    double tail[3]; (*TailPose)[boneid] > tail;
    vtkMath::Subtract(tail, head, *axisY);
}

double vtkSkeleton::GetLength(vtkIdType boneid){
    double ax[3];
    GetYAxis(boneid, &ax);
    return vtkMath::Norm(ax);
}



void vtkSkeleton::RebuildWorldToBoneRestRotation(vtkIdType boneid) // Build world to bone rotation for initialisation
{     
    // Get local axes
    double axisY[3]; GetYAxis(boneid, &axisY); vtkMath::Normalize(axisY);
    double axisZ[3]; (*AxisZ)[boneid] > axisZ; vtkMath::Normalize(axisZ);
    double axisX[3]; vtkMath::Cross(axisY, axisZ, axisX); vtkMath::Normalize(axisX);
  
    // Define local and global coordinate system
    vtkSmartPointer<vtkPoints> globalcoordsyst = vtkSmartPointer<vtkPoints>::New();
    vtkSmartPointer<vtkPoints> localcoordsyst = vtkSmartPointer<vtkPoints>::New();
    globalcoordsyst->InsertNextPoint(GlobalAxisX); globalcoordsyst->InsertNextPoint(GlobalAxisY); globalcoordsyst->InsertNextPoint(GlobalAxisZ);
    localcoordsyst->InsertNextPoint(axisX); localcoordsyst->InsertNextPoint(axisY); localcoordsyst->InsertNextPoint(axisZ);

    // Transform the global coord system to the local one
    vtkSmartPointer<vtkTransform> transform = vtkSmartPointer<vtkTransform>::New();
    vtkSmartPointer<vtkLandmarkTransform> transf = vtkSmartPointer<vtkLandmarkTransform>::New();
    transf->SetSourceLandmarks(globalcoordsyst);
    transf->SetTargetLandmarks(localcoordsyst);
    transf->SetModeToRigidBody(); // rotation, scaling and translation only
    transf->Update();
    transform->Concatenate(transf);
    
    double rotation[4]; double rotangles[3]; 
    transform->GetOrientationWXYZ(rotation);
    vtkQuaterniond quat = vtkQuaterniond(cos(rotation[0]/2 *PI/180), rotation[1] * sin(rotation[0]/2 *PI/180), rotation[2] * sin(rotation[0]/2 *PI/180), rotation[3] * sin(rotation[0]/2 *PI/180));
    
    double ax[3];
    double angle = quat.GetRotationAngleAndAxis(ax);
    
    double array[4]; quat.Get(array);
    WorldToBoneRotation->SetTuple(boneid, array);
    WorldToBoneRestRotation->SetTuple(boneid, array);
}

// Transformation matrix only for anthropometric scaling
void vtkSkeleton::UpdateReferenceToPoseTransform(vtkIdType boneid){
    vtkQuaterniond worldtorest = (*WorldToBoneRestRotation)[boneid]; 
    vtkQuaterniond worldtorest_inv = worldtorest.Inverse();
    
    double rotation_torest[3][3]; worldtorest.ToMatrix3x3(rotation_torest);
    double rotation_torest_inv[3][3]; worldtorest_inv.ToMatrix3x3(rotation_torest_inv);
    
    // Get anthropometric scaling
    double scaling[3];
    AnthroFactor->GetTuple(boneid, scaling);
    double scalingmat[3][3] = {{scaling[0],0,0},{0,scaling[1],0},{0,0,scaling[2]}};
    
    // Combine scaling and rotation
    double transformation[3][3];
    vtkMath::Multiply3x3(scalingmat, rotation_torest_inv, scalingmat);
    vtkMath::Multiply3x3(rotation_torest, scalingmat, transformation);
    
    // Convert 2D array to 1D array and store
    double *transf = (double*)&transformation[0][0];
    (*ReferenceToPoseTransform)[boneid] = transf;
}

void vtkSkeleton::InitRestToPoseRotation(int boneid){
    vtkQuaterniond quat;
    quat.Identity(); double array[4]; quat.Get(array);
    RestToPoseRotation->SetTuple(boneid, array);
    
    double idtransf[9] = {1,0,0,0,1,0,0,0,1};
    ReferenceToPoseTransform->SetTuple(boneid, idtransf);
}

//----------------------------------------------------------------------------
vtkQuaterniond vtkSkeleton::RotationFromReferenceAxis(double referenceAxis[3], double axis[3])
{
  vtkQuaterniond newOrientation;
  // Code greatly inspired by: http://www.fastgraph.com/makegames/3drotation/ .

  // Normalize. This is the unit vector in the "new Z" direction.
  const double epsilon = 1e-6;
  if (vtkMath::Norm(referenceAxis) < epsilon
    || vtkMath::Normalize(axis) < epsilon)
    {
    return newOrientation;
    }

  // The dot product of axis and the referenceAxis gives projection of
  // of axis on referenceAxis.
  double projection = vtkMath::Dot(axis, referenceAxis);

  // First try at making a View Up vector: use World Up.
  double viewUp[3];
  viewUp[0] = referenceAxis[0] - projection*axis[0];
  viewUp[1] = referenceAxis[1] - projection*axis[1];
  viewUp[2] = referenceAxis[2] - projection*axis[2];

  // Check for validity:
  double magnitude = vtkMath::Normalize(viewUp);
  if (magnitude < epsilon)
    {
    // Second try: Use Y axis default  (0,1,0).
    viewUp[0] = -axis[1]*axis[0];
    viewUp[1] = 1-axis[1]*axis[1];
    viewUp[2] = -axis[1]*axis[2];

    // Check for validity:
    magnitude = vtkMath::Normalize(viewUp);

    if (magnitude < epsilon)
      {
      // Final try: Use Z axis default  (0,0,1).
      viewUp[0] = -axis[2]*axis[0];
      viewUp[1] = -axis[2]*axis[1];
      viewUp[2] = 1-axis[2]*axis[2];

      // Check for validity:
      magnitude = vtkMath::Normalize(viewUp);

      if (magnitude < epsilon)
        {
        return newOrientation;
        }
      }
    }

  // Calculate the Right vector. Use cross product of axis and Up.
  double viewRight[3];
  vtkMath::Cross(viewUp, axis, viewRight);
  vtkMath::Normalize(viewRight); //Let's be paranoid about the normalization.

  // Get the rest transform matrix.
  newOrientation.SetRotationAngleAndAxis(acos(projection), viewRight);
  return newOrientation.Normalized();
}

void vtkSkeleton::UpdateParentToOffsetRotation(vtkIdType boneid){
  // recalculate world->offset rotation and offset length
  double offset = 0;
  vtkQuaterniond quat = GetWorldToOffsetRotation(boneid, &offset);
  (*WorldToOffsetRotation)[boneid] = quat;
  
  // recalculate parent->offset rotation
  vtkQuaterniond quat_off = GetParentToOffsetRotation(boneid);
  (*ParentToOffsetRotation)[boneid] = quat_off;
  
  // update offset length
  double fac[3]; (*AnthroFactor)[GetParent(boneid)] > fac;
  if((*BoneOffsetRest)[boneid] > 0 ){(*BoneOffsetScale_length)[boneid] = offset / ((*BoneOffsetRest)[boneid] * fac[1]);}
  (*BoneOffset)[boneid] = (*BoneOffsetRest)[boneid] * fac[1] * (*BoneOffsetScale_length)[boneid];
}

vtkQuaterniond vtkSkeleton::GetWorldToOffsetRotation(vtkIdType boneid, double* offsetlength){
  double head_bone[3]; double tail_parent[3] = {0,0,0};
  (*HeadPose)[boneid] > head_bone;
  if (GetParent(boneid) != -1){ (*TailPose)[GetParent(boneid)] > tail_parent; }
  
  *offsetlength = 0;
  vtkQuaterniond quat; quat.Identity();
  
  if(!CompareVector3(head_bone, tail_parent)){
      double diff[3];
      vtkMath::Subtract(head_bone, tail_parent, diff);
      *offsetlength = vtkMath::Normalize(diff);
      quat = RotationFromReferenceAxis(this->GlobalAxisY, diff); 
  }
  
  return quat;
}



void vtkSkeleton::SetParentToBoneOffset(vtkIdType boneid){
    double offset = 0;
    vtkQuaterniond quat = GetWorldToOffsetRotation(boneid, &offset);
    
    BoneOffsetScale_radius->SetValue(boneid, 1);
    BoneOffsetScale_length->SetValue(boneid, 1);
    BoneOffset->SetValue(boneid, offset);
    BoneOffsetRest->SetValue(boneid, offset);
    double array[4]; quat.Get(array);
    WorldToOffsetRotation->SetTuple(boneid, array);   
    
    double length = GetLength(boneid);
    double scale = 1;
    BoneScale_length->SetTuple(boneid, &scale);
    BoneScale_thick->SetTuple(boneid, &scale);
    BoneLength->SetTuple(boneid, &length);
    BoneLengthRest->SetTuple(boneid, &length);
}

vtkQuaterniond vtkSkeleton::GetParentToBoneRotation(vtkIdType boneid){
    vtkQuaterniond WorldToParent; WorldToParent.Identity();
    if(GetParent(boneid) != -1){ WorldToParent = (*WorldToBoneRotation)[GetParent(boneid)]; }
    vtkQuaterniond WorldToBone = (*WorldToBoneRotation)[boneid];
    
    vtkQuaterniond quat = WorldToParent.Inverse() * WorldToBone;
    quat.Normalize();
    
    return quat;
}

vtkQuaterniond vtkSkeleton::GetParentToOffsetRotation(vtkIdType boneid){
    vtkQuaterniond WorldToParent; WorldToParent.Identity();
    if(GetParent(boneid) != -1){ WorldToParent = (*WorldToBoneRotation)[GetParent(boneid)]; }
    vtkQuaterniond WorldToOffset = (*WorldToOffsetRotation)[boneid];
    
    vtkQuaterniond quat = WorldToParent.Inverse() * WorldToOffset;
    quat.Normalize();    
    
    return quat;
}

void vtkSkeleton::SetParentToBoneRotations(vtkIdType boneid){ 
    vtkQuaterniond quat_bone = GetParentToBoneRotation(boneid);
    double array_bone[4]; quat_bone.Get(array_bone);
    ParentToBoneRotation->SetTuple(boneid, array_bone);
    
    vtkQuaterniond quat_off = GetParentToOffsetRotation(boneid);
    double array_off[4]; quat_off.Get(array_off);
    ParentToOffsetRotation->SetTuple(boneid, array_off);
}

void  vtkSkeleton::ProjectPointOnPlane(double orig[3], double (evec)[3], double point[3], double (*proj)[3]){
    double normal[3]; CopyVector3(evec, normal);
    double d[3]; vtkMath::Subtract(point,orig, d);
    vtkMath::MultiplyScalar(normal, vtkMath::Dot(d, normal));
    vtkMath::Subtract(point, normal, *proj);
}

void vtkSkeleton::SetEulerAngles(vtkIdType boneid)
{
    double angles[3];
    GetEulerAngles(boneid, &angles);
    BoneAngles->SetTuple(boneid, angles);  
    BoneAnglesRest->SetTuple(boneid, angles);
}

void vtkSkeleton::UpdateEulerAngles(vtkIdType boneid)
{
    double angles[3];
    GetEulerAngles(boneid, &angles);
    (*BoneAngles)[boneid] = angles;  
}

void vtkSkeleton::GetEulerAngles(vtkIdType boneid, double (*Angles)[3]){
  
    double relbonepos[3];
    GetRelativeBonePosition(boneid, &relbonepos);
    
    double sagittalpoint[3];
    ProjectPointOnPlane(Origin, this->GlobalAxisZ, relbonepos, &sagittalpoint);

    vtkMath::Normalize(relbonepos);
    (*Angles)[0] = -acos(vtkMath::Dot(relbonepos, this->GlobalAxisZ)) + vtkMath::Pi()/2;// beta
    
    //std::cout << "AxisZ " <<  this->GlobalAxisZ[0] << " "<<  this->GlobalAxisZ[1] << " "<<  this->GlobalAxisZ[2] << std::endl;
    
    //(*Angles)[0] = acos(vtkMath::Norm(sagittalpoint)); // beta
    
    //if (DotProd(&relbonepos, &this->GlobalAxisZ) < 0){(*Angles)[0]=-(*Angles)[0];}
    
    vtkMath::Normalize(sagittalpoint);
    (*Angles)[1] = acos(DotProd(&sagittalpoint, &this->GlobalAxisY)); // alpha
    
    if (DotProd(&sagittalpoint, &this->GlobalAxisX) > 0){(*Angles)[1]=-(*Angles)[1];}
    
   
    // Roll bone over its longitudinal axis
    
    double axis1[3];
    double axis2[3];
    GetLocalRotationAxis(boneid, 0, &axis1);    
    GetLocalRotationAxis(boneid, 1, &axis2);  
    
    // Calculate the residue of the transformations
    vtkQuaterniond QuadBeta;
    QuadBeta.SetRotationAngleAndAxis((*Angles)[0], axis1); // Normalise rotation axis before plugging it into a quaternion!!!
    QuadBeta.Invert();
    QuadBeta.Normalize();
    
    vtkQuaterniond QuadAlfa;
    QuadAlfa.SetRotationAngleAndAxis((*Angles)[1], axis2); // Normalise rotation axis before plugging it into a quaternion!!!
    QuadAlfa.Invert();
    QuadAlfa.Normalize();
    
    vtkQuaterniond ParentToBone = (*ParentToBoneRotation)[boneid];
    
    vtkQuaterniond residue = QuadAlfa * QuadBeta * ParentToBone;
    double residueaxis[3];
    (*Angles)[2] = residue.GetRotationAngleAndAxis(residueaxis);
        
    if (DotProd(&residueaxis, &GlobalAxisY)<0){(*Angles)[2]=-(*Angles)[2]; vtkMath::MultiplyScalar(residueaxis, -1);}
    
    // Can happen that roll angle exceeds -pi, pi range. 
    while((*Angles)[2] < -PI){(*Angles)[2] += 2*PI;}
    while((*Angles)[2] >  PI){(*Angles)[2] -= 2*PI;}
    
    // residueaxis has either to be (0,0,0), (0,1,0) or (0,-1,0)
    if (CompareVector3(residueaxis, GlobalAxisY) && CompareVector3(residueaxis, Origin) ){
        throw std::runtime_error("Was not able to find the kinematical rotation angles!");
    }
}





/////////// CHANGE KINEMATICAL PARAMETERS ////////////////////////////////////////////////


void vtkSkeleton::UpdateParentToBonePoseTransform(vtkIdType boneid, double angles[3]){
    // This function updates the parent to bone transformation and calculates the angles again. It doesnot update the joint positions yet!
    
    // Set the bone in the rest position according to its parent
    vtkQuaterniond quat; quat.Identity();
    
    for(int it=0; it<3; it++){
        int angleidx = (it+2)%3;
        
        double axis[3] = {0,0,0};
        if(angleidx==0){axis[0] = 1;}
        if(angleidx==1){axis[2] = 1;}
        if(angleidx==2){axis[1] = 1;}
        
        //GetLocalRotationAxis(boneid, angleidx, &axis);

        // Get the roll matrix.
        vtkQuaterniond RotQuad;
        RotQuad.SetRotationAngleAndAxis(angles[angleidx] , axis);
        RotQuad.Normalize();
            
        quat = RotQuad * quat;
    }
    
    (*ParentToBoneRotation)[boneid] = quat;
}

void vtkSkeleton::UpdateWorldToBoneRotation(vtkIdType boneid){
    vtkQuaterniond ParentToBone = (*ParentToBoneRotation)[boneid];
    vtkQuaterniond ParentToOffset = (*ParentToOffsetRotation)[boneid];
    vtkQuaterniond WorldToRest = (*WorldToBoneRestRotation)[boneid];
    vtkQuaterniond WorldToParent;
    if(GetParent(boneid) != -1){WorldToParent = (*WorldToBoneRotation)[GetParent(boneid)];}
    else{WorldToParent.Identity();}
    
    // Be carefull on the order of quaternion multiplication !!!
    (*WorldToBoneRotation)[boneid] = WorldToParent * ParentToBone;
    
    // Update REST to POSE rotation
    (*RestToPoseRotation)[boneid] = WorldToParent * ParentToBone * WorldToRest.Inverse();   //ParentToRest.Inverse() * ParentToBone; 
    UpdateReferenceToPoseTransform(boneid);
    
    // Same for the offset
    (*WorldToOffsetRotation)[boneid] =  WorldToParent * ParentToOffset;
}

void vtkSkeleton::MultiplyVectorByQuaternion(vtkQuaterniond rot, double (*vec)[3], double (*res)[3]){
    double rotation[3][3];
    rot.ToMatrix3x3(rotation);
    vtkMath::Multiply3x3(rotation, *vec, *res);
}

void vtkSkeleton::RotateOffset(vtkIdType boneid){
    // Rotation
    vtkQuaterniond rot = (*WorldToOffsetRotation)[boneid];
    double vec[3] = {0, (*BoneOffset)[boneid], 0};
    MultiplyVectorByQuaternion(rot, &vec, &vec);
    
    // Translation
    double pos[3];
    (*TailPose)[GetParent(boneid)] > pos;
    vtkMath::Add(pos, vec, vec);
    
    // Result
    (*HeadPose)[boneid] = vec;
}

void vtkSkeleton::RotateRestOffset(vtkIdType boneid){
    // Rotation
    vtkQuaterniond quat1 = (*WorldToBoneRestRotation)[GetParent(boneid)];
    vtkQuaterniond quat2 = (*ParentToOffsetRotation)[boneid];
    vtkQuaterniond rot = quat1 * quat2; 
    double vec[3] = {0, (*BoneOffset)[boneid], 0};
    MultiplyVectorByQuaternion(rot, &vec, &vec);
    
    // Translation
    double pos[3];
    (*TailRest)[GetParent(boneid)] > pos;
    vtkMath::Add(pos, vec, vec);
    
    // Result
    (*HeadRest)[boneid] = vec;
}

void vtkSkeleton::RotateBone(vtkIdType boneid){
    // Rotation
    vtkQuaterniond rot = (*WorldToBoneRotation)[boneid];
    double vec[3] = {0, (*BoneLength)[boneid], 0};
    MultiplyVectorByQuaternion(rot, &vec, &vec);
    
    double axisZ[3];
    MultiplyVectorByQuaternion(rot, &GlobalAxisZ, &axisZ);
    (*AxisZ)[boneid] = axisZ;
    
    // Translation
    double pos[3];
    (*HeadPose)[boneid] > pos;
    vtkMath::Add(pos, vec, vec);
    
    // Result
    (*TailPose)[boneid] = vec;
}

void vtkSkeleton::RotateRestBone(vtkIdType boneid){
    // Rotation
    vtkQuaterniond rot = (*WorldToBoneRestRotation)[boneid];
    double vec[3] = {0, (*BoneLengthRest)[boneid], 0};
    MultiplyVectorByQuaternion(rot, &vec, &vec);
    
    // Translation
    double pos[3];
    (*HeadRest)[boneid] > pos;
    vtkMath::Add(pos, vec, vec);
    
    // Result
    (*TailRest)[boneid] = vec;
}



void vtkSkeleton::PositionUpdate(vtkIdType boneid){
    // First update the head position, 
    if(GetParent(boneid) != -1){ // the root offset remains constant
        if((*BoneOffset)[boneid] != 0){RotateOffset(boneid);} // in case the bone has an offset with respect to its parent
        else{double pos[3]; (*TailPose)[GetParent(boneid)] > pos; (*HeadPose)[boneid] = pos; } // in case there is no offset, just copy.
    }
    // Secondly update the tail position
    RotateBone(boneid);
}

void vtkSkeleton::RestPositionUpdate(vtkIdType boneid){
    // First update the head position, 
    if(GetParent(boneid) != -1){ // the root offset remains constant
        if((*BoneOffset)[boneid] != 0){RotateRestOffset(boneid);} // in case the bone has an offset with respect to its parent
        else{double pos[3]; (*TailRest)[GetParent(boneid)] > pos; (*HeadRest)[boneid] = pos; } // in case there is no offset, just copy.
    }
    
    // Secondly update the tail position
    RotateRestBone(boneid);
}




void vtkSkeleton::PrintBoneAngles(){
    for(vtkIdType id=0; id<this->GetNumberOfVertices(); id++){
        double triplet[3];
        (*BoneAngles)[id] > triplet;
        std::string name = BoneName_bvh->GetValue(id);
        std::cout << name << " " << triplet[0] << " " << triplet[1] << " " << triplet[2] << "\t" << (*BoneLength)[id] << std::endl;
    }
}

std::string vtkSkeleton::GetBoneName(vtkIdType id){
    return BoneName_bvh->GetValue(id);
}






void vtkSkeleton::GetRotationUnitAxis(vtkIdType boneid, vtkIdType angleid, double (*ax)[3]){
  vtkQuaterniond rot; rot.Identity();
  if(GetParent(boneid) != -1){ rot = (*WorldToBoneRotation)[GetParent(boneid)] ;}
  //else{rot = (*WorldToBoneRotation)[boneid] ;}
  
  double axis[3]; GetLocalRotationAxis(boneid, angleid, &axis); // ie rotation axis in parents local coordinate system
  MultiplyVectorByQuaternion(rot, &axis, ax);
}

void vtkSkeleton::GetRotationCenter(vtkIdType boneid, double (*origin)[3]){
  if(GetParent(boneid) != -1){ double center[3]; (*TailPose)[GetParent(boneid)] > center; CopyVector3(center, *origin);}
  else{GetOrigin(origin);}  
}

void vtkSkeleton::SetRotationAxes(vtkIdType boneid){
  double axisx[3];double axisy[3];double axisz[3];
  GetRotationUnitAxis(boneid, 0, &axisx); 
  GetRotationUnitAxis(boneid, 1, &axisy);
  GetRotationUnitAxis(boneid, 2, &axisz);
  double axes[9] = {axisx[0], axisx[1], axisx[2], axisy[0], axisy[1], axisy[2], axisz[0], axisz[1], axisz[2]};
  (*BoneRotAxes)[boneid] = axes;
}

void vtkSkeleton::InitRotationAxes(vtkIdType boneid){
  double axisx[3];double axisy[3];double axisz[3];
  GetRotationUnitAxis(boneid, 0, &axisx); 
  GetRotationUnitAxis(boneid, 1, &axisy);
  GetRotationUnitAxis(boneid, 2, &axisz);
  double axes[9] = {axisx[0], axisx[1], axisx[2], axisy[0], axisy[1], axisy[2], axisz[0], axisz[1], axisz[2]};
  BoneRotAxes->SetTuple( boneid, axes );
}

vtkSmartPointer<vtkPolyData> GetBasicSkeleton(vtkSmartPointer<vtkSkeleton> skeleton){
    size_t NBones = skeleton->GetNumberOfVertices();
    
    vtkSmartPointer<vtkPolyData> polydata = vtkSmartPointer<vtkPolyData>::New();
    vtkSmartPointer<vtkPoints> pts = vtkSmartPointer<vtkPoints>::New();
    vtkSmartPointer<vtkCellArray> cells = vtkSmartPointer<vtkCellArray>::New();
    
    vtkDataArray* ArrayHead = skeleton->GetVertexData()->GetArray("HeadPose");
    vtkDataArray* ArrayTail = skeleton->GetVertexData()->GetArray("TailPose");
    
    for (vtkIdType boneid=0; boneid<NBones; boneid++){
        double head[3];
        double tail[3];
        ArrayHead->GetTuple(boneid, head);
        ArrayTail->GetTuple(boneid, tail);
        
        pts->InsertNextPoint(head);
        pts->InsertNextPoint(tail);
        
        // Join the two points in one cell
        vtkSmartPointer<vtkLine> line1 = vtkSmartPointer<vtkLine>::New();
        line1->GetPointIds()->SetId(0, 2*boneid); 
        line1->GetPointIds()->SetId(1, 2*boneid+1); 
        
        cells->InsertNextCell(line1);
    }
    
    polydata->SetPoints(pts);
    polydata->SetLines(cells);
    
    return polydata;
}

void vtkSkeleton::ShallowCopy(vtkSkeleton* t){
  
    this->CheckedShallowCopy(t->graph);
    this->AddArmatureArrays();
    
    this->ParentToBoneRotation->SetNumberOfTuples(t->GetNumberOfVertices());
    this->WorldToBoneRestRotation->SetNumberOfTuples(t->GetNumberOfVertices());
    this->ParentToOffsetRotation->SetNumberOfTuples(t->GetNumberOfVertices());
    this->WorldToOffsetRotation->SetNumberOfTuples(t->GetNumberOfVertices());
    this->RestToPoseRotation->SetNumberOfTuples(t->GetNumberOfVertices());
    this->AxisZ->SetNumberOfTuples(t->GetNumberOfVertices());
    this->TailPose->SetNumberOfTuples(t->GetNumberOfVertices());
    this->HeadPose->SetNumberOfTuples(t->GetNumberOfVertices());
    this->TailRest->SetNumberOfTuples(t->GetNumberOfVertices());
    this->HeadRest->SetNumberOfTuples(t->GetNumberOfVertices());
    this->ReferenceToPoseTransform->SetNumberOfTuples(t->GetNumberOfVertices());
    this->AnthroFactor->SetNumberOfTuples(t->GetNumberOfVertices());
    this->BoneAngles->SetNumberOfTuples(t->GetNumberOfVertices());
    this->BoneAnglesRest->SetNumberOfTuples(t->GetNumberOfVertices());
    this->BoneAngleRanges->SetNumberOfTuples(t->GetNumberOfVertices());
    this->BoneRotAxes->SetNumberOfTuples(t->GetNumberOfVertices());
    
    this->BoneScale_length->SetNumberOfTuples(t->GetNumberOfVertices());
    this->BoneScale_thick->SetNumberOfTuples(t->GetNumberOfVertices());
    this->BoneLength->SetNumberOfTuples(t->GetNumberOfVertices());
    this->BoneLengthRest->SetNumberOfTuples(t->GetNumberOfVertices());
    this->BoneOffsetScale_radius->SetNumberOfTuples(t->GetNumberOfVertices());
    this->BoneOffsetScale_length->SetNumberOfTuples(t->GetNumberOfVertices());
    this->BoneOffset->SetNumberOfTuples(t->GetNumberOfVertices());
    this->BoneOffsetRest->SetNumberOfTuples(t->GetNumberOfVertices());
    
    this->BoneName->SetNumberOfTuples(t->GetNumberOfVertices());
    
    
    
    
    this->ParentToBoneRotation = t->GetParentToBoneRotation();
    this->WorldToBoneRestRotation= t->GetWorldToBoneRestRotation() ;
    this->ParentToOffsetRotation= t->GetParentToOffsetRotation() ;
    this->WorldToOffsetRotation= t->GetWorldToOffsetRotation() ;
    this->RestToPoseRotation= t->GetRestToPoseRotation() ;
    
    this->AxisZ=t->GetAxisZ() ;
    this->TailPose=t->GetTailPose()  ;
    this->HeadPose=t->GetHeadPose()  ;
    this->TailRest=t->GetTailRest()  ;
    this->HeadRest=t->GetHeadRest() ;
    this->ReferenceToPoseTransform=t->GetReferenceToPoseTransform()  ;
    this->AnthroFactor=t->GetAnthroFactor()  ;
    this->BoneAngles=t->GetBoneAngles() ;
    this->BoneAnglesRest=t->GetBoneAnglesRest()  ;
    this->BoneAngleRanges=t->GetBoneAngleRanges()  ;
    this->BoneRotAxes=t->GetBoneRotAxes()  ;
    
    this->BoneScale_length=t->GetBoneScale_length()  ;
    this->BoneScale_thick=t->GetBoneScale_thick()  ;
    this->BoneLength=t->GetBoneLength()  ;
    this->BoneLengthRest=t->GetBoneLengthRest()  ;
    this->BoneOffsetScale_radius=t->GetBoneOffsetScale_radius()  ;
    this->BoneOffsetScale_length=t->GetBoneOffsetScale_length() ;
    this->BoneOffset=t->GetBoneOffset()  ;
    this->BoneOffsetRest=t->GetBoneOffsetRest()  ;
    
    this->BoneName=t->GetBoneName() ;
    
  
}


vtkSmartPointer<vtkTransform> vtkSkeleton::GetTransformToParentsLocalFrame(vtkIdType boneid, bool wrtparent, bool getinv){
  
  if(boneid == this->GetNumberOfVertices())
    return NULL;
  
  vtkIdType parentid = boneid;
  double head[3];
  
  if( GetRoot() != boneid && wrtparent){ 
    parentid = GetParent(boneid);
    (*TailRest)[parentid] > head;
  }
  else{
    (*HeadRest)[parentid] > head;
  } 
    
  vtkQuaterniond quat; quat.Identity();
  if(boneid != GetRoot()){ (*WorldToBoneRestRotation)[parentid]; }
  //vtkQuaterniond quat = (*WorldToBoneRestRotation)[parentid]; // using this wil cause the model to drift, because you change the worltobonerest as well
  
  if(getinv){
    quat = quat.Inverse();
  }
  
  double rotation[3][3]; quat.ToMatrix3x3(rotation);
  
  if(getinv){
    vtkMath::Multiply3x3(rotation, head, head);
    vtkMath::MultiplyScalar(head, -1);
  }
  
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
  totaltransformmatrix->SetElement(0,3, head[0]);
  totaltransformmatrix->SetElement(1,3, head[1]);
  totaltransformmatrix->SetElement(2,3, head[2]);

  vtkSmartPointer<vtkTransform> transform = vtkSmartPointer<vtkTransform>::New();
  transform->SetMatrix(totaltransformmatrix);
  return transform;  
}

void vtkSkeleton::ExportSkeletonToTXT(std::string filename){
    std::fstream outputfile(filename.c_str(), ios::out);
    ExportBones(this->GetRoot(), outputfile);
    outputfile.close();
}

void vtkSkeleton::ExportSkeletonToBVH(std::string filename){
    std::fstream outputfile(filename.c_str(), ios::out);
    ExportBonesBVH(this->GetRoot(), outputfile);
    outputfile.close();
}

void vtkSkeleton::ExportBones(vtkIdType parentid, std::fstream& file){
    // get angles, length
    double angles[3];
    (*BoneAngles)[parentid] > angles;
    double scaling[3]; AnthroFactor->GetTuple(parentid, scaling);
    double length = (*BoneLength)[parentid];
    double thicknessfactor = scaling[0];
    double boneoffset = (*BoneOffsetScale_radius)[parentid];
    
    file << std::string(BoneName_bvh->GetValue(parentid)) << " " << angles[0] << " " <<angles[1] << " " << angles[2] << " " << length << " " << thicknessfactor << " " << boneoffset << "\n";
    
    if (this->GetNumberOfChildren(parentid) > 0){ // Preliminar condition, to distinguish radius-MC3 offset from MC3-MC* offsets.    
        vtkSmartPointer<vtkAdjacentVertexIterator> childreniterator = vtkSmartPointer<vtkAdjacentVertexIterator>::New();
        GetChildren(parentid, childreniterator);
        while(childreniterator->HasNext()){
            vtkIdType childid = childreniterator->Next();
            ExportBones(childid, file);
        }
    }
}

void vtkSkeleton::ExportBonesBVH(vtkIdType parentid, std::fstream& file){
    
    double head[3]; double tail[3]; double axis[3]; double rom[6];
    (*HeadPose)[parentid] > head;
    (*TailPose)[parentid] > tail;
    (*AxisZ)[parentid] > axis;
    (*BoneAngleRanges)[parentid] > rom;
    
    std::string pretabs = ""; for(vtkIdType i=0; i< GetLevel(parentid); i++){pretabs += "\t";}
    
    if(parentid == 0){file << "HIERARCHY\nROOT ";}
    else{file << pretabs << "JOINT ";}
    
    file << std::string(BoneName_bvh->GetValue(parentid)) << "\n";
    file << pretabs << "{\n";
    file << pretabs << "\tTAIL " << tail[0] << " " << tail[1] << " " << tail[2] << "\n"; 
    file << pretabs << "\tHEAD " << head[0] << " " << head[1] << " " << head[2] << "\n";
    file << pretabs << "\tROLLAXIS " << axis[0] << " " << axis[1] << " " << axis[2] << "\n"; 
    file << pretabs << "\tRANGE " << rom[0] << " " << rom[1] << " " << rom[2] << " " << rom[3] << " " << rom[4] << " " << rom[5] << "\n"; 
    
    if (this->GetNumberOfChildren(parentid) > 0){ // Preliminar condition, to distinguish radius-MC3 offset from MC3-MC* offsets.    
        vtkSmartPointer<vtkAdjacentVertexIterator> childreniterator = vtkSmartPointer<vtkAdjacentVertexIterator>::New();
        GetChildren(parentid, childreniterator);
        while(childreniterator->HasNext()){
            vtkIdType childid = childreniterator->Next();
            ExportBonesBVH(childid, file);
        }
    }
    
    file << pretabs << "}\n";
}
void vtkSkeleton::SwapBoneRestPose(vtkIdType boneid){
    double headrest[3]; double headpose[3]; double tailrest[3]; double tailpose[3];
    (*HeadPose)[boneid] > headpose; 
    (*HeadRest)[boneid] > headrest;
    (*TailPose)[boneid] > tailpose; 
    (*TailRest)[boneid] > tailrest; 
    vtkQuaterniond resttopose = (*RestToPoseRotation)[boneid];
    
    double scaling[3];
    AnthroFactor->GetTuple(boneid, scaling);
    
    vtkQuaterniond worldtopose = (*WorldToBoneRotation)[boneid];
    vtkQuaterniond worldtorest = (*WorldToBoneRestRotation)[boneid];
    
    // Invert anthropometric scaling matrix
    /*
    double scalingtransf[9];
    (*ReferenceToPoseTransform)[boneid] > scalingtransf;
    double (*scalingmat)[3] = (double (*)[3])&scalingtransf[0];
    double (scalingmat_inv)[3][3];
    vtkMath::Invert3x3(scalingmat, scalingmat_inv);
    double *transf = (double*)&scalingmat_inv[0][0];
    (*ReferenceToPoseTransform)[boneid] = transf;
    */
    
    
    // Swap 
    (*HeadRest)[boneid] = headpose;
    (*HeadPose)[boneid] = headrest;
    (*TailRest)[boneid] = tailpose;
    (*TailPose)[boneid] = tailrest;
    (*RestToPoseRotation)[boneid] = resttopose.Inverse();
    
    (*WorldToBoneRotation)[boneid] = worldtorest;
    (*WorldToBoneRestRotation)[boneid] = worldtopose;
    
    //double scaling_inv[3] = {1., 1., 1.}; // do not normalise scale
    //double scaling_inv[3] = {1., 1./scaling[1], 1.}; //only normalise length scale 
    double scaling_inv[3] = {1./scaling[0], 1./scaling[1], 1./scaling[2]}; // do scale normalisation
    AnthroFactor->SetTuple(boneid, scaling_inv);
    UpdateReferenceToPoseTransform(boneid);
    
    if (this->GetNumberOfChildren(boneid) > 0){
        vtkSmartPointer<vtkAdjacentVertexIterator> childreniterator = vtkSmartPointer<vtkAdjacentVertexIterator>::New();
        GetChildren(boneid, childreniterator);
        while(childreniterator->HasNext()){
            vtkIdType childid = childreniterator->Next();
            SwapBoneRestPose(childid);
        }
    }
}

/*
void vtkSkeleton::ExportBones(vtkIdType parentid, std::fstream& file, std::string prestring){
    std::string prestring_next = prestring + "\t";
    
    std::string type;
    if(prestring == ""){type = "ROOT ";}
    else{type = "JOINT ";}
    
    // get tail, head position
    double head[3]; double tail[3]; double axisz[3]; double rom[8];
    (*TailPose)[parentid] > tail;
    (*HeadPose)[parentid] > head;
    (*AxisZ)[parentid] > axisz;
    (*BoneAngleRanges)[parentid] > rom;
    
    file << prestring << type << std::string(BoneName_bvh->GetValue(parentid)) << "\n";
    file << prestring << "{\n";
    file << prestring_next << "TAIL " << tail[0] << " " << tail[1] << " " <<tail[2] << "\n";
    file << prestring_next << "HEAD " << head[0] << " " << head[1] << " " <<head[2] << "\n";
    file << prestring_next << "ROLLAXIS " << axisz[0] << " " << axisz[1] << " " <<axisz[2] << "\n";
    file << prestring_next << "RANGE " << rom[0] << " " << rom[1] << " " <<rom[2] << " " << rom[3] << " " <<rom[4] << " " << rom[5] << " " <<rom[6] << " " <<rom[7] << "\n";
    
    if (this->GetNumberOfChildren(parentid) > 0){ // Preliminar condition, to distinguish radius-MC3 offset from MC3-MC* offsets.    
        vtkSmartPointer<vtkAdjacentVertexIterator> childreniterator = vtkSmartPointer<vtkAdjacentVertexIterator>::New();
        GetChildren(parentid, childreniterator);
        while(childreniterator->HasNext()){
            vtkIdType childid = childreniterator->Next();
            ExportBones(childid, file, prestring_next);
        }
    }
    file << prestring << "}\n";
}
*/



