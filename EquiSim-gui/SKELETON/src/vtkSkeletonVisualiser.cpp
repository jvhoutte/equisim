#include "SKELETON/vtkSkeletonVisualiser.h"
 
#include "vtkObjectFactory.h"
#include "vtkStreamingDemandDrivenPipeline.h"
#include "vtkInformationVector.h"
#include "vtkInformation.h"
#include "vtkDataObject.h"
#include "vtkSmartPointer.h"
#include "vtkDoubleArray.h"
#include "vtkCellData.h"

#include <thread>

vtkStandardNewMacro(vtkSkeletonVisualiser);

vtkSkeletonVisualiser::vtkSkeletonVisualiser()
{
  this->SetNumberOfInputPorts(1);
}
 
vtkSkeletonVisualiser::~vtkSkeletonVisualiser()
{
}
 
int vtkSkeletonVisualiser::FillInputPortInformation( int port, vtkInformation* info )
{
  if ( port == 0 )
    {
    info->Set(vtkDataObject::DATA_TYPE_NAME(), "vtkSkeleton" );
    return 1;
    }
  return 0;
}
 
 
int vtkSkeletonVisualiser::FillOutputPortInformation(int vtkNotUsed(portNumber), vtkInformation *info)
{
    info->Set(vtkDataObject::DATA_TYPE_NAME(), "vtkPolyData");
    return 1;
}
 
int vtkSkeletonVisualiser::RequestData(vtkInformation *vtkNotUsed(request),
                                          vtkInformationVector **inputVector,
                                          vtkInformationVector *outputVector)
{
    // get the input objects 
    
    this->armature = vtkTree::GetData(inputVector[0],0);
    this->NUMBER_OF_BONES = this->armature->GetNumberOfVertices();
    
    // get the output objects 
    vtkInformation *outInfo = outputVector->GetInformationObject(0);
    this->polyskeleton = vtkPolyData::SafeDownCast(outInfo->Get(vtkDataObject::DATA_OBJECT()));
    
    polyskeleton->SetPoints(GetSkeletonGeometry());
    polyskeleton->SetLines(GetSkeletonTopology());
    polyskeleton->GetCellData()->AddArray(GetSkeletonColors()); 
    polyskeleton->GetCellData()->SetActiveScalars("Colors");
    
    return 1;
}

vtkSmartPointer<vtkPoints> vtkSkeletonVisualiser::GetSkeletonGeometry() {
    vtkSmartPointer<vtkDoubleArray> tailarray = vtkDoubleArray::SafeDownCast( this->armature->GetVertexData()->GetArray("TailPose") );
    vtkSmartPointer<vtkDoubleArray> headarray = vtkDoubleArray::SafeDownCast( this->armature->GetVertexData()->GetArray("HeadPose") );
    vtkSmartPointer<vtkDoubleArray> axisarray = vtkDoubleArray::SafeDownCast( this->armature->GetVertexData()->GetArray("AxisZ") );
    
    vtkSmartPointer<vtkPoints> pts = vtkSmartPointer<vtkPoints>::New();

    for (vtkIdType boneid = 0; boneid < NUMBER_OF_BONES; boneid++) {
        double head[3]; headarray->GetTuple(boneid, head);
        double tail[3]; tailarray->GetTuple(boneid, tail);
        double axisZ[3]; axisarray->GetTuple(boneid, axisZ);

        double diff[3];
        vtkMath::Subtract(tail, head, diff);
        double L = vtkMath::Norm(diff);

        double scale = 8;

        double axisX[3];
        double axisY[3];
        // Z-axis
        vtkMath::Normalize(axisZ);
        vtkMath::MultiplyScalar(axisZ, scale);
        // Y-axis
        vtkMath::Subtract(tail, head, axisY);
        vtkMath::Normalize(axisY);
        vtkMath::MultiplyScalar(axisY, scale);
        // X-axis
        vtkMath::Cross(axisY, axisZ, axisX);
        vtkMath::Normalize(axisX);
        vtkMath::MultiplyScalar(axisX, scale);
        // Translate vectors
        vtkMath::Add(axisX, head, axisX);
        vtkMath::Add(axisY, head, axisY);
        vtkMath::Add(axisZ, head, axisZ);

        pts->InsertNextPoint(head);
        pts->InsertNextPoint(tail);
        pts->InsertNextPoint(axisX);
        pts->InsertNextPoint(axisY);
        pts->InsertNextPoint(axisZ);
    }

    return pts;
}

vtkSmartPointer<vtkUnsignedCharArray> vtkSkeletonVisualiser::GetSkeletonColors() {
    unsigned char black[3] = { 0, 0, 0 };
    unsigned char red[3] = { 255, 0, 0 };
    unsigned char green[3] = { 0, 255, 0 };
    unsigned char blue[3] = { 0, 0, 255 };

    vtkSmartPointer<vtkUnsignedCharArray> colors = vtkSmartPointer<vtkUnsignedCharArray>::New();
    colors->SetNumberOfComponents(3);
    colors->SetName("Colors");
    for (int k = 0; k < NUMBER_OF_BONES; k++) {
        //colors->InsertNextTupleValue(black);
        colors->InsertNextTupleValue(red);
        colors->InsertNextTupleValue(green);
        colors->InsertNextTupleValue(blue);
    }
    return colors;
}
 
vtkSmartPointer<vtkCellArray> vtkSkeletonVisualiser::GetSkeletonTopology() {
    vtkSmartPointer<vtkCellArray> cells = vtkSmartPointer<vtkCellArray>::New();

    for (vtkIdType boneid = 0; boneid < NUMBER_OF_BONES; boneid++) {
        // Join the two points in one cell
        //vtkSmartPointer<vtkLine> line1 = vtkSmartPointer<vtkLine>::New();
        //line1->GetPointIds()->SetId(0, 5 * boneid);
        //line1->GetPointIds()->SetId(1, 5 * boneid + 1);
        vtkSmartPointer<vtkLine> line2 = vtkSmartPointer<vtkLine>::New();
        line2->GetPointIds()->SetId(0, 5 * boneid);
        line2->GetPointIds()->SetId(1, 5 * boneid + 2);
        vtkSmartPointer<vtkLine> line3 = vtkSmartPointer<vtkLine>::New();
        line3->GetPointIds()->SetId(0, 5 * boneid);
        line3->GetPointIds()->SetId(1, 5 * boneid + 3);
        vtkSmartPointer<vtkLine> line4 = vtkSmartPointer<vtkLine>::New();
        line4->GetPointIds()->SetId(0, 5 * boneid);
        line4->GetPointIds()->SetId(1, 5 * boneid + 4);

        //cells->InsertNextCell(line1);
        cells->InsertNextCell(line2);
        cells->InsertNextCell(line3);
        cells->InsertNextCell(line4);
    }
    return cells;
}

void vtkSkeletonVisualiser::PrintSelf(ostream& os, vtkIndent indent)
{
  this->Superclass::PrintSelf(os,indent);
}
 
