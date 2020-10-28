

#include <vtkSTLReader.h>
#include <vtkPolyDataReader.h>
#include <vtkPolyDataNormals.h>
#include <vtkPointData.h>
#include <vtkCompositePolyDataMapper2.h>
#include <vtkCompositeDataGeometryFilter.h>
#include <vtkMultiBlockDataSet.h>
#include <vtkCompositeDataDisplayAttributes.h>
#include <iostream>
#include <exception>
#include <vtkXMLMultiBlockDataWriter.h>

#include <vtkSmartPointer.h>
#include <vtkCellArray.h>
#include <vtkTriangle.h>

#include <vtkPolyDataMapper.h>
#include <vtkActor.h>
#include <vtkRenderWindow.h>
#include <vtkRenderer.h>
#include <vtkRenderWindowInteractor.h>
#include <vtkPolyData.h>
#include <vtkCallbackCommand.h>
#include <vtkCommand.h>
#include <vtkProperty.h>
#include <vtkPolyDataWriter.h>
#include "vtkCamera.h"

#include "vtkLine.h"
#include "vtkCellData.h"
#include "vtkTextActor.h"
#include "vtkTextProperty.h"

/*
// VTK includes
#include <vtkTransformPolyDataFilter.h>
#include <vtkAbstractWidget.h>
#include <vtkCommand.h>
#include <vtkExtractEdges.h>

#include <vtkStdString.h>
#include <vtkLine.h>
#include <vtkLineRepresentation.h>
#include <vtkLineWidget2.h>
#include <vtkCellData.h>
#include <vtkPolyDataMapper.h>
#include <vtkAxesActor.h>
#include <vtkCallbackCommand.h>
#include <vtkCamera.h>
#include <vtkCaptionActor2D.h>
#include <vtkHandleRepresentation.h>
#include <vtkHandleWidget.h>
#include <vtkMath.h>
#include <vtkMatrix3x3.h>
#include <vtkNew.h>
#include <vtkObjectFactory.h>
#include <vtkPointHandleRepresentation3D.h>
#include <vtkProperty.h>
#include <vtkRenderer.h>
#include <vtkRenderWindow.h>
#include <vtkRenderWindowInteractor.h>
#include <vtkSmartPointer.h>
#include <vtkTransform.h>
#include <vtkWidgetCallbackMapper.h>
#include <vtkWidgetEvent.h>
#include <vtkDoubleArray.h>
#include <vtkPolyDataWriter.h>
#include "vtkInformation.h"
#include "vtkInformationVector.h"
#include <vtkTriangle.h>
#include <vtkPointData.h>

#include <vtkTriangleMeshPointNormals.h>

// STD includes
#include <cassert>
*/














