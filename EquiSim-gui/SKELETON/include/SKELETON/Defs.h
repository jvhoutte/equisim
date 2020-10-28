// Bender includes

#include <string>

// VTK includes
#include <vtkLandmarkTransform.h>
#include <vtkTransformPolyDataFilter.h>
#include <vtkExtractEdges.h>
#include <vtkSmartPointer.h>
#include <vtkStdString.h>
#include <vtkLine.h>
#include <vtkPolyLine.h>
#include <vtkWarpVector.h>
#include <vtkCellData.h>
#include <vtkMath.h>
#include <vtkMatrix3x3.h>
#include <vtkNew.h>
#include <vtkSmartPointer.h>
#include <vtkTransform.h>
#include <vtkDoubleArray.h>
#include <vtkFloatArray.h>
#include <vtkStringArray.h>
#include <vtkPolyDataWriter.h>
#include <vtkTriangle.h>
#include <vtkPointData.h>
#include "vtkAppendPolyData.h"
#include "vtkUnsignedCharArray.h"
#include <vtkTypedArray.h>
#include "vtkQuaternion.h"
#include "vtkWeightedTransformFilter.h"
#include "vtkUnstructuredGrid.h"
#include "vtkPoints.h"
#include "vtkPolyDataReader.h"
#include "vtkMinimalStandardRandomSequence.h"
#include "vtkPolyDataAlgorithm.h"
#include "vtkXMLMultiBlockDataReader.h"
#include "vtkTreeDFSIterator.h"
#include "vtkSTLReader.h"
#include "vtkTetra.h"
#include "vtkXMLUnstructuredGridWriter.h"
#include "vtkXMLMultiBlockDataWriter.h"
#include "vtkImplicitPolyDataDistance.h"
#include "vtkKdTreePointLocator.h"
#include "vtkObjectFactory.h"
#include "vtkStreamingDemandDrivenPipeline.h"
#include "vtkInformationVector.h"
#include "vtkInformation.h"
#include "vtkDataObject.h"
#include "vtkSmartPointer.h"
#include "vtkDoubleArray.h"
#include "vtkCellData.h"

// STD includes
#include <cassert>
#include <string>
#include <stdexcept>

// Eigen includes
/*
#include <eigen3/Eigen/Core>
#include <eigen3/Eigen/Eigen>
#include <eigen3/Eigen/QR>
*/
