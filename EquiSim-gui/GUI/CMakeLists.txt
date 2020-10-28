cmake_minimum_required (VERSION 2.8.12)

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${CMAKE_CURRENT_SOURCE_DIR}/../cmake) # main (top) cmake dir
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${CMAKE_CURRENT_SOURCE_DIR}/cmake) # project specific cmake dir
set(CMAKE_CXX_STANDARD 11) # tODO move up to a general cmake config for all sub projects ?

# CMake useful variables
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/bin")
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/lib") 
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/lib")

project (EquiSim)

find_package (Qt5Widgets)

#if(NOT VTK_BINARY_DIR)
#	find_package(VTK REQUIRED)
#	include(${VTK_USE_FILE})
#endif()


Find_package(VTK 7 CONFIG)

if (VTK_VERSION VERSION_LESS "8.90")
  # old system
  Find_package(VTK REQUIRED)
  include(${VTK_USE_FILE})
  # modules are linked via `vtkCommonCore`
  # VTK_DEFINITIONS has autoinit information
else ()
  # modules are linked via `VTK::CommonCore`
  # vtk_module_autoinit is needed
  Find_package(VTK COMPONENTS
  ChartsCore
  CommonColor
  CommonComputationalGeometry
  CommonSystem
  FiltersExtraction
  FiltersFlowPaths
  FiltersGeneric
  FiltersGeometry
  FiltersHybrid
  FiltersModeling
  # FiltersParallelStatistics #  Optional
  FiltersPoints
  FiltersProgrammable
  FiltersTexture
  FiltersVerdict
  GeovisCore
  GUISupportQt
  GUISupportQtOpenGL
  hdf5
  hdf5_hl
  IOExodus
  IOExport
  # IOFFMPEG # Optional
  IOImport
  IOInfovis
  # IOMySQL # Optional
  IOPLY
  IOParallel
  IOParallelXML
  ImagingColor
  ImagingFourier
  ImagingHybrid
  ImagingMath
  ImagingMorphological
  ImagingSources
  ImagingStatistics
  ImagingStencil
  # InfovisBoostGraphAlgorithms # Optional
  InfovisLayout
  InteractionImage
  IOXML
  IOXMLParser
  RenderingAnnotation
  RenderingImage
  RenderingLOD
  RenderingLabel
  RenderingGL2PSOpenGL2
  RenderingVolume
  RenderingVolumeOpenGL2
  TestingRendering
  ViewsContext2D
  ViewsInfovis
)
endif () 

add_dependent_subproject(SKELETON)

enable_testing()
add_subdirectory(include)
add_subdirectory(src)
add_subdirectory(doc)


