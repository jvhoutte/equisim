# @author  Barthélémy von Haller
set(LIB_HEADERS
  "${CMAKE_CURRENT_BINARY_DIR}/../include/SKELETON/Version.h" # Generated header
  ../include/SKELETON/vtkSkeleton.h
  ../include/SKELETON/DeformRigidSurface.h
  ../include/SKELETON/vtkBlockColorSelector.h
  ../include/SKELETON/vtkSkeletonVisualiser.h
  ../include/SKELETON/vtkBVHReader.h
  ../include/SKELETON/vtkStatModelReader.h
  ../include/SKELETON/DLLDefines.h
  ../include/SKELETON/Defs.h
)

set(LIB_SOURCES
  vtkSkeleton.cpp
  DeformRigidSurface.cpp
  vtkBlockColorSelector.cpp
  vtkSkeletonVisualiser.cpp
  vtkBVHReader.cxx
  vtkStatModelReader.cpp
)

# Define two variables in order not to repeat ourselves.
set(LIBRARY_NAME SKELETON)

# Define the library
add_library(${LIBRARY_NAME} SHARED ${LIB_SOURCES})

# Set the build version. It will be used in the name of the lib, with corresponding
# symlinks created. SOVERSION could also be specified for api version. 
set_target_properties(${LIBRARY_NAME} PROPERTIES
  VERSION ${VERSION}  # ${VERSION} was defined in the main CMakeLists.
  FRAMEWORK FALSE
  PUBLIC_HEADER "${LIB_HEADERS}" # specify the public headers
  CUDA_SEPARABLE_COMPILATION ON
)

# Added ${Glue} and ${VTK_LIBRARIES} for QuickView visualisation of itk images
target_link_libraries(${LIBRARY_NAME} ${VTK_LIBRARIES} ) 

# Says how and where to install software
# Targets:
#   * <prefix>/lib/<libraries>
#   * header location after install: <prefix>/include/<project>/*.h
#   * headers can be included by C++ code `#<project>/Bar.hpp>`
install(TARGETS ${LIBRARY_NAME}
  EXPORT ${PROJECT_NAME}Targets            # for downstream dependencies
  ARCHIVE DESTINATION lib COMPONENT libs   # static lib
  RUNTIME DESTINATION bin COMPONENT libs   # binaries
  LIBRARY DESTINATION lib COMPONENT libs   # shared lib
  FRAMEWORK DESTINATION bin COMPONENT libs # for mac
  PUBLIC_HEADER DESTINATION include/${PROJECT_NAME} COMPONENT devel   # headers for mac (note the different component -> different package)
  INCLUDES DESTINATION include             # headers
)

