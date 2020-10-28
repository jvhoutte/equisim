
set (MyProjectLib_src
    ${PROJECT_SOURCE_DIR}/src/mainwindow.cpp 
    ${PROJECT_SOURCE_DIR}/src/CorrelationFileReader.cpp
    )
set (MyProjectLib_hdr 
    ${PROJECT_SOURCE_DIR}/src/mainwindow.h
)

set (MyProjectLib_ui  ${PROJECT_SOURCE_DIR}/src/mainwindow.ui)
set (MyProjectBin_src ${PROJECT_SOURCE_DIR}/src/main.cpp)

qt5_wrap_cpp(MyProjectLib_hdr_moc ${MyProjectLib_hdr})
qt5_wrap_ui (MyProjectLib_ui_moc  ${MyProjectLib_ui})

include_directories (${PROJECT_SOURCE_DIR}/src ${PROJECT_SOURCE_DIR}/include ${PROJECT_SOURCE_DIR}/include)
include_directories (${PROJECT_BINARY_DIR})




add_executable(EquiSim
    main.cpp
    ${MyProjectLib_src}
    ${MyProjectLib_hdr_moc}
    ${MyProjectLib_ui_moc}
    ${MyProjectLib_src}
)

target_link_libraries (EquiSim SKELETON ${VTK_LIBRARIES} Qt5::Widgets )


message(${CMAKE_SOURCE_DIR}/InputData)
message(${CMAKE_BINARY_DIR}/bin)

add_custom_command(TARGET EquiSim POST_BUILD
       COMMAND ${CMAKE_COMMAND} -E copy_directory
       ${CMAKE_SOURCE_DIR}/InputData/ ${CMAKE_BINARY_DIR}/bin)


# Says how and where to install software
# Targets:
#   * <prefix>/lib/<libraries>
#   * header location after install: <prefix>/include/<project>/*.h
#   * headers can be included by C++ code `#<project>/Bar.hpp>`
install(TARGETS ${LIBRARY_NAME} EquiSim
  EXPORT ${PROJECT_NAME}Targets            # for downstream dependencies
  ARCHIVE DESTINATION lib COMPONENT libs   # static lib
  RUNTIME DESTINATION bin COMPONENT libs   # binaries
  LIBRARY DESTINATION lib COMPONENT libs   # shared lib
  FRAMEWORK DESTINATION bin COMPONENT libs # for mac
  PUBLIC_HEADER DESTINATION include/${PROJECT_NAME} COMPONENT devel   # headers for mac (note the different component -> different package)
  INCLUDES DESTINATION include             # headers
)

