CMAKE_MINIMUM_REQUIRED(VERSION 2.6)
PROJECT(diagnsticIndex)

FIND_PACKAGE(statismo REQUIRED)
FIND_PACKAGE(VTK REQUIRED)
include (${VTK_USE_FILE})

include_directories(${statismo_INCLUDE_DIRS})
link_directories( ${VTK_LIBRARY_DIRS} ${statismo_LIBRARY_DIRS})

add_subdirectory( avgAnalysis )

add_executable (diagnosticIndex DiagnosticIndex.cxx)
target_link_libraries (diagnosticIndex ${VTK_LIBRARIES} ${statismo_LIBRARIES})
