cmake_minimum_required(VERSION 2.8)
project(ds)
find_package(VTK REQUIRED)
include(${VTK_USE_FILE})
set(srcs ds.cpp vtkDepthSortPolyData2.cxx)
include_directories(${CMAKE_CURRENT_SOURCE_DIR})
add_executable(ds ${srcs})
target_link_libraries(ds ${VTK_LIBRARIES})
