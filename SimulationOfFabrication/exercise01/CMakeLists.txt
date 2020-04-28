cmake_minimum_required(VERSION 3.10)

project (assignement1)

set(CMAKE_CXX_FLAGS "-std=c++11 -O3")

set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)

add_library(Grid src/Grid.cpp src/GridPoint.cpp)
add_library(Advancer src/Advancer.cpp)

include_directories(include)

add_executable(grid src/grid_exe.cpp)
target_link_libraries(grid Grid)

add_executable(advancer src/advancer_exe.cpp)
target_link_libraries(advancer Grid Advancer)
