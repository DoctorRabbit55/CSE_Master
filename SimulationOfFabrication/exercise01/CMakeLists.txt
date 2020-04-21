cmake_minimum_required(VERSION 3.10)

project (assignement1)

set(CMAKE_CXX_FLAGS "-std=c++11")

set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)

add_library(Grid src/Grid.cpp src/GridPoint.cpp)

include_directories(include)

add_executable(main src/main.cpp)
target_link_libraries(main Grid)
