cmake_minimum_required(VERSION 3.22.1)
project(cmake_test)

find_package(MPI REQUIRED)
# https://stackoverflow.com/questions/23163075/how-to-compile-an-mpi-included-c-program-using-cmake

set(CMAKE_CXX_STANDARD 11)
   
add_executable(cmake_test main.cpp)

target_link_libraries(cmake_test PUBLIC MPI::MPI_CXX)
