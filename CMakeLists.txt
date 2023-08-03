cmake_minimum_required(VERSION 3.22.1)
project(cmake_test)

find_package(MPI REQUIRED)

set(CMAKE_CXX_STANDARD 20)
   
add_executable(Code main.cpp)

target_link_libraries(Code PUBLIC MPI::MPI_CXX)

add_subdirectory(KaGen)
target_link_libraries(Code PUBLIC KaGen::KaGen)
