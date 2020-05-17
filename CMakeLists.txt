cmake_minimum_required(VERSION 3.6)
project(mpi_project)

include_directories(/usr/local/Cellar/open-mpi/4.0.3/include)
set(CMAKE_CXX_STANDARD 11)

add_executable(mpi_project main.cpp)