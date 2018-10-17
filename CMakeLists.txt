cmake_minimum_required(VERSION 3.12)
project(Ray_Tracing)

set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_FLAGS "-fopenmp")

include_directories(../Libs/glm)
include_directories(../Libs/assimp)
include_directories(../Libs)

add_executable(Ray_Tracing main.cpp config.h)