cmake_minimum_required(VERSION 3.10)
project(Ray_Tracing)

set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_FLAGS "-fopenmp")

find_package(glm REQUIRED)
include_directories(${GLM_INCLUDE_DIRS})

include_directories(../Libs/assimp)
include_directories(../Libs)

add_executable(Ray_Tracing main.cpp config.h ray.h material.h triangles.h)