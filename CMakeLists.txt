cmake_minimum_required(VERSION 3.12)
project(Ray_Tracing)

set(CMAKE_CXX_STANDARD 14)

include_directories(../glm)

add_executable(Ray_Tracing main.cpp)

add_executable(Ray_Tracing2 main2.cpp)