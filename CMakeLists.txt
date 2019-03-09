cmake_minimum_required(VERSION 3.10)
project(Ray_Tracing)

set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_FLAGS "-fopenmp")

if (${CMAKE_BUILD_TYPE} STREQUAL  "debug")
    add_compile_definitions(THREADS_TO_RUN=1)
else()
    add_compile_definitions(THREADS_TO_RUN=4)
endif()

find_package(glm REQUIRED)
include_directories(${GLM_INCLUDE_DIRS})

add_executable(Ray_Tracing main.cpp scene.cpp)