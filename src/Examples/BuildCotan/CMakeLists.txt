CMAKE_MINIMUM_REQUIRED(VERSION 3.0)
PROJECT(BuildCotan C CXX)

set(CMAKE_SUPPRESS_REGENERATION true)
set(CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/../../../cmake")
set(SOURCE_DIR "${PROJECT_SOURCE_DIR}")

file(GLOB sources ${SOURCE_DIR}/*.cpp ${SOURCE_DIR}/*.hpp)

add_executable(cotan ${sources})
include_directories(${INCLUDE_DIRS} ${SOURCE_DIR} ${SOURCE_DIR}/../../)
target_link_libraries(cotan Eigen3::Eigen Symbolic -lpthread)

