CMAKE_MINIMUM_REQUIRED(VERSION 3.0)
PROJECT(BaraffWitkinAD C CXX)

set(CMAKE_SUPPRESS_REGENERATION true)
set(CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/../../../cmake")
set(SOURCE_DIR "${PROJECT_SOURCE_DIR}")

file(GLOB sources ${SOURCE_DIR}/*.cpp ${SOURCE_DIR}/*.hpp ${SOURCE_DIR}/*.h)

ADD_EXECUTABLE(baraffwitkinad ${sources})
include_directories(${INCLUDE_DIRS} ${SOURCE_DIR} ${SOURCE_DIR}/../../)
TARGET_LINK_LIBRARIES(baraffwitkinad Eigen3::Eigen Symbolic igl::core)

