CMAKE_MINIMUM_REQUIRED(VERSION 3.0)
PROJECT(Symbolic C CXX)

set(CMAKE_SUPPRESS_REGENERATION true)
set(CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/../cmake")
set(SOURCE_DIR "${PROJECT_SOURCE_DIR}")

file(GLOB sources ${SOURCE_DIR}/*.cpp ${SOURCE_DIR}/*.hpp ${SOURCE_DIR}/*.h	)
include_directories(${INCLUDE_DIRS} ${SOURCE_DIR})

add_library(Symbolic ${sources})
target_compile_features(Symbolic PUBLIC cxx_std_17)
target_link_libraries(Symbolic ${LIBRARIES} Eigen3::Eigen -lpthread -ldl -lm -pthread)


