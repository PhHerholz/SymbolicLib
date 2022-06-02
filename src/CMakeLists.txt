CMAKE_MINIMUM_REQUIRED(VERSION 3.0)
PROJECT(Symbolic C CXX)

set(CMAKE_SUPPRESS_REGENERATION true)
set(CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/../cmake")
set(SOURCE_DIR "${PROJECT_SOURCE_DIR}")

## check if openmp or cuda is available
find_package(OpenMP)
if (OPENMP_FOUND)
    set(CMAKE_FOUND_OMP "true")
else()
    set(CMAKE_FOUND_OMP "false")
endif()

if(DEFINED CMAKE_CUDA_STANDARD)
    set(CMAKE_FOUND_CUDA "true")
else()
    set(CMAKE_FOUND_CUDA "false")
endif()
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/ComputeUnitAvailability.h.in ${CMAKE_CURRENT_SOURCE_DIR}/ComputeUnitAvailability.h)

file(GLOB sources ${SOURCE_DIR}/*.cpp ${SOURCE_DIR}/*.hpp ${SOURCE_DIR}/*.h	)
include_directories(${INCLUDE_DIRS} ${SOURCE_DIR})

add_library(Symbolic ${sources})
target_compile_features(Symbolic PUBLIC cxx_std_17)
target_link_libraries(Symbolic ${LIBRARIES} Eigen3::Eigen -lpthread -ldl -lm -pthread)