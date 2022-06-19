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

find_package(HIP)

if(HIP_FOUND)
    set(CMAKE_FOUND_HIP "true")
else()
    set(CMAKE_FOUND_HIP "false")
endif()

configure_file(${CMAKE_CURRENT_SOURCE_DIR}/ComputeUnitAvailability.h.in ${CMAKE_CURRENT_SOURCE_DIR}/ComputeUnitAvailability.h)

## set various flags for different systems
if (MSVC)
    get_filename_component(MSVC_DIR ${CMAKE_CXX_COMPILER} DIRECTORY)
    file(GLOB_RECURSE VC_VARS "C://Program Files (x86)/Microsoft Visual Studio/" vcvars64.bat  REQUIRED)
    message("msvc_dir: " ${VC_VARS})
    add_compile_options(/std:c++17 /wd4267 /wd4244 /D_CRT_SECURE_NO_WARNINGS)
elseif (APPLE)
    set(OpenMP_CXX_FLAGS "-Xpreprocessor -fopenmp -L /usr/local/opt/libomp/lib -I /usr/local/opt/libomp/include")
    set(OpenMP_CXX_LIB_NAMES "omp")
    set(OpenMP_omp_LIBRARY /usr/local/opt/libomp/lib/libomp.dylib)
endif()

file(GLOB sources ${SOURCE_DIR}/*.cpp ${SOURCE_DIR}/*.hpp ${SOURCE_DIR}/*.h)
include_directories(${INCLUDE_DIRS} ${SOURCE_DIR})

add_library(Symbolic ${sources})
target_compile_features(Symbolic PUBLIC cxx_std_17)

if (MSVC)
    target_link_libraries(Symbolic ${LIBRARIES} Eigen3::Eigen)
else ()
    target_link_libraries(Symbolic ${LIBRARIES} Eigen3::Eigen -lpthread -ldl -lm -pthread)
endif ()

