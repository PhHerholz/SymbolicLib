CMAKE_MINIMUM_REQUIRED(VERSION 3.0)
PROJECT(NeoHookean CXX)

set(CMAKE_SUPPRESS_REGENERATION true)
set(CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/../../../cmake")
set(SOURCE_DIR "${PROJECT_SOURCE_DIR}")

file(GLOB sources ${SOURCE_DIR}/*.cpp ${SOURCE_DIR}/*.hpp)

ADD_EXECUTABLE(NeoHookean ${sources})
include_directories(${INCLUDE_DIRS} ${SOURCE_DIR} ${SOURCE_DIR}/../../)


target_include_directories(NeoHookean PUBLIC /opt/homebrew/opt/libomp/include)
target_link_directories(NeoHookean PUBLIC /opt/homebrew/opt/libomp/lib)
# find_package(OpenMP)

# if(OpenMP_CXX_FOUND)
#     message("OPENMP FOUND")
#     target_link_libraries(NeoHookean PUBLIC OpenMP::OpenMP_CXX)
# endif()

TARGET_LINK_LIBRARIES(NeoHookean Eigen3::Eigen Symbolic -lpthread igl::core -lomp omp)
# target_compile_options(NeoHookean PUBLIC "-L/opt/homebrew/opt/libomp/lib -Wl,-rpath,/opt/intel/lib -lomp -I/opt/homebrew/opt/libomp/include")

