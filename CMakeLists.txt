CMAKE_MINIMUM_REQUIRED(VERSION 3.0)
PROJECT(Symbolic C CXX)

if (MSVC)
  add_compile_options(/std:c++17 /wd4267 /wd4244 /D_CRT_SECURE_NO_WARNINGS)
  include_directories(C:/Users/herholz/Documents/libigl/include)
else()
  add_compile_options(-std=c++17 -O3)
  include_directories(/usr/local/include/)
endif()

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR})
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR})
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_CACHEFILE_DIR ${CMAKE_BINARY_DIR}/build)

find_package (Eigen3 3.3 REQUIRED NO_MODULE)


add_subdirectory(src)

add_subdirectory(src/Examples/BuildCotan)

add_subdirectory(src/Examples/BuildDualCotan)

add_subdirectory(src/Examples/BaraffWitkin)

add_subdirectory(src/Examples/BaraffWitkinAD)

