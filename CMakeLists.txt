if(HAS_PARENT)
  set(SYMBOLIC_TOPLEVEL_PROJECT OFF)
else()
  set(SYMBOLIC_TOPLEVEL_PROJECT ON)
endif()

## set the minimu, version so that some features of cmake can be used
set(REQUIRED_CMAKE_VERSION "3.14.0")
if(SYMBOLIC_TOPLEVEL_PROJECT)
  cmake_minimum_required(VERSION ${REQUIRED_CMAKE_VERSION})
else()
  # Don't use cmake_minimum_required here to avoid implicitly overriding parent policies
  if(${CMAKE_VERSION} VERSION_LESS ${REQUIRED_CMAKE_VERSION})
    message(FATAL_ERROR "CMake required version to build SYMBOLIC_LIB is ${REQUIRED_CMAKE_VERSION}")
  endif()
endif()

PROJECT(Symbolic C CXX)


set(DATA_DIR ${CMAKE_CURRENT_SOURCE_DIR}/data)
add_definitions( -DSYMBOLIC_DATA_DIR=${DATA_DIR} )


## include eigen
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${CMAKE_CURRENT_SOURCE_DIR}/cmake)
include(eigen)

if (MSVC)
  add_compile_options(/std:c++17 /wd4267 /wd4244 /D_CRT_SECURE_NO_WARNINGS)
else()
endif()

if (APPLE)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR})
else()
if (UNIX)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)
else()
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR})
endif()
endif()

set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR})
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_CACHEFILE_DIR ${CMAKE_BINARY_DIR}/build)


add_subdirectory(src)

if(SYMBOLIC_TOPLEVEL_PROJECT)
  add_subdirectory(Examples)
endif()


