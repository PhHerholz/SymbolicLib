include(libigl)
set(DATA_DIR ${CMAKE_CURRENT_SOURCE_DIR}/../data)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/dataPath.hpp.in ${CMAKE_CURRENT_SOURCE_DIR}/dataPath.hpp)

add_subdirectory(BuildCotan)
add_subdirectory(BuildDualCotan)
add_subdirectory(BaraffWitkin)
add_subdirectory(BaraffWitkinAD)
add_subdirectory(Tutorial)
