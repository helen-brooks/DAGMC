set(DRIVERS mesh_unit_test_driver.cc)

set(LINK_LIBS mesh)
set(LINK_LIBS_EXTERN_NAMES)

include_directories(${GTEST_INCLUDE_DIR})
include_directories(${CMAKE_SOURCE_DIR}/src/mesh)
include_directories(${CMAKE_BINARY_DIR}/src/mesh)

dagmc_install_test(libmesh_interface_test  cpp)

dagmc_install_test_file(sphere-withattributes.e)
