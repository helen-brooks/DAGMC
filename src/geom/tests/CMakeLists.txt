set(DRIVERS geom_unit_test_driver.cc)

set(LINK_LIBS geom)
set(LINK_LIBS_EXTERN_NAMES)

include_directories(${GTEST_INCLUDE_DIR})
include_directories(${CMAKE_SOURCE_DIR}/src/utils/tests)
include_directories(${CMAKE_SOURCE_DIR}/src/geom)
include_directories(${CMAKE_BINARY_DIR}/src/geom)

if(LIBMESH)
    dagmc_install_test(topotool_unit_test  cpp)
    dagmc_install_test_file(sphere.e)
endif()
