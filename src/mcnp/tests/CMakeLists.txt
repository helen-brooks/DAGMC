set(DRIVERS
  mcnp5_unit_test_driver.cc
  $<TARGET_OBJECTS:mcnp_funcs>
)

set(LINK_LIBS dagmc pyne_dagmc uwuw)
set(LINK_LIBS_EXTERN_NAMES)

include_directories(${GTEST_INCLUDE_DIR})
include_directories(${CMAKE_SOURCE_DIR}/src/mcnp)

dagmc_install_test(dagmcnp_unit_tests cpp)

dagmc_install_test_file(test_geom_legacy.h5m)
dagmc_install_test_file(test_geom_legacy_comp.h5m)
dagmc_install_test_file(test_reflecting.h5m)
