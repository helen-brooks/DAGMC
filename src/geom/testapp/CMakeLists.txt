set(LINK_LIBS geom)
set(LINK_LIBS_EXTERN_NAMES)

include_directories(${CMAKE_SOURCE_DIR}/src/geom)
include_directories(${CMAKE_BINARY_DIR}/src/geom)

set(SRC_FILES test.cpp)
dagmc_install_exe(geom-test)
