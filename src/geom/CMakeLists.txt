file(GLOB SRC_FILES "*.cpp")
file(GLOB PUB_HEADERS "*.hpp")

# Set libs for linking
set(LINK_LIBS)
set(LINK_LIBS_EXTERN_NAMES libMesh_LIBRARIES)

include_directories(${CMAKE_BINARY_DIR}/src/geom)

dagmc_install_library(geom)

if (BUILD_TESTS)
  add_subdirectory(tests)
endif ()
