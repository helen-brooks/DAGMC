file(GLOB SRC_FILES "*.cpp")
file(GLOB PUB_HEADERS "*.hpp")

# Set libs for linking
if (LIBMESH)
  set(LINK_LIBS mesh utils)
  set(LINK_LIBS_EXTERN_NAMES libMesh_LIBRARIES)
else()
  set(LINK_LIBS mesh)
  set(LINK_LIBS_EXTERN_NAMES )
endif()

include_directories(${CMAKE_SOURCE_DIR}/src/mesh)
include_directories(${CMAKE_BINARY_DIR}/src/mesh)
if (LIBMESH)
  include_directories(${CMAKE_BINARY_DIR}/src/utils)
endif()
include_directories(${CMAKE_BINARY_DIR}/src/geom)

dagmc_install_library(geom)

if (LIBMESH)
add_subdirectory(testapp)
endif()

if (BUILD_TESTS)
  add_subdirectory(tests)
endif ()
