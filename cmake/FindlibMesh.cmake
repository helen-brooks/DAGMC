if (LIBMESH_PREFIX)  
  execute_process(COMMAND ${LIBMESH_PREFIX}/bin/libmesh-config --version
    OUTPUT_VARIABLE libMesh_VERSION OUTPUT_STRIP_TRAILING_WHITESPACE)
  execute_process(COMMAND ${LIBMESH_PREFIX}/bin/libmesh-config --include
    OUTPUT_VARIABLE libMesh_INCLUDE OUTPUT_STRIP_TRAILING_WHITESPACE)
  execute_process(COMMAND ${LIBMESH_PREFIX}/bin/libmesh-config --libs
    OUTPUT_VARIABLE libMesh_LIBRARIES OUTPUT_STRIP_TRAILING_WHITESPACE)
  execute_process(COMMAND ${LIBMESH_PREFIX}/bin/libmesh-config --cppflags --cxxflags
    OUTPUT_VARIABLE libMesh_CPPFLAGS OUTPUT_STRIP_TRAILING_WHITESPACE)

  set(libMesh_INC_DIR "${LIBMESH_PREFIX}/include")

  if (BUILD_STATIC_LIBS)
    set(libMesh_LIBRARIES_STATIC ${libMesh_LIBRARIES})
  endif ()
  if (BUILD_SHARED_LIBS)
    set(libMesh_LIBRARIES_SHARED ${libMesh_LIBRARIES})
  endif ()
  set(LibMesh_LIBRARIES)

  
  message(STATUS "set libMesh version : ${libMesh_VERSION}")
  message(STATUS "set libMesh include : ${libMesh_INCLUDE}")
  message(STATUS "set libMesh static libs : ${libMesh_LIBRARIES_STATIC}")
  message(STATUS "set libMesh shared libs : ${libMesh_LIBRARIES_SHARED}")
  message(STATUS "set libMesh cpp flags : ${libMesh_CPPFLAGS}")

  include_directories(${libMesh_INC_DIR})
  list(APPEND MESH_INCLUDE_DIRS ${libMes_INC_DIR})
  
else ()
  message(FATAL_ERROR "Could not find libMesh. Set -DLIBMESH_PREFIX when using cmake to the top-level installation directory for libMesh (contains bin/libmesh-config)")
endif ()
