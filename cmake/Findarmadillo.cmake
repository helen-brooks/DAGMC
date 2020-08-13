
set(armadillo_LIBDIR "/usr/lib")
if (ARMADILLO_DIR)
  # Allow user to specify the location of their armadillo lib
  set(armadillo_LIBDIR ${ARMADILLO_DIR})
  message(STATUS "User set armadillo library dir ${armadillo_LIBDIR}")
else()
  message(STATUS "Searching for armdillo library in default location ${armadillo_LIBDIR}")
endif()

if(EXISTS ${armadillo_LIBDIR}/libarmadillo.so)

  set(armadillo_LIBRARIES "-Wl,-rpath,${armadillo_LIBDIR} -L${armadillo_LIBDIR} -larmadillo")
  if (BUILD_STATIC_LIBS)
    set(armadillo_LIBRARIES_STATIC ${armadillo_LIBRARIES})
  endif ()
  if (BUILD_SHARED_LIBS)
    set(armadillo_LIBRARIES_SHARED ${armadillo_LIBRARIES})
  endif ()
  set(armadillo_LIBRARIES)
  
  message(STATUS "set armadillo static libs : ${armadillo_LIBRARIES_STATIC}")
  message(STATUS "set armadillo shared libs : ${armadillo_LIBRARIES_SHARED}")
  
else ()    
  message(FATAL_ERROR "Couldn't find armadillo library. Check if -DARMADILLO_DIR was set correctly.")    
endif()

  
