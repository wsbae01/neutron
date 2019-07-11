add_library(DUNE-3DST_BeamMonitoring
  beam.cxx
  beamSpectrum.cxx
  fluxChange.cxx
  fluxRate.cxx
)

include_directories(${INCDIRS})
include_directories(${CMAKE_CURRENT_SOURCE_DIR})

# Requried for fhiclcpp
include_directories(${CMAKE_INSTALL_PREFIX}/include)

install(TARGETS DUNE-3DST_BeamMonitoring DESTINATION lib)

target_link_libraries(DUNE-3DST_BeamMonitoring ${ROOT_LIBS})

if(NOT "${CMAKE_LINK_FLAGS}" STREQUAL "")
  set_target_properties(DUNE-3DST_BeamMonitoring PROPERTIES LINK_FLAGS ${CMAKE_LINK_FLAGS})
endif()

if(NOT "${CMAKE_DEPENDLIB_FLAGS}" STREQUAL "")
  target_link_libraries(DUNE-3DST_BeamMonitoring ${CMAKE_DEPENDLIB_FLAGS})
endif()

SET(MODULETargets ${MODULETargets};DUNE-3DST_BeamMonitoring PARENT_SCOPE)
