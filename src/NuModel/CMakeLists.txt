add_library(DUNE-3DST_NuModel
  analysis_3dst
)

include_directories(${INCDIRS})
include_directories(${CMAKE_CURRENT_SOURCE_DIR})

# Requried for fhiclcpp
include_directories(${CMAKE_INSTALL_PREFIX}/include)

install(TARGETS DUNE-3DST_NuModel DESTINATION lib)

target_link_libraries(DUNE-3DST_NuModel ${ROOT_LIBS})

if(NOT "${CMAKE_LINK_FLAGS}" STREQUAL "")
  set_target_properties(DUNE-3DST_NuModel PROPERTIES LINK_FLAGS ${CMAKE_LINK_FLAGS})
endif()

if(NOT "${CMAKE_DEPENDLIB_FLAGS}" STREQUAL "")
  target_link_libraries(DUNE-3DST_NuModel ${CMAKE_DEPENDLIB_FLAGS})
endif()

SET(MODULETargets ${MODULETargets};DUNE-3DST_NuModel PARENT_SCOPE)
