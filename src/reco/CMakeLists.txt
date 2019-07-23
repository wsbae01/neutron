add_library(DUNE-3DST_reco
  CBond3D.hxx
  CBond3D.cxx
  CCluster3D.cxx
  CCluster3D.hxx
  CCluster3DHits.hxx
  CCreate3DHits.hxx
  CHit2D.cxx
  CHit2D.hxx
  CHit3D.cxx
  CHit3D.hxx
  CInputRead.hxx
  CInputRead_CT.hxx
  CSharedCharge.hxx
  LinearFit.hxx
  TTmplDensityCluster.hxx
  TTmplMinimalSpanningTree.hxx
)

include_directories(${INCDIRS})
include_directories(${CMAKE_CURRENT_SOURCE_DIR})

# Requried for fhiclcpp
include_directories(${CMAKE_INSTALL_PREFIX}/include)

install(TARGETS DUNE-3DST_reco DESTINATION lib)

target_link_libraries(DUNE-3DST_reco ${ROOT_LIBS})

if(NOT "${CMAKE_LINK_FLAGS}" STREQUAL "")
  set_target_properties(DUNE-3DST_reco PROPERTIES LINK_FLAGS ${CMAKE_LINK_FLAGS})
endif()

if(NOT "${CMAKE_DEPENDLIB_FLAGS}" STREQUAL "")
  target_link_libraries(DUNE-3DST_reco ${CMAKE_DEPENDLIB_FLAGS})
endif()

SET(MODULETargets ${MODULETargets};DUNE-3DST_reco PARENT_SCOPE)
