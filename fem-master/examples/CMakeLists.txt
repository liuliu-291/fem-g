set(EXAMPLES
  helmholtz2d/diskScattering
  helmholtz2d/waveguide/longitudinal
  helmholtz2d/waveguide/transverse
  helmholtzFlow/flowPastCylinder
  helmholtzFlow/freefield
  helmholtzFlow/linedWaveguide
#  helmholtzFlow/Nacelle
  helmholtzFlow/waveguide/nonUniform
  helmholtzFlow/waveguide/uniform
  maxwell2d
  maxwell3d
  navier2d/diskScattering
  navier2d/marmousi2
)

foreach(EXAMPLE ${EXAMPLES})
  add_test(NAME "example:${EXAMPLE}:mkdir" COMMAND "${CMAKE_COMMAND}" -E make_directory build WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/${EXAMPLE})
  add_test(NAME "example:${EXAMPLE}:cmake" COMMAND "${CMAKE_COMMAND}" -DCMAKE_PREFIX_PATH=${CMAKE_INSTALL_PREFIX} .. WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/${EXAMPLE}/build)
  add_test(NAME "example:${EXAMPLE}:build" COMMAND "${CMAKE_MAKE_PROGRAM}" WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/${EXAMPLE}/build)
endforeach()
