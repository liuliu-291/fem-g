set(TUTOS
  coupling
  dirichlet
  eigensolve
  electromagnet/magnetostatics
  functions
  globalQuantity/microStrip
  multiRHS
  neumann
  nonLinear
  periodic
  postPro
  shell/polarShell
  waveNewmark
)

foreach(TUTO ${TUTOS})
  add_test(NAME "tuto:${TUTO}:mkdir" COMMAND "${CMAKE_COMMAND}" -E make_directory build WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/${TUTO})
  add_test(NAME "tuto:${TUTO}:cmake" COMMAND "${CMAKE_COMMAND}" -DCMAKE_PREFIX_PATH=${CMAKE_INSTALL_PREFIX} .. WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/${TUTO}/build)
  add_test(NAME "tuto:${TUTO}:build" COMMAND "${CMAKE_MAKE_PROGRAM}" WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/${TUTO}/build)
  add_test(NAME "tuto:${TUTO}:run" COMMAND ./tuto -maxThreads 1 WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/${TUTO}/build)
endforeach()
