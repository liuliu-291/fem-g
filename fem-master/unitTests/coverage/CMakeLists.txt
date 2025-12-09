## GmshFEM - Copyright (C) 2019-2021 Université de Liège
##
## See the LICENSE.txt file for license information. Please report all
## issues on https://gitlab.onelab.info/gmsh/fem/issues

set(SRC
  convergences.cpp
  fields.cpp
  formulations.cpp
  functions.cpp
  geo.cpp
  main.cpp
  posts.cpp
)

file(GLOB HDR RELATIVE ${CMAKE_CURRENT_SOURCE_DIR} *.h)
append_gmshfem_test("${SRC};${HDR}")
