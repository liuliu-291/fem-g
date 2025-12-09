## GmshFEM - Copyright (C) 2019-2022, A. Royer, E. Béchet, C. Geuzaine, Université de Liège
##
## See the LICENSE.txt file for license information. Please report all
## issues on https://gitlab.onelab.info/gmsh/fem/issues

set(SRC
  Helmholtz2D.cpp
  Helmholtz3D.cpp
  Maxwell2D.cpp
  Navier2D.cpp
  Maxwell3D.cpp
)

file(GLOB HDR RELATIVE ${CMAKE_CURRENT_SOURCE_DIR} *.h)
append_gmshfem_src(analytics "${SRC};${HDR}")
