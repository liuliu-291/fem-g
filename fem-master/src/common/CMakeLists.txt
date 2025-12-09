## GmshFEM - Copyright (C) 2019-2022, A. Royer, E. Béchet, C. Geuzaine, Université de Liège
##
## See the LICENSE.txt file for license information. Please report all
## issues on https://gitlab.onelab.info/gmsh/fem/issues

set(SRC
  ArgsManager.cpp
  Color.cpp
  ColorMap.cpp
  CommandLine.cpp
  Exception.cpp
  GmshFem.cpp
  Memory.cpp
  Message.cpp
  OmpInterface.cpp
  Options.cpp
  Timer.cpp
)

file(GLOB HDR RELATIVE ${CMAKE_CURRENT_SOURCE_DIR} *.h)
append_gmshfem_src(common "${SRC};${HDR}")
