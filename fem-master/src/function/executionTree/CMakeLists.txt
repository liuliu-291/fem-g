## GmshFEM - Copyright (C) 2019-2022, A. Royer, E. Béchet, C. Geuzaine, Université de Liège
##
## See the LICENSE.txt file for license information. Please report all
## issues on https://gitlab.onelab.info/gmsh/fem/issues

set(SRC
  ExecutionTree.cpp
  ExecutionTreeIterator.cpp
  fieldNodeFunctions.cpp
  FunctionAllocator.cpp
  OperationsInterface.cpp
)

file(GLOB HDR RELATIVE ${CMAKE_CURRENT_SOURCE_DIR} *.h)
append_gmshfem_src(function/executionTree "${SRC};${HDR}")
