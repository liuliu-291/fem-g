## GmshFEM - Copyright (C) 2019-2022, A. Royer, E. Béchet, C. Geuzaine, Université de Liège
##
## See the LICENSE.txt file for license information. Please report all
## issues on https://gitlab.onelab.info/gmsh/fem/issues

set(SRC
  CompoundField.cpp
  Field.cpp
  FieldInterface.cpp
  FieldModule.cpp
  FieldRoutines.cpp

  FunctionSpace.cpp

  FunctionSpaceConstant.cpp
  FunctionSpaceCurlHCurl.cpp
  FunctionSpaceDivHDiv.cpp
  FunctionSpaceGradH1.cpp
  FunctionSpaceHierarchicalH1.cpp
  FunctionSpaceHierarchicalHCurl.cpp
  FunctionSpaceInterface.cpp
  FunctionSpaceLagrange.cpp
  FunctionSpaceParallel.cpp
  FunctionSpacePerpendicular.cpp
  GlobalQuantity.cpp
)

file(GLOB HDR RELATIVE ${CMAKE_CURRENT_SOURCE_DIR} *.h)
append_gmshfem_src(field "${SRC};${HDR}")
