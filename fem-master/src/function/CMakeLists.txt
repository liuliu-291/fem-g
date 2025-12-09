## GmshFEM - Copyright (C) 2019-2022, A. Royer, E. Béchet, C. Geuzaine, Université de Liège
##
## See the LICENSE.txt file for license information. Please report all
## issues on https://gitlab.onelab.info/gmsh/fem/issues

set(SRC
  GeneralEvaluableObject.cpp

  PiecewiseFunction.cpp	

  ScalarFunction.cpp
  VectorFunction.cpp
  TensorFunction.cpp
  Tensor4Function.cpp

  # scalarNovenaryOperators.cpp
  # vectorNovenaryOperators.cpp
  tensorNovenaryOperators.cpp
  tensor4NovenaryOperators.cpp

  # scalarTernaryOperators.cpp
  vectorTernaryOperators.cpp
  tensorTernaryOperators.cpp

  scalarScalarTypeOperators.cpp
  vectorScalarTypeOperators.cpp
  tensorScalarTypeOperators.cpp

  scalarBinaryOperators.cpp
  vectorBinaryOperators.cpp
  tensorBinaryOperators.cpp
  tensor4BinaryOperators.cpp

  scalarUnaryOperators.cpp
  vectorUnaryOperators.cpp
  tensorUnaryOperators.cpp
  tensor4UnaryOperators.cpp

  scalarNullaryOperators.cpp
  vectorNullaryOperators.cpp
  tensorNullaryOperators.cpp

  scalarFieldOperators.cpp
  vectorFieldOperators.cpp
  tensorFieldOperators.cpp

)

file(GLOB HDR RELATIVE ${CMAKE_CURRENT_SOURCE_DIR} *.h)
append_gmshfem_src(function "${SRC};${HDR}")
