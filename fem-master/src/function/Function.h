// GmshFEM - Copyright (C) 2019-2022, A. Royer, E. Béchet, C. Geuzaine, Université de Liège
//
// See the LICENSE.txt file for license information. Please report all
// issues on https://gitlab.onelab.info/gmsh/fem/issues

#ifndef H_GMSHFEM_FUNCTION
#define H_GMSHFEM_FUNCTION

#include "ExecutionTreeIterator.h"
#include "FieldObject.h"
#include "MathObject.h"
#include "Timer.h"

namespace gmshfem::function
{


  // *****************************
  // Function< T_Scalar, T_Degree >
  // *****************************

  template< class T_Scalar, Degree T_Degree >
  class Function
  {
  };

  template< unsigned int T_Num >
  struct TensorRank {
    static constexpr unsigned int value = 3 * TensorRank< T_Num - 1 >::value;
  };

  template<>
  struct TensorRank< 0 > {
    static constexpr unsigned int value = 1;
  };

  template< class T_Scalar >
  using ScalarFunction = Function< T_Scalar, Degree::Degree0 >;
  template< class T_Scalar >
  using VectorFunction = Function< T_Scalar, Degree::Degree1 >;
  template< class T_Scalar, unsigned int T_Rank = 2 >
  using TensorFunction = Function< T_Scalar, Degree(TensorRank< T_Rank >::value) >;


} // namespace gmshfem::function


#include "PiecewiseFunction.h"
#include "ScalarFunction.h"
#include "Tensor4Function.h"
#include "TensorFunction.h"
#include "VectorFunction.h"
#include "scalarBinaryOperators.h"
#include "scalarFieldOperators.h"
#include "scalarNullaryOperators.h"
#include "scalarScalarTypeOperators.h"
#include "scalarUnaryOperators.h"
#include "tensor4BinaryOperators.h"
#include "tensor4NovenaryOperators.h"
#include "tensor4UnaryOperators.h"
#include "tensorBinaryOperators.h"
#include "tensorFieldOperators.h"
#include "tensorNovenaryOperators.h"
#include "tensorNullaryOperators.h"
#include "tensorScalarTypeOperators.h"
#include "tensorTernaryOperators.h"
#include "tensorUnaryOperators.h"
#include "vectorBinaryOperators.h"
#include "vectorFieldOperators.h"
#include "vectorNullaryOperators.h"
#include "vectorScalarTypeOperators.h"
#include "vectorTernaryOperators.h"
#include "vectorUnaryOperators.h"


#endif // H_GMSHFEM_FUNCTION
