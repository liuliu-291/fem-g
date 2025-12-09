// GmshFEM - Copyright (C) 2019-2022, A. Royer, E. Béchet, C. Geuzaine, Université de Liège
//
// See the LICENSE.txt file for license information. Please report all
// issues on https://gitlab.onelab.info/gmsh/fem/issues

#ifndef H_GMSHFEM_TENSOR4BINARYOPERATORS
#define H_GMSHFEM_TENSOR4BINARYOPERATORS

#include "Function.h"
#include "GeneralEvaluableObject.h"

// ***********************************
// Binary operators for tensor<4> function (not +, -, *, / and %)
// ***********************************

namespace gmshfem::function
{


  template< class T_Scalar >
  Function< T_Scalar, Degree::Degree4 > dyadic(const GeneralEvaluableObject< T_Scalar, Degree::Degree2 > &a, const GeneralEvaluableObject< T_Scalar, Degree::Degree2 > &b);
  
  template< class T_Scalar >
  Function< T_Scalar, Degree::Degree4 > hadamard(const GeneralEvaluableObject< T_Scalar, Degree::Degree4 > &a, const GeneralEvaluableObject< T_Scalar, Degree::Degree4 > &b);


} // namespace gmshfem::function

#endif // H_GMSHFEM_TENSOR4BINARYOPERATORS
