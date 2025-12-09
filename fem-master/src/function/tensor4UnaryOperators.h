// GmshFEM - Copyright (C) 2019-2022, A. Royer, E. Béchet, C. Geuzaine, Université de Liège
//
// See the LICENSE.txt file for license information. Please report all
// issues on https://gitlab.onelab.info/gmsh/fem/issues

#ifndef H_GMSHFEM_TENSOR4UNARYOPERATORS
#define H_GMSHFEM_TENSOR4UNARYOPERATORS

#include "Function.h"
#include "GeneralEvaluableObject.h"

// ***********************************
// Unary operators for tensor4 function
// ***********************************

namespace gmshfem::function
{


  template< class T_Scalar >
  Function< T_Scalar, Degree::Degree4 > operator-(const GeneralEvaluableObject< T_Scalar, Degree::Degree4 > &a);

  template< class T_Scalar >
  Function< T_Scalar, Degree::Degree4 > operator+(const GeneralEvaluableObject< T_Scalar, Degree::Degree4 > &a);

  template< class T_Scalar >
  Function< T_Scalar, Degree::Degree4 > changeOfCoordinates(const GeneralEvaluableObject< T_Scalar, Degree::Degree4 > &a, const GeneralEvaluableObject< scalar::Precision< T_Scalar >, Degree::Degree0 > &x, const GeneralEvaluableObject< scalar::Precision< T_Scalar >, Degree::Degree0 > &y, const GeneralEvaluableObject< scalar::Precision< T_Scalar >, Degree::Degree0 > &z);


} // namespace gmshfem::function

#endif // H_GMSHFEM_TENSOR4UNARYOPERATORS
