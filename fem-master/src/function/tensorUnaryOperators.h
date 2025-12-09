// GmshFEM - Copyright (C) 2019-2022, A. Royer, E. Béchet, C. Geuzaine, Université de Liège
//
// See the LICENSE.txt file for license information. Please report all
// issues on https://gitlab.onelab.info/gmsh/fem/issues

#ifndef H_GMSHFEM_TENSORUNARYOPERATORS
#define H_GMSHFEM_TENSORUNARYOPERATORS

#include "Function.h"
#include "GeneralEvaluableObject.h"

// ***********************************
// Unary operators for tensor function
// ***********************************

namespace gmshfem::function
{


  template< class T_Scalar >
  Function< T_Scalar, Degree::Degree2 > operator-(const GeneralEvaluableObject< T_Scalar, Degree::Degree2 > &a);

  template< class T_Scalar >
  Function< T_Scalar, Degree::Degree2 > operator+(const GeneralEvaluableObject< T_Scalar, Degree::Degree2 > &a);

  template< class T_Scalar >
  Function< T_Scalar, Degree::Degree2 > changeOfCoordinates(const GeneralEvaluableObject< T_Scalar, Degree::Degree2 > &a, const GeneralEvaluableObject< scalar::Precision< T_Scalar >, Degree::Degree0 > &x, const GeneralEvaluableObject< scalar::Precision< T_Scalar >, Degree::Degree0 > &y, const GeneralEvaluableObject< scalar::Precision< T_Scalar >, Degree::Degree0 > &z);

  template< class T_Scalar >
  Function< T_Scalar, Degree::Degree2 > commutator(const GeneralEvaluableObject< T_Scalar, Degree::Degree2 > &a, bool *state);

  template< class T_Scalar >
  Function< T_Scalar, Degree::Degree2 > conj(const GeneralEvaluableObject< T_Scalar, Degree::Degree2 > &a);

  template< class T_Scalar >
  Function< T_Scalar, Degree::Degree2 > inv(const GeneralEvaluableObject< T_Scalar, Degree::Degree2 > &a);


} // namespace gmshfem::function

#endif // H_GMSHFEM_TENSORUNARYOPERATORS
