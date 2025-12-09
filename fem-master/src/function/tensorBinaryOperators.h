// GmshFEM - Copyright (C) 2019-2022, A. Royer, E. Béchet, C. Geuzaine, Université de Liège
//
// See the LICENSE.txt file for license information. Please report all
// issues on https://gitlab.onelab.info/gmsh/fem/issues

#ifndef H_GMSHFEM_TENSORBINARYOPERATORS
#define H_GMSHFEM_TENSORBINARYOPERATORS

#include "Function.h"
#include "GeneralEvaluableObject.h"

// ***********************************
// Binary operators for tensor function (not +, -, *, /)
// ***********************************

namespace gmshfem::function
{


  template< class T_Scalar >
  Function< T_Scalar, Degree::Degree2 > dyadic(const GeneralEvaluableObject< T_Scalar, Degree::Degree1 > &a, const GeneralEvaluableObject< T_Scalar, Degree::Degree1 > &b);

  template< class T_Scalar >
  Function< T_Scalar, Degree::Degree2 > hadamard(const GeneralEvaluableObject< T_Scalar, Degree::Degree2 > &a, const GeneralEvaluableObject< T_Scalar, Degree::Degree2 > &b);


} // namespace gmshfem::function

#endif // H_GMSHFEM_TENSORBINARYOPERATORS
