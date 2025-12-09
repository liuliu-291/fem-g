// GmshFEM - Copyright (C) 2019-2022, A. Royer, E. Béchet, C. Geuzaine, Université de Liège
//
// See the LICENSE.txt file for license information. Please report all
// issues on https://gitlab.onelab.info/gmsh/fem/issues

#ifndef H_GMSHFEM_VECTORBINARYOPERATORS
#define H_GMSHFEM_VECTORBINARYOPERATORS

#include "Function.h"
#include "GeneralEvaluableObject.h"

// ***********************************
// Binary operators for vector function (not +, -, *, / and %)
// ***********************************

namespace gmshfem::function
{


  template< class T_Scalar >
  Function< T_Scalar, Degree::Degree1 > hadamard(const GeneralEvaluableObject< T_Scalar, Degree::Degree1 > &a, const GeneralEvaluableObject< T_Scalar, Degree::Degree1 > &b);


} // namespace gmshfem::function

#endif // H_GMSHFEM_VECTORBINARYOPERATORS
