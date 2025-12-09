// GmshFEM - Copyright (C) 2019-2022, A. Royer, E. Béchet, C. Geuzaine, Université de Liège
//
// See the LICENSE.txt file for license information. Please report all
// issues on https://gitlab.onelab.info/gmsh/fem/issues

#ifndef H_GMSHFEM_SCALARBINARYOPERATORS
#define H_GMSHFEM_SCALARBINARYOPERATORS

#include "Function.h"
#include "GeneralEvaluableObject.h"

// ***********************************
// Binary operators for scalar function (not +, -, *, / and %)
// ***********************************

namespace gmshfem::function
{


  template< class T_PScalar >
  Function< T_PScalar, Degree::Degree0 > atan2(const GeneralEvaluableObject< T_PScalar, Degree::Degree0 > &a, const GeneralEvaluableObject< T_PScalar, Degree::Degree0 > &b);

  template< class T_Scalar >
  Function< T_Scalar, Degree::Degree0 > hadamard(const GeneralEvaluableObject< T_Scalar, Degree::Degree0 > &a, const GeneralEvaluableObject< T_Scalar, Degree::Degree0 > &b);


} // namespace gmshfem::function

#endif // H_GMSHFEM_SCALARBINARYOPERATORS
