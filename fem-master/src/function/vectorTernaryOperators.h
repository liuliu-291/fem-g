// GmshFEM - Copyright (C) 2019-2022, A. Royer, E. Béchet, C. Geuzaine, Université de Liège
//
// See the LICENSE.txt file for license information. Please report all
// issues on https://gitlab.onelab.info/gmsh/fem/issues

#ifndef H_GMSHFEM_VECTORTERNARYOPERATORS
#define H_GMSHFEM_VECTORTERNARYOPERATORS

#include "Function.h"
#include "GeneralEvaluableObject.h"

// *************************************
// Ternary operators for vector function
// *************************************

namespace gmshfem::function
{


  template< class T_Scalar >
  Function< T_Scalar, Degree::Degree1 > vector(const GeneralEvaluableObject< T_Scalar, Degree::Degree0 > &a,
                                               const GeneralEvaluableObject< T_Scalar, Degree::Degree0 > &b,
                                               const GeneralEvaluableObject< T_Scalar, Degree::Degree0 > &c);


} // namespace gmshfem::function

#endif // H_GMSHFEM_VECTORTERNARYOPERATORS
