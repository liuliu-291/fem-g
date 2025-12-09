// GmshFEM - Copyright (C) 2019-2022, A. Royer, E. Béchet, C. Geuzaine, Université de Liège
//
// See the LICENSE.txt file for license information. Please report all
// issues on https://gitlab.onelab.info/gmsh/fem/issues

#ifndef H_GMSHFEM_TENSORNOVENRARYOPERATORS
#define H_GMSHFEM_TENSORNOVENRARYOPERATORS

#include "Function.h"
#include "GeneralEvaluableObject.h"

// *************************************
// Novenary operators for tensor function
// *************************************

namespace gmshfem::function
{


  template< class T_Scalar >
  Function< T_Scalar, Degree::Degree2 > tensor(const GeneralEvaluableObject< T_Scalar, Degree::Degree0 > &a,
                                               const GeneralEvaluableObject< T_Scalar, Degree::Degree0 > &b,
                                               const GeneralEvaluableObject< T_Scalar, Degree::Degree0 > &c,
                                               const GeneralEvaluableObject< T_Scalar, Degree::Degree0 > &d,
                                               const GeneralEvaluableObject< T_Scalar, Degree::Degree0 > &e,
                                               const GeneralEvaluableObject< T_Scalar, Degree::Degree0 > &f,
                                               const GeneralEvaluableObject< T_Scalar, Degree::Degree0 > &g,
                                               const GeneralEvaluableObject< T_Scalar, Degree::Degree0 > &h,
                                               const GeneralEvaluableObject< T_Scalar, Degree::Degree0 > &i);


} // namespace gmshfem::function

#endif // H_GMSHFEM_TENSORNOVENRARYOPERATORS
