// GmshFEM - Copyright (C) 2019-2022, A. Royer, E. Béchet, C. Geuzaine, Université de Liège
//
// See the LICENSE.txt file for license information. Please report all
// issues on https://gitlab.onelab.info/gmsh/fem/issues

#ifndef H_GMSHFEM_TENSOR4NOVENRARYOPERATORS
#define H_GMSHFEM_TENSOR4NOVENRARYOPERATORS

#include "Function.h"
#include "GeneralEvaluableObject.h"

// *************************************
// Novenary operators for tensor<4> function
// *************************************

namespace gmshfem::function
{


  template< class T_Scalar, unsigned int T_Rank >
  Function< T_Scalar, Degree::Degree4 > tensor(const Function< T_Scalar, Degree::Degree2 > &a,
                                               const Function< T_Scalar, Degree::Degree2 > &b,
                                               const Function< T_Scalar, Degree::Degree2 > &c,
                                               const Function< T_Scalar, Degree::Degree2 > &d,
                                               const Function< T_Scalar, Degree::Degree2 > &e,
                                               const Function< T_Scalar, Degree::Degree2 > &f,
                                               const Function< T_Scalar, Degree::Degree2 > &g,
                                               const Function< T_Scalar, Degree::Degree2 > &h,
                                               const Function< T_Scalar, Degree::Degree2 > &i);


} // namespace gmshfem::function

#endif // H_GMSHFEM_TENSOR4NOVENRARYOPERATORS
