// GmshFEM - Copyright (C) 2019-2022, A. Royer, E. Béchet, C. Geuzaine, Université de Liège
//
// See the LICENSE.txt file for license information. Please report all
// issues on https://gitlab.onelab.info/gmsh/fem/issues

#ifndef H_GMSHFEM_VECTORSCALARTYPEOPERATORS
#define H_GMSHFEM_VECTORSCALARTYPEOPERATORS

#include "Function.h"
#include "GeneralEvaluableObject.h"

// ***********************************
// Scalar type operators for vector function
// ***********************************

namespace gmshfem::function
{


  template< class T_Scalar >
  Function< scalar::ComplexPrecision< T_Scalar >, Degree::Degree1 > complex(const GeneralEvaluableObject< T_Scalar, Degree::Degree1 > &a);

  template< class T_Scalar >
  Function< scalar::Precision< T_Scalar >, Degree::Degree1 > realPart(const GeneralEvaluableObject< T_Scalar, Degree::Degree1 > &a);

  template< class T_Scalar >
  Function< scalar::Precision< T_Scalar >, Degree::Degree1 > imagPart(const GeneralEvaluableObject< T_Scalar, Degree::Degree1 > &a);


} // namespace gmshfem::function

#endif // H_GMSHFEM_VECTORSCALARTYPEOPERATORS
