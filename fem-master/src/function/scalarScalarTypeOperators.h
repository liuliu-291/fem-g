// GmshFEM - Copyright (C) 2019-2022, A. Royer, E. Béchet, C. Geuzaine, Université de Liège
//
// See the LICENSE.txt file for license information. Please report all
// issues on https://gitlab.onelab.info/gmsh/fem/issues

#ifndef H_GMSHFEM_SCALARSCALARTYPEOPERATORS
#define H_GMSHFEM_SCALARSCALARTYPEOPERATORS

#include "Function.h"
#include "GeneralEvaluableObject.h"

// ***********************************
// Scalar type operators for scalar function
// ***********************************

namespace gmshfem::function
{


  template< class T_Scalar >
  Function< scalar::ComplexPrecision< T_Scalar >, Degree::Degree0 > complex(const GeneralEvaluableObject< T_Scalar, Degree::Degree0 > &a);

  template< class T_Scalar >
  Function< scalar::Precision< T_Scalar >, Degree::Degree0 > realPart(const GeneralEvaluableObject< T_Scalar, Degree::Degree0 > &a);

  template< class T_Scalar >
  Function< scalar::Precision< T_Scalar >, Degree::Degree0 > imagPart(const GeneralEvaluableObject< T_Scalar, Degree::Degree0 > &a);


} // namespace gmshfem::function

#endif // H_GMSHFEM_SCALARSCALARTYPEOPERATORS
