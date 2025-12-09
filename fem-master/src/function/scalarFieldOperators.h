// GmshFEM - Copyright (C) 2019-2022, A. Royer, E. Béchet, C. Geuzaine, Université de Liège
//
// See the LICENSE.txt file for license information. Please report all
// issues on https://gitlab.onelab.info/gmsh/fem/issues

#ifndef H_GMSHFEM_SCALARFIELDOPERATORS
#define H_GMSHFEM_SCALARFIELDOPERATORS

#include "Function.h"

// ***********************************
// Field operators for scalar function
// ***********************************

namespace gmshfem::function
{


  // ############################
  // Field
  // ############################

  template< class T_Scalar >
  Function< T_Scalar, Degree::Degree0 > div(const field::Field< T_Scalar, field::Form::Form2 > &field);

  template< class T_Scalar >
  Function< T_Scalar, Degree::Degree0 > d(const field::Field< T_Scalar, field::Form::Form2 > &field);

  template< class T_Scalar >
  Function< scalar::ComplexPrecision< T_Scalar >, Degree::Degree0 > eigenfunction(const field::Field< T_Scalar, field::Form::Form0 > &field, const unsigned int eigenTag);

  template< class T_Scalar >
  Function< scalar::ComplexPrecision< T_Scalar >, Degree::Degree0 > eigenfunction(const field::Field< T_Scalar, field::Form::Form3 > &field, const unsigned int eigenTag);


} // namespace gmshfem::function

#endif // H_GMSHFEM_SCALARFIELDOPERATORS
