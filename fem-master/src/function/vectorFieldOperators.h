// GmshFEM - Copyright (C) 2019-2022, A. Royer, E. Béchet, C. Geuzaine, Université de Liège
//
// See the LICENSE.txt file for license information. Please report all
// issues on https://gitlab.onelab.info/gmsh/fem/issues

#ifndef H_GMSHFEM_VECTORFIELDOPERATORS
#define H_GMSHFEM_VECTORFIELDOPERATORS

#include "Function.h"

// ***********************************
// Field operators for vector function
// ***********************************

namespace gmshfem::function
{


  // ############################
  // Field
  // ############################

  template< class T_Scalar >
  Function< T_Scalar, Degree::Degree1 > grad(const field::Field< T_Scalar, field::Form::Form0 > &field);

  template< class T_Scalar >
  Function< T_Scalar, Degree::Degree1 > d(const field::Field< T_Scalar, field::Form::Form0 > &field);

  template< class T_Scalar >
  Function< T_Scalar, Degree::Degree1 > curl(const field::Field< T_Scalar, field::Form::Form1 > &field);

  template< class T_Scalar >
  Function< T_Scalar, Degree::Degree1 > d(const field::Field< T_Scalar, field::Form::Form1 > &field);

  template< class T_Scalar >
  Function< scalar::ComplexPrecision< T_Scalar >, Degree::Degree1 > eigenfunction(const field::Field< T_Scalar, field::Form::Form1 > &field, const unsigned int eigenTag);

  template< class T_Scalar >
  Function< scalar::ComplexPrecision< T_Scalar >, Degree::Degree1 > eigenfunction(const field::Field< T_Scalar, field::Form::Form2 > &field, const unsigned int eigenTag);

  // ############################
  // CompoundField
  // ############################

  template< class T_Scalar, unsigned int T_NumFields >
  Function< T_Scalar, Degree::Degree1 > div(const field::CompoundField< T_Scalar, field::Form::Form2, T_NumFields > &field);

  template< class T_Scalar, unsigned int T_NumFields >
  Function< T_Scalar, Degree::Degree1 > d(const field::CompoundField< T_Scalar, field::Form::Form2, T_NumFields > &field);

  template< class T_Scalar, unsigned int T_NumFields >
  Function< scalar::ComplexPrecision< T_Scalar >, Degree::Degree1 > eigenfunction(const field::CompoundField< T_Scalar, field::Form::Form0, T_NumFields > &field, const unsigned int eigenTag);

  template< class T_Scalar, unsigned int T_NumFields >
  Function< scalar::ComplexPrecision< T_Scalar >, Degree::Degree1 > eigenfunction(const field::CompoundField< T_Scalar, field::Form::Form3, T_NumFields > &field, const unsigned int eigenTag);


} // namespace gmshfem::function

#endif // H_GMSHFEM_VECTORFIELDOPERATORS
