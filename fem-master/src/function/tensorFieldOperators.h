// GmshFEM - Copyright (C) 2019-2022, A. Royer, E. Béchet, C. Geuzaine, Université de Liège
//
// See the LICENSE.txt file for license information. Please report all
// issues on https://gitlab.onelab.info/gmsh/fem/issues

#ifndef H_GMSHFEM_TENSORFIELDOPERATORS
#define H_GMSHFEM_TENSORFIELDOPERATORS

#include "Function.h"

// ***********************************
// Field operators for tensor function
// ***********************************

namespace gmshfem::function
{


  // ############################
  // CompoundField
  // ############################

  template< class T_Scalar, unsigned int T_NumFields >
  Function< T_Scalar, Degree::Degree2 > grad(const field::CompoundField< T_Scalar, field::Form::Form0, T_NumFields > &field);

  template< class T_Scalar, unsigned int T_NumFields >
  Function< T_Scalar, Degree::Degree2 > d(const field::CompoundField< T_Scalar, field::Form::Form0, T_NumFields > &field);

  template< class T_Scalar, unsigned int T_NumFields >
  Function< T_Scalar, Degree::Degree2 > curl(const field::CompoundField< T_Scalar, field::Form::Form1, T_NumFields > &field);

  template< class T_Scalar, unsigned int T_NumFields >
  Function< T_Scalar, Degree::Degree2 > d(const field::CompoundField< T_Scalar, field::Form::Form1, T_NumFields > &field);

  template< class T_Scalar, unsigned int T_NumFields >
  Function< scalar::ComplexPrecision< T_Scalar >, Degree::Degree2 > eigenfunction(const field::CompoundField< T_Scalar, field::Form::Form1, T_NumFields > &field, const unsigned int eigenTag);

  template< class T_Scalar, unsigned int T_NumFields >
  Function< scalar::ComplexPrecision< T_Scalar >, Degree::Degree2 > eigenfunction(const field::CompoundField< T_Scalar, field::Form::Form2, T_NumFields > &field, const unsigned int eigenTag);


} // namespace gmshfem::function

#endif // H_GMSHFEM_TENSORFIELDOPERATORS
