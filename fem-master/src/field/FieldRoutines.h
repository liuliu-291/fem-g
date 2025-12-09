// GmshFEM - Copyright (C) 2019-2022, A. Royer, E. Béchet, C. Geuzaine, Université de Liège
//
// See the LICENSE.txt file for license information. Please report all
// issues on https://gitlab.onelab.info/gmsh/fem/issues

#ifndef H_GMSHFEM_FIELDROUTINES
#define H_GMSHFEM_FIELDROUTINES

#include "FieldInterface.h"
#include "Matrix.h"
#include "scalar.h"

#include <complex>

namespace gmshfem::system
{
  template< class T_Scalar >
  class Matrix;
}

namespace gmshfem::field
{


  // Compute: c = m * a + n * b
  template< class T_Scalar >
  Field< T_Scalar, Form::Form0 > multiplyAdd(const T_Scalar &m, const Field< T_Scalar, Form::Form0 > &a, const T_Scalar &n, const Field< T_Scalar, Form::Form0 > &b);
  // Compute: a = m * a + n * b
  template< class T_Scalar >
  void multiplyAddAssignement(const T_Scalar &m, Field< T_Scalar, Form::Form0 > &a, const T_Scalar &n, const Field< T_Scalar, Form::Form0 > &b);
  // Compute: b = conj(a)
  template< class T_Scalar >
  Field< T_Scalar, Form::Form0 > conjField(const Field< T_Scalar, Form::Form0 > &a);
  template< class T_ComplexScalar >
  Field< std::complex< T_ComplexScalar >, Form::Form0 > conjField(const Field< std::complex< T_ComplexScalar >, Form::Form0 > &a);
  // Compute: b = real(a)
  template< class T_Scalar >
  Field< T_Scalar, Form::Form0 > realField(const Field< T_Scalar, Form::Form0 > &a);
  template< class T_ComplexScalar >
  Field< std::complex< T_ComplexScalar >, Form::Form0 > realField(const Field< std::complex< T_ComplexScalar >, Form::Form0 > &a);
  // Compute: a * M * b
  template< class T_Scalar >
  T_Scalar bilinearProduct(const Field< T_Scalar, Form::Form0 > &a, const Field< T_Scalar, Form::Form0 > &b, const algebra::Matrix< T_Scalar > &M);


} // namespace gmshfem::field


#endif // H_GMSHFEM_FIELDROUTINES
