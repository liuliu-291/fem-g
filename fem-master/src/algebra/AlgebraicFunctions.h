// GmshFEM - Copyright (C) 2019-2022, A. Royer, E. Béchet, C. Geuzaine, Université de Liège
//
// See the LICENSE.txt file for license information. Please report all
// issues on https://gitlab.onelab.info/gmsh/fem/issues

#ifndef H_GMSHFEM_ALGEBRAICFUNCTIONS
#define H_GMSHFEM_ALGEBRAICFUNCTIONS

#include "Matrix.h"
#include "Vector.h"
#include "scalar.h"

namespace gmshfem::algebra
{


  template< class T_Scalar >
  void eigenvalues(Vector< scalar::ComplexPrecision< T_Scalar > > &eigenvalues, std::vector< Vector< scalar::ComplexPrecision< T_Scalar > > > &eigenvectors, Matrix< T_Scalar > &matrix, const bool symmetric = false, const unsigned long long numberOfEigenvalues = 0, const bool saveeigenvectors = true);

  template< class T_Scalar >
  void eigenvalues(Vector< scalar::ComplexPrecision< T_Scalar > > &eigenvalues, Matrix< T_Scalar > &matrix, const bool symmetric = false, const unsigned long long numberOfEigenvalues = 0);

  template< class T_Scalar >
  void eigenvalues(Vector< scalar::ComplexPrecision< T_Scalar > > &eigenvalues, Matrix< T_Scalar > &left, Matrix< T_Scalar > &right, const bool symmetric = false, const unsigned long long numberOfEigenvalues = 0);

  template< class T_Scalar >
  void singularValues(Vector< scalar::Precision< T_Scalar > > &singularValues, Matrix< T_Scalar > &matrix, const bool smallest = false, const unsigned long long numberOfSingularValues = 0);

  template< class T_Scalar >
  T_Scalar bilinear(const Vector< T_Scalar > &x, const Matrix< T_Scalar > &A, const Vector< T_Scalar > &y);

  template< class T_Scalar >
  void linear(Vector< T_Scalar > &b, const Matrix< T_Scalar > &A, const Vector< T_Scalar > &x);

  template< class T_Scalar >
  scalar::Precision< T_Scalar > residual(const Vector< T_Scalar > &b, const Matrix< T_Scalar > &A, const Vector< T_Scalar > &x);


} // namespace gmshfem::algebra

#endif // H_GMSHFEM_ALGEBRAICFUNCTIONS
