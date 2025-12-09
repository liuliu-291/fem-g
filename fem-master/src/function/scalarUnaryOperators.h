// GmshFEM - Copyright (C) 2019-2022, A. Royer, E. Béchet, C. Geuzaine, Université de Liège
//
// See the LICENSE.txt file for license information. Please report all
// issues on https://gitlab.onelab.info/gmsh/fem/issues

#ifndef H_GMSHFEM_SCALARUNARYOPERATORS
#define H_GMSHFEM_SCALARUNARYOPERATORS

#include "Function.h"
#include "GeneralEvaluableObject.h"

// ***********************************
// Unary operators for scalar function
// ***********************************

namespace gmshfem::function
{


  template< class T_Scalar >
  Function< T_Scalar, Degree::Degree0 > operator-(const GeneralEvaluableObject< T_Scalar, Degree::Degree0 > &a);

  template< class T_Scalar >
  Function< T_Scalar, Degree::Degree0 > operator+(const GeneralEvaluableObject< T_Scalar, Degree::Degree0 > &a);

  template< class T_Scalar >
  Function< T_Scalar, Degree::Degree0 > abs(const GeneralEvaluableObject< T_Scalar, Degree::Degree0 > &a);

  template< class T_Scalar >
  Function< T_Scalar, Degree::Degree0 > acos(const GeneralEvaluableObject< T_Scalar, Degree::Degree0 > &a);

  template< class T_Scalar >
  Function< T_Scalar, Degree::Degree0 > acosh(const GeneralEvaluableObject< T_Scalar, Degree::Degree0 > &a);

  template< class T_Scalar >
  Function< T_Scalar, Degree::Degree0 > angularComp(const GeneralEvaluableObject< T_Scalar, Degree::Degree1 > &a);

  template< class T_Scalar >
  Function< T_Scalar, Degree::Degree0 > asin(const GeneralEvaluableObject< T_Scalar, Degree::Degree0 > &a);

  template< class T_Scalar >
  Function< T_Scalar, Degree::Degree0 > asinh(const GeneralEvaluableObject< T_Scalar, Degree::Degree0 > &a);

  template< class T_Scalar >
  Function< T_Scalar, Degree::Degree0 > atan(const GeneralEvaluableObject< T_Scalar, Degree::Degree0 > &a);

  template< class T_Scalar >
  Function< T_Scalar, Degree::Degree0 > atanh(const GeneralEvaluableObject< T_Scalar, Degree::Degree0 > &a);

  template< class T_Scalar >
  Function< T_Scalar, Degree::Degree0 > cbrt(const GeneralEvaluableObject< T_Scalar, Degree::Degree0 > &a);

  template< class T_Scalar >
  Function< T_Scalar, Degree::Degree0 > changeOfCoordinates(const GeneralEvaluableObject< T_Scalar, Degree::Degree0 > &a, const GeneralEvaluableObject< scalar::Precision< T_Scalar >, Degree::Degree0 > &x, const GeneralEvaluableObject< scalar::Precision< T_Scalar >, Degree::Degree0 > &y, const GeneralEvaluableObject< scalar::Precision< T_Scalar >, Degree::Degree0 > &z);

  template< class T_Scalar >
  Function< T_Scalar, Degree::Degree0 > commutator(const GeneralEvaluableObject< T_Scalar, Degree::Degree0 > &a, bool *state);

  template< class T_Scalar >
  Function< T_Scalar, Degree::Degree0 > conj(const GeneralEvaluableObject< T_Scalar, Degree::Degree0 > &a);

  template< class T_Scalar >
  Function< T_Scalar, Degree::Degree0 > cos(const GeneralEvaluableObject< T_Scalar, Degree::Degree0 > &a);

  template< class T_Scalar >
  Function< T_Scalar, Degree::Degree0 > cosh(const GeneralEvaluableObject< T_Scalar, Degree::Degree0 > &a);

  template< class T_Scalar >
  Function< T_Scalar, Degree::Degree0 > cylBesselJ(const scalar::Precision< T_Scalar > v, const GeneralEvaluableObject< T_Scalar, Degree::Degree0 > &a);

  template< class T_Scalar >
  Function< T_Scalar, Degree::Degree0 > cylNeumann(const scalar::Precision< T_Scalar > v, const GeneralEvaluableObject< T_Scalar, Degree::Degree0 > &a);

  template< class T_Scalar >
  Function< T_Scalar, Degree::Degree0 > exp(const GeneralEvaluableObject< T_Scalar, Degree::Degree0 > &a);

  template< class T_Scalar >
  Function< T_Scalar, Degree::Degree0 > heaviside(const GeneralEvaluableObject< T_Scalar, Degree::Degree0 > &a);

  template< class T_Scalar >
  Function< T_Scalar, Degree::Degree0 > imag(const GeneralEvaluableObject< T_Scalar, Degree::Degree0 > &a);

  template< class T_Scalar >
  Function< T_Scalar, Degree::Degree0 > log(const GeneralEvaluableObject< T_Scalar, Degree::Degree0 > &a);

  template< class T_Scalar >
  Function< T_Scalar, Degree::Degree0 > ln(const GeneralEvaluableObject< T_Scalar, Degree::Degree0 > &a);

  template< class T_Scalar >
  Function< T_Scalar, Degree::Degree0 > norm(const GeneralEvaluableObject< T_Scalar, Degree::Degree0 > &a);

  template< class T_Scalar >
  Function< T_Scalar, Degree::Degree0 > norm(const GeneralEvaluableObject< T_Scalar, Degree::Degree1 > &a, const int order = 2);

  template< class T_Scalar >
  Function< T_Scalar, Degree::Degree0 > norm(const GeneralEvaluableObject< T_Scalar, Degree::Degree2 > &a, const int order = 2);

  template< class T_Scalar >
  Function< T_Scalar, Degree::Degree0 > pow(const GeneralEvaluableObject< T_Scalar, Degree::Degree0 > &a, const unsigned int n);

  template< class T_Scalar >
  Function< T_Scalar, Degree::Degree0 > r2dComp(const GeneralEvaluableObject< T_Scalar, Degree::Degree1 > &a);

  template< class T_Scalar >
  Function< T_Scalar, Degree::Degree0 > real(const GeneralEvaluableObject< T_Scalar, Degree::Degree0 > &a);

  template< class T_Scalar >
  Function< T_Scalar, Degree::Degree0 > sin(const GeneralEvaluableObject< T_Scalar, Degree::Degree0 > &a);

  template< class T_Scalar >
  Function< T_Scalar, Degree::Degree0 > sinh(const GeneralEvaluableObject< T_Scalar, Degree::Degree0 > &a);

  template< class T_Scalar >
  Function< T_Scalar, Degree::Degree0 > sqrt(const GeneralEvaluableObject< T_Scalar, Degree::Degree0 > &a);

  template< class T_Scalar >
  Function< T_Scalar, Degree::Degree0 > tan(const GeneralEvaluableObject< T_Scalar, Degree::Degree0 > &a);

  template< class T_Scalar >
  Function< T_Scalar, Degree::Degree0 > tanh(const GeneralEvaluableObject< T_Scalar, Degree::Degree0 > &a);

  template< class T_Scalar >
  Function< T_Scalar, Degree::Degree0 > tr(const GeneralEvaluableObject< T_Scalar, Degree::Degree2 > &a);

  template< class T_Scalar >
  Function< T_Scalar, Degree::Degree0 > xComp(const GeneralEvaluableObject< T_Scalar, Degree::Degree1 > &a);

  template< class T_Scalar >
  Function< T_Scalar, Degree::Degree0 > xxComp(const GeneralEvaluableObject< T_Scalar, Degree::Degree2 > &a);

  template< class T_Scalar >
  Function< T_Scalar, Degree::Degree0 > xyComp(const GeneralEvaluableObject< T_Scalar, Degree::Degree2 > &a);

  template< class T_Scalar >
  Function< T_Scalar, Degree::Degree0 > xzComp(const GeneralEvaluableObject< T_Scalar, Degree::Degree2 > &a);

  template< class T_Scalar >
  Function< T_Scalar, Degree::Degree0 > yComp(const GeneralEvaluableObject< T_Scalar, Degree::Degree1 > &a);

  template< class T_Scalar >
  Function< T_Scalar, Degree::Degree0 > yxComp(const GeneralEvaluableObject< T_Scalar, Degree::Degree2 > &a);

  template< class T_Scalar >
  Function< T_Scalar, Degree::Degree0 > yyComp(const GeneralEvaluableObject< T_Scalar, Degree::Degree2 > &a);

  template< class T_Scalar >
  Function< T_Scalar, Degree::Degree0 > yzComp(const GeneralEvaluableObject< T_Scalar, Degree::Degree2 > &a);

  template< class T_Scalar >
  Function< T_Scalar, Degree::Degree0 > zComp(const GeneralEvaluableObject< T_Scalar, Degree::Degree1 > &a);

  template< class T_Scalar >
  Function< T_Scalar, Degree::Degree0 > zxComp(const GeneralEvaluableObject< T_Scalar, Degree::Degree2 > &a);

  template< class T_Scalar >
  Function< T_Scalar, Degree::Degree0 > zyComp(const GeneralEvaluableObject< T_Scalar, Degree::Degree2 > &a);

  template< class T_Scalar >
  Function< T_Scalar, Degree::Degree0 > zzComp(const GeneralEvaluableObject< T_Scalar, Degree::Degree2 > &a);


} // namespace gmshfem::function

#endif // H_GMSHFEM_SCALARUNARYOPERATORS
