// GmshFEM - Copyright (C) 2019-2022, A. Royer, E. Béchet, C. Geuzaine, Université de Liège
//
// See the LICENSE.txt file for license information. Please report all
// issues on https://gitlab.onelab.info/gmsh/fem/issues

#include "Navier2D.h"

#include "Exception.h"
#include "KahanSum.h"
#include "OmpInterface.h"
#include "gmshfemDefines.h"
#include "instantiate.h"

#include <cmath>
#include <complex>

#ifdef HAVE_BOOST
#include <boost/math/special_functions/bessel.hpp>
#include <boost/math/special_functions/hankel.hpp>
#endif // HAVE_BOOST

namespace gmshfem::analytics::navier2D
{


  // Static usefull functions
#ifdef HAVE_BOOST
  template< class T_Scalar >
  static scalar::Precision< T_Scalar > drBesselJ(const scalar::Precision< T_Scalar > k, const scalar::Precision< T_Scalar > r, unsigned int n)
  {
    if(n == 0) {
      return -k * scalar::Precision< T_Scalar >(boost::math::cyl_bessel_j(1, k * r));
    }
    return scalar::Precision< T_Scalar >(0.5) * k * scalar::Precision< T_Scalar >(boost::math::cyl_bessel_j(n - 1, k * r) - boost::math::cyl_bessel_j(n + 1, k * r));
  }

  template< class T_Scalar >
  static T_Scalar drBesselH1(const scalar::Precision< T_Scalar > k, const scalar::Precision< T_Scalar > r, unsigned int n)
  {
    if(n == 0) {
      return -k * T_Scalar(boost::math::cyl_hankel_1(1, k * r));
    }
    return scalar::Precision< T_Scalar >(0.5) * k * T_Scalar(boost::math::cyl_hankel_1(n - 1, k * r) - boost::math::cyl_hankel_1(n + 1, k * r));
  }

  template< class T_Scalar >
  static T_Scalar drBesselH2(const scalar::Precision< T_Scalar > k, const scalar::Precision< T_Scalar > r, unsigned int n)
  {
    if(n == 0) {
      return -k * T_Scalar(boost::math::cyl_hankel_2(1, k * r));
    }
    return scalar::Precision< T_Scalar >(0.5) * k * T_Scalar(boost::math::cyl_hankel_2(n - 1, k * r) - boost::math::cyl_hankel_2(n + 1, k * r));
  }

  template< class T_Scalar >
  static T_Scalar dr2BesselH1(const scalar::Precision< T_Scalar > k, const scalar::Precision< T_Scalar > r, unsigned int n)
  {
    if(n == 0) {
      // dr2_h1 = -0.5*(k**2) * (besselh(0,1,k*r) - besselh(2,1,k*r) );
      return -scalar::Precision< T_Scalar >(0.5) * k * k * (boost::math::cyl_hankel_1(0, k * r) - boost::math::cyl_hankel_1(2, k * r));
    }
    else if(n == 1) {
      // dr2_h1 = 0.25*(k**2) * (- 3*besselh(1,1,k*b)+besselh(3,1,k*b));
      return scalar::Precision< T_Scalar >(0.25) * k * k * (-scalar::Precision< T_Scalar >(3.) * boost::math::cyl_hankel_1(1, k * r) + boost::math::cyl_hankel_1(3, k * r));
    }
    // dr2_h1 = 0.25*(k**2) * (besselh(n-2,1,k*r)- 2*besselh(n,1,k*r)+besselh(n+2,1,k*r));
    return scalar::Precision< T_Scalar >(0.25) * k * k * (boost::math::cyl_hankel_1(n - 2, k * r) - scalar::Precision< T_Scalar >(2.) * boost::math::cyl_hankel_1(n, k * r) + boost::math::cyl_hankel_1(n + 2, k * r));
  }

  template< class T_Scalar >
  static T_Scalar dr2BesselH2(const scalar::Precision< T_Scalar > k, const scalar::Precision< T_Scalar > r, unsigned int n)
  {
    if(n == 0) {
      // dr2_h2 = -0.5*(k**2) * (besselh(0,2,k*r) - besselh(2,2,k*r) );
      return -scalar::Precision< T_Scalar >(0.5) * k * k * (boost::math::cyl_hankel_2(0, k * r) - boost::math::cyl_hankel_2(2, k * r));
    }
    else if(n == 1) {
      // dr2_h2 = 0.25*(k**2) * (- 3*besselh(1,2,k*b)+besselh(3,2,k*b));
      return scalar::Precision< T_Scalar >(0.25) * k * k * (-scalar::Precision< T_Scalar >(3.) * boost::math::cyl_hankel_2(1, k * r) + boost::math::cyl_hankel_2(3, k * r));
    }
    return scalar::Precision< T_Scalar >(0.25) * k * k * (boost::math::cyl_hankel_2(n - 2, k * r) - scalar::Precision< T_Scalar >(2.) * boost::math::cyl_hankel_2(n, k * r) + boost::math::cyl_hankel_2(n + 2, k * r));
  }

#endif // HAVE_BOOST


  //********************************
  // SoftPWavesScatteringByACylinder
  //********************************

  template< class T_Scalar >
  SoftPWavesScatteringByACylinder< T_Scalar >::SoftPWavesScatteringByACylinder(const scalar::Precision< T_Scalar > kP, const scalar::Precision< T_Scalar > kS, const scalar::Precision< T_Scalar > R0, const scalar::Precision< T_Scalar > x, const scalar::Precision< T_Scalar > y, const unsigned int nbrOfterm, const scalar::Precision< T_Scalar > thetaInc) :
    _kP(kP), _kS(kS), _R0(R0), _x(x), _y(y), _nbrOfterm(nbrOfterm), _thetaInc(thetaInc)
  {
  }

  template< class T_Scalar >
  SoftPWavesScatteringByACylinder< T_Scalar >::SoftPWavesScatteringByACylinder(const SoftPWavesScatteringByACylinder &other) :
    function::AnalyticalOperation< T_Scalar, Degree::Degree1 >(other), _kP(other._kP), _kS(other._kS), _R0(other._R0), _x(other._x), _y(other._y), _nbrOfterm(other._nbrOfterm), _thetaInc(other._thetaInc)
  {
  }

  template< class T_Scalar >
  SoftPWavesScatteringByACylinder< T_Scalar >::~SoftPWavesScatteringByACylinder()
  {
  }

#ifdef HAVE_BOOST
  template< class T_Scalar >
  void SoftPWavesScatteringByACylinder< T_Scalar >::operator()(function::OutputVector< T_Scalar, Degree::Degree1 > &values, const std::vector< scalar::Precision< T_Scalar >, numa::allocator< scalar::Precision< T_Scalar > > > &points) const
  {
    const T_Scalar powI[4] = {T_Scalar(1., 0.), T_Scalar(0., -1.), T_Scalar(-1., 0), T_Scalar(0., 1.)};

    std::vector< T_Scalar > A(_nbrOfterm);
    std::vector< T_Scalar > B(_nbrOfterm);
    for(auto n = 0U; n < _nbrOfterm; ++n) {
      const T_Scalar m11 = drBesselH1< T_Scalar >(_kP, _R0, n); // dr_h_p
      const T_Scalar m12 = (n / _R0) * T_Scalar(boost::math::cyl_hankel_1(n, _kS * _R0)); // (n/a)*besselh(n,1,ks*a)
      const T_Scalar m21 = -(n / _R0) * T_Scalar(boost::math::cyl_hankel_1(n, _kP * _R0)); // -(n/a)*besselh(n,1,kp*a)
      const T_Scalar m22 = -drBesselH1< T_Scalar >(_kS, _R0, n); // -dr_h_s

      const T_Scalar f_1 = -scalar::Precision< T_Scalar >(!n ? 1. : 2.) * powI[n % 4] * drBesselJ< T_Scalar >(_kP, _R0, n); // f_1 =-phi_0*epsilon*((-1i)**n) *dr_j_p;
      const T_Scalar f_2 = scalar::Precision< T_Scalar >(!n ? 1. : 2.) * powI[n % 4] * (n / _R0) * (boost::math::cyl_bessel_j(n, _kP * _R0)); // f_2 = phi_0*epsilon*((-1i)**n) *(n/a)*besselj(n,kp*a);

      A[n] = scalar::Precision< T_Scalar >(1.) / (m11 * m22 - m12 * m21) * (m22 * f_1 - m12 * f_2);
      B[n] = scalar::Precision< T_Scalar >(1.) / (m11 * m22 - m12 * m21) * (-m21 * f_1 + m11 * f_2);
    }

#pragma omp for
    for(auto i = 0ULL; i < values.size(); ++i) {
      const scalar::Precision< T_Scalar > x = points[3 * i] - _x;
      const scalar::Precision< T_Scalar > y = points[3 * i + 1] - _y;
      const scalar::Precision< T_Scalar > r = std::sqrt(x * x + y * y);
      const scalar::Precision< T_Scalar > theta = std::atan2(y, x) + _thetaInc;
      const scalar::Precision< T_Scalar > invr = 1. / r;

      common::KahanSum< T_Scalar > sumP(0.);
      common::KahanSum< T_Scalar > sumQ(0.);
      for(auto n = 0U; n < _nbrOfterm; ++n) {
        // p = p+ ( An*dr_h1(kp,r,n) + Bn*(dble(n)/r)*besselh(int(n),1,ks*r)) *cos(n*theta);
        sumP += (A[n] * drBesselH1< T_Scalar >(_kP, r, n) + B[n] * (n * invr) * boost::math::cyl_hankel_1(n, _kS * r)) * std::cos(n * theta);
        // q = q+ (-An*(dble(n)/r)*besselh(int(n),1,kp*r) - Bn*dr_h1(ks,r,n)) *sin(n*theta);
        sumQ += (-A[n] * (n * invr) * boost::math::cyl_hankel_1(n, _kP * r) - B[n] * drBesselH1< T_Scalar >(_kS, r, n)) * std::sin(n * theta);
      }

      values[i] << std::cos(theta) * sumP.sum() - std::sin(theta) * sumQ.sum(), std::sin(theta) * sumP.sum() + std::cos(theta) * sumQ.sum(), 0.;
    }
  }
#else

  template< class T_Scalar >
  void SoftPWavesScatteringByACylinder< T_Scalar >::operator()(function::OutputVector< T_Scalar, Degree::Degree1 > &values, const std::vector< scalar::Precision< T_Scalar >, numa::allocator< scalar::Precision< T_Scalar > > > &points) const
  {
    throw common::Exception("Cannot use SoftPWavesScatteringByACylinder without Boost");
  }
#endif // HAVE_BOOST


  INSTANTIATE_CLASS(SoftPWavesScatteringByACylinder, 1, TEMPLATE_ARGS(std::complex< double >))


  //********************************
  // SoftSWavesScatteringByACylinder
  //********************************

  template< class T_Scalar >
  SoftSWavesScatteringByACylinder< T_Scalar >::SoftSWavesScatteringByACylinder(const scalar::Precision< T_Scalar > kP, const scalar::Precision< T_Scalar > kS, const scalar::Precision< T_Scalar > R0, const unsigned int nbrOfterm, const scalar::Precision< T_Scalar > thetaInc) :
    _kP(kP), _kS(kS), _R0(R0), _nbrOfterm(nbrOfterm), _thetaInc(thetaInc)
  {
  }

  template< class T_Scalar >
  SoftSWavesScatteringByACylinder< T_Scalar >::SoftSWavesScatteringByACylinder(const SoftSWavesScatteringByACylinder< T_Scalar > &other) :
    function::AnalyticalOperation< T_Scalar, Degree::Degree1 >(other), _kP(other._kP), _kS(other._kS), _R0(other._R0), _nbrOfterm(other._nbrOfterm), _thetaInc(other._thetaInc)
  {
  }

  template< class T_Scalar >
  SoftSWavesScatteringByACylinder< T_Scalar >::~SoftSWavesScatteringByACylinder()
  {
  }

#ifdef HAVE_BOOST
  template< class T_Scalar >
  void SoftSWavesScatteringByACylinder< T_Scalar >::operator()(function::OutputVector< T_Scalar, Degree::Degree1 > &values, const std::vector< scalar::Precision< T_Scalar >, numa::allocator< scalar::Precision< T_Scalar > > > &points) const
  {
    const T_Scalar powI[4] = {T_Scalar(1., 0.), T_Scalar(0., -1.), T_Scalar(-1., 0), T_Scalar(0., 1.)};
    std::vector< T_Scalar > A(_nbrOfterm);
    std::vector< T_Scalar > B(_nbrOfterm);
    for(auto n = 0U; n < _nbrOfterm; ++n) {
      const T_Scalar m11 = drBesselH1< T_Scalar >(_kP, _R0, n); //dr_h_p;
      const T_Scalar m12 = -(n / _R0) * boost::math::cyl_hankel_1(n, _kS * _R0); // -(n/a)*besselh(n,1,ks*a);
      const T_Scalar m21 = (n / _R0) * boost::math::cyl_hankel_1(n, _kP * _R0); //  (n/a)*besselh(n,1,kp*a);
      const T_Scalar m22 = -drBesselH1< T_Scalar >(_kS, _R0, n); //-dr_h_s;

      const T_Scalar f_1 = scalar::Precision< T_Scalar >(!n ? 1. : 2.) * powI[n % 4] * (n / _R0) * (boost::math::cyl_bessel_j(n, _kS * _R0)); // f_1 = phi_0*epsilon*((-1i)**n)*(n/a)*j_s;
      const T_Scalar f_2 = scalar::Precision< T_Scalar >(!n ? 1. : 2.) * powI[n % 4] * drBesselJ< T_Scalar >(_kS, _R0, n); // f_2 = phi_0*epsilon*((-1i)**n) *dr_j_s;

      A[n] = scalar::Precision< T_Scalar >(1.) / (m11 * m22 - m12 * m21) * (m22 * f_1 - m12 * f_2);
      B[n] = scalar::Precision< T_Scalar >(1.) / (m11 * m22 - m12 * m21) * (-m21 * f_1 + m11 * f_2);
    }

#pragma omp for
    for(auto i = 0ULL; i < values.size(); ++i) {
      const scalar::Precision< T_Scalar > x = points[3 * i];
      const scalar::Precision< T_Scalar > y = points[3 * i + 1];
      const scalar::Precision< T_Scalar > r = std::sqrt(x * x + y * y);
      const scalar::Precision< T_Scalar > theta = std::atan2(y, x) + _thetaInc;
      const scalar::Precision< T_Scalar > invr = 1. / r;

      common::KahanSum< T_Scalar > sumP(0.);
      common::KahanSum< T_Scalar > sumQ(0.);
      for(auto n = 0U; n < _nbrOfterm; ++n) {
        // p = p+ ( An*dr_h(1,kp,r,n) - Bn*(dble(n)/r)*besselh(int(n),1,ks*r)) * sin(n*theta);
        sumP += (A[n] * drBesselH1< T_Scalar >(_kP, r, n) - B[n] * (n * invr) * T_Scalar(boost::math::cyl_hankel_1(n, _kS * r))) * std::sin(n * theta);
        // q = q+ (An*(dble(n)/r)*besselh(int(n),1,kp*r) - Bn*dr_h(1,ks,r,n)) * cos(n*theta);
        sumQ += (A[n] * (n * invr) * T_Scalar(boost::math::cyl_hankel_1(n, _kP * r)) - B[n] * drBesselH1< T_Scalar >(_kS, r, n)) * std::cos(n * theta);
      }

      values[i] << std::cos(theta) * sumP.sum() - std::sin(theta) * sumQ.sum(), std::sin(theta) * sumP.sum() + std::cos(theta) * sumQ.sum(), 0.;
    }
  }
#else

  template< class T_Scalar >
  void SoftSWavesScatteringByACylinder< T_Scalar >::operator()(function::OutputVector< T_Scalar, Degree::Degree1 > &values, const std::vector< scalar::Precision< T_Scalar >, numa::allocator< scalar::Precision< T_Scalar > > > &points) const
  {
    throw common::Exception("Cannot use SoftSWavesScatteringByACylinder without Boost");
  }
#endif // HAVE_BOOST

  INSTANTIATE_CLASS(SoftSWavesScatteringByACylinder, 1, TEMPLATE_ARGS(std::complex< double >))


  //*****************************************
  // SoftPWavesScatteringByACylinderWithLOABC
  //*****************************************

  template< class T_Scalar >
  SoftPWavesScatteringByACylinderWithLOABC< T_Scalar >::SoftPWavesScatteringByACylinderWithLOABC(const scalar::Precision< T_Scalar > kP, const scalar::Precision< T_Scalar > kS, const scalar::Precision< T_Scalar > R0, const scalar::Precision< T_Scalar > R1, const scalar::Precision< T_Scalar > lambda, const scalar::Precision< T_Scalar > mu, const unsigned int nbrOfterm, const scalar::Precision< T_Scalar > thetaInc) :
    _kP(kP), _kS(kS), _R0(R0), _R1(R1), _lambda(lambda), _mu(mu), _nbrOfterm(nbrOfterm), _thetaInc(thetaInc)
  {
  }

  template< class T_Scalar >
  SoftPWavesScatteringByACylinderWithLOABC< T_Scalar >::SoftPWavesScatteringByACylinderWithLOABC(const SoftPWavesScatteringByACylinderWithLOABC< T_Scalar > &other) :
    function::AnalyticalOperation< T_Scalar, Degree::Degree1 >(other), _kP(other._kP), _kS(other._kS), _R0(other._R0), _R1(other._R1), _lambda(other._lambda), _mu(other._mu), _nbrOfterm(other._nbrOfterm), _thetaInc(other._thetaInc)
  {
  }

  template< class T_Scalar >
  SoftPWavesScatteringByACylinderWithLOABC< T_Scalar >::~SoftPWavesScatteringByACylinderWithLOABC()
  {
  }

#ifdef HAVE_BOOST
  template< class T_Scalar >
  void SoftPWavesScatteringByACylinderWithLOABC< T_Scalar >::operator()(function::OutputVector< T_Scalar, Degree::Degree1 > &values, const std::vector< scalar::Precision< T_Scalar >, numa::allocator< scalar::Precision< T_Scalar > > > &points) const
  {
    const T_Scalar powI[4] = {T_Scalar(1., 0.), T_Scalar(0., -1.), T_Scalar(-1., 0.), T_Scalar(0., 1.)};
    const scalar::Precision< T_Scalar > c1 = _lambda + 2. * _mu;
    const T_Scalar c2 = T_Scalar(_lambda / _R1, -_kP * c1);
    const T_Scalar c3 = T_Scalar(_mu / _R1, _kS * _mu);

    std::vector< T_Scalar > AB1(_nbrOfterm);
    std::vector< T_Scalar > AB2(_nbrOfterm);
    std::vector< T_Scalar > AB3(_nbrOfterm);
    std::vector< T_Scalar > AB4(_nbrOfterm);
    for(auto n = 0U; n < _nbrOfterm; ++n) {
      const T_Scalar m11 = drBesselH1< T_Scalar >(_kP, _R0, n); //dr_h_p;
      const T_Scalar m12 = drBesselH2< T_Scalar >(_kP, _R0, n); //dr_h_p;
      const T_Scalar m13 = (n / _R0) * boost::math::cyl_hankel_1(n, _kS * _R0); // -(n/a)*besselh(n,1,ks*a);
      const T_Scalar m14 = (n / _R0) * boost::math::cyl_hankel_2(n, _kS * _R0); // -(n/a)*besselh(n,2,ks*a);

      const T_Scalar m21 = -(n / _R0) * boost::math::cyl_hankel_1(n, _kP * _R0); //  (n/a)*besselh(n,1,kp*a);
      const T_Scalar m22 = -(n / _R0) * boost::math::cyl_hankel_2(n, _kP * _R0); //  (n/a)*besselh(n,2,kp*a);
      const T_Scalar m23 = -drBesselH1< T_Scalar >(_kS, _R0, n); //-dr_h_s;
      const T_Scalar m24 = -drBesselH2< T_Scalar >(_kS, _R0, n); //-dr_h_s;

      const T_Scalar m31 = c1 * dr2BesselH1< T_Scalar >(_kP, _R1, n) - (_lambda * n * n / (_R1 * _R1)) * boost::math::cyl_hankel_1(n, _kP * _R1) + c2 * drBesselH1< T_Scalar >(_kP, _R1, n);
      const T_Scalar m32 = c1 * dr2BesselH2< T_Scalar >(_kP, _R1, n) - (_lambda * n * n / (_R1 * _R1)) * boost::math::cyl_hankel_2(n, _kP * _R1) + c2 * drBesselH2< T_Scalar >(_kP, _R1, n);
      const T_Scalar m33 = c1 * (-scalar::Precision< T_Scalar >(n / std::pow(_R1, 2)) * boost::math::cyl_hankel_1(n, _kS * _R1) + (n / _R1) * drBesselH1< T_Scalar >(_kS, _R1, n)) - ((_lambda * n) / _R1) * drBesselH1< T_Scalar >(_kS, _R1, n) + c2 * (n / _R1) * boost::math::cyl_hankel_1(n, _kS * _R1);
      const T_Scalar m34 = c1 * (-(n / (_R1 * _R1)) * boost::math::cyl_hankel_2(n, _kS * _R1) + (n / _R1) * drBesselH2< T_Scalar >(_kS, _R1, n)) - (_lambda * n) / _R1 * drBesselH2< T_Scalar >(_kS, _R1, n) + c2 * (n / _R1) * boost::math::cyl_hankel_2(n, _kS * _R1);

      const T_Scalar m41 = -_mu * (-(n / (_R1 * _R1))) * boost::math::cyl_hankel_1(n, _kP * _R1) - 2 * _mu * (n / _R1) * drBesselH1< T_Scalar >(_kP, _R1, n) + c3 * (n / _R1) * boost::math::cyl_hankel_1(n, _kP * _R1);
      const T_Scalar m42 = -_mu * (-(n / (_R1 * _R1))) * boost::math::cyl_hankel_2(n, _kP * _R1) - 2 * _mu * (n / _R1) * drBesselH2< T_Scalar >(_kP, _R1, n) + c3 * (n / _R1) * boost::math::cyl_hankel_2(n, _kP * _R1);
      const T_Scalar m43 = -(_mu * n * n / (_R1 * _R1)) * boost::math::cyl_hankel_1(n, _kS * _R1) - _mu * dr2BesselH1< T_Scalar >(_kS, _R1, n) + c3 * drBesselH1< T_Scalar >(_kS, _R1, n);
      const T_Scalar m44 = -(_mu * n * n / (_R1 * _R1)) * boost::math::cyl_hankel_2(n, _kS * _R1) - _mu * dr2BesselH2< T_Scalar >(_kS, _R1, n) + c3 * drBesselH2< T_Scalar >(_kS, _R1, n);

      const T_Scalar f1 = -scalar::Precision< T_Scalar >(!n ? 1. : 2.) * powI[n % 4] * drBesselJ< T_Scalar >(_kP, _R0, n);
      const T_Scalar f2 = scalar::Precision< T_Scalar >(!n ? 1. : 2.) * powI[n % 4] * (n / _R0) * (boost::math::cyl_bessel_j(n, _kP * _R0));
      const T_Scalar f3 = 0.;
      const T_Scalar f4 = 0.;

      // Calculate the inverse determinant of the matrix
      const T_Scalar determinantInverse = scalar::Precision< T_Scalar >(1.0) / (m11 * (m22 * (m33 * m44 - m34 * m43) + m23 * (m34 * m42 - m32 * m44) + m24 * (m32 * m43 - m33 * m42)) - m12 * (m21 * (m33 * m44 - m34 * m43) + m23 * (m34 * m41 - m31 * m44) + m24 * (m31 * m43 - m33 * m41)) + m13 * (m21 * (m32 * m44 - m34 * m42) + m22 * (m34 * m41 - m31 * m44) + m24 * (m31 * m42 - m32 * m41)) - m14 * (m21 * (m32 * m43 - m33 * m42) + m22 * (m33 * m41 - m31 * m43) + m23 * (m31 * m42 - m32 * m41)));

      // Calculate the inverse of the matrix B=M^{-1}
      const T_Scalar inverseM11 = determinantInverse * (m22 * (m33 * m44 - m34 * m43) + m23 * (m34 * m42 - m32 * m44) + m24 * (m32 * m43 - m33 * m42));
      const T_Scalar inverseM21 = determinantInverse * (m21 * (m34 * m43 - m33 * m44) + m23 * (m31 * m44 - m34 * m41) + m24 * (m33 * m41 - m31 * m43));
      const T_Scalar inverseM31 = determinantInverse * (m21 * (m32 * m44 - m34 * m42) + m22 * (m34 * m41 - m31 * m44) + m24 * (m31 * m42 - m32 * m41));
      const T_Scalar inverseM41 = determinantInverse * (m21 * (m33 * m42 - m32 * m43) + m22 * (m31 * m43 - m33 * m41) + m23 * (m32 * m41 - m31 * m42));

      const T_Scalar inverseM12 = determinantInverse * (m12 * (m34 * m43 - m33 * m44) + m13 * (m32 * m44 - m34 * m42) + m14 * (m33 * m42 - m32 * m43));
      const T_Scalar inverseM22 = determinantInverse * (m11 * (m33 * m44 - m34 * m43) + m13 * (m34 * m41 - m31 * m44) + m14 * (m31 * m43 - m33 * m41));
      const T_Scalar inverseM32 = determinantInverse * (m11 * (m34 * m42 - m32 * m44) + m12 * (m31 * m44 - m34 * m41) + m14 * (m32 * m41 - m31 * m42));
      const T_Scalar inverseM42 = determinantInverse * (m11 * (m32 * m43 - m33 * m42) + m12 * (m33 * m41 - m31 * m43) + m13 * (m31 * m42 - m32 * m41));

      const T_Scalar inverseM13 = determinantInverse * (m12 * (m23 * m44 - m24 * m43) + m13 * (m24 * m42 - m22 * m44) + m14 * (m22 * m43 - m23 * m42));
      const T_Scalar inverseM23 = determinantInverse * (m11 * (m24 * m43 - m23 * m44) + m13 * (m21 * m44 - m24 * m41) + m14 * (m23 * m41 - m21 * m43));
      const T_Scalar inverseM33 = determinantInverse * (m11 * (m22 * m44 - m24 * m42) + m12 * (m24 * m41 - m21 * m44) + m14 * (m21 * m42 - m22 * m41));
      const T_Scalar inverseM43 = determinantInverse * (m11 * (m23 * m42 - m22 * m43) + m12 * (m21 * m43 - m23 * m41) + m13 * (m22 * m41 - m21 * m42));

      const T_Scalar inverseM14 = determinantInverse * (m12 * (m24 * m33 - m23 * m34) + m13 * (m22 * m34 - m24 * m32) + m14 * (m23 * m32 - m22 * m33));
      const T_Scalar inverseM24 = determinantInverse * (m11 * (m23 * m34 - m24 * m33) + m13 * (m24 * m31 - m21 * m34) + m14 * (m21 * m33 - m23 * m31));
      const T_Scalar inverseM34 = determinantInverse * (m11 * (m24 * m32 - m22 * m34) + m12 * (m21 * m34 - m24 * m31) + m14 * (m22 * m31 - m21 * m32));
      const T_Scalar inverseM44 = determinantInverse * (m11 * (m22 * m33 - m23 * m32) + m12 * (m23 * m31 - m21 * m33) + m13 * (m21 * m32 - m22 * m31));

      // Compute AB = M^{-1}F
      AB1[n] = inverseM11 * f1 + inverseM12 * f2 + inverseM13 * f3 + inverseM14 * f4;
      AB2[n] = inverseM21 * f1 + inverseM22 * f2 + inverseM23 * f3 + inverseM24 * f4;
      AB3[n] = inverseM31 * f1 + inverseM32 * f2 + inverseM33 * f3 + inverseM34 * f4;
      AB4[n] = inverseM41 * f1 + inverseM42 * f2 + inverseM43 * f3 + inverseM44 * f4;
    }

#pragma omp for
    for(auto i = 0ULL; i < values.size(); ++i) {
      const scalar::Precision< T_Scalar > x = points[3 * i];
      const scalar::Precision< T_Scalar > y = points[3 * i + 1];
      const scalar::Precision< T_Scalar > r = std::sqrt(x * x + y * y);
      const scalar::Precision< T_Scalar > theta = std::atan2(y, x) + _thetaInc;
      const scalar::Precision< T_Scalar > invr = 1. / r;

      common::KahanSum< T_Scalar > sumP(0.);
      common::KahanSum< T_Scalar > sumQ(0.);
      for(auto n = 0U; n < _nbrOfterm; ++n) {
        sumP += (AB1[n] * drBesselH1< T_Scalar >(_kP, r, n) + AB2[n] * drBesselH2< T_Scalar >(_kP, r, n) + (n * invr) * AB3[n] * boost::math::cyl_hankel_1(n, _kS * r) + (n * invr) * AB4[n] * boost::math::cyl_hankel_2(n, _kS * r)) * std::cos(n * theta);
        sumQ += (-(n * invr) * AB1[n] * boost::math::cyl_hankel_1(n, _kP * r) - (n * invr) * AB2[n] * boost::math::cyl_hankel_2(n, _kP * r) - AB3[n] * drBesselH1< T_Scalar >(_kS, r, n) - AB4[n] * drBesselH2< T_Scalar >(_kS, r, n)) * std::sin(n * theta);
      }

      values[i] << std::cos(theta) * sumP.sum() - std::sin(theta) * sumQ.sum(), std::sin(theta) * sumP.sum() + std::cos(theta) * sumQ.sum(), 0.;
    }
  }
#else

  template< class T_Scalar >
  void SoftPWavesScatteringByACylinderWithLOABC< T_Scalar >::operator()(function::OutputVector< T_Scalar, Degree::Degree1 > &values, const std::vector< scalar::Precision< T_Scalar >, numa::allocator< scalar::Precision< T_Scalar > > > &points) const
  {
    throw common::Exception("Cannot use SoftPWavesScatteringByACylinderWithLOABC without Boost");
  }
#endif // HAVE_BOOST

  INSTANTIATE_CLASS(SoftPWavesScatteringByACylinderWithLOABC, 1, TEMPLATE_ARGS(std::complex< double >))


  //*****************************************
  // SoftSWavesScatteringByACylinderWithLOABC
  //*****************************************

  template< class T_Scalar >
  SoftSWavesScatteringByACylinderWithLOABC< T_Scalar >::SoftSWavesScatteringByACylinderWithLOABC(const scalar::Precision< T_Scalar > kP, const scalar::Precision< T_Scalar > kS, const scalar::Precision< T_Scalar > R0, const scalar::Precision< T_Scalar > R1, const scalar::Precision< T_Scalar > lambda, const scalar::Precision< T_Scalar > mu, const unsigned int nbrOfterm, const scalar::Precision< T_Scalar > thetaInc) :
    _kP(kP), _kS(kS), _R0(R0), _R1(R1), _lambda(lambda), _mu(mu), _nbrOfterm(nbrOfterm), _thetaInc(thetaInc)
  {
  }

  template< class T_Scalar >
  SoftSWavesScatteringByACylinderWithLOABC< T_Scalar >::SoftSWavesScatteringByACylinderWithLOABC(const SoftSWavesScatteringByACylinderWithLOABC< T_Scalar > &other) :
    function::AnalyticalOperation< T_Scalar, Degree::Degree1 >(other), _kP(other._kP), _kS(other._kS), _R0(other._R0), _R1(other._R1), _lambda(other._lambda), _mu(other._mu), _nbrOfterm(other._nbrOfterm), _thetaInc(other._thetaInc)
  {
  }

  template< class T_Scalar >
  SoftSWavesScatteringByACylinderWithLOABC< T_Scalar >::~SoftSWavesScatteringByACylinderWithLOABC()
  {
  }

#ifdef HAVE_BOOST
  template< class T_Scalar >
  void SoftSWavesScatteringByACylinderWithLOABC< T_Scalar >::operator()(function::OutputVector< T_Scalar, Degree::Degree1 > &values, const std::vector< scalar::Precision< T_Scalar >, numa::allocator< scalar::Precision< T_Scalar > > > &points) const
  {
    const T_Scalar powI[4] = {T_Scalar(1., 0.), T_Scalar(0., -1.), T_Scalar(-1., 0.), T_Scalar(0., 1.)};
    const scalar::Precision< T_Scalar > c1 = _lambda + 2. * _mu;
    const T_Scalar c2 = T_Scalar(_lambda / _R1, -_kP * c1);
    const T_Scalar c3 = T_Scalar(_mu / _R1, _kS * _mu);

    std::vector< T_Scalar > AB1(_nbrOfterm);
    std::vector< T_Scalar > AB2(_nbrOfterm);
    std::vector< T_Scalar > AB3(_nbrOfterm);
    std::vector< T_Scalar > AB4(_nbrOfterm);
    for(auto n = 0U; n < _nbrOfterm; ++n) {
      const T_Scalar m11 = drBesselH1< T_Scalar >(_kP, _R0, n);
      const T_Scalar m12 = drBesselH2< T_Scalar >(_kP, _R0, n);
      const T_Scalar m13 = -(n / _R0) * boost::math::cyl_hankel_1(n, _kS * _R0);
      const T_Scalar m14 = -(n / _R0) * boost::math::cyl_hankel_2(n, _kS * _R0);

      const T_Scalar m21 = (n / _R0) * boost::math::cyl_hankel_1(n, _kP * _R0);
      const T_Scalar m22 = (n / _R0) * boost::math::cyl_hankel_2(n, _kP * _R0);
      const T_Scalar m23 = -drBesselH1< T_Scalar >(_kS, _R0, n);
      const T_Scalar m24 = -drBesselH2< T_Scalar >(_kS, _R0, n);

      const T_Scalar m31 = c1 * dr2BesselH1< T_Scalar >(_kP, _R1, n) - (_lambda * n * n / (_R1 * _R1)) * boost::math::cyl_hankel_1(n, _kP * _R1) + c2 * drBesselH1< T_Scalar >(_kP, _R1, n);
      const T_Scalar m32 = c1 * dr2BesselH2< T_Scalar >(_kP, _R1, n) - (_lambda * n * n / (_R1 * _R1)) * boost::math::cyl_hankel_2(n, _kP * _R1) + c2 * drBesselH2< T_Scalar >(_kP, _R1, n);
      const T_Scalar m33 = -c1 * (-scalar::Precision< T_Scalar >(n / pow(_R1, 2)) * boost::math::cyl_hankel_1(n, _kS * _R1) + (n / _R1) * drBesselH1< T_Scalar >(_kS, _R1, n)) + ((_lambda * n) / _R1) * drBesselH1< T_Scalar >(_kS, _R1, n) - c2 * (n / _R1) * boost::math::cyl_hankel_1(n, _kS * _R1);
      const T_Scalar m34 = -c1 * (-(n / (_R1 * _R1)) * boost::math::cyl_hankel_2(n, _kS * _R1) + (n / _R1) * drBesselH2< T_Scalar >(_kS, _R1, n)) + (_lambda * n) / _R1 * drBesselH2< T_Scalar >(_kS, _R1, n) - c2 * (n / _R1) * boost::math::cyl_hankel_2(n, _kS * _R1);

      const T_Scalar m41 = _mu * (-(n / (_R1 * _R1))) * boost::math::cyl_hankel_1(n, _kP * _R1) + 2 * _mu * (n / _R1) * drBesselH1< T_Scalar >(_kP, _R1, n) - c3 * (n / _R1) * boost::math::cyl_hankel_1(n, _kP * _R1);
      const T_Scalar m42 = _mu * (-(n / (_R1 * _R1))) * boost::math::cyl_hankel_2(n, _kP * _R1) + 2 * _mu * (n / _R1) * drBesselH2< T_Scalar >(_kP, _R1, n) - c3 * (n / _R1) * boost::math::cyl_hankel_2(n, _kP * _R1);
      const T_Scalar m43 = -(_mu * n * n / (_R1 * _R1)) * boost::math::cyl_hankel_1(n, _kS * _R1) - _mu * dr2BesselH1< T_Scalar >(_kS, _R1, n) + c3 * drBesselH1< T_Scalar >(_kS, _R1, n);
      const T_Scalar m44 = -(_mu * n * n / (_R1 * _R1)) * boost::math::cyl_hankel_2(n, _kS * _R1) - _mu * dr2BesselH2< T_Scalar >(_kS, _R1, n) + c3 * drBesselH2< T_Scalar >(_kS, _R1, n);

      const T_Scalar f1 = scalar::Precision< T_Scalar >(!n ? 1. : 2.) * powI[n % 4] * (n / _R0) * (boost::math::cyl_bessel_j(n, _kS * _R0)); // f_1 = phi_0*epsilon*((-1i)**n)*(n/a)*j_s;
      const T_Scalar f2 = scalar::Precision< T_Scalar >(!n ? 1. : 2.) * powI[n % 4] * drBesselJ< T_Scalar >(_kS, _R0, n); // f_2 = phi_0*epsilon*((-1i)**n) *dr_j_s;
      const T_Scalar f3 = 0.;
      const T_Scalar f4 = 0.;

      // Calculate the inverse determinant of the matrix
      const T_Scalar determinantInverse = scalar::Precision< T_Scalar >(1.0) / (m11 * (m22 * (m33 * m44 - m34 * m43) + m23 * (m34 * m42 - m32 * m44) + m24 * (m32 * m43 - m33 * m42)) - m12 * (m21 * (m33 * m44 - m34 * m43) + m23 * (m34 * m41 - m31 * m44) + m24 * (m31 * m43 - m33 * m41)) + m13 * (m21 * (m32 * m44 - m34 * m42) + m22 * (m34 * m41 - m31 * m44) + m24 * (m31 * m42 - m32 * m41)) - m14 * (m21 * (m32 * m43 - m33 * m42) + m22 * (m33 * m41 - m31 * m43) + m23 * (m31 * m42 - m32 * m41)));

      // Calculate the inverse of the matrix B=M^{-1}
      const T_Scalar inverseM11 = determinantInverse * (m22 * (m33 * m44 - m34 * m43) + m23 * (m34 * m42 - m32 * m44) + m24 * (m32 * m43 - m33 * m42));
      const T_Scalar inverseM21 = determinantInverse * (m21 * (m34 * m43 - m33 * m44) + m23 * (m31 * m44 - m34 * m41) + m24 * (m33 * m41 - m31 * m43));
      const T_Scalar inverseM31 = determinantInverse * (m21 * (m32 * m44 - m34 * m42) + m22 * (m34 * m41 - m31 * m44) + m24 * (m31 * m42 - m32 * m41));
      const T_Scalar inverseM41 = determinantInverse * (m21 * (m33 * m42 - m32 * m43) + m22 * (m31 * m43 - m33 * m41) + m23 * (m32 * m41 - m31 * m42));

      const T_Scalar inverseM12 = determinantInverse * (m12 * (m34 * m43 - m33 * m44) + m13 * (m32 * m44 - m34 * m42) + m14 * (m33 * m42 - m32 * m43));
      const T_Scalar inverseM22 = determinantInverse * (m11 * (m33 * m44 - m34 * m43) + m13 * (m34 * m41 - m31 * m44) + m14 * (m31 * m43 - m33 * m41));
      const T_Scalar inverseM32 = determinantInverse * (m11 * (m34 * m42 - m32 * m44) + m12 * (m31 * m44 - m34 * m41) + m14 * (m32 * m41 - m31 * m42));
      const T_Scalar inverseM42 = determinantInverse * (m11 * (m32 * m43 - m33 * m42) + m12 * (m33 * m41 - m31 * m43) + m13 * (m31 * m42 - m32 * m41));

      const T_Scalar inverseM13 = determinantInverse * (m12 * (m23 * m44 - m24 * m43) + m13 * (m24 * m42 - m22 * m44) + m14 * (m22 * m43 - m23 * m42));
      const T_Scalar inverseM23 = determinantInverse * (m11 * (m24 * m43 - m23 * m44) + m13 * (m21 * m44 - m24 * m41) + m14 * (m23 * m41 - m21 * m43));
      const T_Scalar inverseM33 = determinantInverse * (m11 * (m22 * m44 - m24 * m42) + m12 * (m24 * m41 - m21 * m44) + m14 * (m21 * m42 - m22 * m41));
      const T_Scalar inverseM43 = determinantInverse * (m11 * (m23 * m42 - m22 * m43) + m12 * (m21 * m43 - m23 * m41) + m13 * (m22 * m41 - m21 * m42));

      const T_Scalar inverseM14 = determinantInverse * (m12 * (m24 * m33 - m23 * m34) + m13 * (m22 * m34 - m24 * m32) + m14 * (m23 * m32 - m22 * m33));
      const T_Scalar inverseM24 = determinantInverse * (m11 * (m23 * m34 - m24 * m33) + m13 * (m24 * m31 - m21 * m34) + m14 * (m21 * m33 - m23 * m31));
      const T_Scalar inverseM34 = determinantInverse * (m11 * (m24 * m32 - m22 * m34) + m12 * (m21 * m34 - m24 * m31) + m14 * (m22 * m31 - m21 * m32));
      const T_Scalar inverseM44 = determinantInverse * (m11 * (m22 * m33 - m23 * m32) + m12 * (m23 * m31 - m21 * m33) + m13 * (m21 * m32 - m22 * m31));

      // Compute AB = M^{-1}F
      AB1[n] = inverseM11 * f1 + inverseM12 * f2 + inverseM13 * f3 + inverseM14 * f4;
      AB2[n] = inverseM21 * f1 + inverseM22 * f2 + inverseM23 * f3 + inverseM24 * f4;
      AB3[n] = inverseM31 * f1 + inverseM32 * f2 + inverseM33 * f3 + inverseM34 * f4;
      AB4[n] = inverseM41 * f1 + inverseM42 * f2 + inverseM43 * f3 + inverseM44 * f4;
    }

#pragma omp for
    for(auto i = 0ULL; i < values.size(); ++i) {
      const scalar::Precision< T_Scalar > x = points[3 * i];
      const scalar::Precision< T_Scalar > y = points[3 * i + 1];
      const scalar::Precision< T_Scalar > r = std::sqrt(x * x + y * y);
      const scalar::Precision< T_Scalar > theta = std::atan2(y, x) + _thetaInc;
      const scalar::Precision< T_Scalar > invr = 1. / r;

      common::KahanSum< T_Scalar > sumP(0.);
      common::KahanSum< T_Scalar > sumQ(0.);
      for(auto n = 0U; n < _nbrOfterm; ++n) {
        sumP += (AB1[n] * drBesselH1< T_Scalar >(_kP, r, n) + AB2[n] * drBesselH2< T_Scalar >(_kP, r, n) - (n * invr) * AB3[n] * boost::math::cyl_hankel_1(n, _kS * r) - (n * invr) * AB4[n] * boost::math::cyl_hankel_2(n, _kS * r)) * std::sin(n * theta);
        sumQ += ((n * invr) * AB1[n] * boost::math::cyl_hankel_1(n, _kP * r) + (n * invr) * AB2[n] * boost::math::cyl_hankel_2(n, _kP * r) - AB3[n] * drBesselH1< T_Scalar >(_kS, r, n) - AB4[n] * drBesselH2< T_Scalar >(_kS, r, n)) * std::cos(n * theta);
      }

      values[i] << std::cos(theta) * sumP.sum() - std::sin(theta) * sumQ.sum(), std::sin(theta) * sumP.sum() + std::cos(theta) * sumQ.sum(), 0.;
    }
  }
#else

  template< class T_Scalar >
  void SoftSWavesScatteringByACylinderWithLOABC< T_Scalar >::operator()(function::OutputVector< T_Scalar, Degree::Degree1 > &values, const std::vector< scalar::Precision< T_Scalar >, numa::allocator< scalar::Precision< T_Scalar > > > &points) const
  {
    throw common::Exception("Cannot use SoftSWavesScatteringByACylinderWithLOABC without Boost");
  }
#endif // HAVE_BOOST

  INSTANTIATE_CLASS(SoftSWavesScatteringByACylinderWithLOABC, 1, TEMPLATE_ARGS(std::complex< double >))


  //*****************************************
  // SoftPWavesScatteringByACylinderWithHOABC
  //*****************************************

  template< class T_Scalar >
  SoftPWavesScatteringByACylinderWithHOABC< T_Scalar >::SoftPWavesScatteringByACylinderWithHOABC(const scalar::Precision< T_Scalar > f, const scalar::Precision< T_Scalar > R0, const scalar::Precision< T_Scalar > R1, const unsigned int nbrOfterm, const scalar::Precision< T_Scalar > lambda, const scalar::Precision< T_Scalar > mu, const scalar::Precision< T_Scalar > rho, const unsigned int padeOrder, const scalar::Precision< T_Scalar > padeAngle, const scalar::Precision< T_Scalar > thetaInc) :
    _f(f), _R0(R0), _R1(R1), _nbrOfterm(nbrOfterm), _lambda(lambda), _mu(mu), _rho(rho), _padeOrder(padeOrder), _padeAngle(padeAngle), _thetaInc(thetaInc)
  {
  }

  template< class T_Scalar >
  SoftPWavesScatteringByACylinderWithHOABC< T_Scalar >::SoftPWavesScatteringByACylinderWithHOABC(const SoftPWavesScatteringByACylinderWithHOABC< T_Scalar > &other) :
    function::AnalyticalOperation< T_Scalar, Degree::Degree1 >(other), _f(other._f), _R0(other._R0), _R1(other._R1), _nbrOfterm(other._nbrOfterm), _lambda(other._lambda), _mu(other._mu), _rho(other._rho), _padeOrder(other._padeOrder), _padeAngle(other._padeAngle), _thetaInc(other._thetaInc)
  {
  }

  template< class T_Scalar >
  SoftPWavesScatteringByACylinderWithHOABC< T_Scalar >::~SoftPWavesScatteringByACylinderWithHOABC()
  {
  }

#ifdef HAVE_BOOST
  template< class T_Scalar >
  void SoftPWavesScatteringByACylinderWithHOABC< T_Scalar >::operator()(function::OutputVector< T_Scalar, Degree::Degree1 > &values, const std::vector< scalar::Precision< T_Scalar >, numa::allocator< scalar::Precision< T_Scalar > > > &points) const
  {
    const T_Scalar powI[4] = {T_Scalar(1., 0.), T_Scalar(0., -1.), T_Scalar(-1., 0.), T_Scalar(0., 1.)};
    const T_Scalar im = T_Scalar(0., 1.);
    const double pi = 3.14159265358979323846264338327950288;

    const double w = 2 * pi * _f;
    const double cP = std::sqrt((_lambda + 2 * _mu) / _rho);
    const double cS = std::sqrt(_mu / _rho);
    const double kP = w / cP;
    const double kS = w / cS;

    T_Scalar padeD[_padeOrder];
    T_Scalar padeC[_padeOrder];
    T_Scalar padeR[_padeOrder];
    T_Scalar padeS[_padeOrder];

    for(auto i = 0U; i < _padeOrder; ++i) {
      padeD[i] = 1.0 + pow(std::tan((pi / (2.0 * (double)_padeOrder)) * (0.5 + (double)i)), 2);
      padeC[i] = padeD[i] / (double)_padeOrder;
      padeR[i] = padeC[i] * T_Scalar(std::cos(_padeAngle / 2.0), std::sin(_padeAngle / 2.0));
      padeS[i] = 1.0 + T_Scalar(std::cos(_padeAngle), std::sin(_padeAngle)) * (-1.0 + padeD[i]);
    }

    const T_Scalar kPEps = kP + im * 0.39 * pow(kP, 1. / 3) * pow(_R1, -2. / 3);
    const T_Scalar kSEps = kS + im * 0.39 * pow(kS, 1. / 3) * pow(_R1, -2. / 3);
    const scalar::Precision< T_Scalar > d1 = _lambda + 2. * _mu;

    std::vector< T_Scalar > AB1(_nbrOfterm);
    std::vector< T_Scalar > AB2(_nbrOfterm);
    std::vector< T_Scalar > AB3(_nbrOfterm);
    std::vector< T_Scalar > AB4(_nbrOfterm);
    for(auto n = 0U; n < _nbrOfterm; ++n) {

      T_Scalar padeXiP = 0.;
      T_Scalar padeXiS = 0.;

      for(auto i = 0U; i < _padeOrder; ++i) {
        padeXiP += padeR[i] / (padeS[i] - pow(n * 1. / (kPEps * _R1), 2));
        padeXiS += padeR[i] / (padeS[i] - pow(n * 1. / (kSEps * _R1), 2));
      }
      padeXiP /= kPEps;
      padeXiS /= kSEps;

      const T_Scalar den = T_Scalar(1. + pow(n / _R1, 2) * padeXiP * padeXiS);
      const T_Scalar d2 = T_Scalar(_rho * pow(w, 2) * padeXiP * padeXiS / den);
      const T_Scalar d3 = T_Scalar(_lambda / _R1 - (im * _rho * pow(w, 2) * padeXiP) / den);
      const T_Scalar d4 = T_Scalar(-_mu / _R1 - (im * _rho * pow(w, 2) * padeXiS) / den);

      const T_Scalar m11 = drBesselH1< T_Scalar >(kP, _R0, n); //dr_h_p;
      const T_Scalar m12 = drBesselH2< T_Scalar >(kP, _R0, n); //dr_h_p;
      const T_Scalar m13 = (n / _R0) * boost::math::cyl_hankel_1(n, kS * _R0); // -(n/a)*besselh(n,1,ks*a);
      const T_Scalar m14 = (n / _R0) * boost::math::cyl_hankel_2(n, kS * _R0); // -(n/a)*besselh(n,2,ks*a);

      const T_Scalar m21 = -(n / _R0) * boost::math::cyl_hankel_1(n, kP * _R0); //  (n/a)*besselh(n,1,kp*a);
      const T_Scalar m22 = -(n / _R0) * boost::math::cyl_hankel_2(n, kP * _R0); //  (n/a)*besselh(n,2,kp*a);
      const T_Scalar m23 = -drBesselH1< T_Scalar >(kS, _R0, n); //-dr_h_s;
      const T_Scalar m24 = -drBesselH2< T_Scalar >(kS, _R0, n); //-dr_h_s;

      const T_Scalar m31 = d1 * dr2BesselH1< T_Scalar >(kP, _R1, n) - (d1 - d2) * pow(n / _R1, 2) * boost::math::cyl_hankel_1(n, kP * _R1) + d3 * drBesselH1< T_Scalar >(kP, _R1, n);
      const T_Scalar m32 = d1 * dr2BesselH2< T_Scalar >(kP, _R1, n) - (d1 - d2) * pow(n / _R1, 2) * boost::math::cyl_hankel_2(n, kP * _R1) + d3 * drBesselH2< T_Scalar >(kP, _R1, n);
      const T_Scalar m33 = d1 * (-scalar::Precision< T_Scalar >(n / std::pow(_R1, 2)) * boost::math::cyl_hankel_1(n, kS * _R1) + (n / _R1) * drBesselH1< T_Scalar >(kS, _R1, n)) - (d1 - d2) * (n / _R1) * drBesselH1< T_Scalar >(kS, _R1, n) + d3 * (n / _R1) * boost::math::cyl_hankel_1(n, kS * _R1);
      const T_Scalar m34 = d1 * (-scalar::Precision< T_Scalar >(n / std::pow(_R1, 2)) * boost::math::cyl_hankel_2(n, kS * _R1) + (n / _R1) * drBesselH2< T_Scalar >(kS, _R1, n)) - (d1 - d2) * (n / _R1) * drBesselH2< T_Scalar >(kS, _R1, n) + d3 * (n / _R1) * boost::math::cyl_hankel_2(n, kS * _R1);

      const T_Scalar m41 = -_mu * (-(n / (_R1 * _R1)) * boost::math::cyl_hankel_1(n, kP * _R1) + (n / _R1) * drBesselH1< T_Scalar >(kP, _R1, n)) - (d2 - _mu) * (n / _R1) * drBesselH1< T_Scalar >(kP, _R1, n) - d4 * (n / _R1) * boost::math::cyl_hankel_1(n, kP * _R1);
      const T_Scalar m42 = -_mu * (-(n / (_R1 * _R1)) * boost::math::cyl_hankel_2(n, kP * _R1) + (n / _R1) * drBesselH2< T_Scalar >(kP, _R1, n)) - (d2 - _mu) * (n / _R1) * drBesselH2< T_Scalar >(kP, _R1, n) - d4 * (n / _R1) * boost::math::cyl_hankel_2(n, kP * _R1);
      const T_Scalar m43 = -_mu * dr2BesselH1< T_Scalar >(kS, _R1, n) - (d2 - _mu) * pow(n / _R1, 2) * boost::math::cyl_hankel_1(n, kS * _R1) - d4 * drBesselH1< T_Scalar >(kS, _R1, n);
      const T_Scalar m44 = -_mu * dr2BesselH2< T_Scalar >(kS, _R1, n) - (d2 - _mu) * pow(n / _R1, 2) * boost::math::cyl_hankel_2(n, kS * _R1) - d4 * drBesselH2< T_Scalar >(kS, _R1, n);

      const T_Scalar f1 = -scalar::Precision< T_Scalar >(!n ? 1. : 2.) * powI[n % 4] * drBesselJ< T_Scalar >(kP, _R0, n);
      const T_Scalar f2 = scalar::Precision< T_Scalar >(!n ? 1. : 2.) * powI[n % 4] * (n / _R0) * (boost::math::cyl_bessel_j(n, kP * _R0));
      const T_Scalar f3 = 0.;
      const T_Scalar f4 = 0.;

      // Calculate the inverse determinant of the matrix
      const T_Scalar determinantInverse = scalar::Precision< T_Scalar >(1.0) / (m11 * (m22 * (m33 * m44 - m34 * m43) + m23 * (m34 * m42 - m32 * m44) + m24 * (m32 * m43 - m33 * m42)) - m12 * (m21 * (m33 * m44 - m34 * m43) + m23 * (m34 * m41 - m31 * m44) + m24 * (m31 * m43 - m33 * m41)) + m13 * (m21 * (m32 * m44 - m34 * m42) + m22 * (m34 * m41 - m31 * m44) + m24 * (m31 * m42 - m32 * m41)) - m14 * (m21 * (m32 * m43 - m33 * m42) + m22 * (m33 * m41 - m31 * m43) + m23 * (m31 * m42 - m32 * m41)));

      // Calculate the inverse of the matrix B=M^{-1}
      const T_Scalar inverseM11 = determinantInverse * (m22 * (m33 * m44 - m34 * m43) + m23 * (m34 * m42 - m32 * m44) + m24 * (m32 * m43 - m33 * m42));
      const T_Scalar inverseM21 = determinantInverse * (m21 * (m34 * m43 - m33 * m44) + m23 * (m31 * m44 - m34 * m41) + m24 * (m33 * m41 - m31 * m43));
      const T_Scalar inverseM31 = determinantInverse * (m21 * (m32 * m44 - m34 * m42) + m22 * (m34 * m41 - m31 * m44) + m24 * (m31 * m42 - m32 * m41));
      const T_Scalar inverseM41 = determinantInverse * (m21 * (m33 * m42 - m32 * m43) + m22 * (m31 * m43 - m33 * m41) + m23 * (m32 * m41 - m31 * m42));

      const T_Scalar inverseM12 = determinantInverse * (m12 * (m34 * m43 - m33 * m44) + m13 * (m32 * m44 - m34 * m42) + m14 * (m33 * m42 - m32 * m43));
      const T_Scalar inverseM22 = determinantInverse * (m11 * (m33 * m44 - m34 * m43) + m13 * (m34 * m41 - m31 * m44) + m14 * (m31 * m43 - m33 * m41));
      const T_Scalar inverseM32 = determinantInverse * (m11 * (m34 * m42 - m32 * m44) + m12 * (m31 * m44 - m34 * m41) + m14 * (m32 * m41 - m31 * m42));
      const T_Scalar inverseM42 = determinantInverse * (m11 * (m32 * m43 - m33 * m42) + m12 * (m33 * m41 - m31 * m43) + m13 * (m31 * m42 - m32 * m41));

      const T_Scalar inverseM13 = determinantInverse * (m12 * (m23 * m44 - m24 * m43) + m13 * (m24 * m42 - m22 * m44) + m14 * (m22 * m43 - m23 * m42));
      const T_Scalar inverseM23 = determinantInverse * (m11 * (m24 * m43 - m23 * m44) + m13 * (m21 * m44 - m24 * m41) + m14 * (m23 * m41 - m21 * m43));
      const T_Scalar inverseM33 = determinantInverse * (m11 * (m22 * m44 - m24 * m42) + m12 * (m24 * m41 - m21 * m44) + m14 * (m21 * m42 - m22 * m41));
      const T_Scalar inverseM43 = determinantInverse * (m11 * (m23 * m42 - m22 * m43) + m12 * (m21 * m43 - m23 * m41) + m13 * (m22 * m41 - m21 * m42));

      const T_Scalar inverseM14 = determinantInverse * (m12 * (m24 * m33 - m23 * m34) + m13 * (m22 * m34 - m24 * m32) + m14 * (m23 * m32 - m22 * m33));
      const T_Scalar inverseM24 = determinantInverse * (m11 * (m23 * m34 - m24 * m33) + m13 * (m24 * m31 - m21 * m34) + m14 * (m21 * m33 - m23 * m31));
      const T_Scalar inverseM34 = determinantInverse * (m11 * (m24 * m32 - m22 * m34) + m12 * (m21 * m34 - m24 * m31) + m14 * (m22 * m31 - m21 * m32));
      const T_Scalar inverseM44 = determinantInverse * (m11 * (m22 * m33 - m23 * m32) + m12 * (m23 * m31 - m21 * m33) + m13 * (m21 * m32 - m22 * m31));

      // Compute AB = M^{-1}F
      AB1[n] = inverseM11 * f1 + inverseM12 * f2 + inverseM13 * f3 + inverseM14 * f4;
      AB2[n] = inverseM21 * f1 + inverseM22 * f2 + inverseM23 * f3 + inverseM24 * f4;
      AB3[n] = inverseM31 * f1 + inverseM32 * f2 + inverseM33 * f3 + inverseM34 * f4;
      AB4[n] = inverseM41 * f1 + inverseM42 * f2 + inverseM43 * f3 + inverseM44 * f4;
    }

#pragma omp for
    for(auto i = 0ULL; i < values.size(); ++i) {
      const scalar::Precision< T_Scalar > x = points[3 * i];
      const scalar::Precision< T_Scalar > y = points[3 * i + 1];
      const scalar::Precision< T_Scalar > r = std::sqrt(x * x + y * y);
      const scalar::Precision< T_Scalar > theta = std::atan2(y, x) + _thetaInc;
      const scalar::Precision< T_Scalar > invr = 1. / r;

      common::KahanSum< T_Scalar > sumP(0.);
      common::KahanSum< T_Scalar > sumQ(0.);
      for(auto n = 0U; n < _nbrOfterm; ++n) {
        sumP += (AB1[n] * drBesselH1< T_Scalar >(kP, r, n) + AB2[n] * drBesselH2< T_Scalar >(kP, r, n) + (n * invr) * AB3[n] * boost::math::cyl_hankel_1(n, kS * r) + (n * invr) * AB4[n] * boost::math::cyl_hankel_2(n, kS * r)) * std::cos(n * theta);
        sumQ += (-(n * invr) * AB1[n] * boost::math::cyl_hankel_1(n, kP * r) - (n * invr) * AB2[n] * boost::math::cyl_hankel_2(n, kP * r) - AB3[n] * drBesselH1< T_Scalar >(kS, r, n) - AB4[n] * drBesselH2< T_Scalar >(kS, r, n)) * std::sin(n * theta);
      }

      values[i] << std::cos(theta) * sumP.sum() - std::sin(theta) * sumQ.sum(), std::sin(theta) * sumP.sum() + std::cos(theta) * sumQ.sum(), 0.;
    }
  }
#else

  template< class T_Scalar >
  void SoftPWavesScatteringByACylinderWithHOABC< T_Scalar >::operator()(function::OutputVector< T_Scalar, Degree::Degree1 > &values, const std::vector< scalar::Precision< T_Scalar >, numa::allocator< scalar::Precision< T_Scalar > > > &points) const
  {
    throw common::Exception("Cannot use SoftPWavesScatteringByACylinderWithHOABC without Boost");
  }
#endif // HAVE_BOOST

  INSTANTIATE_CLASS(SoftPWavesScatteringByACylinderWithHOABC, 1, TEMPLATE_ARGS(std::complex< double >))


  //*****************************************
  // SoftSWavesScatteringByACylinderWithHOABC
  //*****************************************

  template< class T_Scalar >
  SoftSWavesScatteringByACylinderWithHOABC< T_Scalar >::SoftSWavesScatteringByACylinderWithHOABC(const scalar::Precision< T_Scalar > f, const scalar::Precision< T_Scalar > R0, const scalar::Precision< T_Scalar > R1, const unsigned int nbrOfterm, const scalar::Precision< T_Scalar > lambda, const scalar::Precision< T_Scalar > mu, const scalar::Precision< T_Scalar > rho, const unsigned int padeOrder, const scalar::Precision< T_Scalar > padeAngle, const scalar::Precision< T_Scalar > thetaInc) :
    _f(f), _R0(R0), _R1(R1), _nbrOfterm(nbrOfterm), _lambda(lambda), _mu(mu), _rho(rho), _padeOrder(padeOrder), _padeAngle(padeAngle), _thetaInc(thetaInc)
  {
  }

  template< class T_Scalar >
  SoftSWavesScatteringByACylinderWithHOABC< T_Scalar >::SoftSWavesScatteringByACylinderWithHOABC(const SoftSWavesScatteringByACylinderWithHOABC< T_Scalar > &other) :
    function::AnalyticalOperation< T_Scalar, Degree::Degree1 >(other), _f(other._f), _R0(other._R0), _R1(other._R1), _nbrOfterm(other._nbrOfterm), _lambda(other._lambda), _mu(other._mu), _rho(other._rho), _padeOrder(other._padeOrder), _padeAngle(other._padeAngle), _thetaInc(other._thetaInc)
  {
  }

  template< class T_Scalar >
  SoftSWavesScatteringByACylinderWithHOABC< T_Scalar >::~SoftSWavesScatteringByACylinderWithHOABC()
  {
  }

#ifdef HAVE_BOOST
  template< class T_Scalar >
  void SoftSWavesScatteringByACylinderWithHOABC< T_Scalar >::operator()(function::OutputVector< T_Scalar, Degree::Degree1 > &values, const std::vector< scalar::Precision< T_Scalar >, numa::allocator< scalar::Precision< T_Scalar > > > &points) const
  {
    const T_Scalar powI[4] = {T_Scalar(1., 0.), T_Scalar(0., -1.), T_Scalar(-1., 0.), T_Scalar(0., 1.)};
    const T_Scalar im = T_Scalar(0., 1.);
    const double pi = 3.14159265358979323846264338327950288;

    const double w = 2 * pi * _f;
    const double cP = std::sqrt((_lambda + 2 * _mu) / _rho);
    const double cS = std::sqrt(_mu / _rho);
    const double kP = w / cP;
    const double kS = w / cS;

    T_Scalar padeD[_padeOrder];
    T_Scalar padeC[_padeOrder];
    T_Scalar padeR[_padeOrder];
    T_Scalar padeS[_padeOrder];

    for(auto i = 0U; i < _padeOrder; ++i) {
      padeD[i] = 1.0 + pow(std::tan((pi / (2.0 * (double)_padeOrder)) * (0.5 + (double)i)), 2);
      padeC[i] = padeD[i] / (double)_padeOrder;
      padeR[i] = padeC[i] * T_Scalar(std::cos(_padeAngle / 2.0), std::sin(_padeAngle / 2.0));
      padeS[i] = 1.0 + T_Scalar(std::cos(_padeAngle), std::sin(_padeAngle)) * (-1.0 + padeD[i]);
    }

    const T_Scalar kPEps = kP + im * 0.39 * pow(kP, 1. / 3) * pow(_R1, -2. / 3);
    const T_Scalar kSEps = kS + im * 0.39 * pow(kS, 1. / 3) * pow(_R1, -2. / 3);
    const scalar::Precision< T_Scalar > d1 = _lambda + 2. * _mu;

    std::vector< T_Scalar > AB1(_nbrOfterm);
    std::vector< T_Scalar > AB2(_nbrOfterm);
    std::vector< T_Scalar > AB3(_nbrOfterm);
    std::vector< T_Scalar > AB4(_nbrOfterm);
    for(auto n = 0U; n < _nbrOfterm; ++n) {

      T_Scalar padeXiP = 0.;
      T_Scalar padeXiS = 0.;

      for(auto i = 0U; i < _padeOrder; ++i) {
        padeXiP += padeR[i] / (padeS[i] - pow(n * 1. / (kPEps * _R1), 2));
        padeXiS += padeR[i] / (padeS[i] - pow(n * 1. / (kSEps * _R1), 2));
      }
      padeXiP /= kPEps;
      padeXiS /= kSEps;

      const T_Scalar den = T_Scalar(1. + pow(n / _R1, 2) * padeXiP * padeXiS);
      const T_Scalar d2 = T_Scalar(_rho * pow(w, 2) * padeXiP * padeXiS / den);
      const T_Scalar d3 = T_Scalar(_lambda / _R1 - (im * _rho * pow(w, 2) * padeXiP) / den);
      const T_Scalar d4 = T_Scalar(-_mu / _R1 - (im * _rho * pow(w, 2) * padeXiS) / den);

      const T_Scalar m11 = drBesselH1< T_Scalar >(kP, _R0, n);
      const T_Scalar m12 = drBesselH2< T_Scalar >(kP, _R0, n);
      const T_Scalar m13 = -(n / _R0) * boost::math::cyl_hankel_1(n, kS * _R0);
      const T_Scalar m14 = -(n / _R0) * boost::math::cyl_hankel_2(n, kS * _R0);

      const T_Scalar m21 = (n / _R0) * boost::math::cyl_hankel_1(n, kP * _R0);
      const T_Scalar m22 = (n / _R0) * boost::math::cyl_hankel_2(n, kP * _R0);
      const T_Scalar m23 = -drBesselH1< T_Scalar >(kS, _R0, n);
      const T_Scalar m24 = -drBesselH2< T_Scalar >(kS, _R0, n);

      const T_Scalar m31 = d1 * dr2BesselH1< T_Scalar >(kP, _R1, n) - (d1 - d2) * pow(n / _R1, 2) * boost::math::cyl_hankel_1(n, kP * _R1) + d3 * drBesselH1< T_Scalar >(kP, _R1, n);
      const T_Scalar m32 = d1 * dr2BesselH2< T_Scalar >(kP, _R1, n) - (d1 - d2) * pow(n / _R1, 2) * boost::math::cyl_hankel_2(n, kP * _R1) + d3 * drBesselH2< T_Scalar >(kP, _R1, n);
      const T_Scalar m33 = -d1 * (-scalar::Precision< T_Scalar >(n / std::pow(_R1, 2)) * boost::math::cyl_hankel_1(n, kS * _R1) + (n / _R1) * drBesselH1< T_Scalar >(kS, _R1, n)) + (d1 - d2) * (n / _R1) * drBesselH1< T_Scalar >(kS, _R1, n) - d3 * (n / _R1) * boost::math::cyl_hankel_1(n, kS * _R1);
      const T_Scalar m34 = -d1 * (-scalar::Precision< T_Scalar >(n / std::pow(_R1, 2)) * boost::math::cyl_hankel_2(n, kS * _R1) + (n / _R1) * drBesselH2< T_Scalar >(kS, _R1, n)) + (d1 - d2) * (n / _R1) * drBesselH2< T_Scalar >(kS, _R1, n) - d3 * (n / _R1) * boost::math::cyl_hankel_2(n, kS * _R1);

      const T_Scalar m41 = _mu * (-(n / (_R1 * _R1)) * boost::math::cyl_hankel_1(n, kP * _R1) + (n / _R1) * drBesselH1< T_Scalar >(kP, _R1, n)) + (d2 - _mu) * (n / _R1) * drBesselH1< T_Scalar >(kP, _R1, n) + d4 * (n / _R1) * boost::math::cyl_hankel_1(n, kP * _R1);
      const T_Scalar m42 = _mu * (-(n / (_R1 * _R1)) * boost::math::cyl_hankel_2(n, kP * _R1) + (n / _R1) * drBesselH2< T_Scalar >(kP, _R1, n)) + (d2 - _mu) * (n / _R1) * drBesselH2< T_Scalar >(kP, _R1, n) + d4 * (n / _R1) * boost::math::cyl_hankel_2(n, kP * _R1);
      const T_Scalar m43 = -_mu * dr2BesselH1< T_Scalar >(kS, _R1, n) - (d2 - _mu) * pow(n / _R1, 2) * boost::math::cyl_hankel_1(n, kS * _R1) - d4 * drBesselH1< T_Scalar >(kS, _R1, n);
      const T_Scalar m44 = -_mu * dr2BesselH2< T_Scalar >(kS, _R1, n) - (d2 - _mu) * pow(n / _R1, 2) * boost::math::cyl_hankel_2(n, kS * _R1) - d4 * drBesselH2< T_Scalar >(kS, _R1, n);

      const T_Scalar f1 = scalar::Precision< T_Scalar >(!n ? 1. : 2.) * powI[n % 4] * (n / _R0) * (boost::math::cyl_bessel_j(n, kS * _R0));
      const T_Scalar f2 = scalar::Precision< T_Scalar >(!n ? 1. : 2.) * powI[n % 4] * drBesselJ< T_Scalar >(kS, _R0, n);
      const T_Scalar f3 = 0.;
      const T_Scalar f4 = 0.;

      // Calculate the inverse determinant of the matrix
      const T_Scalar determinantInverse = scalar::Precision< T_Scalar >(1.0) / (m11 * (m22 * (m33 * m44 - m34 * m43) + m23 * (m34 * m42 - m32 * m44) + m24 * (m32 * m43 - m33 * m42)) - m12 * (m21 * (m33 * m44 - m34 * m43) + m23 * (m34 * m41 - m31 * m44) + m24 * (m31 * m43 - m33 * m41)) + m13 * (m21 * (m32 * m44 - m34 * m42) + m22 * (m34 * m41 - m31 * m44) + m24 * (m31 * m42 - m32 * m41)) - m14 * (m21 * (m32 * m43 - m33 * m42) + m22 * (m33 * m41 - m31 * m43) + m23 * (m31 * m42 - m32 * m41)));

      // Calculate the inverse of the matrix B=M^{-1}
      const T_Scalar inverseM11 = determinantInverse * (m22 * (m33 * m44 - m34 * m43) + m23 * (m34 * m42 - m32 * m44) + m24 * (m32 * m43 - m33 * m42));
      const T_Scalar inverseM21 = determinantInverse * (m21 * (m34 * m43 - m33 * m44) + m23 * (m31 * m44 - m34 * m41) + m24 * (m33 * m41 - m31 * m43));
      const T_Scalar inverseM31 = determinantInverse * (m21 * (m32 * m44 - m34 * m42) + m22 * (m34 * m41 - m31 * m44) + m24 * (m31 * m42 - m32 * m41));
      const T_Scalar inverseM41 = determinantInverse * (m21 * (m33 * m42 - m32 * m43) + m22 * (m31 * m43 - m33 * m41) + m23 * (m32 * m41 - m31 * m42));

      const T_Scalar inverseM12 = determinantInverse * (m12 * (m34 * m43 - m33 * m44) + m13 * (m32 * m44 - m34 * m42) + m14 * (m33 * m42 - m32 * m43));
      const T_Scalar inverseM22 = determinantInverse * (m11 * (m33 * m44 - m34 * m43) + m13 * (m34 * m41 - m31 * m44) + m14 * (m31 * m43 - m33 * m41));
      const T_Scalar inverseM32 = determinantInverse * (m11 * (m34 * m42 - m32 * m44) + m12 * (m31 * m44 - m34 * m41) + m14 * (m32 * m41 - m31 * m42));
      const T_Scalar inverseM42 = determinantInverse * (m11 * (m32 * m43 - m33 * m42) + m12 * (m33 * m41 - m31 * m43) + m13 * (m31 * m42 - m32 * m41));

      const T_Scalar inverseM13 = determinantInverse * (m12 * (m23 * m44 - m24 * m43) + m13 * (m24 * m42 - m22 * m44) + m14 * (m22 * m43 - m23 * m42));
      const T_Scalar inverseM23 = determinantInverse * (m11 * (m24 * m43 - m23 * m44) + m13 * (m21 * m44 - m24 * m41) + m14 * (m23 * m41 - m21 * m43));
      const T_Scalar inverseM33 = determinantInverse * (m11 * (m22 * m44 - m24 * m42) + m12 * (m24 * m41 - m21 * m44) + m14 * (m21 * m42 - m22 * m41));
      const T_Scalar inverseM43 = determinantInverse * (m11 * (m23 * m42 - m22 * m43) + m12 * (m21 * m43 - m23 * m41) + m13 * (m22 * m41 - m21 * m42));

      const T_Scalar inverseM14 = determinantInverse * (m12 * (m24 * m33 - m23 * m34) + m13 * (m22 * m34 - m24 * m32) + m14 * (m23 * m32 - m22 * m33));
      const T_Scalar inverseM24 = determinantInverse * (m11 * (m23 * m34 - m24 * m33) + m13 * (m24 * m31 - m21 * m34) + m14 * (m21 * m33 - m23 * m31));
      const T_Scalar inverseM34 = determinantInverse * (m11 * (m24 * m32 - m22 * m34) + m12 * (m21 * m34 - m24 * m31) + m14 * (m22 * m31 - m21 * m32));
      const T_Scalar inverseM44 = determinantInverse * (m11 * (m22 * m33 - m23 * m32) + m12 * (m23 * m31 - m21 * m33) + m13 * (m21 * m32 - m22 * m31));

      // Compute AB = M^{-1}F
      AB1[n] = inverseM11 * f1 + inverseM12 * f2 + inverseM13 * f3 + inverseM14 * f4;
      AB2[n] = inverseM21 * f1 + inverseM22 * f2 + inverseM23 * f3 + inverseM24 * f4;
      AB3[n] = inverseM31 * f1 + inverseM32 * f2 + inverseM33 * f3 + inverseM34 * f4;
      AB4[n] = inverseM41 * f1 + inverseM42 * f2 + inverseM43 * f3 + inverseM44 * f4;
    }

#pragma omp for
    for(auto i = 0ULL; i < values.size(); ++i) {
      const scalar::Precision< T_Scalar > x = points[3 * i];
      const scalar::Precision< T_Scalar > y = points[3 * i + 1];
      const scalar::Precision< T_Scalar > r = std::sqrt(x * x + y * y);
      const scalar::Precision< T_Scalar > theta = std::atan2(y, x) + _thetaInc;
      const scalar::Precision< T_Scalar > invr = 1. / r;

      common::KahanSum< T_Scalar > sumP(0.);
      common::KahanSum< T_Scalar > sumQ(0.);
      for(auto n = 0U; n < _nbrOfterm; ++n) {
        sumP += (AB1[n] * drBesselH1< T_Scalar >(kP, r, n) + AB2[n] * drBesselH2< T_Scalar >(kP, r, n) - (n * invr) * AB3[n] * boost::math::cyl_hankel_1(n, kS * r) - (n * invr) * AB4[n] * boost::math::cyl_hankel_2(n, kS * r)) * std::sin(n * theta);
        sumQ += ((n * invr) * AB1[n] * boost::math::cyl_hankel_1(n, kP * r) + (n * invr) * AB2[n] * boost::math::cyl_hankel_2(n, kP * r) - AB3[n] * drBesselH1< T_Scalar >(kS, r, n) - AB4[n] * drBesselH2< T_Scalar >(kS, r, n)) * std::cos(n * theta);
      }

      values[i] << std::cos(theta) * sumP.sum() - std::sin(theta) * sumQ.sum(), std::sin(theta) * sumP.sum() + std::cos(theta) * sumQ.sum(), 0.;
    }
  }
#else

  template< class T_Scalar >
  void SoftSWavesScatteringByACylinderWithHOABC< T_Scalar >::operator()(function::OutputVector< T_Scalar, Degree::Degree1 > &values, const std::vector< scalar::Precision< T_Scalar >, numa::allocator< scalar::Precision< T_Scalar > > > &points) const
  {
    throw common::Exception("Cannot use SoftSWavesScatteringByACylinderWithHOABC without Boost");
  }
#endif // HAVE_BOOST

  INSTANTIATE_CLASS(SoftSWavesScatteringByACylinderWithHOABC, 1, TEMPLATE_ARGS(std::complex< double >))


} // namespace gmshfem::analytics::navier2D
