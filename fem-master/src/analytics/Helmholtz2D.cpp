// GmshFEM - Copyright (C) 2019-2022, A. Royer, E. Béchet, C. Geuzaine, Université de Liège
//
// See the LICENSE.txt file for license information. Please report all
// issues on https://gitlab.onelab.info/gmsh/fem/issues

#include "Helmholtz2D.h"

#include "Exception.h"
#include "KahanSum.h"
#include "Message.h"
#include "OmpInterface.h"
#include "gmshfemDefines.h"
#include "instantiate.h"

#include <cmath>
#include <complex>

#ifdef HAVE_BOOST
#include <boost/math/special_functions/bessel.hpp>
#include <boost/math/special_functions/hankel.hpp>
#endif // HAVE_BOOST

#ifdef HAVE_BESSEL
#include "Bessel.h"
#endif // HAVE_BOOST

namespace gmshfem::analytics::helmholtz2D
{


  //**************************
  // ScatteringByASoftCylinder
  //**************************

  template< class T_Scalar >
  ScatteringByASoftCylinder< T_Scalar >::ScatteringByASoftCylinder(const ScatteringByASoftCylinder &other) :
    function::AnalyticalOperation< T_Scalar, Degree::Degree0 >(other), _k(other._k), _R0(other._R0), _x(other._x), _y(other._y), _nbrOfterm(other._nbrOfterm), _thetaInc(other._thetaInc)
  {
  }

  template< class T_Scalar >
  ScatteringByASoftCylinder< T_Scalar >::ScatteringByASoftCylinder(const scalar::Precision< T_Scalar > k, const scalar::Precision< T_Scalar > R0, const scalar::Precision< T_Scalar > x, const scalar::Precision< T_Scalar > y, const unsigned int nbrOfterm, const scalar::Precision< T_Scalar > thetaInc) :
    _k(k), _R0(R0), _x(x), _y(y), _nbrOfterm(nbrOfterm), _thetaInc(thetaInc)
  {
  }

  template< class T_Scalar >
  ScatteringByASoftCylinder< T_Scalar >::~ScatteringByASoftCylinder()
  {
  }

#ifdef HAVE_BOOST
  template< class T_Scalar >
  void ScatteringByASoftCylinder< T_Scalar >::operator()(function::OutputVector< T_Scalar, Degree::Degree0 > &values, const std::vector< scalar::Precision< T_Scalar >, numa::allocator< scalar::Precision< T_Scalar > > > &points) const
  {
    const scalar::Precision< T_Scalar > kR0 = _k * _R0;
    const T_Scalar powI[4] = {1., T_Scalar(0., 1.), -1., -T_Scalar(0., 1.)};

    std::vector< T_Scalar > precomputed(_nbrOfterm);
    for(auto n = 0U; n < _nbrOfterm; ++n) {
      const scalar::Precision< T_Scalar > JnkR0 = boost::math::cyl_bessel_j(n, kR0);
      const T_Scalar invHn1kR0 = T_Scalar(1., 0.) / T_Scalar(JnkR0, boost::math::cyl_neumann(n, kR0));
      precomputed[n] = -(!n ? scalar::Precision< T_Scalar >(1.) : scalar::Precision< T_Scalar >(2.)) * JnkR0 * invHn1kR0 * powI[n % 4];
    }

#pragma omp for
    for(auto i = 0ULL; i < values.size(); ++i) {
      const scalar::Precision< T_Scalar > x = points[3 * i] - _x;
      const scalar::Precision< T_Scalar > y = points[3 * i + 1] - _y;
      const scalar::Precision< T_Scalar > r = std::sqrt(x * x + y * y);
      const scalar::Precision< T_Scalar > theta = std::atan2(y, x) + _thetaInc;
      const scalar::Precision< T_Scalar > kr = _k * r;
      common::KahanSum< T_Scalar > sum(0.);
      for(auto n = 0U; n < _nbrOfterm; ++n) {
        sum += std::cos(n * theta) * T_Scalar(boost::math::cyl_hankel_1(n, kr)) * precomputed[n];
      }
      values[i] = sum.sum();
    }
  }
#else
  template< class T_Scalar >
  void ScatteringByASoftCylinder< T_Scalar >::operator()(function::OutputVector< T_Scalar, Degree::Degree0 > &values, const std::vector< scalar::Precision< T_Scalar >, numa::allocator< scalar::Precision< T_Scalar > > > &points) const
  {
    throw common::Exception("Cannot use ScatteringByASoftCylinder without Boost");
  }
#endif // HAVE_BOOST

  INSTANTIATE_CLASS(ScatteringByASoftCylinder, 2, TEMPLATE_ARGS(std::complex< double >, std::complex< float >))


  //**************************
  // ScatteringByAHardCylinder
  //**************************

  template< class T_Scalar >
  ScatteringByAHardCylinder< T_Scalar >::ScatteringByAHardCylinder(const ScatteringByAHardCylinder &other) :
    function::AnalyticalOperation< T_Scalar, Degree::Degree0 >(other), _k(other._k), _R0(other._R0), _x(other._x), _y(other._y), _nbrOfterm(other._nbrOfterm), _thetaInc(other._thetaInc)
  {
  }

  template< class T_Scalar >
  ScatteringByAHardCylinder< T_Scalar >::ScatteringByAHardCylinder(const scalar::Precision< T_Scalar > k, const scalar::Precision< T_Scalar > R0, const scalar::Precision< T_Scalar > x, const scalar::Precision< T_Scalar > y, const unsigned int nbrOfterm, const scalar::Precision< T_Scalar > thetaInc) :
    _k(k), _R0(R0), _x(x), _y(y), _nbrOfterm(nbrOfterm), _thetaInc(thetaInc)
  {
  }

  template< class T_Scalar >
  ScatteringByAHardCylinder< T_Scalar >::~ScatteringByAHardCylinder()
  {
  }

#ifdef HAVE_BOOST
  template< class T_Scalar >
  void ScatteringByAHardCylinder< T_Scalar >::operator()(function::OutputVector< T_Scalar, Degree::Degree0 > &values, const std::vector< scalar::Precision< T_Scalar >, numa::allocator< scalar::Precision< T_Scalar > > > &points) const
  {
    const scalar::Precision< T_Scalar > kR0 = _k * _R0;
    const T_Scalar powI[4] = {1., T_Scalar(0., 1.), -1., -T_Scalar(0., 1.)};

    std::vector< T_Scalar > precomputed(_nbrOfterm);
    for(auto n = 0U; n < _nbrOfterm; ++n) {
      const T_Scalar dHn1kR0 = (!n ? -T_Scalar(boost::math::cyl_hankel_1(1, kR0)) : T_Scalar(boost::math::cyl_hankel_1(n - 1, kR0)) - scalar::Precision< T_Scalar >(n) / kR0 * T_Scalar(boost::math::cyl_hankel_1(n, kR0)));
      const T_Scalar invdHn1kR0 = T_Scalar(1., 0.) / dHn1kR0;
      precomputed[n] = -(!n ? scalar::Precision< T_Scalar >(1.) : scalar::Precision< T_Scalar >(2.)) * std::real(dHn1kR0) * invdHn1kR0 * powI[n % 4];
    }

#pragma omp for
    for(auto i = 0ULL; i < values.size(); ++i) {
      const scalar::Precision< T_Scalar > x = points[3 * i] - _x;
      const scalar::Precision< T_Scalar > y = points[3 * i + 1] - _y;
      const scalar::Precision< T_Scalar > r = std::sqrt(x * x + y * y);
      const scalar::Precision< T_Scalar > theta = std::atan2(y, x) - _thetaInc;
      const scalar::Precision< T_Scalar > kr = _k * r;
      common::KahanSum< T_Scalar > sum(0.);
      for(auto n = 0U; n < _nbrOfterm; ++n) {
        sum += std::cos(n * theta) * T_Scalar(boost::math::cyl_hankel_1(n, kr)) * precomputed[n];
      }
      values[i] = sum.sum();
    }
  }
#else
  template< class T_Scalar >
  void ScatteringByAHardCylinder< T_Scalar >::operator()(function::OutputVector< T_Scalar, Degree::Degree0 > &values, const std::vector< scalar::Precision< T_Scalar >, numa::allocator< scalar::Precision< T_Scalar > > > &points) const
  {
    throw common::Exception("Cannot use ScatteringByAHardCylinder without Boost");
  }
#endif // HAVE_BOOST

  INSTANTIATE_CLASS(ScatteringByAHardCylinder, 2, TEMPLATE_ARGS(std::complex< double >, std::complex< float >))


  //*****************************
  // dr_ScatteringByAHardCylinder
  //*****************************

  template< class T_Scalar >
  dr_ScatteringByAHardCylinder< T_Scalar >::dr_ScatteringByAHardCylinder(const dr_ScatteringByAHardCylinder &other) :
    function::AnalyticalOperation< T_Scalar, Degree::Degree0 >(other), _k(other._k), _R0(other._R0), _x(other._x), _y(other._y), _nbrOfterm(other._nbrOfterm), _thetaInc(other._thetaInc)
  {
  }

  template< class T_Scalar >
  dr_ScatteringByAHardCylinder< T_Scalar >::dr_ScatteringByAHardCylinder(const scalar::Precision< T_Scalar > k, const scalar::Precision< T_Scalar > R0, const scalar::Precision< T_Scalar > x, const scalar::Precision< T_Scalar > y, const unsigned int nbrOfterm, const scalar::Precision< T_Scalar > thetaInc) :
    _k(k), _R0(R0), _x(x), _y(y), _nbrOfterm(nbrOfterm), _thetaInc(thetaInc)
  {
  }

  template< class T_Scalar >
  dr_ScatteringByAHardCylinder< T_Scalar >::~dr_ScatteringByAHardCylinder()
  {
  }

#ifdef HAVE_BOOST
  template< class T_Scalar >
  void dr_ScatteringByAHardCylinder< T_Scalar >::operator()(function::OutputVector< T_Scalar, Degree::Degree0 > &values, const std::vector< scalar::Precision< T_Scalar >, numa::allocator< scalar::Precision< T_Scalar > > > &points) const
  {
    const scalar::Precision< T_Scalar > kR0 = _k * _R0;
    const T_Scalar powI[4] = {1., T_Scalar(0., 1.), -1., -T_Scalar(0., 1.)};

    std::vector< T_Scalar > precomputed(_nbrOfterm);
    for(auto n = 0U; n < _nbrOfterm; ++n) {
      const T_Scalar dHn1kR0 = (!n ? -T_Scalar(boost::math::cyl_hankel_1(1, kR0)) : T_Scalar(boost::math::cyl_hankel_1(n - 1, kR0)) - scalar::Precision< T_Scalar >(n) / kR0 * T_Scalar(boost::math::cyl_hankel_1(n, kR0)));
      const T_Scalar invdHn1kR0 = T_Scalar(1., 0.) / dHn1kR0;
      precomputed[n] = -(!n ? scalar::Precision< T_Scalar >(1.) : scalar::Precision< T_Scalar >(2.)) * std::real(dHn1kR0) * invdHn1kR0 * powI[n % 4] * _k;
    }

#pragma omp for
    for(auto i = 0ULL; i < values.size(); ++i) {
      const scalar::Precision< T_Scalar > x = points[3 * i] - _x;
      const scalar::Precision< T_Scalar > y = points[3 * i + 1] - _y;
      const scalar::Precision< T_Scalar > r = std::sqrt(x * x + y * y);
      const scalar::Precision< T_Scalar > theta = std::atan2(y, x) - _thetaInc;
      const scalar::Precision< T_Scalar > kr = _k * r;
      common::KahanSum< T_Scalar > sum(0.);
      for(auto n = 0U; n < _nbrOfterm; ++n) {
        const T_Scalar dHn1kr = (!n ? -T_Scalar(boost::math::cyl_hankel_1(1, kr)) : T_Scalar(boost::math::cyl_hankel_1(n - 1, kr)) - scalar::Precision< T_Scalar >(n) / kr * T_Scalar(boost::math::cyl_hankel_1(n, kr)));
        sum += std::cos(n * theta) * dHn1kr * precomputed[n];
      }
      values[i] = sum.sum();
    }
  }
#else
  template< class T_Scalar >
  void dr_ScatteringByAHardCylinder< T_Scalar >::operator()(function::OutputVector< T_Scalar, Degree::Degree0 > &values, const std::vector< scalar::Precision< T_Scalar >, numa::allocator< scalar::Precision< T_Scalar > > > &points) const
  {
    throw common::Exception("Cannot use dr_ScatteringByAHardCylinder without Boost");
  }
#endif // HAVE_BOOST

  INSTANTIATE_CLASS(dr_ScatteringByAHardCylinder, 2, TEMPLATE_ARGS(std::complex< double >, std::complex< float >))


  //**************************
  // DuctModeSolution
  //**************************

  template< class T_Scalar >
  DuctModeSolution< T_Scalar >::DuctModeSolution(const DuctModeSolution &other) :
    function::AnalyticalOperation< T_Scalar, Degree::Degree0 >(other), _k(other._k), _M(other._M), _H(other._H), _x(other._x), _y(other._y), _n(other._n), _BC(other._BC)
  {
  }

  template< class T_Scalar >
  DuctModeSolution< T_Scalar >::DuctModeSolution(const scalar::Precision< T_Scalar > k, const scalar::Precision< T_Scalar > M, const scalar::Precision< T_Scalar > H, const scalar::Precision< T_Scalar > x, const scalar::Precision< T_Scalar > y, const unsigned int n, const bool BC) :
    _k(k), _M(M), _H(H), _x(x), _y(y), _n(n), _BC(BC)
  {
  }

  template< class T_Scalar >
  DuctModeSolution< T_Scalar >::~DuctModeSolution()
  {
  }

  template< class T_Scalar >
  void DuctModeSolution< T_Scalar >::operator()(function::OutputVector< T_Scalar, Degree::Degree0 > &values, const std::vector< scalar::Precision< T_Scalar >, numa::allocator< scalar::Precision< T_Scalar > > > &points) const
  {
    const T_Scalar im = T_Scalar(0., 1.);
    const double pi = 3.14159265358979323846264338327950288;
    const scalar::Precision< T_Scalar > ky = (_n * pi / _H);
    const scalar::Precision< T_Scalar > ArgSqrt = _k * _k - (1. - _M * _M) * ky * ky;
    T_Scalar kx;
    if(std::real(ArgSqrt) >= 0.) {
      kx = T_Scalar((-_M * _k + std::sqrt(ArgSqrt)) / (1. - _M * _M), 0.);
    }
    else {
      kx = T_Scalar(-_M * _k / (1. - _M * _M), -std::sqrt(std::abs(ArgSqrt)) / (1. - _M * _M));
    }
#pragma omp for
    for(auto i = 0ULL; i < values.size(); ++i) {
      const scalar::Precision< T_Scalar > x = points[3 * i] - _x;
      const scalar::Precision< T_Scalar > y = points[3 * i + 1] - _y;
      if(_BC == 0) {
        values[i] = std::cos(ky * y) * std::exp(-im * kx * x);
      }
      else {
        values[i] = std::sin(ky * y) * std::exp(-im * kx * x);
      }
    }
  }

  INSTANTIATE_CLASS(DuctModeSolution, 2, TEMPLATE_ARGS(std::complex< double >, std::complex< float >))

  //**************************
  // DuctModeSolutionAiry
  //**************************

  template< class T_Scalar >
  DuctModeSolutionAiry< T_Scalar >::DuctModeSolutionAiry(const DuctModeSolutionAiry &other) :
    function::AnalyticalOperation< T_Scalar, Degree::Degree0 >(other), _k(other._k), _H(other._H), _a(other._a), _b(other._b), _x(other._x), _y(other._y), _n(other._n), _BC(other._BC)
  {
  }

  template< class T_Scalar >
  DuctModeSolutionAiry< T_Scalar >::DuctModeSolutionAiry(const scalar::Precision< T_Scalar > k, const scalar::Precision< T_Scalar > H, const scalar::Precision< T_Scalar > a, const scalar::Precision< T_Scalar > b, const scalar::Precision< T_Scalar > x, const scalar::Precision< T_Scalar > y, const unsigned int n, const bool BC) :
    _k(k), _H(H), _a(a), _b(b), _x(x), _y(y), _n(n), _BC(BC)
  {
  }

  template< class T_Scalar >
  DuctModeSolutionAiry< T_Scalar >::~DuctModeSolutionAiry()
  {
  }

#ifdef HAVE_BESSEL
  template< class T_Scalar >
  void DuctModeSolutionAiry< T_Scalar >::operator()(function::OutputVector< T_Scalar, Degree::Degree0 > &values, const std::vector< scalar::Precision< T_Scalar >, numa::allocator< scalar::Precision< T_Scalar > > > &points) const
  {
    const double pi = 3.14159265358979323846264338327950288;
    double valr, vali;
    const scalar::Precision< T_Scalar > ky = (_n * pi / _H);
    const scalar::Precision< T_Scalar > den = std::pow(_a * _k * _k, 2. / 3.);
#pragma omp for
    for(auto i = 0ULL; i < values.size(); ++i) {
      const scalar::Precision< T_Scalar > x = points[3 * i] - _x;
      const scalar::Precision< T_Scalar > y = points[3 * i + 1] - _y;
      // solution for exp(i*w*t) convention !
      double z = (ky * ky - _k * _k * (_a * x + _b)) / den;
      int err = AiryComplex(std::cos(2. * pi / 3.) * z, -std::sin(2. * pi / 3.) * z, 0, &valr, &vali);
      if(err != 0) {
        msg::warning << "Issue with complex airy function, error output" << err << msg::endl;
      }

      if(_BC == 0) {
        values[i] = std::cos(ky * y) * T_Scalar(valr, vali);
      }
      else {
        values[i] = std::sin(ky * y) * T_Scalar(valr, vali);
      }
    }
  }
#else
  template< class T_Scalar >
  void DuctModeSolutionAiry< T_Scalar >::operator()(function::OutputVector< T_Scalar, Degree::Degree0 > &values, const std::vector< scalar::Precision< T_Scalar >, numa::allocator< scalar::Precision< T_Scalar > > > &points) const
  {
    throw common::Exception("Cannot use DuctModeSolutionAiry without Bessel functions");
  }
#endif // HAVE_BOOST

  INSTANTIATE_CLASS(DuctModeSolutionAiry, 2, TEMPLATE_ARGS(std::complex< double >, std::complex< float >))


  //**************************
  // DuctModeSolutionMulti
  //**************************

  template< class T_Scalar >
  DuctModeSolutionMulti< T_Scalar >::DuctModeSolutionMulti(const DuctModeSolutionMulti &other) :
    function::AnalyticalOperation< T_Scalar, Degree::Degree0 >(other), _k(other._k), _M(other._M), _H(other._H), _A(other._A), _n(other._n), _x(other._x), _y(other._y), _BC(other._BC)
  {
  }

  template< class T_Scalar >
  DuctModeSolutionMulti< T_Scalar >::DuctModeSolutionMulti(const scalar::Precision< T_Scalar > k, const scalar::Precision< T_Scalar > M, const scalar::Precision< T_Scalar > H, const std::vector< scalar::Precision< T_Scalar > > A, const std::vector< unsigned int > n, const scalar::Precision< T_Scalar > x, const scalar::Precision< T_Scalar > y, const bool BC) :
    _k(k), _M(M), _H(H), _A(A), _n(n), _x(x), _y(y), _BC(BC)
  {
  }

  template< class T_Scalar >
  DuctModeSolutionMulti< T_Scalar >::~DuctModeSolutionMulti()
  {
  }

  template< class T_Scalar >
  void DuctModeSolutionMulti< T_Scalar >::operator()(function::OutputVector< T_Scalar, Degree::Degree0 > &values, const std::vector< scalar::Precision< T_Scalar >, numa::allocator< scalar::Precision< T_Scalar > > > &points) const
  {
    const T_Scalar im = T_Scalar(0., 1.);
    const double pi = 3.14159265358979323846264338327950288;

    unsigned int N_modes = _n.size();
    if(N_modes != _A.size()) {
      msg::error << "The amplitudes and indices do not have the same size" << msg::endl;
    }

    std::vector< T_Scalar > kx(N_modes);
    std::vector< scalar::Precision< T_Scalar > > ky(N_modes);
    std::vector< scalar::Precision< T_Scalar > > ArgSqrt(N_modes);
    // compute each wavenumbers via the dispersion relation
    for(auto j = 0U; j < N_modes; ++j) {
      ky[j] = _n[j] * pi / _H;
      ArgSqrt[j] = _k * _k - (1. - _M * _M) * (ky[j] * ky[j]);
      if(std::real(ArgSqrt[j]) >= 0.) {
        kx[j] = T_Scalar((-_M * _k + std::sqrt(ArgSqrt[j])) / (1. - _M * _M), 0.);
      }
      else {
        kx[j] = T_Scalar(-_M * _k / (1. - _M * _M), -std::sqrt(std::abs(ArgSqrt[j])) / (1. - _M * _M));
      }
    }

#pragma omp for
    for(auto i = 0ULL; i < values.size(); ++i) {
      const scalar::Precision< T_Scalar > x = points[3 * i] - _x;
      const scalar::Precision< T_Scalar > y = points[3 * i + 1] - _y;
      common::KahanSum< T_Scalar > sum(0.);
      for(auto j = 0U; j < N_modes; ++j) {
        if(_BC == 0) {
          sum += _A[j] * std::cos(ky[j] * y) * std::exp(-im * kx[j] * x);
        }
        else {
          sum += _A[j] * std::sin(ky[j] * y) * std::exp(-im * kx[j] * x);
        }
      }
      values[i] = sum.sum();
    }
  }

  INSTANTIATE_CLASS(DuctModeSolutionMulti, 2, TEMPLATE_ARGS(std::complex< double >, std::complex< float >))


  //**************************
  // MonopoleFreeField
  //**************************

  template< class T_Scalar >
  MonopoleFreeField< T_Scalar >::MonopoleFreeField(const MonopoleFreeField &other) :
    function::AnalyticalOperation< T_Scalar, Degree::Degree0 >(other), _k(other._k), _M(other._M), _theta(other._theta), _A(other._A), _xs(other._xs), _ys(other._ys), _x(other._x), _y(other._y)
  {
  }

  template< class T_Scalar >
  MonopoleFreeField< T_Scalar >::MonopoleFreeField(const scalar::Precision< T_Scalar > k, const scalar::Precision< T_Scalar > M, const scalar::Precision< T_Scalar > theta, const scalar::Precision< T_Scalar > A, const scalar::Precision< T_Scalar > xs, const scalar::Precision< T_Scalar > ys, const scalar::Precision< T_Scalar > x, const scalar::Precision< T_Scalar > y) :
    _k(k), _M(M), _theta(theta), _A(A), _xs(xs), _ys(ys), _x(x), _y(y)
  {
  }

  template< class T_Scalar >
  MonopoleFreeField< T_Scalar >::~MonopoleFreeField()
  {
  }

#ifdef HAVE_BOOST
  template< class T_Scalar >
  void MonopoleFreeField< T_Scalar >::operator()(function::OutputVector< T_Scalar, Degree::Degree0 > &values, const std::vector< scalar::Precision< T_Scalar >, numa::allocator< scalar::Precision< T_Scalar > > > &points) const
  {
    const T_Scalar im = T_Scalar(0., 1.);
    const scalar::Precision< T_Scalar > Mx = _M * std::cos(_theta);
    const scalar::Precision< T_Scalar > My = _M * std::sin(_theta);
    const scalar::Precision< T_Scalar > beta = std::sqrt(1 - _M * _M);
    const scalar::Precision< T_Scalar > khat = _k / beta;
    const scalar::Precision< T_Scalar > Alphax = 1 + Mx * Mx / (beta * (1 + beta));
    const scalar::Precision< T_Scalar > Alphay = 1 + My * My / (beta * (1 + beta));
    const scalar::Precision< T_Scalar > K = Mx * My / (beta * (1 + beta));

#pragma omp for
    for(auto i = 0ULL; i < values.size(); ++i) {
      const scalar::Precision< T_Scalar > xhat = Alphax * ((points[3 * i] - _x) - _xs) + K * ((points[3 * i + 1] - _y) - _ys);
      const scalar::Precision< T_Scalar > yhat = K * ((points[3 * i] - _x) - _xs) + Alphay * ((points[3 * i + 1] - _y) - _ys);
      const scalar::Precision< T_Scalar > rhat = std::sqrt(xhat * xhat + yhat * yhat);
      const scalar::Precision< T_Scalar > krhat = khat * rhat;
      // add 1e-16 to avoid evaluation of cyl_hankel_2 at the source node (xs,ys)
      // FIXME do not evaluate the field at the source node when plotting the solution in gmsh
      values[i] = -im * _A / (4 * beta) * T_Scalar(boost::math::cyl_hankel_2(0, krhat + 1e-16)) * std::exp(im * khat * Mx * xhat) * std::exp(im * khat * My * yhat);
    }
  }
#else
  template< class T_Scalar >
  void MonopoleFreeField< T_Scalar >::operator()(function::OutputVector< T_Scalar, Degree::Degree0 > &values, const std::vector< scalar::Precision< T_Scalar >, numa::allocator< scalar::Precision< T_Scalar > > > &points) const
  {
    throw common::Exception("Cannot use MonopoleFreeField without Boost");
  }
#endif // HAVE_BOOST

  INSTANTIATE_CLASS(MonopoleFreeField, 2, TEMPLATE_ARGS(std::complex< double >, std::complex< float >))


} // namespace gmshfem::analytics::helmholtz2D
