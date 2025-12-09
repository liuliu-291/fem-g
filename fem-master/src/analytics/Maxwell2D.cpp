// GmshFEM - Copyright (C) 2019-2022, A. Royer, E. Béchet, C. Geuzaine, Université de Liège
//
// See the LICENSE.txt file for license information. Please report all
// issues on https://gitlab.onelab.info/gmsh/fem/issues

#include "Maxwell2D.h"

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

namespace gmshfem::analytics::maxwell2D
{


  //***************************
  // ScatteringByAPECCylinderTE
  //***************************

  template< class T_Scalar >
  ScatteringByAPECCylinderTE< T_Scalar >::ScatteringByAPECCylinderTE(const ScatteringByAPECCylinderTE &other) :
    function::AnalyticalOperation< T_Scalar, Degree::Degree1 >(other), _helmholtz(other._helmholtz), _values()
  {
  }

  template< class T_Scalar >
  ScatteringByAPECCylinderTE< T_Scalar >::ScatteringByAPECCylinderTE(const scalar::Precision< T_Scalar > k, const scalar::Precision< T_Scalar > R0, const scalar::Precision< T_Scalar > x, const scalar::Precision< T_Scalar > y, const unsigned int nbrOfterm, const scalar::Precision< T_Scalar > thetaInc) :
    _helmholtz(k, R0, x, y, nbrOfterm, thetaInc), _values()
  {
  }

  template< class T_Scalar >
  ScatteringByAPECCylinderTE< T_Scalar >::~ScatteringByAPECCylinderTE()
  {
  }

#ifdef HAVE_BOOST
  template< class T_Scalar >
  void ScatteringByAPECCylinderTE< T_Scalar >::operator()(function::OutputVector< T_Scalar, Degree::Degree1 > &values, const std::vector< scalar::Precision< T_Scalar >, numa::allocator< scalar::Precision< T_Scalar > > > &points) const
  {
#pragma omp single
    _values.resize(values.size());

    _helmholtz(_values, points);

#pragma omp for
    for(auto i = 0ULL; i < values.size(); ++i) {
      values[i] = typename MathObject< T_Scalar, Degree::Degree1 >::Object(0., 0., _values[i]);
    }
  }
#else
  template< class T_Scalar >
  void ScatteringByAPECCylinderTE< T_Scalar >::operator()(function::OutputVector< T_Scalar, Degree::Degree1 > &values, const std::vector< scalar::Precision< T_Scalar >, numa::allocator< scalar::Precision< T_Scalar > > > &points) const
  {
    throw common::Exception("Cannot use ScatteringByAPECCylinderTE without Boost");
  }
#endif // HAVE_BOOST

  INSTANTIATE_CLASS(ScatteringByAPECCylinderTE, 2, TEMPLATE_ARGS(std::complex< double >, std::complex< float >))


  //***************************
  // ScatteringByAPECCylinderTM
  //***************************

  template< class T_Scalar >
  ScatteringByAPECCylinderTM< T_Scalar >::ScatteringByAPECCylinderTM(const ScatteringByAPECCylinderTM &other) :
    function::AnalyticalOperation< T_Scalar, Degree::Degree1 >(other), _k(other._k), _R0(other._R0), _x(other._x), _y(other._y), _nbrOfterm(other._nbrOfterm), _thetaInc(other._thetaInc)
  {
  }

  template< class T_Scalar >
  ScatteringByAPECCylinderTM< T_Scalar >::ScatteringByAPECCylinderTM(const scalar::Precision< T_Scalar > k, const scalar::Precision< T_Scalar > R0, const scalar::Precision< T_Scalar > x, const scalar::Precision< T_Scalar > y, const unsigned int nbrOfterm, const scalar::Precision< T_Scalar > thetaInc) :
    _k(k), _R0(R0), _x(x), _y(y), _nbrOfterm(nbrOfterm), _thetaInc(thetaInc)
  {
  }

  template< class T_Scalar >
  ScatteringByAPECCylinderTM< T_Scalar >::~ScatteringByAPECCylinderTM()
  {
  }

#ifdef HAVE_BOOST
  template< class T_Scalar >
  void ScatteringByAPECCylinderTM< T_Scalar >::operator()(function::OutputVector< T_Scalar, Degree::Degree1 > &values, const std::vector< scalar::Precision< T_Scalar >, numa::allocator< scalar::Precision< T_Scalar > > > &points) const
  {
    const scalar::Precision< T_Scalar > kR0 = _k * _R0;
    const T_Scalar powI[4] = {1., T_Scalar(0., 1.), -1., -T_Scalar(0., 1.)};

    std::vector< T_Scalar > precomputed(_nbrOfterm);
    for(auto n = 0U; n < _nbrOfterm; ++n) {
      const T_Scalar dHn1kR0 = (!n ? -T_Scalar(boost::math::cyl_hankel_1(1, kR0)) : T_Scalar(boost::math::cyl_hankel_1(n - 1, kR0)) - scalar::Precision< T_Scalar >(n) / kR0 * T_Scalar(boost::math::cyl_hankel_1(n, kR0)));
      const T_Scalar invdHn1kR0 = T_Scalar(1., 0.) / dHn1kR0;
      precomputed[n] = (!n ? scalar::Precision< T_Scalar >(1.) : scalar::Precision< T_Scalar >(2.)) * std::real(dHn1kR0) * invdHn1kR0 * (!n ? powI[3] : powI[(n - 1) % 4]) / _k;
    }

#pragma omp for
    for(auto i = 0ULL; i < values.size(); ++i) {
      const scalar::Precision< T_Scalar > x = points[3 * i] - _x;
      const scalar::Precision< T_Scalar > y = points[3 * i + 1] - _y;
      const scalar::Precision< T_Scalar > r = std::sqrt(x * x + y * y);
      const scalar::Precision< T_Scalar > th = std::atan2(y, x);
      const scalar::Precision< T_Scalar > theta = th - _thetaInc;
      const scalar::Precision< T_Scalar > kr = _k * r;
      const scalar::Precision< T_Scalar > invr = scalar::Precision< T_Scalar >(1.) / r;
      common::KahanSum< T_Scalar > sumR(0.); // radius (r)
      common::KahanSum< T_Scalar > sumT(0.); // angular (theta)
      for(auto n = 0U; n < _nbrOfterm; ++n) {
        sumR += n * invr * precomputed[n] * T_Scalar(boost::math::cyl_hankel_1(n, kr)) * std::sin(n * theta);
        sumT += _k * precomputed[n] * (!n ? -T_Scalar(boost::math::cyl_hankel_1(1, kr)) : T_Scalar(boost::math::cyl_hankel_1(n - 1, kr)) - scalar::Precision< T_Scalar >(n) / kr * T_Scalar(boost::math::cyl_hankel_1(n, kr))) * std::cos(n * theta);
      }
      const T_Scalar totR = -sumR.sum();
      const T_Scalar totT = -sumT.sum();
      values[i] = typename MathObject< T_Scalar, Degree::Degree1 >::Object(std::cos(th) * totR - std::sin(th) * totT, std::sin(th) * totR + std::cos(th) * totT, 0.);
    }
  }
#else
  template< class T_Scalar >
  void ScatteringByAPECCylinderTM< T_Scalar >::operator()(function::OutputVector< T_Scalar, Degree::Degree1 > &values, const std::vector< scalar::Precision< T_Scalar >, numa::allocator< scalar::Precision< T_Scalar > > > &points) const
  {
    throw common::Exception("Cannot use ScatteringByAPECCylinderTM without Boost");
  }
#endif // HAVE_BOOST

  INSTANTIATE_CLASS(ScatteringByAPECCylinderTM, 2, TEMPLATE_ARGS(std::complex< double >, std::complex< float >))


  //***************************
  // ScatteringByAPMCCylinderTE
  //***************************

  template< class T_Scalar >
  ScatteringByAPMCCylinderTE< T_Scalar >::ScatteringByAPMCCylinderTE(const ScatteringByAPMCCylinderTE &other) :
    function::AnalyticalOperation< T_Scalar, Degree::Degree1 >(other), _helmholtz(other._helmholtz), _values()
  {
  }

  template< class T_Scalar >
  ScatteringByAPMCCylinderTE< T_Scalar >::ScatteringByAPMCCylinderTE(const scalar::Precision< T_Scalar > k, const scalar::Precision< T_Scalar > R0, const scalar::Precision< T_Scalar > x, const scalar::Precision< T_Scalar > y, const unsigned int nbrOfterm, const scalar::Precision< T_Scalar > thetaInc) :
    _helmholtz(k, R0, x, y, nbrOfterm, thetaInc), _values()
  {
  }

  template< class T_Scalar >
  ScatteringByAPMCCylinderTE< T_Scalar >::~ScatteringByAPMCCylinderTE()
  {
  }

#ifdef HAVE_BOOST
  template< class T_Scalar >
  void ScatteringByAPMCCylinderTE< T_Scalar >::operator()(function::OutputVector< T_Scalar, Degree::Degree1 > &values, const std::vector< scalar::Precision< T_Scalar >, numa::allocator< scalar::Precision< T_Scalar > > > &points) const
  {
#pragma omp single
    _values.resize(values.size());

    _helmholtz(_values, points);

#pragma omp for
    for(auto i = 0ULL; i < values.size(); ++i) {
      values[i] = typename MathObject< T_Scalar, Degree::Degree1 >::Object(0., 0., _values[i]);
    }
  }
#else
  template< class T_Scalar >
  void ScatteringByAPMCCylinderTE< T_Scalar >::operator()(function::OutputVector< T_Scalar, Degree::Degree1 > &values, const std::vector< scalar::Precision< T_Scalar >, numa::allocator< scalar::Precision< T_Scalar > > > &points) const
  {
    throw common::Exception("Cannot use ScatteringByAPMCCylinderTE without Boost");
  }
#endif // HAVE_BOOST

  INSTANTIATE_CLASS(ScatteringByAPMCCylinderTE, 2, TEMPLATE_ARGS(std::complex< double >, std::complex< float >))


  //******************************
  // dr_ScatteringByAPMCCylinderTE
  //******************************

  template< class T_Scalar >
  dr_ScatteringByAPMCCylinderTE< T_Scalar >::dr_ScatteringByAPMCCylinderTE(const dr_ScatteringByAPMCCylinderTE &other) :
    function::AnalyticalOperation< T_Scalar, Degree::Degree1 >(other), _helmholtz(other._helmholtz), _values()
  {
  }

  template< class T_Scalar >
  dr_ScatteringByAPMCCylinderTE< T_Scalar >::dr_ScatteringByAPMCCylinderTE(const scalar::Precision< T_Scalar > k, const scalar::Precision< T_Scalar > R0, const scalar::Precision< T_Scalar > x, const scalar::Precision< T_Scalar > y, const unsigned int nbrOfterm, const scalar::Precision< T_Scalar > thetaInc) :
    _helmholtz(k, R0, x, y, nbrOfterm, thetaInc), _values()
  {
  }

  template< class T_Scalar >
  dr_ScatteringByAPMCCylinderTE< T_Scalar >::~dr_ScatteringByAPMCCylinderTE()
  {
  }

#ifdef HAVE_BOOST
  template< class T_Scalar >
  void dr_ScatteringByAPMCCylinderTE< T_Scalar >::operator()(function::OutputVector< T_Scalar, Degree::Degree1 > &values, const std::vector< scalar::Precision< T_Scalar >, numa::allocator< scalar::Precision< T_Scalar > > > &points) const
  {
#pragma omp single
    _values.resize(values.size());

    _helmholtz(_values, points);

#pragma omp for
    for(auto i = 0ULL; i < values.size(); ++i) {
      values[i] = typename MathObject< T_Scalar, Degree::Degree1 >::Object(0., 0., -_values[i]);
    }
  }
#else
  template< class T_Scalar >
  void dr_ScatteringByAPMCCylinderTE< T_Scalar >::operator()(function::OutputVector< T_Scalar, Degree::Degree1 > &values, const std::vector< scalar::Precision< T_Scalar >, numa::allocator< scalar::Precision< T_Scalar > > > &points) const
  {
    throw common::Exception("Cannot use dr_ScatteringByAPMCCylinderTE without Boost");
  }
#endif // HAVE_BOOST

  INSTANTIATE_CLASS(dr_ScatteringByAPMCCylinderTE, 2, TEMPLATE_ARGS(std::complex< double >, std::complex< float >))


} // namespace gmshfem::analytics::maxwell2D
