// GmshFEM - Copyright (C) 2019-2022, A. Royer, E. Béchet, C. Geuzaine, Université de Liège
//
// See the LICENSE.txt file for license information. Please report all
// issues on https://gitlab.onelab.info/gmsh/fem/issues

#include "Helmholtz3D.h"

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
#include <boost/math/special_functions/legendre.hpp>
#endif // HAVE_BOOST

namespace gmshfem::analytics::helmholtz3D
{


  //************************
  // ScatteringByASoftSphere
  //************************

  template< class T_Scalar >
  ScatteringByASoftSphere< T_Scalar >::ScatteringByASoftSphere(const scalar::Precision< T_Scalar > k, const scalar::Precision< T_Scalar > R0, const scalar::Precision< T_Scalar > x, const scalar::Precision< T_Scalar > y, const scalar::Precision< T_Scalar > z, const unsigned int nbrOfterm, const scalar::Precision< T_Scalar > thetaInc, const scalar::Precision< T_Scalar > phiInc) :
    _k(k), _R0(R0), _x(x), _y(y), _z(z), _nbrOfterm(nbrOfterm), _thetaInc(thetaInc), _phiInc(phiInc)
  {
  }

  template< class T_Scalar >
  ScatteringByASoftSphere< T_Scalar >::ScatteringByASoftSphere(const ScatteringByASoftSphere &other) :
    function::AnalyticalOperation< T_Scalar, Degree::Degree0 >(other), _k(other._k), _R0(other._R0), _x(other._x), _y(other._y), _z(other._z), _nbrOfterm(other._nbrOfterm), _thetaInc(other._thetaInc), _phiInc(other._phiInc)
  {
  }

  template< class T_Scalar >
  ScatteringByASoftSphere< T_Scalar >::~ScatteringByASoftSphere()
  {
  }

#ifdef HAVE_BOOST
  template< class T_Scalar >
  void ScatteringByASoftSphere< T_Scalar >::operator()(function::OutputVector< T_Scalar, Degree::Degree0 > &values, const std::vector< scalar::Precision< T_Scalar >, numa::allocator< scalar::Precision< T_Scalar > > > &points) const
  {
    const scalar::Precision< T_Scalar > kR0 = _k * _R0;
    const T_Scalar powI[4] = {1., T_Scalar(0., 1.), -1., -T_Scalar(0., 1.)};

    std::vector< T_Scalar > precomputed(_nbrOfterm);
    for(auto n = 0U; n < _nbrOfterm; ++n) {
      const scalar::Precision< T_Scalar > jnkR0 = boost::math::sph_bessel(n, kR0);
      const T_Scalar invhn1kR0 = T_Scalar(1., 0.) / T_Scalar(jnkR0, boost::math::sph_neumann(n, kR0));
      precomputed[n] = -powI[n % 4] * T_Scalar((2 * n + 1), 0.) * jnkR0 * invhn1kR0;
    }

    const scalar::Precision< T_Scalar > sinPhiInc = std::sin(_phiInc);
    const scalar::Precision< T_Scalar > cosPhiInc = std::cos(_phiInc);

    const scalar::Precision< T_Scalar > sinThetaInc = std::sin(_thetaInc);
    const scalar::Precision< T_Scalar > cosThetaInc = std::cos(_thetaInc);

#pragma omp for
    for(auto i = 0ULL; i < values.size(); ++i) {
      const scalar::Precision< T_Scalar > x = points[3 * i] - _x;
      const scalar::Precision< T_Scalar > y = points[3 * i + 1] - _y;
      const scalar::Precision< T_Scalar > z = points[3 * i + 2] - _z;
      const scalar::Precision< T_Scalar > r = std::sqrt(x * x + y * y + z * z);

      const scalar::Precision< T_Scalar > dx = cosThetaInc * sinPhiInc;
      const scalar::Precision< T_Scalar > dy = sinThetaInc * sinPhiInc;
      const scalar::Precision< T_Scalar > dz = cosPhiInc;
      const scalar::Precision< T_Scalar > dr = std::sqrt(dx * dx + dy * dy + dz * dz);

      const scalar::Precision< T_Scalar > kx = (x * dx + y * dy + z * dz) / (r * dr);

      const scalar::Precision< T_Scalar > kr = _k * r;
      common::KahanSum< T_Scalar > sum(0.);
      for(auto n = 0U; n < _nbrOfterm; ++n) {
        sum += T_Scalar(boost::math::legendre_p(n, kx)) * T_Scalar(boost::math::sph_hankel_1(n, kr)) * precomputed[n];
      }
      values[i] = sum.sum();
    }
  }
#else
  template< class T_Scalar >
  void ScatteringByASoftSphere< T_Scalar >::operator()(function::OutputVector< T_Scalar, Degree::Degree0 > &values, const std::vector< scalar::Precision< T_Scalar >, numa::allocator< scalar::Precision< T_Scalar > > > &points) const
  {
    throw common::Exception("Cannot use ScatteringByASoftSphere without Boost");
  }
#endif // HAVE_BOOST

  INSTANTIATE_CLASS(ScatteringByASoftSphere, 2, TEMPLATE_ARGS(std::complex< double >, std::complex< float >))


  //************************
  // ScatteringByAHardSphere
  //************************

  template< class T_Scalar >
  ScatteringByAHardSphere< T_Scalar >::ScatteringByAHardSphere(const scalar::Precision< T_Scalar > k, const scalar::Precision< T_Scalar > R0, const scalar::Precision< T_Scalar > x, const scalar::Precision< T_Scalar > y, const scalar::Precision< T_Scalar > z, const unsigned int nbrOfterm, const scalar::Precision< T_Scalar > thetaInc, const scalar::Precision< T_Scalar > phiInc) :
    _k(k), _R0(R0), _x(x), _y(y), _z(z), _nbrOfterm(nbrOfterm), _thetaInc(thetaInc), _phiInc(phiInc)
  {
  }

  template< class T_Scalar >
  ScatteringByAHardSphere< T_Scalar >::ScatteringByAHardSphere(const ScatteringByAHardSphere &other) :
    function::AnalyticalOperation< T_Scalar, Degree::Degree0 >(other), _k(other._k), _R0(other._R0), _x(other._x), _y(other._y), _z(other._z), _nbrOfterm(other._nbrOfterm), _thetaInc(other._thetaInc), _phiInc(other._phiInc)
  {
  }

  template< class T_Scalar >
  ScatteringByAHardSphere< T_Scalar >::~ScatteringByAHardSphere()
  {
  }

#ifdef HAVE_BOOST
  template< class T_Scalar >
  void ScatteringByAHardSphere< T_Scalar >::operator()(function::OutputVector< T_Scalar, Degree::Degree0 > &values, const std::vector< scalar::Precision< T_Scalar >, numa::allocator< scalar::Precision< T_Scalar > > > &points) const
  {
    const scalar::Precision< T_Scalar > kR0 = _k * _R0;
    const T_Scalar powI[4] = {1., T_Scalar(0., 1.), -1., -T_Scalar(0., 1.)};

    std::vector< T_Scalar > precomputed(_nbrOfterm);
    for(auto n = 0U; n < _nbrOfterm; ++n) {
      const T_Scalar dhn1kR0 = (!n ? -T_Scalar(boost::math::sph_hankel_1(1, kR0)) : T_Scalar(boost::math::sph_hankel_1(n - 1, kR0)) - scalar::Precision< T_Scalar >(n + 1) / kR0 * T_Scalar(boost::math::sph_hankel_1(n, kR0)));
      const T_Scalar invdhn1kR0 = T_Scalar(1., 0.) / dhn1kR0;
      precomputed[n] = -powI[n % 4] * T_Scalar((2 * n + 1), 0.) * std::real(dhn1kR0) * invdhn1kR0;
    }

    const scalar::Precision< T_Scalar > sinPhiInc = std::sin(_phiInc);
    const scalar::Precision< T_Scalar > cosPhiInc = std::cos(_phiInc);

    const scalar::Precision< T_Scalar > sinThetaInc = std::sin(_thetaInc);
    const scalar::Precision< T_Scalar > cosThetaInc = std::cos(_thetaInc);

#pragma omp for
    for(auto i = 0ULL; i < values.size(); ++i) {
      const scalar::Precision< T_Scalar > x = points[3 * i] - _x;
      const scalar::Precision< T_Scalar > y = points[3 * i + 1] - _y;
      const scalar::Precision< T_Scalar > z = points[3 * i + 2] - _z;
      const scalar::Precision< T_Scalar > r = std::sqrt(x * x + y * y + z * z);

      const scalar::Precision< T_Scalar > dx = cosThetaInc * sinPhiInc;
      const scalar::Precision< T_Scalar > dy = sinThetaInc * sinPhiInc;
      const scalar::Precision< T_Scalar > dz = cosPhiInc;
      const scalar::Precision< T_Scalar > dr = std::sqrt(dx * dx + dy * dy + dz * dz);

      const scalar::Precision< T_Scalar > kx = (x * dx + y * dy + z * dz) / (r * dr);

      const scalar::Precision< T_Scalar > kr = _k * r;
      common::KahanSum< T_Scalar > sum(0.);
      for(auto n = 0U; n < _nbrOfterm; ++n) {
        sum += T_Scalar(boost::math::legendre_p(n, kx)) * T_Scalar(boost::math::sph_hankel_1(n, kr)) * precomputed[n];
      }
      values[i] = sum.sum();
    }
  }
#else
  template< class T_Scalar >
  void ScatteringByAHardSphere< T_Scalar >::operator()(function::OutputVector< T_Scalar, Degree::Degree0 > &values, const std::vector< scalar::Precision< T_Scalar >, numa::allocator< scalar::Precision< T_Scalar > > > &points) const
  {
    throw common::Exception("Cannot use ScatteringByAHardSphere without Boost");
  }
#endif // HAVE_BOOST

  INSTANTIATE_CLASS(ScatteringByAHardSphere, 2, TEMPLATE_ARGS(std::complex< double >, std::complex< float >))


  //***************************
  // dr_ScatteringByAHardSphere
  //***************************

  template< class T_Scalar >
  dr_ScatteringByAHardSphere< T_Scalar >::dr_ScatteringByAHardSphere(const scalar::Precision< T_Scalar > k, const scalar::Precision< T_Scalar > R0, const scalar::Precision< T_Scalar > x, const scalar::Precision< T_Scalar > y, const scalar::Precision< T_Scalar > z, const unsigned int nbrOfterm, const scalar::Precision< T_Scalar > thetaInc, const scalar::Precision< T_Scalar > phiInc) :
    _k(k), _R0(R0), _x(x), _y(y), _z(z), _nbrOfterm(nbrOfterm), _thetaInc(thetaInc), _phiInc(phiInc)
  {
  }

  template< class T_Scalar >
  dr_ScatteringByAHardSphere< T_Scalar >::dr_ScatteringByAHardSphere(const dr_ScatteringByAHardSphere &other) :
    function::AnalyticalOperation< T_Scalar, Degree::Degree0 >(other), _k(other._k), _R0(other._R0), _x(other._x), _y(other._y), _z(other._z), _nbrOfterm(other._nbrOfterm), _thetaInc(other._thetaInc), _phiInc(other._phiInc)
  {
  }

  template< class T_Scalar >
  dr_ScatteringByAHardSphere< T_Scalar >::~dr_ScatteringByAHardSphere()
  {
  }

#ifdef HAVE_BOOST
  template< class T_Scalar >
  void dr_ScatteringByAHardSphere< T_Scalar >::operator()(function::OutputVector< T_Scalar, Degree::Degree0 > &values, const std::vector< scalar::Precision< T_Scalar >, numa::allocator< scalar::Precision< T_Scalar > > > &points) const
  {
    const scalar::Precision< T_Scalar > kR0 = _k * _R0;
    const T_Scalar powI[4] = {1., T_Scalar(0., 1.), -1., -T_Scalar(0., 1.)};

    std::vector< T_Scalar > precomputed(_nbrOfterm);
    for(auto n = 0U; n < _nbrOfterm; ++n) {
      const T_Scalar dhn1kR0 = (!n ? -T_Scalar(boost::math::sph_hankel_1(1, kR0)) : T_Scalar(boost::math::sph_hankel_1(n - 1, kR0)) - scalar::Precision< T_Scalar >(n + 1) / kR0 * T_Scalar(boost::math::sph_hankel_1(n, kR0)));
      const T_Scalar invdhn1kR0 = T_Scalar(1., 0.) / dhn1kR0;
      precomputed[n] = -powI[n % 4] * T_Scalar((2 * n + 1), 0.) * std::real(dhn1kR0) * invdhn1kR0 * _k;
    }

    const scalar::Precision< T_Scalar > sinPhiInc = std::sin(_phiInc);
    const scalar::Precision< T_Scalar > cosPhiInc = std::cos(_phiInc);

    const scalar::Precision< T_Scalar > sinThetaInc = std::sin(_thetaInc);
    const scalar::Precision< T_Scalar > cosThetaInc = std::cos(_thetaInc);

#pragma omp for
    for(auto i = 0ULL; i < values.size(); ++i) {
      const scalar::Precision< T_Scalar > x = points[3 * i] - _x;
      const scalar::Precision< T_Scalar > y = points[3 * i + 1] - _y;
      const scalar::Precision< T_Scalar > z = points[3 * i + 2] - _z;
      const scalar::Precision< T_Scalar > r = std::sqrt(x * x + y * y + z * z);

      const scalar::Precision< T_Scalar > dx = cosThetaInc * sinPhiInc;
      const scalar::Precision< T_Scalar > dy = sinThetaInc * sinPhiInc;
      const scalar::Precision< T_Scalar > dz = cosPhiInc;
      const scalar::Precision< T_Scalar > dr = std::sqrt(dx * dx + dy * dy + dz * dz);

      const scalar::Precision< T_Scalar > kx = (x * dx + y * dy + z * dz) / (r * dr);

      const scalar::Precision< T_Scalar > kr = _k * r;
      common::KahanSum< T_Scalar > sum(0.);
      for(auto n = 0U; n < _nbrOfterm; ++n) {
        const T_Scalar dhn1kr = (!n ? -T_Scalar(boost::math::sph_hankel_1(1, kr)) : T_Scalar(boost::math::sph_hankel_1(n - 1, kr)) - scalar::Precision< T_Scalar >(n + 1) / kr * T_Scalar(boost::math::sph_hankel_1(n, kr)));
        sum += T_Scalar(boost::math::legendre_p(n, kx)) * dhn1kr * precomputed[n];
      }
      values[i] = sum.sum();
    }
  }
#else
  template< class T_Scalar >
  void dr_ScatteringByAHardSphere< T_Scalar >::operator()(function::OutputVector< T_Scalar, Degree::Degree0 > &values, const std::vector< scalar::Precision< T_Scalar >, numa::allocator< scalar::Precision< T_Scalar > > > &points) const
  {
    throw common::Exception("Cannot use dr_ScatteringByAHardSphere without Boost");
  }
#endif // HAVE_BOOST

  INSTANTIATE_CLASS(dr_ScatteringByAHardSphere, 2, TEMPLATE_ARGS(std::complex< double >, std::complex< float >))


  //***************************
  // transmissionScatteringByASphere
  //***************************

  template< class T_Scalar >
  transmissionScatteringByASphere< T_Scalar >::transmissionScatteringByASphere(const scalar::Precision< T_Scalar > kint, const scalar::Precision< T_Scalar > kout, const scalar::Precision< T_Scalar > R0, const scalar::Precision< T_Scalar > rhoint, const scalar::Precision< T_Scalar > rhoout, const scalar::Precision< T_Scalar > x, const scalar::Precision< T_Scalar > y, const scalar::Precision< T_Scalar > z, const unsigned int nbrOfterm) :
    _kint(kint), _kout(kout), _R0(R0), _rhoint(rhoint), _rhoout(rhoout), _x(x), _y(y), _z(z), _nbrOfterm(nbrOfterm)
  {
  }

  template< class T_Scalar >
  transmissionScatteringByASphere< T_Scalar >::transmissionScatteringByASphere(const transmissionScatteringByASphere &other) :
    function::AnalyticalOperation< T_Scalar, Degree::Degree0 >(other), _kint(other._kint), _kout(other._kout), _R0(other._R0), _rhoint(other._rhoint), _rhoout(other._rhoout), _x(other._x), _y(other._y), _z(other._z), _nbrOfterm(other._nbrOfterm)
  {
  }

  template< class T_Scalar >
  transmissionScatteringByASphere< T_Scalar >::~transmissionScatteringByASphere()
  {
  }

#ifdef HAVE_BOOST
  template< class T_Scalar >
  void transmissionScatteringByASphere< T_Scalar >::operator()(function::OutputVector< T_Scalar, Degree::Degree0 > &values, const std::vector< scalar::Precision< T_Scalar >, numa::allocator< scalar::Precision< T_Scalar > > > &points) const
  {
    const scalar::Precision< T_Scalar > kintR0 = _kint * _R0;
    const scalar::Precision< T_Scalar > koutR0 = _kout * _R0;
    const T_Scalar powI[4] = {1., T_Scalar(0., 1.), -1., -T_Scalar(0., 1.)};
    const scalar::Precision< T_Scalar > rho = _rhoout / _rhoint;
    const scalar::Precision< T_Scalar > rhokint = rho * _kint;
    std::vector< T_Scalar > precomputed(_nbrOfterm);
    for(auto n = 0U; n < _nbrOfterm; ++n) {
      const T_Scalar dhn1koutR0 = (!n ? -T_Scalar(boost::math::sph_hankel_1(1, koutR0)) : T_Scalar(boost::math::sph_hankel_1(n - 1, koutR0)) - scalar::Precision< T_Scalar >(n + 1) / koutR0 * T_Scalar(boost::math::sph_hankel_1(n, koutR0)));
      const T_Scalar djn1koutR0 = (!n ? -T_Scalar(boost::math::sph_bessel(1, koutR0)) : T_Scalar(boost::math::sph_bessel(n - 1, koutR0)) - scalar::Precision< T_Scalar >(n + 1) / koutR0 * T_Scalar(boost::math::sph_bessel(n, koutR0)));
      const T_Scalar djn1kintR0 = (!n ? -T_Scalar(boost::math::sph_bessel(1, kintR0)) : T_Scalar(boost::math::sph_bessel(n - 1, kintR0)) - scalar::Precision< T_Scalar >(n + 1) / kintR0 * T_Scalar(boost::math::sph_bessel(n, kintR0)));
      const T_Scalar hn1koutR0 = T_Scalar(boost::math::sph_hankel_1(n, koutR0));
      const T_Scalar jn1koutR0 = T_Scalar(boost::math::sph_bessel(n, koutR0));
      const T_Scalar jn1kintR0 = T_Scalar(boost::math::sph_bessel(n, kintR0));
      precomputed[n] = powI[n % 4] * T_Scalar((2 * n + 1), 0.) * _kout * (djn1koutR0 * hn1koutR0 - dhn1koutR0 * jn1koutR0) / (rhokint * djn1kintR0 * hn1koutR0 - _kout * dhn1koutR0 * jn1kintR0);
    }

#pragma omp for
    for(auto i = 0ULL; i < values.size(); ++i) {
      const scalar::Precision< T_Scalar > x = points[3 * i] - _x;
      const scalar::Precision< T_Scalar > y = points[3 * i + 1] - _y;
      const scalar::Precision< T_Scalar > z = points[3 * i + 2] - _z;
      const scalar::Precision< T_Scalar > r2_2 = x * x + y * y;
      const scalar::Precision< T_Scalar > r3_2 = r2_2 + z * z;
      const scalar::Precision< T_Scalar > r = std::sqrt(r3_2);

      const scalar::Precision< T_Scalar > cosPhi = z / r; // std::cos(std::acos(z / r))

      const scalar::Precision< T_Scalar > kintr = _kint * r;
      common::KahanSum< T_Scalar > sum(0.);

      for(auto n = 0U; n < _nbrOfterm; ++n) {
        sum += T_Scalar(boost::math::sph_bessel(n, kintr)) * T_Scalar(boost::math::legendre_p(n, cosPhi)) * precomputed[n];
      }
      values[i] = sum.sum();
    }
  }
#else
  template< class T_Scalar >
  void transmissionScatteringByASphere< T_Scalar >::operator()(function::OutputVector< T_Scalar, Degree::Degree0 > &values, const std::vector< scalar::Precision< T_Scalar >, numa::allocator< scalar::Precision< T_Scalar > > > &points) const
  {
    throw common::Exception("Cannot use transmissionScatteringByASphere without Boost");
  }
#endif // HAVE_BOOST

  INSTANTIATE_CLASS(transmissionScatteringByASphere, 2, TEMPLATE_ARGS(std::complex< double >, std::complex< float >))


} // namespace gmshfem::analytics::helmholtz3D
