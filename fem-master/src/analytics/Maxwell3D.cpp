// GmshFEM - Copyright (C) 2019-2022, A. Royer, E. Béchet, C. Geuzaine, Université de Liège
//
// See the LICENSE.txt file for license information. Please report all
// issues on https://gitlab.onelab.info/gmsh/fem/issues

#include "Maxwell3D.h"

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


namespace gmshfem::analytics::maxwell3D
{


  //***************************
  // Diffraction by a conducting sphere, J=H^n
  //***************************

  template< class T_Scalar >
  JElectricSurfaceCurrent< T_Scalar >::JElectricSurfaceCurrent(const scalar::Precision< T_Scalar > k_int, const scalar::Precision< T_Scalar > k_out, const scalar::Precision< std::complex< T_Scalar > > k2_int, const scalar::Precision< std::complex< T_Scalar > > k2_out, const scalar::Precision< T_Scalar > R0, const scalar::Precision< T_Scalar > x, const scalar::Precision< T_Scalar > y, const scalar::Precision< T_Scalar > z, const unsigned int nbrOfterm) :
    _k_int(k_int), _k_out(k_out), _k2_int(k2_int), _k2_out(k2_out), _R0(R0), _x(x), _y(y), _z(z), _nbrOfterm(nbrOfterm)
  {
  }

  template< class T_Scalar >
  JElectricSurfaceCurrent< T_Scalar >::JElectricSurfaceCurrent(const JElectricSurfaceCurrent &other) :
    function::AnalyticalOperation< T_Scalar, Degree::Degree1 >(other), _k_int(other._k_int), _k_out(other._k_out), _k2_int(other._k2_int), _k2_out(other._k2_out), _R0(other._R0), _x(other._x), _y(other._y), _z(other._z), _nbrOfterm(other._nbrOfterm)
  {
  }

  template< class T_Scalar >
  JElectricSurfaceCurrent< T_Scalar >::~JElectricSurfaceCurrent()
  {
  }

#ifdef HAVE_BOOST
  template< class T_Scalar >
  void JElectricSurfaceCurrent< T_Scalar >::operator()(function::OutputVector< T_Scalar, Degree::Degree1 > &values, const std::vector< scalar::Precision< T_Scalar >, numa::allocator< scalar::Precision< T_Scalar > > > &points) const
  {
    const scalar::Precision< T_Scalar > kintR0 = _k_int * _R0;
    const scalar::Precision< T_Scalar > koutR0 = _k_out * _R0;
    const T_Scalar powI[4] = {1., T_Scalar(0., 1.), -1., -T_Scalar(0., 1.)};
    std::vector< T_Scalar > precomputedBle(_nbrOfterm);
    std::vector< T_Scalar > precomputedBlm(_nbrOfterm);
    for(auto n = 1U; n <= _nbrOfterm; ++n) {
      const T_Scalar jn1koutR0 = T_Scalar(boost::math::sph_bessel(n, koutR0));
      const T_Scalar jn1kintR0 = T_Scalar(boost::math::sph_bessel(n, kintR0));
      const T_Scalar koutR0jn1koutR0 = koutR0 * jn1koutR0;
      const T_Scalar kintR0jn1kintR0 = kintR0 * jn1kintR0;
      const T_Scalar dkoutR0jn1koutR0 = koutR0 * (T_Scalar(boost::math::sph_bessel(n - 1, koutR0)) - scalar::Precision< T_Scalar >(n + 1) / koutR0 * jn1koutR0) + jn1koutR0;
      const T_Scalar dkintR0jn1kintR0 = kintR0 * (T_Scalar(boost::math::sph_bessel(n - 1, kintR0)) - scalar::Precision< T_Scalar >(n + 1) / kintR0 * jn1kintR0) + jn1kintR0;

      const T_Scalar hn1koutR0 = T_Scalar(boost::math::sph_hankel_1(n, koutR0));
      const T_Scalar koutR0hn1koutR0 = koutR0 * hn1koutR0;
      const T_Scalar dkoutR0hn1koutR0 = hn1koutR0 + koutR0 * (T_Scalar(boost::math::sph_hankel_1(n - 1, koutR0)) - scalar::Precision< T_Scalar >(n + 1) / koutR0 * hn1koutR0);
      precomputedBle[n - 1] = powI[(n + 1) % 4] * T_Scalar((2. * n + 1) / (n * (n + 1)), 0.) * (_k2_out * _k_int * dkoutR0jn1koutR0 * kintR0jn1kintR0 - _k2_int * _k_out * koutR0jn1koutR0 * dkintR0jn1kintR0) / (_k2_out * _k_int * dkoutR0hn1koutR0 * kintR0jn1kintR0 - _k2_int * _k_out * koutR0hn1koutR0 * dkintR0jn1kintR0);
      precomputedBlm[n - 1] = powI[(n + 1) % 4] * T_Scalar((2. * n + 1) / (n * (n + 1)), 0.) * (_k2_out * _k_int * koutR0jn1koutR0 * dkintR0jn1kintR0 - _k2_int * _k_out * dkoutR0jn1koutR0 * kintR0jn1kintR0) / (_k2_out * _k_int * koutR0hn1koutR0 * dkintR0jn1kintR0 - _k2_int * _k_out * dkoutR0hn1koutR0 * kintR0jn1kintR0);
    }

#pragma omp for
    for(auto i = 0ULL; i < values.size(); ++i) {
      const scalar::Precision< T_Scalar > x = points[3 * i] - _x;
      const scalar::Precision< T_Scalar > y = points[3 * i + 1] - _y;
      const scalar::Precision< T_Scalar > z = points[3 * i + 2] - _z;
      const scalar::Precision< T_Scalar > r = std::sqrt(x * x + y * y + z * z);
      const scalar::Precision< T_Scalar > cosPhi = z / r;
      const scalar::Precision< T_Scalar > cosPhiSquare = cosPhi * cosPhi;
      const scalar::Precision< T_Scalar > sinPhi = std::sin(std::acos(z / r));
      const scalar::Precision< T_Scalar > teta = std::atan2(y, x);
      const scalar::Precision< T_Scalar > cosTeta = std::cos(teta);
      const scalar::Precision< T_Scalar > sinTeta = std::sin(teta);
      const scalar::Precision< T_Scalar > koutr = _k_out * r;
      common::KahanSum< T_Scalar > sumJphi(0.);
      common::KahanSum< T_Scalar > sumJteta(0.);
      for(auto n = 1U; n <= _nbrOfterm; ++n) {
        const T_Scalar hn1koutr = T_Scalar(boost::math::sph_hankel_1(n, koutr));
        const T_Scalar koutrhn1koutr = koutr * hn1koutr;
        const T_Scalar dkoutrhn1koutr = hn1koutr + koutr * (T_Scalar(boost::math::sph_hankel_1(n - 1, koutr)) - scalar::Precision< T_Scalar >(n + 1) / koutr * hn1koutr);
        T_Scalar dlegendre;
        if(cosPhi != 1 && cosPhi != -1) {
          dlegendre = T_Scalar(n, 0.) * (-cosPhi * T_Scalar(boost::math::legendre_p(n, 0, cosPhi)) + T_Scalar(boost::math::legendre_p(n - 1, 0, cosPhi))) / (1 - cosPhiSquare);
        }
        else if(cosPhi == 1) {
          dlegendre = T_Scalar(n * (n + 1) / 2., 0.);
        }
        else {
          if(n % 2 == 0) {
            dlegendre = T_Scalar(-n * (n + 1) / 2., 0.);
          }
          else {
            dlegendre = T_Scalar(n * (n + 1) / 2., 0.);
          }
        }
        const T_Scalar Legendre2n = T_Scalar(boost::math::legendre_p(n, 2, cosPhi));
        sumJphi += -precomputedBle[n - 1] * koutrhn1koutr * (cosPhi * dlegendre - Legendre2n) + T_Scalar(0., 1.) * precomputedBlm[n - 1] * dkoutrhn1koutr * dlegendre;
        sumJteta += precomputedBle[n - 1] * koutrhn1koutr * dlegendre - T_Scalar(0., 1.) * precomputedBlm[n - 1] * dkoutrhn1koutr * (cosPhi * dlegendre - Legendre2n);
      }

      const T_Scalar Jphi_sc = sumJphi.sum() * cosTeta / (r * _k2_out);
      const T_Scalar Jteta_sc = sumJteta.sum() * sinTeta / (r * _k2_out);
      const T_Scalar Hy_inc = T_Scalar(0., 1.) * _k_out / _k2_out * std::exp(T_Scalar(0., 1.) * _k_out * z);
      const T_Scalar Jphi_inc = Hy_inc * cosTeta;
      const T_Scalar Jteta_inc = -Hy_inc * cosPhi * sinTeta;
      values[i] = typename MathObject< T_Scalar, Degree::Degree1 >::Object(cosPhi * cosTeta * (Jphi_inc + Jphi_sc) - sinTeta * (Jteta_sc + Jteta_inc), cosPhi * sinTeta * (Jphi_inc + Jphi_sc) + cosTeta * (Jteta_sc + Jteta_inc), -sinPhi * (Jphi_inc + Jphi_sc));
    }
  }
#else
  template< class T_Scalar >
  void JElectricSurfaceCurrent< T_Scalar >::operator()(function::OutputVector< T_Scalar, Degree::Degree1 > &values, const std::vector< scalar::Precision< T_Scalar >, numa::allocator< scalar::Precision< T_Scalar > > > &points) const
  {
    throw common::Exception("Cannot use JElectricSurfaceCurrent without Boost");
  }
#endif // HAVE_BOOST

  INSTANTIATE_CLASS(JElectricSurfaceCurrent, 2, TEMPLATE_ARGS(std::complex< double >, std::complex< float >))


  //***************************
  // Diffraction by a conducting sphere, M=E^n
  //***************************

  template< class T_Scalar >
  MMagneticSurfaceCurrent< T_Scalar >::MMagneticSurfaceCurrent(const scalar::Precision< T_Scalar > k_int, const scalar::Precision< T_Scalar > k_out, const scalar::Precision< std::complex< T_Scalar > > k2_int, const scalar::Precision< std::complex< T_Scalar > > k2_out, const scalar::Precision< T_Scalar > R0, const scalar::Precision< T_Scalar > x, const scalar::Precision< T_Scalar > y, const scalar::Precision< T_Scalar > z, const unsigned int nbrOfterm) :
    _k_int(k_int), _k_out(k_out), _k2_int(k2_int), _k2_out(k2_out), _R0(R0), _x(x), _y(y), _z(z), _nbrOfterm(nbrOfterm)
  {
  }

  template< class T_Scalar >
  MMagneticSurfaceCurrent< T_Scalar >::MMagneticSurfaceCurrent(const MMagneticSurfaceCurrent &other) :
    function::AnalyticalOperation< T_Scalar, Degree::Degree1 >(other), _k_int(other._k_int), _k_out(other._k_out), _k2_int(other._k2_int), _k2_out(other._k2_out), _R0(other._R0), _x(other._x), _y(other._y), _z(other._z), _nbrOfterm(other._nbrOfterm)
  {
  }

  template< class T_Scalar >
  MMagneticSurfaceCurrent< T_Scalar >::~MMagneticSurfaceCurrent()
  {
  }

#ifdef HAVE_BOOST
  template< class T_Scalar >
  void MMagneticSurfaceCurrent< T_Scalar >::operator()(function::OutputVector< T_Scalar, Degree::Degree1 > &values, const std::vector< scalar::Precision< T_Scalar >, numa::allocator< scalar::Precision< T_Scalar > > > &points) const
  {
    const scalar::Precision< T_Scalar > kintR0 = _k_int * _R0;
    const scalar::Precision< T_Scalar > koutR0 = _k_out * _R0;
    const T_Scalar powI[4] = {1., T_Scalar(0., 1.), -1., -T_Scalar(0., 1.)};
    std::vector< T_Scalar > precomputedBle(_nbrOfterm);
    std::vector< T_Scalar > precomputedBlm(_nbrOfterm);
    for(auto n = 1U; n <= _nbrOfterm; ++n) {
      const T_Scalar jn1koutR0 = T_Scalar(boost::math::sph_bessel(n, koutR0));
      const T_Scalar jn1kintR0 = T_Scalar(boost::math::sph_bessel(n, kintR0));
      const T_Scalar koutR0jn1koutR0 = koutR0 * jn1koutR0;
      const T_Scalar kintR0jn1kintR0 = kintR0 * jn1kintR0;
      const T_Scalar dkoutR0jn1koutR0 = koutR0 * (T_Scalar(boost::math::sph_bessel(n - 1, koutR0)) - scalar::Precision< T_Scalar >(n + 1) / koutR0 * jn1koutR0) + jn1koutR0;
      const T_Scalar dkintR0jn1kintR0 = kintR0 * (T_Scalar(boost::math::sph_bessel(n - 1, kintR0)) - scalar::Precision< T_Scalar >(n + 1) / kintR0 * jn1kintR0) + jn1kintR0;

      const T_Scalar hn1koutR0 = T_Scalar(boost::math::sph_hankel_1(n, koutR0));
      const T_Scalar koutR0hn1koutR0 = koutR0 * hn1koutR0;
      const T_Scalar dkoutR0hn1koutR0 = hn1koutR0 + koutR0 * (T_Scalar(boost::math::sph_hankel_1(n - 1, koutR0)) - scalar::Precision< T_Scalar >(n + 1) / koutR0 * hn1koutR0);
      precomputedBle[n - 1] = powI[(n + 1) % 4] * T_Scalar((2. * n + 1) / (n * (n + 1)), 0.) * (_k2_out * _k_int * dkoutR0jn1koutR0 * kintR0jn1kintR0 - _k2_int * _k_out * koutR0jn1koutR0 * dkintR0jn1kintR0) / (_k2_out * _k_int * dkoutR0hn1koutR0 * kintR0jn1kintR0 - _k2_int * _k_out * koutR0hn1koutR0 * dkintR0jn1kintR0);
      precomputedBlm[n - 1] = powI[(n + 1) % 4] * T_Scalar((2. * n + 1) / (n * (n + 1)), 0.) * (_k2_out * _k_int * koutR0jn1koutR0 * dkintR0jn1kintR0 - _k2_int * _k_out * dkoutR0jn1koutR0 * kintR0jn1kintR0) / (_k2_out * _k_int * koutR0hn1koutR0 * dkintR0jn1kintR0 - _k2_int * _k_out * dkoutR0hn1koutR0 * kintR0jn1kintR0);
    }

#pragma omp for
    for(auto i = 0ULL; i < values.size(); ++i) {
      const scalar::Precision< T_Scalar > x = points[3 * i] - _x;
      const scalar::Precision< T_Scalar > y = points[3 * i + 1] - _y;
      const scalar::Precision< T_Scalar > z = points[3 * i + 2] - _z;
      const scalar::Precision< T_Scalar > r = std::sqrt(x * x + y * y + z * z);
      const scalar::Precision< T_Scalar > cosPhi = z / r;
      const scalar::Precision< T_Scalar > cosPhiSquare = cosPhi * cosPhi;
      const scalar::Precision< T_Scalar > sinPhi = std::sin(std::acos(z / r));
      const scalar::Precision< T_Scalar > teta = std::atan2(y, x);
      const scalar::Precision< T_Scalar > cosTeta = std::cos(teta);
      const scalar::Precision< T_Scalar > sinTeta = std::sin(teta);
      const scalar::Precision< T_Scalar > koutr = _k_out * r;
      common::KahanSum< T_Scalar > sumMphi(0.);
      common::KahanSum< T_Scalar > sumMteta(0.);
      for(auto n = 1U; n <= _nbrOfterm; ++n) {
        const T_Scalar hn1koutr = T_Scalar(boost::math::sph_hankel_1(n, koutr));
        const T_Scalar koutrhn1koutr = koutr * hn1koutr;
        const T_Scalar dkoutrhn1koutr = hn1koutr + koutr * (T_Scalar(boost::math::sph_hankel_1(n - 1, koutr)) - scalar::Precision< T_Scalar >(n + 1) / koutr * hn1koutr);
        T_Scalar dlegendre;
        if(cosPhi != 1 && cosPhi != -1) {
          dlegendre = T_Scalar(n, 0.) * (-cosPhi * T_Scalar(boost::math::legendre_p(n, 0, cosPhi)) + T_Scalar(boost::math::legendre_p(n - 1, 0, cosPhi))) / (1 - cosPhiSquare);
        }
        else if(cosPhi == 1) {
          dlegendre = T_Scalar(n * (n + 1) / 2., 0.);
        }
        else {
          if(n % 2 == 0) {
            dlegendre = T_Scalar(-n * (n + 1) / 2., 0.);
          }
          else {
            dlegendre = T_Scalar(n * (n + 1) / 2., 0.);
          }
        }
        const T_Scalar Legendre2n = T_Scalar(boost::math::legendre_p(n, 2, cosPhi));
        sumMphi += precomputedBle[n - 1] * dkoutrhn1koutr * dlegendre + T_Scalar(0., 1.) * precomputedBlm[n - 1] * koutrhn1koutr * (cosPhi * dlegendre - Legendre2n);
        sumMteta += precomputedBle[n - 1] * dkoutrhn1koutr * (cosPhi * dlegendre - Legendre2n) + T_Scalar(0., 1.) * precomputedBlm[n - 1] * koutrhn1koutr * dlegendre;
      }

      const T_Scalar Mphi_sc = -sumMphi.sum() * sinTeta / (r * _k_out);
      const T_Scalar Mteta_sc = -sumMteta.sum() * cosTeta / (r * _k_out);
      const T_Scalar Ex_inc = std::exp(T_Scalar(0., 1.) * _k_out * z);
      const T_Scalar Mphi_inc = -Ex_inc * sinTeta;
      const T_Scalar Mteta_inc = -Ex_inc * cosPhi * cosTeta;

      values[i] = typename MathObject< T_Scalar, Degree::Degree1 >::Object(cosPhi * cosTeta * (Mphi_inc + Mphi_sc) - sinTeta * (Mteta_sc + Mteta_inc), cosPhi * sinTeta * (Mphi_inc + Mphi_sc) + cosTeta * (Mteta_sc + Mteta_inc), -sinPhi * (Mphi_inc + Mphi_sc));
    }
  }
#else
  template< class T_Scalar >
  void MMagneticSurfaceCurrent< T_Scalar >::operator()(function::OutputVector< T_Scalar, Degree::Degree1 > &values, const std::vector< scalar::Precision< T_Scalar >, numa::allocator< scalar::Precision< T_Scalar > > > &points) const
  {
    throw common::Exception("Cannot use MMagneticSurfaceCurrent without Boost");
  }
#endif // HAVE_BOOST

  INSTANTIATE_CLASS(MMagneticSurfaceCurrent, 2, TEMPLATE_ARGS(std::complex< double >, std::complex< float >))


} // namespace gmshfem::analytics::maxwell3D
