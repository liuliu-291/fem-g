// GmshFEM - Copyright (C) 2019-2022, A. Royer, E. Béchet, C. Geuzaine, Université de Liège
//
// See the LICENSE.txt file for license information. Please report all
// issues on https://gitlab.onelab.info/gmsh/fem/issues

#include "AxisymmetryShell.h"

#include "Exception.h"
#include "Message.h"
#include "OmpInterface.h"
#include "instantiate.h"
#include "scalar.h"

namespace gmshfem::domain
{


  template< class T_PScalar >
  AxisymmetryShell< T_PScalar >::AxisymmetryShell(const T_PScalar Rmin, const T_PScalar Rmax, const T_PScalar centerY, const unsigned int exponant) :
    _Rmin(Rmin), _Rmax(Rmax), _centerY(centerY), _exponant(exponant)
  {
  }

  template< class T_PScalar >
  AxisymmetryShell< T_PScalar >::~AxisymmetryShell()
  {
  }

  template< class T_PScalar >
  bool AxisymmetryShell< T_PScalar >::needGaussCoordinates() const
  {
    return true;
  }

  template< class T_PScalar >
  void AxisymmetryShell< T_PScalar >::apply(const std::vector< T_PScalar, numa::allocator< T_PScalar > > &points, std::vector< T_PScalar, numa::allocator< T_PScalar > > &determinants, std::vector< T_PScalar, numa::allocator< T_PScalar > > &jacobians) const
  {
    const unsigned long long nbrValues = points.size() / 3;

#pragma omp for
    for(auto i = 0ULL; i < nbrValues; ++i) {
      T_PScalar R = std::sqrt(points[3 * i + 0] * points[3 * i + 0] + (points[3 * i + 1] - _centerY) * (points[3 * i + 1] - _centerY));
      T_PScalar S = points[3 * i + 0];

      if(R < _Rmin - 10. * scalar::Epsilon< T_PScalar >::value || R > _Rmax + 10. * scalar::Epsilon< T_PScalar >::value) {
        throw common::Exception("Wrong paramameters for 'AxisymmetryShell': the integration point radius is " + std::to_string(R) + " whereas the shell boundaries are (" + std::to_string(_Rmin) + ", " + std::to_string(_Rmax) + ")");
      }

      if(R >= _Rmax) {
        R = _Rmax - 10 * scalar::Epsilon< T_PScalar >::value;
      }

      if(S == 0.) {
        S = scalar::Epsilon< T_PScalar >::value;
#pragma omp critical
        msg::debug << "A zero axisymmetric shell coordinate is shifted by eps0(" << S << ")" << msg::endl;
      }

      T_PScalar dRdx[2];
      dRdx[0] = points[3 * i + 0] / R;
      dRdx[1] = (points[3 * i + 1] - _centerY) / R;

      const T_PScalar f = std::pow((_Rmin * (_Rmax - _Rmin)) / (R * (_Rmax - R)), _exponant);
      const T_PScalar theta = _exponant * (_Rmax - 2. * R) / (_Rmax - R);

      T_PScalar jac[2][2];
      jac[0][0] = f * (1.0 - theta * dRdx[0] * dRdx[0]);
      jac[0][1] = f * (-theta * dRdx[0] * dRdx[1]);

      jac[1][0] = f * (-theta * dRdx[1] * dRdx[0]);
      jac[1][1] = f * (1.0 - theta * dRdx[1] * dRdx[1]);

      const T_PScalar det = f * f * f * (1.0 - theta);

      if(!jacobians.empty()) {
        T_PScalar newJac[9];
        for(auto k = 0; k < 2; ++k) {
          for(auto l = 0; l < 2; ++l) {
            newJac[k * 3 + l] = jacobians[i * 9 + k * 3 + 0] * jac[0][l] + jacobians[i * 9 + k * 3 + 1] * jac[1][l];
          }
        }
        newJac[0 * 3 + 2] = newJac[1 * 3 + 2] = newJac[2 * 3 + 0] = newJac[2 * 3 + 1] = 0.;
        newJac[2 * 3 + 2] = S * f;
        for(auto k = 0; k < 9; ++k) {
          jacobians[9 * i + k] = newJac[k];
        }
      }
      if(!determinants.empty()) {
        determinants[i] *= S * det;
      }
    }
  }

  template< class T_PScalar >
  AxisymmetryShell< T_PScalar > *AxisymmetryShell< T_PScalar >::copy() const
  {
    return new AxisymmetryShell< T_PScalar >(_Rmin, _Rmax, _centerY, _exponant);
  }

  INSTANTIATE_CLASS(AxisymmetryShell, 2, TEMPLATE_ARGS(double, float))


} // namespace gmshfem::domain
