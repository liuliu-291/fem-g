// GmshFEM - Copyright (C) 2019-2022, A. Royer, E. Béchet, C. Geuzaine, Université de Liège
//
// See the LICENSE.txt file for license information. Please report all
// issues on https://gitlab.onelab.info/gmsh/fem/issues

#include "CylindricalShell.h"

#include "Exception.h"
#include "Message.h"
#include "OmpInterface.h"
#include "instantiate.h"
#include "scalar.h"

#include <cmath>

namespace gmshfem::domain
{


  template< class T_PScalar, unsigned int T_axis >
  CylindricalShell< T_PScalar, T_axis >::CylindricalShell(const T_PScalar Rmin, const T_PScalar Rmax, const T_PScalar centerX, const T_PScalar centerY, const T_PScalar centerZ, const unsigned int exponant) :
    _Rmin(Rmin), _Rmax(Rmax), _center{centerX, centerY, centerZ}, _exponant(exponant)
  {
  }

  template< class T_PScalar, unsigned int T_axis >
  CylindricalShell< T_PScalar, T_axis >::~CylindricalShell()
  {
  }

  template< class T_PScalar, unsigned int T_axis >
  bool CylindricalShell< T_PScalar, T_axis >::needGaussCoordinates() const
  {
    return true;
  }

  template< class T_PScalar, unsigned int T_axis >
  void CylindricalShell< T_PScalar, T_axis >::apply(const std::vector< T_PScalar, numa::allocator< T_PScalar > > &points, std::vector< T_PScalar, numa::allocator< T_PScalar > > &determinants, std::vector< T_PScalar, numa::allocator< T_PScalar > > &jacobians) const
  {
    const unsigned long long nbrValues = points.size() / 3;

#pragma omp for
    for(auto i = 0ULL; i < nbrValues; ++i) {
      T_PScalar R2 = 0.;
      for(auto j = 0; j < 3; ++j) {
        if(T_axis == j) continue;
        R2 += (points[3 * i + j] - _center[j]) * (points[3 * i + j] - _center[j]);
      }
      T_PScalar R = std::sqrt(R2);

      if(R < _Rmin - 10. * scalar::Epsilon< T_PScalar >::value || R > _Rmax + 10. * scalar::Epsilon< T_PScalar >::value) {
        throw common::Exception("Wrong paramameters for 'CylindricalShell' the Gauss point radius is " + std::to_string(R) + " whereas the shell boundaries are (" + std::to_string(_Rmin) + ", " + std::to_string(_Rmax) + ")");
      }

      if(R >= _Rmax) {
        R = _Rmax - 10 * scalar::Epsilon< T_PScalar >::value;
      }

      T_PScalar dRdx[3];
      for(auto j = 0; j < 3; ++j) {
        if(T_axis == j) {
          dRdx[j] = 0.;
        }
        else {
          dRdx[j] = (points[3 * i + j] - _center[j]) / R;
        }
      }
      const T_PScalar f = std::pow((_Rmin * (_Rmax - _Rmin)) / (R * (_Rmax - R)), _exponant);
      const T_PScalar theta = _exponant * (_Rmax - 2. * R) / (_Rmax - R);

      T_PScalar jac[3][3];
      jac[0][0] = f * (1.0 - theta * dRdx[0] * dRdx[0]);
      jac[0][1] = f * (-theta * dRdx[0] * dRdx[1]);
      jac[0][2] = f * (-theta * dRdx[0] * dRdx[2]);

      jac[1][0] = f * (-theta * dRdx[1] * dRdx[0]);
      jac[1][1] = f * (1.0 - theta * dRdx[1] * dRdx[1]);
      jac[1][2] = f * (-theta * dRdx[1] * dRdx[2]);

      jac[2][0] = f * (-theta * dRdx[2] * dRdx[0]);
      jac[2][1] = f * (-theta * dRdx[2] * dRdx[1]);
      jac[2][2] = f * (1.0 - theta * dRdx[2] * dRdx[2]);

      jac[T_axis][T_axis] = 1.;

      const T_PScalar det = f * f * (1.0 - theta);

      if(!jacobians.empty()) {
        T_PScalar newJac[9];
        for(auto k = 0; k < 3; ++k) {
          for(auto l = 0; l < 3; ++l) {
            newJac[k * 3 + l] = jacobians[i * 9 + k * 3 + 0] * jac[0][l] + jacobians[i * 9 + k * 3 + 1] * jac[1][l] + jacobians[i * 9 + k * 3 + 2] * jac[2][l];
          }
        }
        for(auto k = 0; k < 9; ++k) {
          jacobians[9 * i + k] = newJac[k];
        }
      }
      if(!determinants.empty()) {
        determinants[i] *= det;
      }
    }
  }

  template< class T_PScalar, unsigned int T_axis >
  CylindricalShell< T_PScalar, T_axis > *CylindricalShell< T_PScalar, T_axis >::copy() const
  {
    return new CylindricalShell< T_PScalar, T_axis >(_Rmin, _Rmax, _center[0], _center[1], _center[2], _exponant);
  }

  INSTANTIATE_CLASS_2(CylindricalShell, 2, 3, TEMPLATE_ARGS(double, float), TEMPLATE_ARGS(0, 1, 2))


} // namespace gmshfem::domain
