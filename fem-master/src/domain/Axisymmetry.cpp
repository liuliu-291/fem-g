// GmshFEM - Copyright (C) 2019-2022, A. Royer, E. Béchet, C. Geuzaine, Université de Liège
//
// See the LICENSE.txt file for license information. Please report all
// issues on https://gitlab.onelab.info/gmsh/fem/issues

#include "Axisymmetry.h"

#include "Message.h"
#include "OmpInterface.h"
#include "instantiate.h"
#include "scalar.h"

namespace gmshfem::domain
{


  template< class T_PScalar >
  Axisymmetry< T_PScalar >::Axisymmetry()
  {
  }

  template< class T_PScalar >
  Axisymmetry< T_PScalar >::~Axisymmetry()
  {
  }

  template< class T_PScalar >
  bool Axisymmetry< T_PScalar >::needGaussCoordinates() const
  {
    return true;
  }

  template< class T_PScalar >
  void Axisymmetry< T_PScalar >::apply(const std::vector< T_PScalar, numa::allocator< T_PScalar > > &points, std::vector< T_PScalar, numa::allocator< T_PScalar > > &determinants, std::vector< T_PScalar, numa::allocator< T_PScalar > > &jacobians) const
  {
    const unsigned long long nbrValues = points.size() / 3;

#pragma omp for
    for(auto i = 0ULL; i < nbrValues; ++i) {
      T_PScalar R = points[3 * i + 0];

      if(R == 0.) {
        R = scalar::Epsilon< T_PScalar >::value;
#pragma omp critical
        msg::debug << "A zero axisymmetric shell coordinate is shifted by eps0 (" << R << ")" << msg::endl;
      }

      if(!jacobians.empty()) {
        jacobians[9 * i + 8] = R;
      }
      if(!determinants.empty()) {
        determinants[i] *= R;
      }
    }
  }

  template< class T_PScalar >
  Axisymmetry< T_PScalar > *Axisymmetry< T_PScalar >::copy() const
  {
    return new Axisymmetry< T_PScalar >();
  }

  INSTANTIATE_CLASS(Axisymmetry, 2, TEMPLATE_ARGS(double, float))


} // namespace gmshfem::domain
