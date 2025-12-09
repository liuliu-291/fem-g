// GmshFEM - Copyright (C) 2019-2022, A. Royer, E. Béchet, C. Geuzaine, Université de Liège
//
// See the LICENSE.txt file for license information. Please report all
// issues on https://gitlab.onelab.info/gmsh/fem/issues

#ifndef H_GMSHFEM_AXISYMMETRY
#define H_GMSHFEM_AXISYMMETRY

#include "JacobiansModificator.h"

namespace gmshfem::domain
{


  // Mesh must be in the plane (0, ex, ey) with axis along ey
  template< class T_PScalar >
  class Axisymmetry final : public JacobiansModificator
  {
   public:
    Axisymmetry();
    virtual ~Axisymmetry();

    virtual bool needGaussCoordinates() const override;

    void apply(const std::vector< T_PScalar, numa::allocator< T_PScalar > > &points, std::vector< T_PScalar, numa::allocator< T_PScalar > > &determinants, std::vector< T_PScalar, numa::allocator< T_PScalar > > &jacobians) const override;
    virtual Axisymmetry< T_PScalar > *copy() const override;
  };


} // namespace gmshfem::domain

#endif // H_GMSHFEM_AXISYMMETRY
