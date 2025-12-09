// GmshFEM - Copyright (C) 2019-2022, A. Royer, E. Béchet, C. Geuzaine, Université de Liège
//
// See the LICENSE.txt file for license information. Please report all
// issues on https://gitlab.onelab.info/gmsh/fem/issues

#ifndef H_GMSHFEM_AXISYMMETRYSHELL
#define H_GMSHFEM_AXISYMMETRYSHELL

#include "JacobiansModificator.h"

namespace gmshfem::domain
{


  // Mesh must be in the plane (0, ex, ey) with axis along ey
  template< class T_PScalar >
  class AxisymmetryShell final : public JacobiansModificator
  {
   private:
    const T_PScalar _Rmin;
    const T_PScalar _Rmax;
    const T_PScalar _centerY;
    const unsigned int _exponant;

   public:
    AxisymmetryShell(const T_PScalar Rmin, const T_PScalar Rmax, const T_PScalar centerY = 0., const unsigned int exponant = 1);
    virtual ~AxisymmetryShell();

    virtual bool needGaussCoordinates() const override;

    void apply(const std::vector< T_PScalar, numa::allocator< T_PScalar > > &points, std::vector< T_PScalar, numa::allocator< T_PScalar > > &determinants, std::vector< T_PScalar, numa::allocator< T_PScalar > > &jacobians) const override;
    virtual AxisymmetryShell< T_PScalar > *copy() const override;
  };


} // namespace gmshfem::domain

#endif // H_GMSHFEM_AXISYMMETRYSHELL
