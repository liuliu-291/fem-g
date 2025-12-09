// GmshFEM - Copyright (C) 2019-2022, A. Royer, E. Béchet, C. Geuzaine, Université de Liège
//
// See the LICENSE.txt file for license information. Please report all
// issues on https://gitlab.onelab.info/gmsh/fem/issues

#ifndef H_GMSHFEM_CYLINDRICALSHELL
#define H_GMSHFEM_CYLINDRICALSHELL

#include "JacobiansModificator.h"

namespace gmshfem::domain
{


  // T_axis = 0 for X, 1 for Y, 2 for Z
  template< class T_PScalar, unsigned int T_axis >
  class CylindricalShell final : JacobiansModificator
  {
   private:
    const T_PScalar _Rmin;
    const T_PScalar _Rmax;
    const T_PScalar _center[3];
    const unsigned int _exponant;

   public:
    CylindricalShell(const T_PScalar Rmin, const T_PScalar Rmax, const T_PScalar centerX = 0., const T_PScalar centerY = 0., const T_PScalar centerZ = 0., const unsigned int exponant = 1);
    virtual ~CylindricalShell();

    virtual bool needGaussCoordinates() const override;

    void apply(const std::vector< T_PScalar, numa::allocator< T_PScalar > > &points, std::vector< T_PScalar, numa::allocator< T_PScalar > > &determinants, std::vector< T_PScalar, numa::allocator< T_PScalar > > &jacobians) const override;
    virtual CylindricalShell< T_PScalar, T_axis > *copy() const override;
  };


  template< class T_PScalar >
  using CylindricalShellX = CylindricalShell< T_PScalar, 0 >;
  template< class T_PScalar >
  using CylindricalShellY = CylindricalShell< T_PScalar, 1 >;
  template< class T_PScalar >
  using CylindricalShellZ = CylindricalShell< T_PScalar, 2 >;


} // namespace gmshfem::domain

#endif // H_GMSHFEM_CYLINDRICALSHELL
