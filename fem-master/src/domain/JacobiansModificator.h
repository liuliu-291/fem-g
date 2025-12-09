// GmshFEM - Copyright (C) 2019-2022, A. Royer, E. Béchet, C. Geuzaine, Université de Liège
//
// See the LICENSE.txt file for license information. Please report all
// issues on https://gitlab.onelab.info/gmsh/fem/issues

#ifndef H_GMSHFEM_JACOBIANSMODIFICATOR
#define H_GMSHFEM_JACOBIANSMODIFICATOR

#include "numa.h"

#include <vector>

namespace gmshfem::domain
{


  class JacobiansModificator
  {
   public:
    JacobiansModificator();
    virtual ~JacobiansModificator();

    virtual bool needGaussCoordinates() const = 0;

    virtual void apply(const std::vector< double, numa::allocator< double > > &points, std::vector< double, numa::allocator< double > > &determinants, std::vector< double, numa::allocator< double > > &jacobians) const;
    virtual void apply(const std::vector< float, numa::allocator< float > > &points, std::vector< float, numa::allocator< float > > &determinants, std::vector< float, numa::allocator< float > > &jacobians) const;

    virtual JacobiansModificator *copy() const = 0;
  };


} // namespace gmshfem::domain

#include "Axisymmetry.h"
#include "AxisymmetryShell.h"
#include "CylindricalShell.h"
#include "PolarShell.h"
#include "SphericalShell.h"

#endif // H_GMSHFEM_JACOBIANSMODIFICATOR
