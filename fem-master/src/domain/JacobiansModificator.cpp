// GmshFEM - Copyright (C) 2019-2022, A. Royer, E. Béchet, C. Geuzaine, Université de Liège
//
// See the LICENSE.txt file for license information. Please report all
// issues on https://gitlab.onelab.info/gmsh/fem/issues

#include "JacobiansModificator.h"

#include "Exception.h"

namespace gmshfem::domain
{


  JacobiansModificator::JacobiansModificator()
  {
  }

  JacobiansModificator::~JacobiansModificator()
  {
  }

  void JacobiansModificator::apply(const std::vector< double, numa::allocator< double > > &points, std::vector< double, numa::allocator< double > > &determinants, std::vector< double, numa::allocator< double > > &jacobians) const
  {
    throw common::Exception("Cannot apply jacobians modificator. Unknown operations for type 'double'");
  }

  void JacobiansModificator::apply(const std::vector< float, numa::allocator< float > > &points, std::vector< float, numa::allocator< float > > &determinants, std::vector< float, numa::allocator< float > > &jacobians) const
  {
    throw common::Exception("Cannot apply jacobians modificator. Unknown operations for type 'float'");
  }


} // namespace gmshfem::domain
