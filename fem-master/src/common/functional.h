// GmshFEM - Copyright (C) 2019-2022, A. Royer, E. Béchet, C. Geuzaine, Université de Liège
//
// See the LICENSE.txt file for license information. Please report all
// issues on https://gitlab.onelab.info/gmsh/fem/issues

#ifndef H_GMSHFEM_FUNCTIONAL
#define H_GMSHFEM_FUNCTIONAL

#include "Dof.h"
#include "FunctionSpaceInterface.h"

namespace gmshfem::common
{


  // hash

  template< class T_Class >
  struct hash {
  };

  template<>
  struct hash< gmshfem::dofs::Dof * > {
    inline size_t operator()(const gmshfem::dofs::Dof *dof) const noexcept
    {
      return dof->hash();
    }
  };

  template<>
  struct hash< const gmshfem::dofs::Dof * > {
    inline size_t operator()(const gmshfem::dofs::Dof *dof) const noexcept
    {
      return dof->hash();
    }
  };

  // equal_to

  template< class T_Class >
  struct equal_to {
  };

  template<>
  struct equal_to< gmshfem::dofs::Dof * > {
    inline bool operator()(const gmshfem::dofs::Dof *x, const gmshfem::dofs::Dof *y) const noexcept
    {
      return *x == *y;
    }
  };

  template<>
  struct equal_to< const gmshfem::dofs::Dof * > {
    inline bool operator()(const gmshfem::dofs::Dof *x, const gmshfem::dofs::Dof *y) const noexcept
    {
      return *x == *y;
    }
  };

  template< class T_PScalar >
  struct equal_to< gmshfem::field::FunctionSpaceInterface< T_PScalar > * > {
    inline bool operator()(const gmshfem::field::FunctionSpaceInterface< T_PScalar > *x, const gmshfem::field::FunctionSpaceInterface< T_PScalar > *y) const noexcept
    {
      return *x == *y;
    }
  };

  // less

  template< class T_Class >
  struct less {
  };

  template< class T_PScalar >
  struct less< gmshfem::field::FunctionSpaceInterface< T_PScalar > * > {
    inline bool operator()(const gmshfem::field::FunctionSpaceInterface< T_PScalar > *x, const gmshfem::field::FunctionSpaceInterface< T_PScalar > *y) const noexcept
    {
      return *x < *y;
    }
  };

  template<>
  struct less< gmshfem::domain::Domain > {
    inline bool operator()(const gmshfem::domain::Domain &x, const gmshfem::domain::Domain &y) const noexcept
    {
      return x.hash() < y.hash();
    }
  };


} // namespace gmshfem::common

#endif // H_GMSHFEM_FUNCTIONAL
