// GmshFEM - Copyright (C) 2019-2022, A. Royer, E. Béchet, C. Geuzaine, Université de Liège
//
// See the LICENSE.txt file for license information. Please report all
// issues on https://gitlab.onelab.info/gmsh/fem/issues

#include "ElementBucket.h"

#include "instantiate.h"

namespace gmshfem::problem
{


  template< class T_PScalar >
  ElementBucket< T_PScalar >::ElementBucket() :
    _determinants(), _jacobians()
  {
  }

  template< class T_PScalar >
  ElementBucket< T_PScalar >::~ElementBucket()
  {
  }

  template< class T_PScalar >
  void ElementBucket< T_PScalar >::clear()
  {
    _determinants.clear();
    _jacobians.clear();
  }

  template< class T_PScalar >
  void ElementBucket< T_PScalar >::determinants(const std::string &gauss, std::vector< T_PScalar, numa::allocator< T_PScalar > > &determinants)
  {
    _determinants[gauss] = std::move(determinants);
  }

  template< class T_PScalar >
  void ElementBucket< T_PScalar >::jacobians(const std::string &gauss, std::vector< T_PScalar, numa::allocator< T_PScalar > > &jacobians)
  {
    _jacobians[gauss] = std::move(jacobians);
  }

  template< class T_PScalar >
  const std::vector< T_PScalar, numa::allocator< T_PScalar > > *ElementBucket< T_PScalar >::determinants(const std::string &gauss) const
  {
    auto it = _determinants.find(gauss);
    if(it != _determinants.end()) {
      return &(it->second);
    }
    return nullptr;
  }

  template< class T_PScalar >
  const std::vector< T_PScalar, numa::allocator< T_PScalar > > *ElementBucket< T_PScalar >::jacobians(const std::string &gauss) const
  {
    auto it = _jacobians.find(gauss);
    if(it != _jacobians.end()) {
      return &(it->second);
    }
    return nullptr;
  }

  template< class T_PScalar >
  common::Memory ElementBucket< T_PScalar >::memory() const
  {
    common::Memory memory;
    for(auto it = _determinants.begin(); it != _determinants.end(); ++it) {
      memory += it->second.size() * sizeof(T_PScalar);
    }

    for(auto it = _jacobians.begin(); it != _jacobians.end(); ++it) {
      memory += it->second.size() * sizeof(T_PScalar);
    }

    return memory;
  }

  INSTANTIATE_CLASS(ElementBucket, 2, TEMPLATE_ARGS(double, float))


} // namespace gmshfem::problem
