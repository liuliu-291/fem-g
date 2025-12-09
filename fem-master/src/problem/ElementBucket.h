// GmshFEM - Copyright (C) 2019-2022, A. Royer, E. Béchet, C. Geuzaine, Université de Liège
//
// See the LICENSE.txt file for license information. Please report all
// issues on https://gitlab.onelab.info/gmsh/fem/issues

#ifndef H_GMSHFEM_ELEMENTBUCKET
#define H_GMSHFEM_ELEMENTBUCKET

#include "Memory.h"
#include "numa.h"

#include <string>
#include <unordered_map>
#include <vector>

namespace gmshfem::problem
{


  template< class T_PScalar >
  class ElementBucket
  {
   private:
    std::unordered_map< std::string, std::vector< T_PScalar, numa::allocator< T_PScalar > > > _determinants;
    std::unordered_map< std::string, std::vector< T_PScalar, numa::allocator< T_PScalar > > > _jacobians;

   public:
    ElementBucket();
    ~ElementBucket();

    void clear();

    ElementBucket(const ElementBucket &other) = delete;
    ElementBucket &operator=(const ElementBucket &other) = delete;

    void determinants(const std::string &gauss, std::vector< T_PScalar, numa::allocator< T_PScalar > > &determinants);
    void jacobians(const std::string &gauss, std::vector< T_PScalar, numa::allocator< T_PScalar > > &jacobians);

    const std::vector< T_PScalar, numa::allocator< T_PScalar > > *determinants(const std::string &gauss) const;
    const std::vector< T_PScalar, numa::allocator< T_PScalar > > *jacobians(const std::string &gauss) const;

    common::Memory memory() const;
  };


} // namespace gmshfem::problem

#endif // H_GMSHFEM_ELEMENTBUCKET
