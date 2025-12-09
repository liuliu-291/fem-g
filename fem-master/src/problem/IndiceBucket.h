// GmshFEM - Copyright (C) 2019-2022, A. Royer, E. Béchet, C. Geuzaine, Université de Liège
//
// See the LICENSE.txt file for license information. Please report all
// issues on https://gitlab.onelab.info/gmsh/fem/issues

#ifndef H_GMSHFEM_INDICEBUCKET
#define H_GMSHFEM_INDICEBUCKET

#include "Memory.h"

#include <unordered_map>
#include <vector>

namespace gmshfem::problem
{


  class IndiceBucket
  {
   private:
    std::unordered_map< unsigned int, std::vector< std::pair< unsigned long long, int > > > _indices;

   public:
    IndiceBucket();
    ~IndiceBucket();

    void clear();

    IndiceBucket(const IndiceBucket &other) = delete;
    IndiceBucket &operator=(const IndiceBucket &other) = delete;

    void indices(const unsigned int field, std::vector< unsigned long long > &numIndices, std::vector< int > &typeIndices);

    const std::vector< std::pair< unsigned long long, int > > *indices(const unsigned int field) const;

    bool have(const unsigned int field) const;

    common::Memory memory() const;
  };


} // namespace gmshfem::problem

#endif // H_GMSHFEM_INDICEBUCKET
