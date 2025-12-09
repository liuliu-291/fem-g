// GmshFEM - Copyright (C) 2019-2022, A. Royer, E. Béchet, C. Geuzaine, Université de Liège
//
// See the LICENSE.txt file for license information. Please report all
// issues on https://gitlab.onelab.info/gmsh/fem/issues

#include "IndiceBucket.h"

namespace gmshfem::problem
{


  IndiceBucket::IndiceBucket() :
    _indices()
  {
  }

  IndiceBucket::~IndiceBucket()
  {
  }

  void IndiceBucket::clear()
  {
    _indices.clear();
  }

  void IndiceBucket::indices(const unsigned int field, std::vector< unsigned long long > &numIndices, std::vector< int > &typeIndices)
  {
    _indices[field].resize(numIndices.size());
    for(auto i = 0ULL; i < numIndices.size(); ++i) {
      _indices[field][i] = std::pair< unsigned long long, int >(numIndices[i], typeIndices[i]);
    }
  }

  const std::vector< std::pair< unsigned long long, int > > *IndiceBucket::indices(const unsigned int field) const
  {
    auto it = _indices.find(field);
    if(it != _indices.end()) {
      return &(it->second);
    }
    return nullptr;
  }

  bool IndiceBucket::have(const unsigned int field) const
  {
    auto it = _indices.find(field);
    if(it == _indices.end()) return false;
    return true;
  }

  common::Memory IndiceBucket::memory() const
  {
    common::Memory memory;
    for(auto it = _indices.begin(); it != _indices.end(); ++it) {
      memory += it->second.size() * sizeof(std::pair< unsigned int, unsigned int >);
    }

    return memory;
  }


} // namespace gmshfem::problem
