// GmshFEM - Copyright (C) 2019-2022, A. Royer, E. Béchet, C. Geuzaine, Université de Liège
//
// See the LICENSE.txt file for license information. Please report all
// issues on https://gitlab.onelab.info/gmsh/fem/issues

#include "Domain.h"

#include "Exception.h"
#include "Message.h"
#include "OmpInterface.h"
#include "Post.h"
#include "ScalarFunction.h"

#include <algorithm>
#include <gmsh.h>
#include <unordered_set>

namespace gmshfem::domain
{


  Domain::Domain() :
    GeometricObject()
  {
  }

  Domain::Domain(const int dim, const int tag) :
    GeometricObject()
  {
    std::vector< int > nums;
    try {
      gmsh::model::getEntitiesForPhysicalGroup(dim, tag, nums);
    }
    catch(...) {
    }

    for(auto i = 0ULL; i < nums.size(); ++i) {
      _addEntity(dim, nums[i]);
    }
  }

  Domain::Domain(const std::pair< int, int > &dimTag) :
    GeometricObject()
  {
    *this = Domain(dimTag.first, dimTag.second);
  }

  Domain::Domain(const int dim, const int tag, const JacobiansModificator &jacMod) :
    GeometricObject()
  {
    std::vector< int > nums;
    try {
      gmsh::model::getEntitiesForPhysicalGroup(dim, tag, nums);
    }
    catch(...) {
    }

    for(auto i = 0ULL; i < nums.size(); ++i) {
      _addEntity(dim, nums[i], jacMod.copy());
    }
  }

  Domain::Domain(const std::pair< int, int > &dimTag, const JacobiansModificator &jacMod) :
    GeometricObject()
  {
    *this = Domain(dimTag.first, dimTag.second, jacMod);
  }

  Domain::Domain(const std::string &name) :
    GeometricObject()
  {
    std::vector< std::pair< int, int > > dimTags;
    gmsh::model::getEntitiesForPhysicalName(name, dimTags);
    if(dimTags.size() != 0) {
      for(auto dt : dimTags) {
        _addEntity(dt.first, dt.second);
      }
    }
    else {
      msg::warning << "No entity found for domain '" << name << "'" << msg::endl;
    }
  }

  Domain::Domain(const std::string &name, const JacobiansModificator &jacMod) :
    GeometricObject()
  {
    std::vector< std::pair< int, int > > dimTags;
    gmsh::model::getEntitiesForPhysicalName(name, dimTags);
    if(dimTags.size() != 0) {
      for(auto dt : dimTags) {
        _addEntity(dt.first, dt.second, jacMod.copy());
      }
    }
    else {
      msg::warning << "No entity found for domain '" << name << "'" << msg::endl;
    }
  }

  Domain::Domain(const Domain &other) :
    GeometricObject(other)
  {
  }

  Domain::Domain(Domain &&other) :
    GeometricObject(other)
  {
  }

  Domain::~Domain()
  {
  }

  void Domain::addEntity(const int dimEntity, const int tagEntity)
  {
    _entities.insert(std::make_pair(dimEntity, tagEntity));
  }

  Domain &Domain::operator=(const Domain &other)
  {
    _entities = other._entities;
    for(auto it = other._jacobiansModificators.begin(); it != other._jacobiansModificators.end(); ++it) {
      _jacobiansModificators.insert(std::pair< std::pair< int, int >, JacobiansModificator * >(it->first, it->second->copy()));
    }

    return *this;
  }

  Domain &Domain::operator=(Domain &&other)
  {
    _entities = std::move(other._entities);
    _jacobiansModificators = std::move(other._jacobiansModificators);

    return *this;
  }

  bool Domain::operator==(const Domain &other) const
  {
    if(_entities != other._entities) {
      return false;
    }

    return true;
  }

  bool Domain::operator!=(const Domain &other) const
  {
    return !((*this) == other);
  }

  Domain &Domain::operator|=(const Domain &other)
  {
    for(auto it = other._entities.begin(); it != other._entities.end(); ++it) {
      _entities.insert(*it);
    }

    for(auto it = other._jacobiansModificators.begin(); it != other._jacobiansModificators.end(); ++it) {
      _jacobiansModificators.insert(std::pair< std::pair< int, int >, JacobiansModificator * >(it->first, it->second->copy()));
    }

    return *this;
  }

  Domain Domain::operator|(const Domain &other) const
  {
    Domain domain(*this);
    domain |= other;

    return domain;
  }

  Domain &Domain::operator&=(const Domain &other)
  {
    std::set< std::pair< int, int > > newEntities;
    for(auto it = other._entities.begin(); it != other._entities.end(); ++it) {
      auto itFind = _entities.find(*it);
      if(itFind != _entities.end()) {
        newEntities.insert(*it);
      }
    }
    _entities = std::move(newEntities);

    for(auto it = other._jacobiansModificators.begin(); it != other._jacobiansModificators.end(); ++it) {
      auto itFind = _entities.find(it->first);
      if(itFind != _entities.end()) {
        _jacobiansModificators.insert(std::pair< std::pair< int, int >, JacobiansModificator * >(it->first, it->second->copy()));
      }
    }

    return *this;
  }

  Domain Domain::operator&(const Domain &other) const
  {
    Domain domain(*this);
    domain &= other;

    return domain;
  }

  Domain &Domain::operator^=(const Domain &other)
  {
    std::set< std::pair< int, int > > newEntities;
    for(auto it = other._entities.begin(); it != other._entities.end(); ++it) {
      auto itFind = _entities.find(*it);
      if(itFind == _entities.end()) {
        newEntities.insert(*it);
      }
    }
    for(auto it = _entities.begin(); it != _entities.end(); ++it) {
      auto itFind = other._entities.find(*it);
      if(itFind == other._entities.end()) {
        newEntities.insert(*it);
      }
    }
    _entities = std::move(newEntities);

    for(auto it = other._jacobiansModificators.begin(); it != other._jacobiansModificators.end(); ++it) {
      auto itFind = _entities.find(it->first);
      if(itFind != _entities.end()) {
        _jacobiansModificators.insert(std::pair< std::pair< int, int >, JacobiansModificator * >(it->first, it->second->copy()));
      }
    }

    return *this;
  }

  Domain Domain::operator^(const Domain &other) const
  {
    Domain domain(*this);
    domain ^= other;

    return domain;
  }

  Domain Domain::operator~() const
  {
    Domain domain;
    std::vector< std::pair< int, int > > entities;
    gmsh::model::getEntities(entities);
    for(auto i = 0ULL; i < entities.size(); ++i) {
      auto it = _entities.find(entities[i]);
      if(it == _entities.end()) {
        domain._entities.insert(entities[i]);
      }
    }

    return domain;
  }

  Domain Domain::getBoundary(const bool combined) const
  {
    std::vector< std::pair< int, int > > dimTags;
    dimTags.reserve(_entities.size());
    for(auto it = _entities.begin(); it != _entities.end(); ++it) {
      dimTags.push_back(*it);
    }
    std::vector< std::pair< int, int > > outDimTags;
    gmsh::model::getBoundary(dimTags, outDimTags, combined, false, false);

    Domain boundary;
    for(auto i = 0ULL; i < outDimTags.size(); ++i) {
      boundary._entities.insert(outDimTags[i]);
    }
    return boundary;
  }

  SkinLayer Domain::getSkinLayer(const Domain &tool) const
  {
    unsigned int objectDim = maxDim();
    unsigned int toolDim = tool.maxDim();
    if(objectDim != minDim()) {
      msg::error << "Unable to get the skin layer of a domain made of entities of different dimension" << msg::endl;
      return SkinLayer();
    }
    if(toolDim != tool.minDim()) {
      msg::error << "Unable to get the skin layer along a tool made of entities of different dimension" << msg::endl;
      return SkinLayer();
    }

    std::unordered_set< std::size_t > toolNodes;

    std::vector< std::size_t > nodeTags;
    std::vector< double > coord;
    std::vector< double > parametricCoord;
    for(auto it = tool.cbegin(); it != tool.cend(); ++it) {
      gmsh::model::mesh::getNodes(nodeTags, coord, parametricCoord, it->first, it->second, true, false);
      toolNodes.reserve(toolNodes.size() + nodeTags.size());
      for(unsigned int i = 0; i < nodeTags.size(); ++i) {
        toolNodes.insert(nodeTags[i]);
      }
    }

    // create a new discrete entity
    int newEntityTag = gmsh::model::addDiscreteEntity(objectDim);

    std::vector< std::size_t > layerElementTags;
    std::vector< std::size_t > layerNodeTags;

    std::vector< int > elementTypes;
    std::vector< std::size_t > elementTags;
    for(auto it = cbegin(); it != cend(); ++it) {
      gmsh::model::mesh::getElementTypes(elementTypes, it->first, it->second);
      for(unsigned int iType = 0; iType < elementTypes.size(); ++iType) {
#pragma omp parallel num_threads(omp::getMaxThreads())
        {
          const unsigned int numThreads = omp::getNumThreads();
          const unsigned int myThreadID = omp::getThreadNum();

#pragma omp single
          gmsh::model::mesh::preallocateElementsByType(elementTypes[iType], true, true, elementTags, nodeTags, it->second);
          gmsh::model::mesh::getElementsByType(elementTypes[iType], elementTags, nodeTags, it->second, myThreadID, numThreads);
          const unsigned int nodeByElement = nodeTags.size() / elementTags.size();
#pragma omp for
          for(unsigned int i = 0; i < elementTags.size(); ++i) {
            for(unsigned int j = 0; j < nodeByElement; ++j) {
              auto itFind = toolNodes.find(nodeTags[nodeByElement * i + j]);
              if(itFind != toolNodes.end()) {
#pragma omp critical
                {
                  layerElementTags.push_back(elementTags[i]);
                  layerNodeTags.reserve(layerNodeTags.size() + nodeByElement);
                  for(unsigned int k = 0; k < nodeByElement; ++k) {
                    layerNodeTags.push_back(nodeTags[nodeByElement * i + k]);
                  }
                }
                break;
              }
            }
          }
        }

        // add to the new entity
        gmsh::model::mesh::addElementsByType(newEntityTag, elementTypes[iType], layerElementTags, layerNodeTags);
        layerElementTags.clear();
        layerNodeTags.clear();
      }
    }

    return SkinLayer(objectDim, newEntityTag, *this, tool, _jacobiansModificators);
  }

  bool Domain::isIncludedInto(const Domain &other) const
  {
    for(auto it = _entities.begin(); it != _entities.end(); ++it) {
      auto itFind = other._entities.find(*it);
      if(itFind == other._entities.end()) {
        return false;
      }
    }

    return true;
  }

  bool Domain::included(const Domain &other) const
  {
    return other.isIncludedInto(*this);
  }

  void Domain::printDebug() const
  {
    msg::debug << "Domain has:" << msg::endl;
    msg::debug << " * elementary entities (" << _entities.size() << "):" << msg::endl;
    for(auto it = _entities.begin(); it != _entities.end(); ++it) {
      msg::debug << "  - " << it->first << " / " << it->second << msg::endl;
    }
  }

  void Domain::saveDebug(const std::string &name) const
  {
    post::save(function::ScalarFunction< double >(1.), *this, "Domain debug" + (name == "" ? "" : " : " + name), "pos");
  }


} // namespace gmshfem::domain
