// GmshFEM - Copyright (C) 2019-2022, A. Royer, E. Béchet, C. Geuzaine, Université de Liège
//
// See the LICENSE.txt file for license information. Please report all
// issues on https://gitlab.onelab.info/gmsh/fem/issues

#include "GeometricObject.h"

#include "Exception.h"

namespace gmshfem::domain
{


  void GeometricObject::_addEntity(const int dim, const int tag)
  {
    _entities.insert(std::make_pair(dim, tag));
  }

  void GeometricObject::_addEntity(const int dim, const int tag, JacobiansModificator *jacMod)
  {
    _entities.insert(std::make_pair(dim, tag));
    _jacobiansModificators.insert(std::pair< std::pair< int, int >, JacobiansModificator * >(std::make_pair(dim, tag), jacMod));
  }

  GeometricObject::GeometricObject() :
    _entities(), _jacobiansModificators()
  {
  }

  GeometricObject::GeometricObject(const GeometricObject &other) :
    _entities(other._entities), _jacobiansModificators()
  {
    for(auto it = other._jacobiansModificators.begin(); it != other._jacobiansModificators.end(); ++it) {
      _jacobiansModificators.insert(std::pair< std::pair< int, int >, JacobiansModificator * >(it->first, it->second->copy()));
    }
  }

  GeometricObject::GeometricObject(GeometricObject &&other) :
    _entities(std::move(other._entities)), _jacobiansModificators(std::move(other._jacobiansModificators))
  {
  }

  GeometricObject::~GeometricObject()
  {
    for(auto it = _jacobiansModificators.begin(); it != _jacobiansModificators.end(); ++it) {
      delete it->second;
    }
  }

  bool GeometricObject::isEmpty() const
  {
    return (_entities.size() == 0);
  }

  void GeometricObject::clear()
  {
    _entities.clear();
  }

  std::set< std::pair< int, int > >::const_iterator GeometricObject::cbegin() const
  {
    return _entities.cbegin();
  }

  std::set< std::pair< int, int > >::const_iterator GeometricObject::cend() const
  {
    return _entities.cend();
  }

  unsigned int GeometricObject::numberOfEntities() const
  {
    return _entities.size();
  }

  bool GeometricObject::have(const std::pair< int, int > &entity) const
  {
    auto it = _entities.find(entity);
    if(it == _entities.end()) {
      return false;
    }
    return true;
  }

  unsigned int GeometricObject::maxDim() const
  {
    unsigned int dim = 0;
    for(auto it = _entities.begin(); it != _entities.end(); ++it) {
      if(dim < static_cast< unsigned int >(it->first)) {
        dim = it->first;
      }
    }
    return dim;
  }

  unsigned int GeometricObject::minDim() const
  {
    unsigned int dim = 3;
    for(auto it = _entities.begin(); it != _entities.end(); ++it) {
      if(dim > static_cast< unsigned int >(it->first)) {
        dim = it->first;
      }
    }
    return dim;
  }

  size_t GeometricObject::hash() const
  {
    size_t seed = 0;
    for(auto it = _entities.begin(); it != _entities.end(); ++it) {
      seed ^= std::hash< unsigned long long >{}(it->first) + 0x9e3779b9;
      seed ^= std::hash< unsigned long long >{}(it->second) + 0x9e3779b9 + (seed << 6) + (seed >> 2);
    }
    return seed;
  }

  bool GeometricObject::haveJacobiansModificators(const std::pair< int, int > &dimTag) const
  {
    auto it = _jacobiansModificators.find(dimTag);
    if(it == _jacobiansModificators.end()) {
      return false;
    }
    return true;
  }

  bool GeometricObject::needGaussCoordinates() const
  {
    for(auto it = _jacobiansModificators.begin(); it != _jacobiansModificators.end(); ++it) {
      if(it->second->needGaussCoordinates()) {
        return true;
      }
    }
    return false;
  }

  void GeometricObject::applyJacobiansModificator(const std::vector< double, numa::allocator< double > > &points, std::vector< double, numa::allocator< double > > &determinants, std::vector< double, numa::allocator< double > > &jacobians, const std::pair< int, int > &dimTag) const
  {
    auto it = _jacobiansModificators.find(dimTag);
    if(it == _jacobiansModificators.end()) {
      throw common::Exception("Cannot find Jacobian modificator associated to entity (" + std::to_string(dimTag.first) + ", " + std::to_string(dimTag.second) + ")");
    }
    it->second->apply(points, determinants, jacobians);
  }

  void GeometricObject::applyJacobiansModificator(const std::vector< float, numa::allocator< float > > &points, std::vector< float, numa::allocator< float > > &determinants, std::vector< float, numa::allocator< float > > &jacobians, const std::pair< int, int > &dimTag) const
  {
    auto it = _jacobiansModificators.find(dimTag);
    if(it == _jacobiansModificators.end()) {
      throw common::Exception("Cannot find Jacobian modificator associated to entity (" + std::to_string(dimTag.first) + ", " + std::to_string(dimTag.second) + ")");
    }
    it->second->apply(points, determinants, jacobians);
  }


} // namespace gmshfem::domain
