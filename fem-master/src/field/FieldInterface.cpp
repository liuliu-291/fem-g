// GmshFEM - Copyright (C) 2019-2022, A. Royer, E. Béchet, C. Geuzaine, Université de Liège
//
// See the LICENSE.txt file for license information. Please report all
// issues on https://gitlab.onelab.info/gmsh/fem/issues

#include "FieldInterface.h"

#include "Dof.h"
#include "Exception.h"
#include "IndiceBucket.h"
#include "Message.h"
#include "OmpInterface.h"
#include "Options.h"
#include "gmshfemDefines.h"
#include "instantiate.h"

#include <complex>
#include <gmsh.h>
#include <map>

namespace gmshfem::field
{


  //
  // class DofPool
  //

  template< class T_Scalar >
  FieldInterface< T_Scalar >::DofPool::DofPool() :
    _unknownDofPool(), _fixedDofPool(), _linkedDofPool(), _unknownDofPos(GMSHFEM_DOF_MEMORY_POOL_SIZE), _fixedDofPos(GMSHFEM_DOF_MEMORY_POOL_SIZE), _linkedDofPos(GMSHFEM_DOF_MEMORY_POOL_SIZE)
  {
  }

  template< class T_Scalar >
  FieldInterface< T_Scalar >::DofPool::~DofPool()
  {
    clear();
  }

  template< class T_Scalar >
  void FieldInterface< T_Scalar >::DofPool::clear()
  {
    clearUnknown();
    clearFixed();
    clearLinked();
  }

  template< class T_Scalar >
  void FieldInterface< T_Scalar >::DofPool::clearUnknown()
  {
    for(auto i = 0ULL; i < _unknownDofPool.size(); ++i) {
      ::operator delete(_unknownDofPool[i]);
    }
    _unknownDofPool.clear();
    _unknownDofPos = GMSHFEM_DOF_MEMORY_POOL_SIZE;
  }

  template< class T_Scalar >
  void FieldInterface< T_Scalar >::DofPool::clearFixed()
  {
    for(auto i = 0ULL; i < _fixedDofPool.size(); ++i) {
      ::operator delete(_fixedDofPool[i]);
    }
    _fixedDofPool.clear();
    _fixedDofPos = GMSHFEM_DOF_MEMORY_POOL_SIZE;
  }

  template< class T_Scalar >
  void FieldInterface< T_Scalar >::DofPool::clearLinked()
  {
    for(auto i = 0ULL; i < _linkedDofPool.size(); ++i) {
      ::operator delete(_linkedDofPool[i]);
    }
    _linkedDofPool.clear();
    _linkedDofPos = GMSHFEM_DOF_MEMORY_POOL_SIZE;
  }

  template< class T_Scalar >
  dofs::UnknownDof *FieldInterface< T_Scalar >::DofPool::getNextUnknownDofMemoryPlace()
  {
    if(_unknownDofPos == GMSHFEM_DOF_MEMORY_POOL_SIZE) {
      _unknownDofPool.push_back((dofs::UnknownDof *)::operator new(GMSHFEM_DOF_MEMORY_POOL_SIZE * sizeof(dofs::UnknownDof)));
      _unknownDofPos = 0;
    }
    return &_unknownDofPool.back()[_unknownDofPos++];
  }

  template< class T_Scalar >
  dofs::FixedDof *FieldInterface< T_Scalar >::DofPool::getNextFixedDofMemoryPlace()
  {
    if(_fixedDofPos == GMSHFEM_DOF_MEMORY_POOL_SIZE) {
      _fixedDofPool.push_back((dofs::FixedDof *)::operator new(GMSHFEM_DOF_MEMORY_POOL_SIZE * sizeof(dofs::FixedDof)));
      _fixedDofPos = 0;
    }
    return &_fixedDofPool.back()[_fixedDofPos++];
  }

  template< class T_Scalar >
  dofs::LinkedDof *FieldInterface< T_Scalar >::DofPool::getNextLinkedDofMemoryPlace()
  {
    if(_linkedDofPos == GMSHFEM_DOF_MEMORY_POOL_SIZE) {
      _linkedDofPool.push_back((dofs::LinkedDof *)::operator new(GMSHFEM_DOF_MEMORY_POOL_SIZE * sizeof(dofs::LinkedDof)));
      _linkedDofPos = 0;
    }
    return &_linkedDofPool.back()[_linkedDofPos++];
  }

  template< class T_Scalar >
  void FieldInterface< T_Scalar >::DofPool::lastUnknownDofIsUnused()
  {
    --_unknownDofPos;
  }

  template< class T_Scalar >
  void FieldInterface< T_Scalar >::DofPool::lastFixedDofIsUnused()
  {
    --_fixedDofPos;
  }

  template< class T_Scalar >
  void FieldInterface< T_Scalar >::DofPool::lastLinkedDofIsUnused()
  {
    --_linkedDofPos;
  }

  template< class T_Scalar >
  common::Memory FieldInterface< T_Scalar >::DofPool::memory() const
  {
    return common::Memory((_unknownDofPool.size() * sizeof(dofs::UnknownDof) + _fixedDofPool.size() * sizeof(dofs::FixedDof) + _linkedDofPool.size() * sizeof(dofs::LinkedDof)) * GMSHFEM_DOF_MEMORY_POOL_SIZE);
  }

  //
  // class FieldInterface
  //

  static unsigned int s_nbrField = 0;
  static std::vector< bool > s_listOfField;

  template< class T_Scalar >
  void FieldInterface< T_Scalar >::s_incrementNbrOfFields()
  {
    s_listOfField.resize(++s_nbrField);
    s_listOfField[s_nbrField - 1] = true;
  }

  template< class T_Scalar >
  FieldInterface< T_Scalar >::FieldInterface() :
    _functionSpace(nullptr), _name(""), _model(""), _tag(s_nbrField++), _domain(), _values(), _pool(), _globalQuantities(), _modules(), _periodicConstraints(), _orderedDofs()
  {
    gmsh::model::getCurrent(_model);
    s_listOfField.resize(s_nbrField);
    s_listOfField[_tag] = true;
  }

  template< class T_Scalar >
  FieldInterface< T_Scalar >::FieldInterface(const std::string &name, const domain::Domain &domain, const std::string &model) :
    _functionSpace(nullptr), _name(name), _model(model), _tag(s_nbrField++), _domain(domain), _values(), _pool(), _globalQuantities(), _modules(), _periodicConstraints(), _orderedDofs()
  {
    if(_model == "") {
      gmsh::model::getCurrent(_model);
    }
    s_listOfField.resize(s_nbrField);
    s_listOfField[_tag] = true;
  }

  template< class T_Scalar >
  FieldInterface< T_Scalar >::FieldInterface(const FieldInterface< T_Scalar > &other) :
    _functionSpace(nullptr), _name(other._name), _model(other._model), _tag(s_nbrField++), _domain(other._domain), _values(), _pool(), _globalQuantities(), _modules(), _periodicConstraints(), _orderedDofs()
  {
    s_listOfField.resize(s_nbrField);
    s_listOfField[_tag] = true;
    _copy(other);
  }

  template< class T_Scalar >
  void FieldInterface< T_Scalar >::_copy(const FieldInterface< T_Scalar > &other)
  {
    _pool.clear();
    _functionSpace = (other._functionSpace == nullptr ? nullptr : other._functionSpace->copy());
    _domain = other._domain;

    _pool.clear();
    _values.clear();
    bool hasLinkedDof = false;
    for(auto it = other._values.begin(); it != other._values.end(); ++it) {
      dofs::Dof *dof = nullptr;
      if(it->first->assembleType() == dofs::AssembleType::Fixed) {
        dof = new(_pool.getNextFixedDofMemoryPlace()) dofs::FixedDof(it->first->numType() - GMSHFEM_DOF_FIELD_OFFSET * other._tag + GMSHFEM_DOF_FIELD_OFFSET * _tag, it->first->entity());
      }
      else if(it->first->assembleType() == dofs::AssembleType::Unknown) {
        dof = new(_pool.getNextUnknownDofMemoryPlace()) dofs::UnknownDof(it->first->numType() - GMSHFEM_DOF_FIELD_OFFSET * other._tag + GMSHFEM_DOF_FIELD_OFFSET * _tag, it->first->entity());
      }
      else if(it->first->assembleType() == dofs::AssembleType::UnknownBubble) {
        dof = new(_pool.getNextUnknownDofMemoryPlace()) dofs::UnknownBubbleDof(it->first->numType() - GMSHFEM_DOF_FIELD_OFFSET * other._tag + GMSHFEM_DOF_FIELD_OFFSET * _tag, it->first->entity());
      }
      else if(it->first->assembleType() == dofs::AssembleType::UnknownGlobal) {
        dof = new(_pool.getNextUnknownDofMemoryPlace()) dofs::UnknownGlobalDof(it->first->numType() - GMSHFEM_DOF_FIELD_OFFSET * other._tag + GMSHFEM_DOF_FIELD_OFFSET * _tag, it->first->entity());
      }
      else if(it->first->assembleType() == dofs::AssembleType::Linked || it->first->assembleType() == dofs::AssembleType::LinkedBubble) {
        hasLinkedDof = true;
      }
      dof->coordinates(it->first->x(), it->first->y(), it->first->z());
      dof->numDof(it->first->numDof());
      _values.insert(hashmap::pair< dofs::Dof *, T_Scalar >(dof, it->second));
    }
    if(hasLinkedDof) {
      for(auto it = other._values.begin(); it != other._values.end(); ++it) {
        dofs::LinkedDof *dof = nullptr;
        if(it->first->assembleType() == dofs::AssembleType::Linked) {
          const dofs::Dof *master = searchDof(static_cast< const dofs::LinkedDof * >(it->first)->master()->numType() - GMSHFEM_DOF_FIELD_OFFSET * other._tag + GMSHFEM_DOF_FIELD_OFFSET * _tag, static_cast< const dofs::LinkedDof * >(it->first)->master()->entity());
          if(master != nullptr) {
            dof = new(_pool.getNextLinkedDofMemoryPlace()) dofs::LinkedDof(it->first->numType() - GMSHFEM_DOF_FIELD_OFFSET * other._tag + GMSHFEM_DOF_FIELD_OFFSET * _tag, it->first->entity(), master);
          }
          else {
            msg::error << "Unable to find the master dof of a linked dof" << msg::endl;
          }
        }
        else if(it->first->assembleType() == dofs::AssembleType::LinkedBubble) {
          const dofs::Dof *master = searchDof(static_cast< const dofs::LinkedDof * >(it->first)->master()->numType() - GMSHFEM_DOF_FIELD_OFFSET * other._tag + GMSHFEM_DOF_FIELD_OFFSET * _tag, static_cast< const dofs::LinkedDof * >(it->first)->master()->entity());
          if(master != nullptr) {
            dof = new(_pool.getNextLinkedDofMemoryPlace()) dofs::LinkedBubbleDof(it->first->numType() - GMSHFEM_DOF_FIELD_OFFSET * other._tag + GMSHFEM_DOF_FIELD_OFFSET * _tag, it->first->entity(), master);
          }
          else {
            msg::error << "Unable to find the master dof of a linked dof" << msg::endl;
          }
        }
        dof->coordinates(it->first->x(), it->first->y(), it->first->z());
        dof->numDof(it->first->numDof());
        _values.insert(hashmap::pair< dofs::Dof *, T_Scalar >(dof, it->second));
      }
    }
    for(auto i = 0ULL; i < other._modules.size(); ++i) {
      _modules.push_back(other._modules[i]->copy());
    }
    _periodicConstraints.clear();
    for(auto it = other._periodicConstraints.begin(); it != other._periodicConstraints.end(); ++it) {
      PeriodicConstraint< T_Scalar > constraint(it->_link, it->_bondageCoefficient);
      _periodicConstraints.push_back(constraint);
    }
  }

  template< class T_Scalar >
  FieldInterface< T_Scalar >::~FieldInterface()
  {
    clear();
    s_listOfField[_tag] = false;
  }

  template< class T_Scalar >
  void FieldInterface< T_Scalar >::clear()
  {
    if(_functionSpace) {
      delete _functionSpace;
      _functionSpace = nullptr;
    }
    _values.clear();
    _pool.clear();
    for(auto i = 0ULL; i < _modules.size(); ++i) {
      delete _modules[i];
    }
    _periodicConstraints.clear();
    if(!_globalQuantities.empty()) {
      FieldInterface< T_Scalar > *dual = nullptr;
      for(auto i = 0ULL; i < _globalQuantities.size(); ++i) {
        if(globalValueIsStillValid(_globalQuantities[i].first)) {
          dual = _globalQuantities[0].second->getAssociatedDualField();
          _globalQuantities[i].second->setAssociatedPrimalField(nullptr);
          _globalQuantities[i].second->setAssociatedDualField(nullptr);
        }
      }
      if(dual != nullptr) {
        delete dual;
      }
      _globalQuantities.clear();
    }
  }

  template< class T_Scalar >
  void FieldInterface< T_Scalar >::clearDofs()
  {
    _values.clear();
    _pool.clear();
  }

  template< class T_Scalar >
  typename hashmap::unordered_map< dofs::Dof *, T_Scalar, common::hash< dofs::Dof * >, common::equal_to< dofs::Dof * > >::iterator FieldInterface< T_Scalar >::begin()
  {
    return _values.begin();
  }

  template< class T_Scalar >
  typename hashmap::unordered_map< dofs::Dof *, T_Scalar, common::hash< dofs::Dof * >, common::equal_to< dofs::Dof * > >::iterator FieldInterface< T_Scalar >::end()
  {
    return _values.end();
  }

  template< class T_Scalar >
  typename hashmap::unordered_map< dofs::Dof *, T_Scalar, common::hash< dofs::Dof * >, common::equal_to< dofs::Dof * > >::const_iterator FieldInterface< T_Scalar >::begin() const
  {
    return _values.cbegin();
  }

  template< class T_Scalar >
  typename hashmap::unordered_map< dofs::Dof *, T_Scalar, common::hash< dofs::Dof * >, common::equal_to< dofs::Dof * > >::const_iterator FieldInterface< T_Scalar >::end() const
  {
    return _values.cend();
  }

  template< class T_Scalar >
  typename hashmap::unordered_map< dofs::Dof *, T_Scalar, common::hash< dofs::Dof * >, common::equal_to< dofs::Dof * > >::iterator FieldInterface< T_Scalar >::findDof(dofs::Dof *dof)
  {
    return _values.find(dof);
  }

  template< class T_Scalar >
  typename std::vector< std::pair< unsigned int, GlobalQuantity< T_Scalar > * > >::const_iterator FieldInterface< T_Scalar >::firstGlobalQuantity() const
  {
    return _globalQuantities.cbegin();
  }

  template< class T_Scalar >
  typename std::vector< std::pair< unsigned int, GlobalQuantity< T_Scalar > * > >::const_iterator FieldInterface< T_Scalar >::lastGlobalQuantity() const
  {
    return _globalQuantities.cend();
  }

  template< class T_Scalar >
  void FieldInterface< T_Scalar >::addPeriodicConstraint(const domain::PeriodicLink &link, const T_Scalar &coefficient)
  {
    PeriodicConstraint< T_Scalar > constraint(link, coefficient);
    _periodicConstraints.push_back(constraint);
  }

  template< class T_Scalar >
  std::vector< PeriodicConstraint< T_Scalar > > FieldInterface< T_Scalar >::getPeriodicConstraints() const
  {
    return _periodicConstraints;
  }

  template< class T_Scalar >
  std::string FieldInterface< T_Scalar >::name() const
  {
    return _name;
  }

  template< class T_Scalar >
  std::string FieldInterface< T_Scalar >::model() const
  {
    return _model;
  }

  template< class T_Scalar >
  void FieldInterface< T_Scalar >::name(const std::string &name)
  {
    _name = name;
  }

  template< class T_Scalar >
  unsigned int FieldInterface< T_Scalar >::tag() const
  {
    return _tag;
  }

  template< class T_Scalar >
  domain::Domain FieldInterface< T_Scalar >::domain() const
  {
    return _domain;
  }

  template< class T_Scalar >
  void FieldInterface< T_Scalar >::domain(const domain::Domain &domain)
  {
    _domain = domain;
    _checkDomain();
  }

  template< class T_Scalar >
  common::Memory FieldInterface< T_Scalar >::memory() const
  {
    common::Memory mem = _pool.memory() + (_values.size() * sizeof(hashmap::pair< dofs::Dof *, T_Scalar >));
    for(auto i = 0ULL; i < _modules.size(); ++i) {
      mem += _modules[i]->memory();
    }
    for(auto it = _orderedDofs.begin(); it != _orderedDofs.end(); ++it) {
      mem += it->second.size() * sizeof(dofs::Dof *);
    }
    return mem;
  }

  template< class T_Scalar >
  FieldModule< T_Scalar > *FieldInterface< T_Scalar >::getModule(const std::string &name) const
  {
    for(auto i = 0ULL; i < _modules.size(); ++i) {
      if(_modules[i]->name() == name) {
        return _modules[i];
      }
    }
    return nullptr;
  }

  template< class T_Scalar >
  void FieldInterface< T_Scalar >::setModule(FieldModule< T_Scalar > *module)
  {
    _modules.push_back(module);
  }

  template< class T_Scalar >
  unsigned long long FieldInterface< T_Scalar >::numberOfDofs() const
  {
    return _values.size();
  }

  template< class T_Scalar >
  unsigned long long FieldInterface< T_Scalar >::numberOfUnknownDofs() const
  {
    unsigned long long count = 0;
    for(auto it = _values.begin(); it != _values.end(); ++it) {
      if(it->first->type() == dofs::Type::Unknown) {
        count++;
      }
    }
    return count;
  }

  template< class T_Scalar >
  unsigned long long FieldInterface< T_Scalar >::numberOfFixedDofs() const
  {
    unsigned long long count = 0;
    for(auto it = _values.begin(); it != _values.end(); ++it) {
      if(it->first->type() == dofs::Type::Fixed) {
        count++;
      }
    }
    return count;
  }

  template< class T_Scalar >
  unsigned long long FieldInterface< T_Scalar >::numberOfLinkedDofs() const
  {
    unsigned long long count = 0;
    for(auto it = _values.begin(); it != _values.end(); ++it) {
      if(it->first->type() == dofs::Type::Linked) {
        count++;
      }
    }
    return count;
  }

  template< class T_Scalar >
  void FieldInterface< T_Scalar >::reserve(const unsigned long long size)
  {
    _values.reserve(size);
  }

  template< class T_Scalar >
  bool FieldInterface< T_Scalar >::setDof(dofs::UnknownDof *dof)
  {
    auto out = _values.insert(hashmap::pair< dofs::Dof *, T_Scalar >(dof, 0.));
    if(!out.second) {
      _pool.lastUnknownDofIsUnused();
    }
    return out.second;
  }

  template< class T_Scalar >
  bool FieldInterface< T_Scalar >::setDof(dofs::LinkedDof *dof, const T_Scalar &value)
  {
    auto out = _values.insert(hashmap::pair< dofs::Dof *, T_Scalar >(dof, value));
    if(!out.second) {
      _pool.lastLinkedDofIsUnused();
    }
    return out.second;
  }

  template< class T_Scalar >
  bool FieldInterface< T_Scalar >::setValue(dofs::FixedDof *dof, const T_Scalar &value)
  {
    auto out = _values.insert(hashmap::pair< dofs::Dof *, T_Scalar >(dof, value));
    if(!out.second) {
      _pool.lastFixedDofIsUnused();
    }
    return out.second;
  }

  template< class T_Scalar >
  T_Scalar FieldInterface< T_Scalar >::getValue(const dofs::Dof *dof) const
  {
    auto out = _values.find(const_cast< dofs::Dof * >(dof));
    if(out != _values.end()) {
      if(out->first->type() == dofs::Type::Linked) {
        return out->second * getValue(static_cast< const dofs::LinkedDof * >(out->first)->master());
      }
      return out->second;
    }
    return 0.;
  }

  template< class T_Scalar >
  const dofs::Dof *FieldInterface< T_Scalar >::searchDof(const int numType, const unsigned long long entity, const unsigned int multiplicityIndex) const
  {
    dofs::SearchDof dof(numType + GMSHFEM_DOF_FIELD_OFFSET * (_tag + multiplicityIndex), entity);
    auto it = _values.find(&dof);
    if(it != _values.end()) {
      return it->first;
    }
    return nullptr;
  }

  template< class T_Scalar >
  dofs::UnknownDof *FieldInterface< T_Scalar >::getNextUnknownDofMemoryPlace()
  {
    return _pool.getNextUnknownDofMemoryPlace();
  }

  template< class T_Scalar >
  dofs::FixedDof *FieldInterface< T_Scalar >::getNextFixedDofMemoryPlace()
  {
    return _pool.getNextFixedDofMemoryPlace();
  }

  template< class T_Scalar >
  dofs::LinkedDof *FieldInterface< T_Scalar >::getNextLinkedDofMemoryPlace()
  {
    return _pool.getNextLinkedDofMemoryPlace();
  }

  template< class T_Scalar >
  void FieldInterface< T_Scalar >::setValuesToZero()
  {
    for(auto it = _values.begin(); it != _values.end(); ++it) {
      it->second = 0.;
    }
  }

  // Internal functions
  template< class T_Scalar >
  void FieldInterface< T_Scalar >::assignValues(const std::vector< T_Scalar > &values)
  {
    for(auto it = _values.begin(); it != _values.end(); ++it) {
      if(it->first->type() == dofs::Type::Unknown) {
        it->second = values[it->first->numDof() - 1];
      }
    }
  }

  template< class T_Scalar >
  void FieldInterface< T_Scalar >::assignValueTo(const unsigned long long numDof, const dofs::Type &type, const T_Scalar &value)
  {
    for(auto it = _values.begin(); it != _values.end(); ++it) {
      if(it->first->type() == type) {
        if(it->first->numDof() - 1 == numDof) {
          it->second = value;
        }
      }
    }
  }

  template< class T_Scalar >
  void FieldInterface< T_Scalar >::getValues(const std::vector< int > &typeKeys, const std::vector< unsigned long long > &entityKeys, std::vector< T_Scalar > &values, const unsigned long long begin, const unsigned long long end) const
  {
    const unsigned int multi = multiplicity();
    for(auto i = begin; i < end; ++i) {
      for(unsigned int j = 0; j < multi; ++j) {
        dofs::SearchDof dof(typeKeys[i] + GMSHFEM_DOF_FIELD_OFFSET * (_tag + j), entityKeys[i]);
        auto it = _values.find(&dof);
        if(it != _values.end()) {
          if(it->first->type() == dofs::Type::Linked) { // slave dof
            values[i * multi + j] = it->second * getValue(static_cast< const dofs::LinkedDof * >(it->first)->master());
          }
          else {
            values[i * multi + j] = it->second;
          }
        }
        else {
          values[i * multi + j] = 0.;
        }
      }
    }
  }

  template< class T_Scalar >
  void FieldInterface< T_Scalar >::getUnknownValues(std::vector< T_Scalar > &values)
  {
    values.reserve(_values.size());
    for(auto it = _values.begin(); it != _values.end(); ++it) {
      if(it->first->type() == dofs::Type::Unknown) {
        values[it->first->numDof() - 1] = it->second;
      }
    }
  }

  template< class T_Scalar >
  void FieldInterface< T_Scalar >::getFixedValues(std::vector< T_Scalar > &values) const
  {
    for(auto it = _values.begin(); it != _values.end(); ++it) {
      if(it->first->type() == dofs::Type::Fixed) {
        values[it->first->numDof() - 1] = it->second;
      }
    }
  }

  template< class T_Scalar >
  void FieldInterface< T_Scalar >::getLinkedValues(std::vector< T_Scalar > &values) const
  {
    for(auto it = _values.begin(); it != _values.end(); ++it) {
      if(it->first->type() == dofs::Type::Linked) {
        values[it->first->numDof() - 1] = it->second;
      }
    }
  }

  template< class T_Scalar >
  void FieldInterface< T_Scalar >::getLinkedIndices(std::vector< unsigned long long > &indices) const
  {
    for(auto it = _values.begin(); it != _values.end(); ++it) {
      if(it->first->type() == dofs::Type::Linked) {
        indices[it->first->numDof() - 1] = static_cast< const dofs::LinkedDof * >(it->first)->master()->numDof() - 1;
      }
    }
  }

  // External functions
  template< class T_Scalar >
  void FieldInterface< T_Scalar >::getAllVector(algebra::Vector< T_Scalar > &vector, const dofs::RawOrder order) const
  {
    if(order == dofs::RawOrder::Hash) {
      std::vector< T_Scalar > values;
      values.reserve(numberOfDofs());
      for(auto it = _values.begin(); it != _values.end(); ++it) {
        values.push_back(it->second);
      }
      vector = std::move(values);
    }
    else {
      auto less_num = [](const dofs::Dof *lhs, const dofs::Dof *rhs) { return lhs->less_num(*rhs); };
      auto less_pair = [](const dofs::Dof *lhs, const dofs::Dof *rhs) { return lhs->less_pair(*rhs); };
      std::map< dofs::Dof *, T_Scalar, decltype(order == dofs::RawOrder::System ? less_num : less_pair) > sortedMap(order == dofs::RawOrder::System ? less_num : less_pair);
      for(auto it = _values.begin(); it != _values.end(); ++it) {
        sortedMap.insert(std::make_pair(it->first, it->second));
      }

      std::vector< T_Scalar > values;
      values.reserve(sortedMap.size());
      for(auto it = sortedMap.begin(); it != sortedMap.end(); ++it) {
        values.push_back(it->second);
      }
      vector = std::move(values);
    }
  }

  template< class T_Scalar >
  void FieldInterface< T_Scalar >::setAllVector(const algebra::Vector< T_Scalar > &vector, const dofs::RawOrder order)
  {
    const std::vector< T_Scalar > &valuesVector = vector.getStdVector();
    if(order == dofs::RawOrder::Hash) {
      unsigned long long i = 0;
      for(auto it = _values.begin(); it != _values.end(); ++it) {
        it->second = valuesVector[i];
        ++i;
      }
    }
    else {
      auto less_num = [](const dofs::Dof *lhs, const dofs::Dof *rhs) { return lhs->less_num(*rhs); };
      auto less_pair = [](const dofs::Dof *lhs, const dofs::Dof *rhs) { return lhs->less_pair(*rhs); };
      std::map< dofs::Dof *, T_Scalar *, decltype(order == dofs::RawOrder::System ? less_num : less_pair) > sortedMap(order == dofs::RawOrder::System ? less_num : less_pair);
      for(auto it = _values.begin(); it != _values.end(); ++it) {
        sortedMap.insert(std::make_pair(it->first, &it->second));
      }

      unsigned long long i = 0;
      for(auto it = sortedMap.begin(); it != sortedMap.end(); ++it) {
        *it->second = valuesVector[i];
        ++i;
      }
    }
  }

  template< class T_Scalar >
  void FieldInterface< T_Scalar >::getAllDofVector(std::vector< dofs::RawDof > &vector) const
  {
    std::vector< dofs::RawDof > dofs;
    dofs.reserve(_values.size());
    for(auto it = _values.begin(); it != _values.end(); ++it) {
      dofs.push_back(it->first->getRawDof());
    }
    vector = std::move(dofs);
  }

  template< class T_Scalar >
  void FieldInterface< T_Scalar >::setAllDofVector(const std::vector< dofs::RawDof > &vector, const bool append)
  {
    if(!append) {
      _pool.clear();
      _values.clear();
    }
    for(auto i = 0ULL; i < vector.size(); ++i) {
      dofs::Dof *dof;
      if(vector[i].assembleType == dofs::AssembleType::Fixed) {
        dof = new(_pool.getNextFixedDofMemoryPlace()) dofs::FixedDof(vector[i].numType, vector[i].entity);
        dof->coordinates(vector[i].coordinates[0], vector[i].coordinates[1], vector[i].coordinates[2]);
        dof->numDof(vector[i].numDof);
      }
      else if(vector[i].assembleType == dofs::AssembleType::Unknown) {
        dof = new(_pool.getNextUnknownDofMemoryPlace()) dofs::UnknownDof(vector[i].numType, vector[i].entity);
        dof->coordinates(vector[i].coordinates[0], vector[i].coordinates[1], vector[i].coordinates[2]);
        dof->numDof(vector[i].numDof);
      }
      else if(vector[i].assembleType == dofs::AssembleType::UnknownBubble) {
        dof = new(_pool.getNextUnknownDofMemoryPlace()) dofs::UnknownBubbleDof(vector[i].numType, vector[i].entity);
        dof->coordinates(vector[i].coordinates[0], vector[i].coordinates[1], vector[i].coordinates[2]);
        dof->numDof(vector[i].numDof);
      }
      else if(vector[i].assembleType == dofs::AssembleType::UnknownGlobal) {
        dof = new(_pool.getNextUnknownDofMemoryPlace()) dofs::UnknownGlobalDof(vector[i].numType, vector[i].entity);
        dof->coordinates(vector[i].coordinates[0], vector[i].coordinates[1], vector[i].coordinates[2]);
        dof->numDof(vector[i].numDof);
      }
      else {
        msg::warning << "Trying to set a dof with unkown assembleType: ignoring" << msg::endl;
        continue;
      }
      auto out = _values.insert(hashmap::pair< dofs::Dof *, T_Scalar >(dof, 0.));
      if(!out.second) {
        if(vector[i].assembleType == 1)
          _pool.lastFixedDofIsUnused();
        else
          _pool.lastUnknownDofIsUnused();
      }
    }
  }

  template< class T_Scalar >
  void FieldInterface< T_Scalar >::getUnknownVector(algebra::Vector< T_Scalar > &vector, const dofs::RawOrder order) const
  {
    if(order == dofs::RawOrder::Hash) {
      std::vector< T_Scalar > values;
      values.reserve(numberOfUnknownDofs());
      for(auto it = _values.begin(); it != _values.end(); ++it) {
        if(it->first->type() == dofs::Type::Unknown) {
          values.push_back(it->second);
        }
      }
      vector = std::move(values);
    }
    else {
      auto less_num = [](const dofs::Dof *lhs, const dofs::Dof *rhs) { return lhs->less_num(*rhs); };
      auto less_pair = [](const dofs::Dof *lhs, const dofs::Dof *rhs) { return lhs->less_pair(*rhs); };
      std::map< dofs::Dof *, T_Scalar, decltype(order == dofs::RawOrder::System ? less_num : less_pair) > sortedMap(order == dofs::RawOrder::System ? less_num : less_pair);
      for(auto it = _values.begin(); it != _values.end(); ++it) {
        if(it->first->type() == dofs::Type::Unknown) {
          sortedMap.insert(std::make_pair(it->first, it->second));
        }
      }

      std::vector< T_Scalar > values;
      values.reserve(sortedMap.size());
      for(auto it = sortedMap.begin(); it != sortedMap.end(); ++it) {
        values.push_back(it->second);
      }
      vector = std::move(values);
    }
  }

  template< class T_Scalar >
  void FieldInterface< T_Scalar >::setUnknownVector(const algebra::Vector< T_Scalar > &vector, const dofs::RawOrder order)
  {
    const std::vector< T_Scalar > &valuesVector = vector.getStdVector();
    if(order == dofs::RawOrder::Hash) {
      unsigned long long i = 0;
      for(auto it = _values.begin(); it != _values.end(); ++it) {
        if(it->first->type() == dofs::Type::Unknown) {
          it->second = valuesVector[i];
          ++i;
        }
      }
    }
    else {
      auto less_num = [](const dofs::Dof *lhs, const dofs::Dof *rhs) { return lhs->less_num(*rhs); };
      auto less_pair = [](const dofs::Dof *lhs, const dofs::Dof *rhs) { return lhs->less_pair(*rhs); };
      std::map< dofs::Dof *, T_Scalar *, decltype(order == dofs::RawOrder::System ? less_num : less_pair) > sortedMap(order == dofs::RawOrder::System ? less_num : less_pair);
      for(auto it = _values.begin(); it != _values.end(); ++it) {
        if(it->first->type() == dofs::Type::Unknown) {
          sortedMap.insert(std::make_pair(it->first, &it->second));
        }
      }

      unsigned long long i = 0;
      for(auto it = sortedMap.begin(); it != sortedMap.end(); ++it) {
        *it->second = valuesVector[i];
        ++i;
      }
    }
  }

  template< class T_Scalar >
  void FieldInterface< T_Scalar >::getUnknownDofVector(std::vector< dofs::RawDof > &vector) const
  {
    std::vector< dofs::RawDof > dofs;
    dofs.reserve(_values.size());
    for(auto it = _values.begin(); it != _values.end(); ++it) {
      if(it->first->type() == dofs::Type::Unknown) {
        dofs.push_back(it->first->getRawDof());
      }
    }
    vector = std::move(dofs);
  }

  template< class T_Scalar >
  void FieldInterface< T_Scalar >::setUnknownDofVector(const std::vector< dofs::RawDof > &vector, const bool append)
  {
    if(!append) {
      _pool.clear();
      _values.clear();
    }
    for(auto i = 0ULL; i < vector.size(); ++i) {
      dofs::UnknownDof *dof;
      if(vector[i].assembleType == dofs::AssembleType::Unknown) {
        dof = new(_pool.getNextUnknownDofMemoryPlace()) dofs::UnknownDof(vector[i].numType, vector[i].entity);
        dof->coordinates(vector[i].coordinates[0], vector[i].coordinates[1], vector[i].coordinates[2]);
        dof->numDof(vector[i].numDof);
      }
      else if(vector[i].assembleType == dofs::AssembleType::UnknownBubble) {
        dof = new(_pool.getNextUnknownDofMemoryPlace()) dofs::UnknownBubbleDof(vector[i].numType, vector[i].entity);
        dof->coordinates(vector[i].coordinates[0], vector[i].coordinates[1], vector[i].coordinates[2]);
        dof->numDof(vector[i].numDof);
      }
      else if(vector[i].assembleType == dofs::AssembleType::UnknownGlobal) {
        dof = new(_pool.getNextUnknownDofMemoryPlace()) dofs::UnknownGlobalDof(vector[i].numType, vector[i].entity);
        dof->coordinates(vector[i].coordinates[0], vector[i].coordinates[1], vector[i].coordinates[2]);
        dof->numDof(vector[i].numDof);
      }
      else {
        msg::warning << "A non unknown dof is set to field '" << _name << "' with 'seUnknownDofVector': ignoring" << msg::endl;
        continue;
      }

      auto out = _values.insert(hashmap::pair< dofs::Dof *, T_Scalar >(dof, 0.));
      if(!out.second) {
        _pool.lastUnknownDofIsUnused();
      }
    }
  }

  template< class T_Scalar >
  void FieldInterface< T_Scalar >::getFixedVector(algebra::Vector< T_Scalar > &vector, const dofs::RawOrder order) const
  {
    if(order == dofs::RawOrder::Hash) {
      std::vector< T_Scalar > values;
      values.reserve(numberOfFixedDofs());
      for(auto it = _values.begin(); it != _values.end(); ++it) {
        if(it->first->type() == dofs::Type::Fixed) {
          values.push_back(it->second);
        }
      }
      vector = std::move(values);
    }
    else {
      auto less_num = [](const dofs::Dof *lhs, const dofs::Dof *rhs) { return lhs->less_num(*rhs); };
      auto less_pair = [](const dofs::Dof *lhs, const dofs::Dof *rhs) { return lhs->less_pair(*rhs); };
      std::map< dofs::Dof *, T_Scalar, decltype(order == dofs::RawOrder::System ? less_num : less_pair) > sortedMap(order == dofs::RawOrder::System ? less_num : less_pair);
      for(auto it = _values.begin(); it != _values.end(); ++it) {
        if(it->first->type() == dofs::Type::Fixed) {
          sortedMap.insert(std::make_pair(it->first, it->second));
        }
      }

      std::vector< T_Scalar > values;
      values.reserve(sortedMap.size());
      for(auto it = sortedMap.begin(); it != sortedMap.end(); ++it) {
        values.push_back(it->second);
      }
      vector = std::move(values);
    }
  }

  template< class T_Scalar >
  void FieldInterface< T_Scalar >::setFixedVector(const algebra::Vector< T_Scalar > &vector, const dofs::RawOrder order)
  {
    const std::vector< T_Scalar > &valuesVector = vector.getStdVector();
    if(order == dofs::RawOrder::Hash) {
      unsigned long long i = 0;
      for(auto it = _values.begin(); it != _values.end(); ++it) {
        if(it->first->type() == dofs::Type::Fixed) {
          it->second = valuesVector[i];
          ++i;
        }
      }
    }
    else {
      auto less_num = [](const dofs::Dof *lhs, const dofs::Dof *rhs) { return lhs->less_num(*rhs); };
      auto less_pair = [](const dofs::Dof *lhs, const dofs::Dof *rhs) { return lhs->less_pair(*rhs); };
      std::map< dofs::Dof *, T_Scalar *, decltype(order == dofs::RawOrder::System ? less_num : less_pair) > sortedMap(order == dofs::RawOrder::System ? less_num : less_pair);
      for(auto it = _values.begin(); it != _values.end(); ++it) {
        if(it->first->type() == dofs::Type::Fixed) {
          sortedMap.insert(std::make_pair(it->first, &it->second));
        }
      }

      unsigned long long i = 0;
      for(auto it = sortedMap.begin(); it != sortedMap.end(); ++it) {
        *it->second = valuesVector[i];
        ++i;
      }
    }
  }

  template< class T_Scalar >
  void FieldInterface< T_Scalar >::getFixedDofVector(std::vector< dofs::RawDof > &vector) const
  {
    std::vector< dofs::RawDof > dofs;
    dofs.reserve(_values.size());
    for(auto it = _values.begin(); it != _values.end(); ++it) {
      if(it->first->type() == dofs::Type::Fixed) {
        dofs.push_back(it->first->getRawDof());
      }
    }
    vector = std::move(dofs);
  }

  template< class T_Scalar >
  void FieldInterface< T_Scalar >::setFixedDofVector(const std::vector< dofs::RawDof > &vector, const bool append)
  {
    if(!append) {
      _pool.clear();
      _values.clear();
    }
    for(auto i = 0ULL; i < vector.size(); ++i) {
      dofs::FixedDof *dof;
      if(vector[i].assembleType == dofs::AssembleType::Fixed) {
        dof = new(_pool.getNextFixedDofMemoryPlace()) dofs::FixedDof(vector[i].numType, vector[i].entity);
        dof->coordinates(vector[i].coordinates[0], vector[i].coordinates[1], vector[i].coordinates[2]);
        dof->numDof(vector[i].numDof);
      }
      else {
        msg::warning << "A non fixed dof is set to field '" << _name << "' with 'setFixedDofVector': ignoring" << msg::endl;
        continue;
      }

      auto out = _values.insert(hashmap::pair< dofs::Dof *, T_Scalar >(dof, 0.));
      if(!out.second) {
        _pool.lastFixedDofIsUnused();
      }
    }
  }

  template< class T_Scalar >
  void FieldInterface< T_Scalar >::getUnknownIndices(std::vector< unsigned long long > &indices) const
  {
    indices.reserve(numberOfDofs());
    for(auto it = _values.begin(); it != _values.end(); ++it) {
      if(it->first->type() == dofs::Type::Unknown) {
        indices.push_back(it->first->numDof() - 1);
      }
    }
  }

  template< class T_Scalar >
  void FieldInterface< T_Scalar >::getIndices(const std::vector< int > &typeKeys, const std::vector< unsigned long long > &entityKeys, std::vector< unsigned long long > &numIndices, std::vector< int > &typeIndices, const unsigned long long begin, const unsigned long long end) const
  {
    const unsigned int multi = multiplicity();
    for(auto i = begin; i < end; ++i) {
      for(unsigned int j = 0; j < multi; ++j) {
        dofs::SearchDof dof(typeKeys[i] + GMSHFEM_DOF_FIELD_OFFSET * (_tag + j), entityKeys[i]);
        auto it = _values.find(&dof);
        if(it != _values.end()) {
          numIndices[i * multi + j] = it->first->numDof() - 1;
          typeIndices[i * multi + j] = it->first->assembleType();
        }
        else {
          numIndices[i * multi + j] = 0;
          typeIndices[i * multi + j] = 0;
        }
      }
    }
  }

  template< class T_Scalar >
  void FieldInterface< T_Scalar >::getIndices(const std::vector< int > &typeKeys, const std::vector< unsigned long long > &entityKeys, std::vector< unsigned long long > &numIndices, std::vector< int > &typeIndices, std::vector< dofs::Dof * > &dofs, const unsigned long long begin, const unsigned long long end) const
  {
    const unsigned int multi = multiplicity();
    for(auto i = begin; i < end; ++i) {
      for(unsigned int j = 0; j < multi; ++j) {
        dofs::SearchDof dof(typeKeys[i] + GMSHFEM_DOF_FIELD_OFFSET * (_tag + j), entityKeys[i]);
        auto it = _values.find(&dof);
        if(it != _values.end()) {
          numIndices[i * multi + j] = it->first->numDof() - 1;
          typeIndices[i * multi + j] = it->first->assembleType();
          dofs[i * multi + j] = it->first;
        }
        else {
          numIndices[i * multi + j] = 0;
          typeIndices[i * multi + j] = 0;
          dofs[i * multi + j] = nullptr;
        }
      }
    }
  }

  template< class T_Scalar >
  void FieldInterface< T_Scalar >::invalidateOrderedDofCache()
  {
    _orderedDofs.clear();
    for(auto i = 0ULL; i < _globalQuantities.size(); ++i) {
      if(globalValueIsStillValid(_globalQuantities[i].first)) {
        _globalQuantities[i].second->getAssociatedDualField()->invalidateOrderedDofCache();
      }
    }
  }

  template< class T_Scalar >
  unsigned int FieldInterface< T_Scalar >::fillIndices(problem::IndiceBucket &indices, const std::pair< int, int > &entity, const int elementType) const
  {
    int nbrKeysByElements = _functionSpace->getNumberOfKeysByElement(elementType);
    if(!indices.have(_tag)) {
      std::string modelName;
      gmsh::model::getCurrent(modelName);
      auto it = _orderedDofs.find(std::pair(modelName, std::pair(elementType, entity)));
      if(it != _orderedDofs.end()) {
        std::vector< unsigned long long > numIndices(it->second.size());
        std::vector< int > typeIndices(it->second.size());

#pragma omp parallel for num_threads(omp::getMaxThreads())
        for(auto i = 0ULL; i < it->second.size(); ++i) {
          if(it->second[i] != nullptr) {
            numIndices[i] = it->second[i]->numDof() - 1;
            typeIndices[i] = it->second[i]->assembleType();
          }
          else {
            numIndices[i] = 0;
            typeIndices[i] = 0;
          }
        }
        indices.indices(_tag, numIndices, typeIndices);
      }
      else {
        std::vector< int > typeKeys;
        std::vector< unsigned long long > entityKeys;
        std::vector< scalar::Precision< T_Scalar > > coord;

        _functionSpace->getKeys(false, typeKeys, entityKeys, coord, elementType, entity);

        std::vector< unsigned long long > *numIndices = nullptr;
        std::vector< int > *typeIndices = nullptr;
        std::vector< dofs::Dof * > dofs;

        if(multiplicity() == 1) {
          numIndices = &entityKeys;
          typeIndices = &typeKeys;
        }
        else {
          numIndices = new std::vector< unsigned long long >(entityKeys.size() * multiplicity());
          typeIndices = new std::vector< int >(typeKeys.size() * multiplicity());
        }
        dofs.resize(typeKeys.size() * multiplicity());

#pragma omp parallel num_threads(omp::getMaxThreads())
        {
          const unsigned int numThreads = omp::getNumThreads();
          const unsigned int myThreadID = omp::getThreadNum();
          const unsigned long long begin = (myThreadID * typeKeys.size()) / numThreads;
          const unsigned long long end = ((myThreadID + 1) * typeKeys.size()) / numThreads;

          getIndices(typeKeys, entityKeys, *numIndices, *typeIndices, dofs, begin, end);
        }
        _orderedDofs.insert(std::pair(std::pair(modelName, std::pair(elementType, entity)), std::move(dofs)));
        indices.indices(_tag, *numIndices, *typeIndices);

        if(multiplicity() != 1) {
          delete numIndices;
          delete typeIndices;
        }
      }
    }
    return nbrKeysByElements;
  }

  template< class T_Scalar >
  void FieldInterface< T_Scalar >::fillIndices(problem::IndiceBucket &indices, const std::vector< int > &typeKeys, const std::vector< unsigned long long > &entityKeys) const
  {
    if(!indices.have(_tag)) {
      std::vector< unsigned long long > numIndices(multiplicity() * typeKeys.size());
      std::vector< int > typeIndices(multiplicity() * typeKeys.size());
#pragma omp parallel num_threads(omp::getMaxThreads())
      {
        const unsigned int numThreads = omp::getNumThreads();
        const unsigned int myThreadID = omp::getThreadNum();
        const unsigned long long begin = (myThreadID * typeKeys.size()) / numThreads;
        const unsigned long long end = ((myThreadID + 1) * typeKeys.size()) / numThreads;

        getIndices(typeKeys, entityKeys, numIndices, typeIndices, begin, end);
      }
      indices.indices(_tag, numIndices, typeIndices);
    }
  }

  template< class T_Scalar >
  void FieldInterface< T_Scalar >::assignGlobalQuantity(GlobalQuantity< T_Scalar > &globalQuantity)
  {
    if((globalQuantity.domain() & this->_domain).isEmpty() == false) {
      throw common::Exception("Domains of the field and of the global quantity overlap");
    }

    FieldInterface< T_Scalar > *dual = nullptr;
    if(!_globalQuantities.empty()) {
      for(auto i = 0ULL; i < _globalQuantities.size(); ++i) {
        if(globalValueIsStillValid(_globalQuantities[i].first)) {
          dual = _globalQuantities[i].second->getAssociatedDualField();
          dual->domain(dual->domain() | globalQuantity.domain());
        }
      }
    }
    else {
      switch(form()) {
      case Form::Form0:
        dual = new Field< T_Scalar, Form::Form0 >("Dual of " + name(), globalQuantity.domain(), *static_cast< FunctionSpace< scalar::Precision< T_Scalar >, field::Form::Form0 > * >(getFunctionSpace()), model());
        break;
      case Form::Form1:
        dual = new Field< T_Scalar, Form::Form1 >("Dual of " + name(), globalQuantity.domain(), *static_cast< FunctionSpace< scalar::Precision< T_Scalar >, field::Form::Form1 > * >(getFunctionSpace()), model());
        break;
      case Form::Form2:
        dual = new Field< T_Scalar, Form::Form2 >("Dual of " + name(), globalQuantity.domain(), *static_cast< FunctionSpace< scalar::Precision< T_Scalar >, field::Form::Form2 > * >(getFunctionSpace()), model());
        break;
      case Form::Form3:
        dual = new Field< T_Scalar, Form::Form3 >("Dual of " + name(), globalQuantity.domain(), *static_cast< FunctionSpace< scalar::Precision< T_Scalar >, field::Form::Form3 > * >(getFunctionSpace()), model());
        break;
      default:
        break;
      }
    }
    globalQuantity.setAssociatedPrimalField(this);
    globalQuantity.setAssociatedDualField(dual);
    _globalQuantities.push_back(std::make_pair(globalQuantity.tag(), &globalQuantity));
  }

  INSTANTIATE_CLASS(FieldInterface, 4, TEMPLATE_ARGS(std::complex< double >, std::complex< float >, double, float))

  bool fieldIsStillValid(const unsigned int tag)
  {
    if(tag >= s_listOfField.size()) {
      return false;
    }

    return s_listOfField[tag];
  }


} // namespace gmshfem::field
