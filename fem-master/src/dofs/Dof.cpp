// GmshFEM - Copyright (C) 2019-2022, A. Royer, E. Béchet, C. Geuzaine, Université de Liège
//
// See the LICENSE.txt file for license information. Please report all
// issues on https://gitlab.onelab.info/gmsh/fem/issues

#include "Dof.h"

namespace gmshfem::dofs
{


  //
  //  class Dof
  //

  Dof::Dof(const unsigned long long numType, const unsigned long long entity) :
    _numType(numType), _entity(entity), _numDof(0) /*, _coordinates()*/
  {
  }

  Dof::~Dof()
  {
  }

  unsigned long long Dof::numType() const
  {
    return _numType;
  }

  unsigned long long Dof::entity() const
  {
    return _entity;
  }

  void Dof::coordinates(const float &x, const float &y, const float &z)
  {
    _coordinates[0] = x;
    _coordinates[1] = y;
    _coordinates[2] = z;
  }

  float Dof::x() const
  {
    return _coordinates[0];
  }

  float Dof::y() const
  {
    return _coordinates[1];
  }

  float Dof::z() const
  {
    return _coordinates[2];
  }

  void Dof::numDof(const unsigned long long numDof)
  {
    _numDof = numDof;
  }

  unsigned long long Dof::numDof() const
  {
    return _numDof;
  }

  RawDof Dof::getRawDof() const
  {
    RawDof raw{_numType, _entity, _numDof, assembleType(), {_coordinates[0], _coordinates[1], _coordinates[2]}};
    return raw;
  }

  //
  //  class SearchDof : public Dof
  //

  SearchDof::SearchDof(const unsigned long long numType, const unsigned long long entity) :
    Dof(numType, entity)
  {
  }

  SearchDof::~SearchDof()
  {
  }

  Type SearchDof::type() const
  {
    return Type::Search;
  }

  AssembleType SearchDof::assembleType() const
  {
    return AssembleType::Undefined;
  }

  //
  //  class FixedDof : public Dof
  //

  FixedDof::FixedDof(const unsigned long long numType, const unsigned long long entity) :
    Dof(numType, entity)
  {
  }

  FixedDof::~FixedDof()
  {
  }

  Type FixedDof::type() const
  {
    return Type::Fixed;
  }

  AssembleType FixedDof::assembleType() const
  {
    return AssembleType::Fixed;
  }

  //
  //  class UnknownDof : public Dof
  //

  UnknownDof::UnknownDof(const unsigned long long numType, const unsigned long long entity) :
    Dof(numType, entity)
  {
  }

  UnknownDof::~UnknownDof()
  {
  }

  Type UnknownDof::type() const
  {
    return Type::Unknown;
  }

  AssembleType UnknownDof::assembleType() const
  {
    return AssembleType::Unknown;
  }

  bool UnknownDof::isBubble() const
  {
    return false;
  }

  //
  //  class UnknownBubbleDof : public UnknownDof
  //

  UnknownBubbleDof::UnknownBubbleDof(const unsigned long long numType, const unsigned long long entity) :
    UnknownDof(numType, entity)
  {
  }

  UnknownBubbleDof::~UnknownBubbleDof()
  {
  }

  AssembleType UnknownBubbleDof::assembleType() const
  {
    return AssembleType::UnknownBubble;
  }

  bool UnknownBubbleDof::isBubble() const
  {
    return true;
  }

  //
  //  class UnknownGlobalDof : public UnknownDof
  //

  UnknownGlobalDof::UnknownGlobalDof(const unsigned long long numType, const unsigned long long entity) :
    UnknownDof(numType, entity)
  {
  }

  UnknownGlobalDof::~UnknownGlobalDof()
  {
  }

  AssembleType UnknownGlobalDof::assembleType() const
  {
    return AssembleType::UnknownGlobal;
  }

  bool UnknownGlobalDof::isBubble() const
  {
    return false;
  }

  //
  //  class LinkedDof : public Dof
  //

  LinkedDof::LinkedDof(const unsigned long long numType, const unsigned long long entity, const Dof *master) :
    Dof(numType, entity), _master(master)
  {
  }

  LinkedDof::~LinkedDof()
  {
  }

  Type LinkedDof::type() const
  {
    return Type::Linked;
  }

  AssembleType LinkedDof::assembleType() const
  {
    return AssembleType::Linked;
  }

  bool LinkedDof::isBubble() const
  {
    return false;
  }

  const Dof *LinkedDof::master() const
  {
    return _master;
  }

  //
  //  class LinkedBubbleDof : public LinkedDof
  //

  LinkedBubbleDof::LinkedBubbleDof(const unsigned long long numType, const unsigned long long entity, const Dof *master) :
    LinkedDof(numType, entity, master)
  {
  }

  LinkedBubbleDof::~LinkedBubbleDof()
  {
  }

  AssembleType LinkedBubbleDof::assembleType() const
  {
    return AssembleType::LinkedBubble;
  }

  bool LinkedBubbleDof::isBubble() const
  {
    return true;
  }

  const Dof *LinkedBubbleDof::master() const
  {
    return _master;
  }


} // namespace gmshfem::dofs
