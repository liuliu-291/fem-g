// GmshFEM - Copyright (C) 2019-2022, A. Royer, E. Béchet, C. Geuzaine, Université de Liège
//
// See the LICENSE.txt file for license information. Please report all
// issues on https://gitlab.onelab.info/gmsh/fem/issues

#ifndef H_GMSHFEM_DOF
#define H_GMSHFEM_DOF

#include "Message.h"
#include "scalar.h"

#include <complex>
#include <functional>

namespace gmshfem::dofs
{


  enum class Type {
    Search = -1,
    NotFound = 0,
    Unknown = 1,
    Fixed = 2,
    Linked = 3
  };

  enum AssembleType {
    Undefined = 0,
    Fixed = 1,
    Unknown = 2,
    UnknownBubble = 3,
    UnknownGlobal = 4,
    Linked = 5,
    LinkedBubble = 6
  };

  enum class RawOrder {
    System,
    Pair,
    Hash
  };

  // The structure used to externalize the dof dictionary stored in fields.
  struct RawDof {
    unsigned long long numType; // 8 bytes
    unsigned long long entity; // 8 bytes
    unsigned long long numDof; // 8 bytes
    unsigned int assembleType; // 4 bytes
    float coordinates[3]; // 3*4 bytes
  };

  // The structure used to externalize the dof definition.
  struct RawDofKey {
    unsigned long long numType; // 8 bytes
    unsigned long long entity; // 8 bytes
  };

  class Dof
  {
   private:
    // Dof type
    const unsigned long long _numType;
    const unsigned long long _entity;
    unsigned long long _numDof; //in [1;+inf]
    float _coordinates[3];

   public:
    Dof(const unsigned long long numType, const unsigned long long entity);
    virtual ~Dof();

    unsigned long long numType() const;
    unsigned long long entity() const;

    void coordinates(const float &x, const float &y, const float &z);
    float x() const;
    float y() const;
    float z() const;

    virtual Type type() const = 0;
    virtual AssembleType assembleType() const = 0;

    void numDof(const unsigned long long numDof);
    unsigned long long numDof() const;

    RawDof getRawDof() const;

    inline size_t hash() const
    {
      size_t seed = std::hash< unsigned long long >{}(_entity) + 0x9e3779b9;
      seed ^= std::hash< unsigned long long >{}(_numType) + 0x9e3779b9 + (seed << 6) + (seed >> 2);
      return seed;
    }

    inline bool less_pair(const Dof &b) const
    {
      if(_entity < b._entity) {
        return true;
      }
      else if(_entity == b._entity) {
        return _numType < b._numType;
      }
      return false;
    }

    inline bool less_num(const Dof &b) const
    {
      if(assembleType() < b.assembleType()) {
        return true;
      }
      else if(assembleType() == b.assembleType()) {
        return _numDof < b._numDof;
      }
      return false;
    }

    inline bool less_hash(const Dof &b) const
    {
      if(hash() < b.hash()) {
        return true;
      }
      return false;
    }

    inline bool operator==(const Dof &b) const
    {
      if(_entity != b._entity) return false;
      if(_numType != b._numType) return false;
      return true;
    }

    inline bool operator!=(const Dof &b) const
    {
      if(_entity == b._entity && _numType == b._numType) return false;
      return true;
    }
  };

  class SearchDof final : public Dof
  {
   public:
    SearchDof(const unsigned long long numType, const unsigned long long entity);
    ~SearchDof() final;

    Type type() const override;
    AssembleType assembleType() const override;
  };

  class FixedDof final : public Dof
  {
   public:
    FixedDof(const unsigned long long numType, const unsigned long long entity);
    ~FixedDof() final;

    Type type() const override;
    AssembleType assembleType() const override;
  };

  class UnknownDof : public Dof
  {
   public:
    UnknownDof(const unsigned long long numType, const unsigned long long entity);
    virtual ~UnknownDof();

    virtual Type type() const override;
    virtual AssembleType assembleType() const override;

    virtual bool isBubble() const;
  };

  class UnknownBubbleDof final : public UnknownDof
  {
   public:
    UnknownBubbleDof(const unsigned long long numType, const unsigned long long entity);
    ~UnknownBubbleDof();

    AssembleType assembleType() const override;
    virtual bool isBubble() const override;
  };

  class UnknownGlobalDof final : public UnknownDof
  {
   public:
    UnknownGlobalDof(const unsigned long long numType, const unsigned long long entity);
    ~UnknownGlobalDof();

    AssembleType assembleType() const override;
    virtual bool isBubble() const override;
  };

  class LinkedDof : public Dof
  {
   protected:
    const Dof *_master;

   public:
    LinkedDof(const unsigned long long numType, const unsigned long long entity, const Dof *master);
    virtual ~LinkedDof();

    virtual Type type() const override;
    virtual AssembleType assembleType() const override;

    virtual bool isBubble() const;
    virtual const Dof *master() const;
  };

  class LinkedBubbleDof final : public LinkedDof
  {
   public:
    LinkedBubbleDof(const unsigned long long numType, const unsigned long long entity, const Dof *master);
    ~LinkedBubbleDof();

    AssembleType assembleType() const override;

    virtual bool isBubble() const override;
    virtual const Dof *master() const override;
  };


} // namespace gmshfem::dofs

#endif // H_GMSHFEM_DOF
