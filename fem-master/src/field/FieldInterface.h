// GmshFEM - Copyright (C) 2019-2022, A. Royer, E. Béchet, C. Geuzaine, Université de Liège
//
// See the LICENSE.txt file for license information. Please report all
// issues on https://gitlab.onelab.info/gmsh/fem/issues

#ifndef H_GMSHFEM_FIELDINTERFACE
#define H_GMSHFEM_FIELDINTERFACE

#include "Domain.h"
#include "FieldModule.h"
#include "FieldObject.h"
#include "Function.h"
#include "FunctionSpaceInterface.h"
#include "GlobalQuantity.h"
#include "MathObject.h"
#include "Memory.h"
#include "Vector.h"
#include "functional.h"

#include <list>
#include <set>
#include <string>
#include <vector>

#ifdef HAVE_ROBINHOOD
#include <robin_hood.h>
namespace hashmap = robin_hood;
#else
#include <unordered_map>
namespace hashmap = std;
#endif // HAVE_ROBINHOOD

namespace gmshfem::dofs
{
  class Dof;
  class FixedDof;
  class UnknownDof;
  class LinkedDof;
} // namespace gmshfem::dofs

namespace gmshfem::problem
{
  class IndiceBucket;
}

namespace gmshfem::field
{


  template< class T_Scalar, Degree T_Degree >
  struct Constraint {
    domain::Domain _domain;
    function::Function< T_Scalar, T_Degree > _function;

    Constraint(const domain::Domain &domain, const function::Function< T_Scalar, T_Degree > &function) :
      _domain(domain), _function(function)
    {
    }
  };

  template< class T_Scalar >
  struct PeriodicConstraint {
    domain::PeriodicLink _link;
    T_Scalar _bondageCoefficient;

    PeriodicConstraint(const domain::PeriodicLink &link, const T_Scalar &bondageCoefficient) :
      _link(link), _bondageCoefficient(bondageCoefficient)
    {
    }
  };

  template< class T_Scalar >
  class FieldInterface
  {
   private:
    class DofPool
    {
     private:
      std::vector< dofs::UnknownDof * > _unknownDofPool;
      std::vector< dofs::FixedDof * > _fixedDofPool;
      std::vector< dofs::LinkedDof * > _linkedDofPool;
      unsigned long long _unknownDofPos;
      unsigned long long _fixedDofPos;
      unsigned long long _linkedDofPos;

     public:
      DofPool();
      ~DofPool();

      void clear();
      void clearUnknown();
      void clearFixed();
      void clearLinked();

      dofs::UnknownDof *getNextUnknownDofMemoryPlace();
      dofs::FixedDof *getNextFixedDofMemoryPlace();
      dofs::LinkedDof *getNextLinkedDofMemoryPlace();

      void lastUnknownDofIsUnused();
      void lastFixedDofIsUnused();
      void lastLinkedDofIsUnused();

      common::Memory memory() const;
    };

   protected:
    field::FunctionSpaceInterface< scalar::Precision< T_Scalar > > *_functionSpace;
    std::string _name;
    std::string _model;
    const unsigned int _tag;
    domain::Domain _domain;
    hashmap::unordered_map< dofs::Dof *, T_Scalar, common::hash< dofs::Dof * >, common::equal_to< dofs::Dof * > > _values;
    DofPool _pool;
    std::vector< std::pair< unsigned int, GlobalQuantity< T_Scalar > * > > _globalQuantities;
    std::vector< FieldModule< T_Scalar > * > _modules;
    std::vector< PeriodicConstraint< T_Scalar > > _periodicConstraints;
    mutable std::map< std::pair< std::string, std::pair< int, std::pair< int, int > > >, std::vector< dofs::Dof * > > _orderedDofs;

    virtual void _copy(const FieldInterface< T_Scalar > &other);
    virtual bool _checkDomain() const = 0;

    static void s_incrementNbrOfFields();

   public:
    FieldInterface();
    FieldInterface(const std::string &name, const domain::Domain &domain, const std::string &model);
    FieldInterface(const FieldInterface< T_Scalar > &other);
    virtual ~FieldInterface();

    virtual void clear();
    void clearDofs();

    typename hashmap::unordered_map< dofs::Dof *, T_Scalar, common::hash< dofs::Dof * >, common::equal_to< dofs::Dof * > >::iterator begin();
    typename hashmap::unordered_map< dofs::Dof *, T_Scalar, common::hash< dofs::Dof * >, common::equal_to< dofs::Dof * > >::iterator end();
    typename hashmap::unordered_map< dofs::Dof *, T_Scalar, common::hash< dofs::Dof * >, common::equal_to< dofs::Dof * > >::const_iterator begin() const;
    typename hashmap::unordered_map< dofs::Dof *, T_Scalar, common::hash< dofs::Dof * >, common::equal_to< dofs::Dof * > >::const_iterator end() const;

    typename hashmap::unordered_map< dofs::Dof *, T_Scalar, common::hash< dofs::Dof * >, common::equal_to< dofs::Dof * > >::iterator findDof(dofs::Dof *dof);

    typename std::vector< std::pair< unsigned int, GlobalQuantity< T_Scalar > * > >::const_iterator firstGlobalQuantity() const;
    typename std::vector< std::pair< unsigned int, GlobalQuantity< T_Scalar > * > >::const_iterator lastGlobalQuantity() const;

    // slave = coefficient * master
    void addPeriodicConstraint(const domain::PeriodicLink &link, const T_Scalar &coefficient);
    std::vector< PeriodicConstraint< T_Scalar > > getPeriodicConstraints() const;

    std::string name() const;
    std::string model() const;
    void name(const std::string &name);
    virtual field::Form form() const = 0;
    virtual unsigned int multiplicity() const = 0;
    unsigned int tag() const;
    domain::Domain domain() const;
    void domain(const domain::Domain &domain);
    common::Memory memory() const;

    FieldModule< T_Scalar > *getModule(const std::string &name) const;
    void setModule(FieldModule< T_Scalar > *module);
    virtual field::FunctionSpaceInterface< scalar::Precision< T_Scalar > > *getFunctionSpace() const = 0;

    unsigned long long numberOfDofs() const;
    unsigned long long numberOfUnknownDofs() const;
    unsigned long long numberOfFixedDofs() const;
    unsigned long long numberOfLinkedDofs() const;
    void reserve(const unsigned long long size);

    // Manage memory via DofPool
    bool setDof(dofs::UnknownDof *dof);
    bool setDof(dofs::LinkedDof *dof, const T_Scalar &value);
    bool setValue(dofs::FixedDof *dof, const T_Scalar &value);
    T_Scalar getValue(const dofs::Dof *dof) const;
    const dofs::Dof *searchDof(const int numType, const unsigned long long entity, const unsigned int multiplicityIndex = 0) const;
    dofs::UnknownDof *getNextUnknownDofMemoryPlace();
    dofs::FixedDof *getNextFixedDofMemoryPlace();
    dofs::LinkedDof *getNextLinkedDofMemoryPlace();

    void setValuesToZero();

    // Internal functions
    virtual void assignValues(const std::vector< T_Scalar > &values);
    virtual void assignValueTo(const unsigned long long numDof, const dofs::Type &type, const T_Scalar &value);
    virtual void getValues(const std::vector< int > &typeKeys, const std::vector< unsigned long long > &entityKeys, std::vector< T_Scalar > &values, const unsigned long long begin, const unsigned long long end) const;
    virtual void getUnknownValues(std::vector< T_Scalar > &values);
    virtual void getFixedValues(std::vector< T_Scalar > &values) const;
    virtual void getLinkedValues(std::vector< T_Scalar > &values) const;
    virtual void getLinkedIndices(std::vector< unsigned long long > &indices) const;

    // External functions
    void getAllVector(algebra::Vector< T_Scalar > &vector, const dofs::RawOrder order = dofs::RawOrder::System) const;
    void setAllVector(const algebra::Vector< T_Scalar > &vector, const dofs::RawOrder order = dofs::RawOrder::System);

    void getAllDofVector(std::vector< dofs::RawDof > &vector) const;
    void setAllDofVector(const std::vector< dofs::RawDof > &vector, const bool append = false);

    void getUnknownVector(algebra::Vector< T_Scalar > &vector, const dofs::RawOrder order = dofs::RawOrder::System) const;
    void setUnknownVector(const algebra::Vector< T_Scalar > &vector, const dofs::RawOrder order = dofs::RawOrder::System);

    void getUnknownDofVector(std::vector< dofs::RawDof > &vector) const;
    void setUnknownDofVector(const std::vector< dofs::RawDof > &vector, const bool append = false);

    void getFixedVector(algebra::Vector< T_Scalar > &vector, const dofs::RawOrder order = dofs::RawOrder::System) const;
    void setFixedVector(const algebra::Vector< T_Scalar > &vector, const dofs::RawOrder order = dofs::RawOrder::System);

    void getFixedDofVector(std::vector< dofs::RawDof > &vector) const;
    void setFixedDofVector(const std::vector< dofs::RawDof > &vector, const bool append = false);


    void getUnknownIndices(std::vector< unsigned long long > &indices) const;
    void getIndices(const std::vector< int > &typeKeys, const std::vector< unsigned long long > &entityKeys, std::vector< unsigned long long > &numIndices, std::vector< int > &typeIndices, const unsigned long long begin, const unsigned long long end) const;
    void getIndices(const std::vector< int > &typeKeys, const std::vector< unsigned long long > &entityKeys, std::vector< unsigned long long > &numIndices, std::vector< int > &typeIndices, std::vector< dofs::Dof * > &dofs, const unsigned long long begin, const unsigned long long end) const;

    void invalidateOrderedDofCache();
    unsigned int fillIndices(problem::IndiceBucket &indices, const std::pair< int, int > &entity, const int elementType) const;
    void fillIndices(problem::IndiceBucket &indices, const std::vector< int > &typeKeys, const std::vector< unsigned long long > &entityKeys) const;

    void assignGlobalQuantity(GlobalQuantity< T_Scalar > &globalQuantity);
  };

  bool fieldIsStillValid(const unsigned int tag);


} // namespace gmshfem::field


#include "CompoundField.h"
#include "Field.h"

#endif // H_GMSHFEM_FIELDINTERFACE
