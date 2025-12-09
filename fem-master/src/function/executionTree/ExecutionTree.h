// GmshFEM - Copyright (C) 2019-2022, A. Royer, E. Béchet, C. Geuzaine, Université de Liège
//
// See the LICENSE.txt file for license information. Please report all
// issues on https://gitlab.onelab.info/gmsh/fem/issues

#ifndef H_GMSHFEM_EXECUTIONTREE
#define H_GMSHFEM_EXECUTIONTREE

#include "Domain.h"
#include "ExecutionTreeIterator.h"
#include "FunctionAllocator.h"
#include "HTMLTree.h"
#include "Leaves.h"
#include "MathObject.h"
#include "NodeObject.h"
#include "numa.h"

#include <vector>

namespace gmshfem::function
{


  class ExecutionTreeInterface
  {
   protected:
    const unsigned int _tag;

    virtual unsigned int numberOfLeaves() const = 0;
    virtual const ExecutionTreeInterface *getLeaf(unsigned int n) const = 0;

    friend class ExecutionTreeIterator;

    ExecutionTreeInterface(const ExecutionTreeInterface &other);

   public:
    ExecutionTreeInterface();
    ExecutionTreeInterface(const unsigned int tag);
    virtual ~ExecutionTreeInterface();

    virtual Arity arity() const = 0;
    virtual NodeType nodeType() const = 0;
    virtual bool isConstant(const std::pair< int, int > &entity) const = 0;
    virtual Degree degree() const = 0;

    virtual unsigned int tag() const;

    virtual ExecutionTreeIterator begin() const;
    virtual ExecutionTreeIterator end() const;

    virtual ExecutionTreeInterface *copy() const = 0;

    virtual void exportTree(const domain::Domain &domain, common::HTMLNode *node) const = 0;
    virtual common::Memory peakMemoryByGaussPoint(const domain::Domain &domain) const = 0;
    virtual bool operator==(const ExecutionTreeInterface &other) const = 0;
  };

  template< class T_Scalar >
  class ExecutionTree : public ExecutionTreeInterface
  {
   protected:
    virtual unsigned int numberOfLeaves() const override = 0;
    virtual const ExecutionTreeInterface *getLeaf(unsigned int n) const override = 0;

   public:
    ExecutionTree() :
      ExecutionTreeInterface() {}
    ExecutionTree(const unsigned int tag) :
      ExecutionTreeInterface(tag) {}
    virtual ~ExecutionTree() {}

    virtual Arity arity() const override = 0;
    virtual NodeType nodeType() const override = 0;
    virtual bool isConstant(const std::pair< int, int > &entity) const override = 0;
    virtual Degree degree() const override = 0;

    virtual void evaluate(OutputVector< T_Scalar, Degree::Degree0 > &values, const std::vector< scalar::Precision< T_Scalar >, numa::allocator< scalar::Precision< T_Scalar > > > &points, const std::vector< scalar::Precision< T_Scalar > > &gaussPoints, const int elementType, const std::pair< int, int > &entity) const
    {
      throw common::Exception("Unable to evalute an expression");
    }
    virtual void evaluate(OutputVector< T_Scalar, Degree::Degree1 > &values, const std::vector< scalar::Precision< T_Scalar >, numa::allocator< scalar::Precision< T_Scalar > > > &points, const std::vector< scalar::Precision< T_Scalar > > &gaussPoints, const int elementType, const std::pair< int, int > &entity) const
    {
      throw common::Exception("Unable to evalute an expression");
    }
    virtual void evaluate(OutputVector< T_Scalar, Degree::Degree2 > &values, const std::vector< scalar::Precision< T_Scalar >, numa::allocator< scalar::Precision< T_Scalar > > > &points, const std::vector< scalar::Precision< T_Scalar > > &gaussPoints, const int elementType, const std::pair< int, int > &entity) const
    {
      throw common::Exception("Unable to evalute an expression");
    }
    virtual void evaluate(OutputVector< T_Scalar, Degree::Degree4 > &values, const std::vector< scalar::Precision< T_Scalar >, numa::allocator< scalar::Precision< T_Scalar > > > &points, const std::vector< scalar::Precision< T_Scalar > > &gaussPoints, const int elementType, const std::pair< int, int > &entity) const
    {
      throw common::Exception("Unable to evalute an expression");
    }

    void evaluate(OutputVector< T_Scalar, Degree::Degree0 > &values, const scalar::Precision< T_Scalar > x, const scalar::Precision< T_Scalar > y, const scalar::Precision< T_Scalar > z, const std::pair< int, int > &entity) const
    {
      std::vector< scalar::Precision< T_Scalar >, numa::allocator< scalar::Precision< T_Scalar > > > points{x, y, z};
      std::vector< scalar::Precision< T_Scalar > > gaussPoints;
      evaluate(values, points, gaussPoints, 0, entity);
    }
    void evaluate(OutputVector< T_Scalar, Degree::Degree1 > &values, const scalar::Precision< T_Scalar > x, const scalar::Precision< T_Scalar > y, const scalar::Precision< T_Scalar > z, const std::pair< int, int > &entity) const
    {
      std::vector< scalar::Precision< T_Scalar >, numa::allocator< scalar::Precision< T_Scalar > > > points{x, y, z};
      std::vector< scalar::Precision< T_Scalar > > gaussPoints;
      evaluate(values, points, gaussPoints, 0, entity);
    }
    void evaluate(OutputVector< T_Scalar, Degree::Degree2 > &values, const scalar::Precision< T_Scalar > x, const scalar::Precision< T_Scalar > y, const scalar::Precision< T_Scalar > z, const std::pair< int, int > &entity) const
    {
      std::vector< scalar::Precision< T_Scalar >, numa::allocator< scalar::Precision< T_Scalar > > > points{x, y, z};
      std::vector< scalar::Precision< T_Scalar > > gaussPoints;
      evaluate(values, points, gaussPoints, 0, entity);
    }
    void evaluate(OutputVector< T_Scalar, Degree::Degree4 > &values, const scalar::Precision< T_Scalar > x, const scalar::Precision< T_Scalar > y, const scalar::Precision< T_Scalar > z, const std::pair< int, int > &entity) const
    {
      std::vector< scalar::Precision< T_Scalar >, numa::allocator< scalar::Precision< T_Scalar > > > points{x, y, z};
      std::vector< scalar::Precision< T_Scalar > > gaussPoints;
      evaluate(values, points, gaussPoints, 0, entity);
    }

    void evaluate(typename MathObject< T_Scalar, Degree::Degree0 >::Object &value, const scalar::Precision< T_Scalar > x, const scalar::Precision< T_Scalar > y, const scalar::Precision< T_Scalar > z, const std::pair< int, int > &entity) const
    {
      OutputVector< T_Scalar, Degree::Degree0 > values;
      std::vector< scalar::Precision< T_Scalar >, numa::allocator< scalar::Precision< T_Scalar > > > points{x, y, z};
      const std::vector< scalar::Precision< T_Scalar > > gaussPoints;
      evaluate(values, points, gaussPoints, 0, entity);
      value = values[0];
    }
    void evaluate(typename MathObject< T_Scalar, Degree::Degree1 >::Object &value, const scalar::Precision< T_Scalar > x, const scalar::Precision< T_Scalar > y, const scalar::Precision< T_Scalar > z, const std::pair< int, int > &entity) const
    {
      OutputVector< T_Scalar, Degree::Degree1 > values;
      std::vector< scalar::Precision< T_Scalar >, numa::allocator< scalar::Precision< T_Scalar > > > points{x, y, z};
      const std::vector< scalar::Precision< T_Scalar > > gaussPoints;
      evaluate(values, points, gaussPoints, 0, entity);
      value = values[0];
    }
    void evaluate(typename MathObject< T_Scalar, Degree::Degree2 >::Object &value, const scalar::Precision< T_Scalar > x, const scalar::Precision< T_Scalar > y, const scalar::Precision< T_Scalar > z, const std::pair< int, int > &entity) const
    {
      OutputVector< T_Scalar, Degree::Degree2 > values;
      std::vector< scalar::Precision< T_Scalar >, numa::allocator< scalar::Precision< T_Scalar > > > points{x, y, z};
      const std::vector< scalar::Precision< T_Scalar > > gaussPoints;
      evaluate(values, points, gaussPoints, 0, entity);
      value = values[0];
    }
    void evaluate(typename MathObject< T_Scalar, Degree::Degree4 >::Object &value, const scalar::Precision< T_Scalar > x, const scalar::Precision< T_Scalar > y, const scalar::Precision< T_Scalar > z, const std::pair< int, int > &entity) const
    {
      OutputVector< T_Scalar, Degree::Degree4 > values;
      std::vector< scalar::Precision< T_Scalar >, numa::allocator< scalar::Precision< T_Scalar > > > points{x, y, z};
      const std::vector< scalar::Precision< T_Scalar > > gaussPoints;
      evaluate(values, points, gaussPoints, 0, entity);
      value = values[0];
    }

    virtual ExecutionTree< T_Scalar > *copy() const override = 0;
  };

  template< class T_Scalar, Degree T_Degree >
  class ExecutionTreeWithDegree : public ExecutionTree< T_Scalar >
  {
   protected:
    virtual unsigned int numberOfLeaves() const override = 0;
    virtual const ExecutionTreeInterface *getLeaf(unsigned int n) const override = 0;

   public:
    ExecutionTreeWithDegree() :
      ExecutionTree< T_Scalar >() {}
    ExecutionTreeWithDegree(const unsigned int tag) :
      ExecutionTree< T_Scalar >(tag) {}
    virtual ~ExecutionTreeWithDegree() {}

    virtual Arity arity() const override = 0;
    virtual NodeType nodeType() const override = 0;
    virtual bool isConstant(const std::pair< int, int > &entity) const override = 0;
    virtual Degree degree() const override
    {
      return T_Degree;
    }

    virtual ExecutionTreeWithDegree< T_Scalar, T_Degree > *copy() const override = 0;
  };

  template< class T_ScalarOut, class T_ScalarIn, Degree T_Degree, Degree... TT_LeavesDegrees >
  class ExecutionTreeWithDegreeAndScalarModification : public ExecutionTreeWithDegree< T_ScalarOut, T_Degree >
  {
   protected:
    Leaves< T_ScalarIn, TT_LeavesDegrees... > _leaf;

    virtual unsigned int numberOfLeaves() const override
    {
      return _leaf.size();
    }

    virtual const ExecutionTreeInterface *getLeaf(unsigned int n) const override
    {
      return _leaf.get(n);
    }

   public:
    ExecutionTreeWithDegreeAndScalarModification() :
      ExecutionTreeWithDegree< T_ScalarOut, T_Degree >() {}
    ExecutionTreeWithDegreeAndScalarModification(const unsigned int tag) :
      ExecutionTreeWithDegree< T_ScalarOut, T_Degree >(tag) {}
    virtual ~ExecutionTreeWithDegreeAndScalarModification() {}

    Arity arity() const final
    {
      switch(_leaf.size()) {
      case 0:
        return Arity::Nullary;
        break;
      case 1:
        return Arity::Unary;
        break;
      case 2:
        return Arity::Binary;
        break;
      case 3:
        return Arity::Ternary;
        break;
      case 9:
        return Arity::Novenary;
        break;
      default:
        break;
      }
      return Arity::Multary;
    }
    virtual NodeType nodeType() const override = 0;
    bool isConstant(const std::pair< int, int > &entity) const override
    {
      return _leaf.get(0)->isConstant(entity);
    }

    virtual ExecutionTreeWithDegreeAndScalarModification< T_ScalarOut, T_ScalarIn, T_Degree, TT_LeavesDegrees... > *copy() const override = 0;
  };

  template< class T_Scalar, Degree T_Degree, Degree... TT_LeavesDegrees >
  class ExecutionTreeWithDegreeAndLeaves : public ExecutionTreeWithDegree< T_Scalar, T_Degree >
  {
   protected:
    Leaves< T_Scalar, TT_LeavesDegrees... > _leaves;

    virtual unsigned int numberOfLeaves() const override
    {
      return _leaves.size();
    }

    virtual const ExecutionTree< T_Scalar > *getLeaf(unsigned int n) const override
    {
      return _leaves.get(n);
    }

   public:
    ExecutionTreeWithDegreeAndLeaves() :
      ExecutionTreeWithDegree< T_Scalar, T_Degree >() {}
    ExecutionTreeWithDegreeAndLeaves(const unsigned int tag) :
      ExecutionTreeWithDegree< T_Scalar, T_Degree >(tag) {}
    virtual ~ExecutionTreeWithDegreeAndLeaves() {}

    Arity arity() const final
    {
      switch(_leaves.size()) {
      case 0:
        return Arity::Nullary;
        break;
      case 1:
        return Arity::Unary;
        break;
      case 2:
        return Arity::Binary;
        break;
      case 3:
        return Arity::Ternary;
        break;
      case 9:
        return Arity::Novenary;
        break;
      default:
        break;
      }
      return Arity::Multary;
    }
    virtual NodeType nodeType() const override = 0;
    bool isConstant(const std::pair< int, int > &entity) const override
    {
      for(auto i = 0ULL; i < _leaves.size(); ++i) {
        if(!_leaves.get(i)->isConstant(entity)) {
          return false;
        }
      }
      return true;
    }

    virtual ExecutionTreeWithDegreeAndLeaves< T_Scalar, T_Degree, TT_LeavesDegrees... > *copy() const override = 0;
  };


} // namespace gmshfem::function


#endif // H_GMSHFEM_EXECUTIONTREE
