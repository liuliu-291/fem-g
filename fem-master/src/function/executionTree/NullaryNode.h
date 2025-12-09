// GmshFEM - Copyright (C) 2019-2022, A. Royer, E. Béchet, C. Geuzaine, Université de Liège
//
// See the LICENSE.txt file for license information. Please report all
// issues on https://gitlab.onelab.info/gmsh/fem/issues

#ifndef H_GMSHFEM_NULLARYNODE
#define H_GMSHFEM_NULLARYNODE

#include "ExecutionTree.h"
#include "OmpInterface.h"
#include "OperationsInterface.h"

namespace gmshfem::function
{


  template< class T_Scalar, Degree T_Degree >
  class NullaryOperation : public OperationsInterface< T_Scalar, T_Degree >
  {
   public:
    typedef T_Scalar D_Scalar;
    constexpr static Degree D_Degree = T_Degree;

    NullaryOperation() :
      OperationsInterface< T_Scalar, T_Degree >()
    {
    }

    NullaryOperation(const NullaryOperation &other) :
      OperationsInterface< T_Scalar, T_Degree >(other)
    {
    }

    virtual bool isConstant() const = 0;

    virtual bool operator==(const NullaryOperation &other) const
    {
      return true;
    }
  };

  template< class T_Op >
  class NullaryNode final : public ExecutionTreeWithDegreeAndLeaves< typename T_Op::D_Scalar, T_Op::D_Degree >
  {
   private:
    T_Op _op;

    NullaryNode(const NullaryNode &other) :
      ExecutionTreeWithDegreeAndLeaves< typename T_Op::D_Scalar, T_Op::D_Degree >(other), _op(other._op)
    {
    }

   public:
    NullaryNode() :
      ExecutionTreeWithDegreeAndLeaves< typename T_Op::D_Scalar, T_Op::D_Degree >(), _op()
    {
    }

    NullaryNode(const T_Op &other) :
      ExecutionTreeWithDegreeAndLeaves< typename T_Op::D_Scalar, T_Op::D_Degree >(), _op(other)
    {
    }

    ~NullaryNode() {}

    virtual NodeType nodeType() const override
    {
      return NodeType::Nullary;
    }

    virtual bool isConstant(const std::pair< int, int > &entity) const override
    {
      return _op.isConstant();
    }

    virtual void evaluate(OutputVector< typename T_Op::D_Scalar, T_Op::D_Degree > &values, const std::vector< scalar::Precision< typename T_Op::D_Scalar >, numa::allocator< scalar::Precision< typename T_Op::D_Scalar > > > &points, const std::vector< scalar::Precision< typename T_Op::D_Scalar > > &gaussPoints, const int elementType, const std::pair< int, int > &entity) const override
    {
      if(values.size() < points.size() / 3) {
#pragma omp barrier
        FunctionAllocator< typename MathObject< typename T_Op::D_Scalar, T_Op::D_Degree >::Object >::initalizeNUMAMemory(points.size() / 3);
#pragma omp master
        values.resize(points.size() / 3);
#pragma omp barrier
      }
      _op(values, points, gaussPoints, elementType, entity);
#pragma omp barrier
    }

    virtual NullaryNode< T_Op > *copy() const override
    {
      return new NullaryNode< T_Op >(*this);
    }

    virtual void exportTree(const domain::Domain &domain, common::HTMLNode *node) const override
    {
      node->name = _op.name();
      node->scalar = scalar::Name< typename T_Op::D_Scalar >::name;
      node->degree = T_Op::D_Degree;
      node->nodeType = this->nodeType();
      node->constant = true;
      for(auto it = domain.cbegin(); it != domain.cend(); ++it) {
        node->constant &= this->isConstant(*it);
      }
    }

    virtual common::Memory peakMemoryByGaussPoint(const domain::Domain &domain) const override
    {
      return common::Memory(0);
    }

    virtual bool operator==(const ExecutionTreeInterface &other) const override
    {
      if(this->tag() == other.tag()) {
        return true;
      }

      const std::type_info &ext1 = typeid(*this);
      const std::type_info &ext2 = typeid(other);
      if(ext1.hash_code() == ext2.hash_code()) {
        if(this->_op == static_cast< const T_Op & >(static_cast< const NullaryNode< T_Op > & >(other)._op)) {
          return true;
        }
      }

      return false;
    }
  };


} // namespace gmshfem::function


#endif // H_GMSHFEM_NULLARYNODE
