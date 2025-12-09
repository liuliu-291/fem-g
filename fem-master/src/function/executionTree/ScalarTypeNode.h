// GmshFEM - Copyright (C) 2019-2022, A. Royer, E. Béchet, C. Geuzaine, Université de Liège
//
// See the LICENSE.txt file for license information. Please report all
// issues on https://gitlab.onelab.info/gmsh/fem/issues

#ifndef H_GMSHFEM_SCALARTYPENODE
#define H_GMSHFEM_SCALARTYPENODE

#include "ExecutionTree.h"
#include "OmpInterface.h"
#include "OperationsInterface.h"

namespace gmshfem::function
{


  template< class T_ScalarOut, class T_ScalarIn, Degree T_Degree >
  class ScalarTypeOperation : public OperationsInterface< T_ScalarOut, T_Degree >
  {
   public:
    typedef T_ScalarOut D_ScalarOut;
    typedef T_ScalarIn D_ScalarIn;
    constexpr static Degree D_Degree = T_Degree;

    ScalarTypeOperation() :
      OperationsInterface< D_ScalarOut, T_Degree >()
    {
    }

    ScalarTypeOperation(const ScalarTypeOperation &other) :
      OperationsInterface< D_ScalarOut, T_Degree >(other)
    {
    }

    virtual bool operator==(const ScalarTypeOperation &other) const
    {
      return true;
    }
  };

  template< class T_Op >
  class ScalarTypeNode final : public ExecutionTreeWithDegreeAndScalarModification< typename T_Op::D_ScalarOut, typename T_Op::D_ScalarIn, T_Op::D_Degree, T_Op::D_Degree >
  {
   private:
    T_Op _op;
    mutable OutputVector< typename T_Op::D_ScalarIn, T_Op::D_Degree > _tmpA;

    ScalarTypeNode(const ScalarTypeNode &other) :
      ExecutionTreeWithDegreeAndScalarModification< typename T_Op::D_ScalarOut, typename T_Op::D_ScalarIn, T_Op::D_Degree, T_Op::D_Degree >(other), _op(other._op)
    {
      this->_leaf.set(0, other._leaf.get(0)->copy());
    }

   public:
    ScalarTypeNode(const function::ExecutionTreeWithDegree< typename T_Op::D_ScalarIn, T_Op::D_Degree > *const leaf0) :
      ExecutionTreeWithDegreeAndScalarModification< typename T_Op::D_ScalarOut, typename T_Op::D_ScalarIn, T_Op::D_Degree, T_Op::D_Degree >(), _op()
    {
      this->_leaf.set(0, leaf0);
    }

    ScalarTypeNode(const T_Op &other, const function::ExecutionTreeWithDegree< typename T_Op::D_ScalarIn, T_Op::D_Degree > *const leaf0) :
      ExecutionTreeWithDegreeAndScalarModification< typename T_Op::D_ScalarOut, typename T_Op::D_ScalarIn, T_Op::D_Degree, T_Op::D_Degree >(), _op(other)
    {
      this->_leaf.set(0, leaf0);
    }

    virtual NodeType nodeType() const override
    {
      return NodeType::ScalarType;
    }

    ~ScalarTypeNode()
    {
    }

    virtual void evaluate(OutputVector< typename T_Op::D_ScalarOut, T_Op::D_Degree > &values, const std::vector< scalar::Precision< typename T_Op::D_ScalarIn >, numa::allocator< scalar::Precision< typename T_Op::D_ScalarIn > > > &points, const std::vector< scalar::Precision< typename T_Op::D_ScalarIn > > &gaussPoints, const int elementType, const std::pair< int, int > &entity) const override
    {
      this->_leaf.get(0)->evaluate(_tmpA, points, gaussPoints, elementType, entity);
      if(values.size() < points.size() / 3) {
#pragma omp barrier
        FunctionAllocator< typename MathObject< typename T_Op::D_ScalarOut, T_Op::D_Degree >::Object >::initalizeNUMAMemory(points.size() / 3);
#pragma omp master
        values.resize(points.size() / 3);
#pragma omp barrier
      }

      _op(values, _tmpA, points, gaussPoints, elementType, entity);
#pragma omp barrier
#pragma omp single
      {
        _tmpA.clear();
        _tmpA.shrink_to_fit();
      }
    }

    virtual ScalarTypeNode< T_Op > *copy() const override
    {
      return new ScalarTypeNode< T_Op >(*this);
    }

    virtual void exportTree(const domain::Domain &domain, common::HTMLNode *node) const override
    {
      node->name = _op.name();
      node->scalar = scalar::Name< typename T_Op::D_ScalarOut >::name;
      bool isConst = true;
      for(auto it = domain.cbegin(); it != domain.cend(); ++it) {
        if(!this->_leaf.get(0)->isConstant(*it)) {
          isConst = false;
          break;
        }
      }
      node->localMemory = isConst ? common::Memory(0) : common::Memory(sizeof(typename MathObject< typename T_Op::D_ScalarIn, T_Op::D_Degree >::Object));
      node->degree = T_Op::D_Degree;
      node->nodeType = this->nodeType();
      node->constant = true;
      for(auto it = domain.cbegin(); it != domain.cend(); ++it) {
        node->constant &= this->isConstant(*it);
      }
      node->leaves.resize(1);
      this->_leaf.get(0)->exportTree(domain, &(node->leaves[0]));
    }

    virtual common::Memory peakMemoryByGaussPoint(const domain::Domain &domain) const override
    {
      bool isConst = true;
      for(auto it = domain.cbegin(); it != domain.cend(); ++it) {
        if(!this->_leaf.get(0)->isConstant(*it)) {
          isConst = false;
          break;
        }
      }
      return (isConst ? common::Memory(0) : common::Memory(sizeof(typename MathObject< typename T_Op::D_ScalarIn, T_Op::D_Degree >::Object))) + this->_leaf.get(0)->peakMemoryByGaussPoint(domain);
    }

    virtual bool operator==(const ExecutionTreeInterface &other) const override
    {
      if(this->tag() == other.tag()) {
        return true;
      }

      const std::type_info &ext1 = typeid(*this);
      const std::type_info &ext2 = typeid(other);
      if(ext1.hash_code() == ext2.hash_code()) {
        if(this->_op == static_cast< const T_Op & >(static_cast< const ScalarTypeNode< T_Op > & >(other)._op)) {
          if(*(this->_leaf.get(0)) == *static_cast< const ScalarTypeNode< T_Op > & >(other)._leaf.get(0)) {
            return true;
          }
        }
      }

      return false;
    }
  };


} // namespace gmshfem::function


#endif // H_GMSHFEM_SCALARTYPENODE
