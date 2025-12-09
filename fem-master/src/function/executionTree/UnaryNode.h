// GmshFEM - Copyright (C) 2019-2022, A. Royer, E. Béchet, C. Geuzaine, Université de Liège
//
// See the LICENSE.txt file for license information. Please report all
// issues on https://gitlab.onelab.info/gmsh/fem/issues

#ifndef H_GMSHFEM_UNARYNODE
#define H_GMSHFEM_UNARYNODE

#include "ExecutionTree.h"
#include "OmpInterface.h"
#include "OperationsInterface.h"

namespace gmshfem::function
{


  template< class T_Scalar, Degree T_Degree, Degree T_DegreeA >
  class UnaryOperation : public OperationsInterface< T_Scalar, T_Degree >
  {
   public:
    typedef T_Scalar D_Scalar;
    constexpr static Degree D_Degree = T_Degree;
    constexpr static Degree D_DegreeA = T_DegreeA;

    UnaryOperation() :
      OperationsInterface< T_Scalar, T_Degree >()
    {
    }

    UnaryOperation(const UnaryOperation &other) :
      OperationsInterface< T_Scalar, T_Degree >(other)
    {
    }

    virtual bool operator==(const UnaryOperation &other) const
    {
      return true;
    }
  };

  template< class T_Op >
  class UnaryNode final : public ExecutionTreeWithDegreeAndLeaves< typename T_Op::D_Scalar, T_Op::D_Degree, T_Op::D_DegreeA >
  {
   private:
    const std::vector< scalar::Precision< typename T_Op::D_Scalar >, numa::allocator< scalar::Precision< typename T_Op::D_Scalar > > > _zeros{0., 0., 0.};
    T_Op _op;
    mutable OutputVector< typename T_Op::D_Scalar, T_Op::D_DegreeA > _tmpA;

    UnaryNode(const UnaryNode &other) :
      ExecutionTreeWithDegreeAndLeaves< typename T_Op::D_Scalar, T_Op::D_Degree, T_Op::D_DegreeA >(other), _op(other._op)
    {
      this->_leaves.set(0, other._leaves.get(0)->copy());
    }

   public:
    UnaryNode(const function::ExecutionTreeWithDegree< typename T_Op::D_Scalar, T_Op::D_DegreeA > *const leaf0) :
      ExecutionTreeWithDegreeAndLeaves< typename T_Op::D_Scalar, T_Op::D_Degree, T_Op::D_DegreeA >(), _op()
    {
      this->_leaves.set(0, leaf0);
    }

    UnaryNode(const T_Op &other, const function::ExecutionTreeWithDegree< typename T_Op::D_Scalar, T_Op::D_DegreeA > *const leaf0) :
      ExecutionTreeWithDegreeAndLeaves< typename T_Op::D_Scalar, T_Op::D_Degree, T_Op::D_DegreeA >(), _op(other)
    {
      this->_leaves.set(0, leaf0);
    }

    ~UnaryNode()
    {
    }

    virtual NodeType nodeType() const override
    {
      return NodeType::Unary;
    }

    // EfficientEvaluate

    struct EfficientEvaluate {
      void operator()(const UnaryNode< T_Op > &enclose, OutputVector< typename T_Op::D_Scalar, T_Op::D_Degree > &values, const std::vector< scalar::Precision< typename T_Op::D_Scalar >, numa::allocator< scalar::Precision< typename T_Op::D_Scalar > > > &points, const std::vector< scalar::Precision< typename T_Op::D_Scalar > > &gaussPoints, const int elementType, const std::pair< int, int > &entity) const
      {
        bool isConstant = enclose._leaves.get(0)->isConstant(entity);
        OutputVector< typename T_Op::D_Scalar, T_Op::D_DegreeA > *A = nullptr;

        if(enclose._op.canUseSameVectorsForOutputAndInputs()) {
          if constexpr(T_Op::D_Degree == T_Op::D_DegreeA) {
            A = &values;
          }
        }

        if(A == nullptr) {
          A = &enclose._tmpA;
        }

        enclose._leaves.get(0)->evaluate(*A, isConstant ? enclose._zeros : points, gaussPoints, elementType, entity);
#pragma omp barrier
        if(isConstant) {
          enclose._op(values, InputVector< typename T_Op::D_Scalar, T_Op::D_DegreeA, std::true_type >(*A), points, gaussPoints, elementType, entity);
        }
        else {
          enclose._op(values, InputVector< typename T_Op::D_Scalar, T_Op::D_DegreeA, std::false_type >(*A), points, gaussPoints, elementType, entity);
        }
      }

      common::Memory localMemory(const UnaryNode< T_Op > &enclose, const domain::Domain &domain)
      {
        bool isConstant = true;
        for(auto it = domain.cbegin(); it != domain.cend(); ++it) {
          if(!enclose._leaves.get(0)->isConstant(*it)) {
            isConstant = false;
          }
        }

        if(enclose._op.canUseSameVectorsForOutputAndInputs()) {
          if constexpr(T_Op::D_Degree == T_Op::D_DegreeA) {
            return common::Memory(0);
          }
        }

        return isConstant ? common::Memory(0) : common::Memory(sizeof(typename MathObject< typename T_Op::D_Scalar, T_Op::D_DegreeA >::Object));
      }

      common::Memory peakMemory(const UnaryNode< T_Op > &enclose, const domain::Domain &domain)
      {
        bool isConstant = true;
        for(auto it = domain.cbegin(); it != domain.cend(); ++it) {
          if(!enclose._leaves.get(0)->isConstant(*it)) {
            isConstant = false;
          }
        }

        if(enclose._op.canUseSameVectorsForOutputAndInputs()) {
          if constexpr(T_Op::D_Degree == T_Op::D_DegreeA) {
            return enclose._leaves.get(0)->peakMemoryByGaussPoint(domain);
          }
        }

        return (isConstant ? common::Memory(0) : common::Memory(sizeof(typename MathObject< typename T_Op::D_Scalar, T_Op::D_DegreeA >::Object))) + enclose._leaves.get(0)->peakMemoryByGaussPoint(domain);
      }
    };

    virtual void evaluate(OutputVector< typename T_Op::D_Scalar, T_Op::D_Degree > &values, const std::vector< scalar::Precision< typename T_Op::D_Scalar >, numa::allocator< scalar::Precision< typename T_Op::D_Scalar > > > &points, const std::vector< scalar::Precision< typename T_Op::D_Scalar > > &gaussPoints, const int elementType, const std::pair< int, int > &entity) const override
    {
      if(values.size() < points.size() / 3) {
#pragma omp barrier
        FunctionAllocator< typename MathObject< typename T_Op::D_Scalar, T_Op::D_Degree >::Object >::initalizeNUMAMemory(points.size() / 3);
#pragma omp master
        values.resize(points.size() / 3);
#pragma omp barrier
      }
      EfficientEvaluate evaluator;
      evaluator(*this, values, points, gaussPoints, elementType, entity);
#pragma omp barrier
#pragma omp single
      {
        _tmpA.clear();
        _tmpA.shrink_to_fit();
      }
    }

    virtual UnaryNode< T_Op > *copy() const override
    {
      return new UnaryNode< T_Op >(*this);
    }

    virtual void exportTree(const domain::Domain &domain, common::HTMLNode *node) const override
    {
      node->name = _op.name();
      node->scalar = scalar::Name< typename T_Op::D_Scalar >::name;
      EfficientEvaluate evaluator;
      node->localMemory = evaluator.localMemory(*this, domain);
      node->degree = T_Op::D_Degree;
      node->nodeType = this->nodeType();
      node->constant = true;
      node->swappable = _op.canUseSameVectorsForOutputAndInputs();
      for(auto it = domain.cbegin(); it != domain.cend(); ++it) {
        node->constant &= this->isConstant(*it);
      }
      node->leaves.resize(this->numberOfLeaves());
      for(unsigned int i = 0; i < this->numberOfLeaves(); ++i) {
        this->_leaves.get(i)->exportTree(domain, &(node->leaves[i]));
      }
    }

    virtual common::Memory peakMemoryByGaussPoint(const domain::Domain &domain) const override
    {
      EfficientEvaluate evaluator;
      return evaluator.peakMemory(*this, domain);
    }

    virtual bool operator==(const ExecutionTreeInterface &other) const override
    {
      if(this->tag() == other.tag()) {
        return true;
      }

      if(auto cast = dynamic_cast< const UnaryNode< T_Op > * >(&other)) {
        if(_op == cast->_op) {
          for(int i = 0; i < 1; ++i) {
            if(!(*(this->_leaves.get(i)) == *(cast->_leaves.get(i)))) {
              return false;
            }
          }
        }
        else {
          return false;
        }
      }
      else {
        return false;
      }

      return true;
    }
  };


} // namespace gmshfem::function


#endif // H_GMSHFEM_UNARYNODE
