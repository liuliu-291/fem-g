// GmshFEM - Copyright (C) 2019-2022, A. Royer, E. Béchet, C. Geuzaine, Université de Liège
//
// See the LICENSE.txt file for license information. Please report all
// issues on https://gitlab.onelab.info/gmsh/fem/issues

#ifndef H_GMSHFEM_ANALYTICALNODE
#define H_GMSHFEM_ANALYTICALNODE

#include "ExecutionTree.h"
#include "OmpInterface.h"
#include "OperationsInterface.h"

namespace gmshfem::function
{


  template< class T_Scalar, Degree T_Degree >
  class AnalyticalOperation : public OperationsInterface< T_Scalar, T_Degree >
  {
   public:
    typedef T_Scalar D_Scalar;
    constexpr static Degree D_Degree = T_Degree;

    AnalyticalOperation() :
      OperationsInterface< T_Scalar, T_Degree >()
    {
    }

    AnalyticalOperation(const AnalyticalOperation &other) :
      OperationsInterface< T_Scalar, T_Degree >(other)
    {
    }

    virtual void operator()(OutputVector< T_Scalar, T_Degree > &values, const std::vector< scalar::Precision< T_Scalar >, numa::allocator< scalar::Precision< T_Scalar > > > &points) const = 0;
  };

  template< class T_Op >
  class AnalyticalNode final : public ExecutionTreeWithDegreeAndLeaves< typename T_Op::D_Scalar, T_Op::D_Degree >
  {
   private:
    T_Op _op;

    AnalyticalNode(const AnalyticalNode &other) :
      ExecutionTreeWithDegreeAndLeaves< typename T_Op::D_Scalar, T_Op::D_Degree >(other), _op(other._op)
    {
    }

   public:
    AnalyticalNode() :
      ExecutionTreeWithDegreeAndLeaves< typename T_Op::D_Scalar, T_Op::D_Degree >(), _op()
    {
    }

    template< class... T_Params >
    AnalyticalNode(T_Params... params) :
      ExecutionTreeWithDegreeAndLeaves< typename T_Op::D_Scalar, T_Op::D_Degree >(), _op(params...)
    {
    }

    AnalyticalNode(const T_Op &op) :
      ExecutionTreeWithDegreeAndLeaves< typename T_Op::D_Scalar, T_Op::D_Degree >(), _op(op)
    {
    }

    ~AnalyticalNode()
    {
    }

    virtual NodeType nodeType() const override
    {
      return NodeType::Analytical;
    }

    virtual bool isConstant(const std::pair< int, int > &entity) const override
    {
      return false; // AnalyticalNode are always non-constant
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
      _op(values, points);
#pragma omp barrier
    }

    virtual AnalyticalNode< T_Op > *copy() const override
    {
      return new AnalyticalNode< T_Op >(*this);
    }

    virtual void exportTree(const domain::Domain &domain, common::HTMLNode *node) const override
    {
      node->name = "analytic function";
      node->scalar = scalar::Name< typename T_Op::D_Scalar >::name;
      node->degree = T_Op::D_Degree;
      node->nodeType = this->nodeType();
      node->constant = false;
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
      return false;
    }
  };


} // namespace gmshfem::function


#endif // H_GMSHFEM_ANALYTICALNODE
