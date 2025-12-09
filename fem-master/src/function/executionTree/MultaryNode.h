// GmshFEM - Copyright (C) 2019-2022, A. Royer, E. Béchet, C. Geuzaine, Université de Liège
//
// See the LICENSE.txt file for license information. Please report all
// issues on https://gitlab.onelab.info/gmsh/fem/issues

#ifndef H_GMSHFEM_MULTARYNODE
#define H_GMSHFEM_MULTARYNODE

#include "Domain.h"
#include "ExecutionTree.h"
#include "OmpInterface.h"
#include "functional.h"

namespace gmshfem::function
{


  template< class T_Scalar, Degree T_Degree >
  class MultaryNode final : public ExecutionTreeWithDegreeAndLeaves< T_Scalar, T_Degree >
  {
   private:
    mutable std::map< domain::Domain, const function::ExecutionTreeWithDegree< T_Scalar, T_Degree > *, common::less< gmshfem::domain::Domain > > _leaves;
    mutable typename std::map< domain::Domain, const function::ExecutionTreeWithDegree< T_Scalar, T_Degree > *, common::less< gmshfem::domain::Domain > >::iterator _it;

    MultaryNode(const MultaryNode &other) :
      ExecutionTreeWithDegreeAndLeaves< T_Scalar, T_Degree >(other)
    {
      for(auto it = other._leaves.begin(); it != other._leaves.end(); ++it) {
        _leaves.insert(std::make_pair(it->first, it->second->copy()));
      }
    }

   public:
    MultaryNode() {}
    ~MultaryNode()
    {
      for(auto it = _leaves.begin(); it != _leaves.end(); ++it) {
        delete it->second;
      }
    }

    virtual NodeType nodeType() const override
    {
      return NodeType::Multary;
    }

    virtual bool isConstant(const std::pair< int, int > &entity) const override
    {
      if(entity.first == -1 && entity.second == -1) {
        return false; // piecewise functions can be locally constant but they are never globally constant.
      }
      for(auto it = _leaves.begin(); it != _leaves.end(); ++it) {
        if(it->first.have(entity)) {
          return it->second->isConstant(entity);
        }
      }
      return true;
    }

    virtual void evaluate(OutputVector< T_Scalar, T_Degree > &values, const std::vector< scalar::Precision< T_Scalar >, numa::allocator< scalar::Precision< T_Scalar > > > &points, const std::vector< scalar::Precision< T_Scalar > > &gaussPoints, const int elementType, const std::pair< int, int > &entity) const override
    {
      if(values.size() < points.size() / 3) {
#pragma omp barrier
        FunctionAllocator< typename MathObject< T_Scalar, T_Degree >::Object >::initalizeNUMAMemory(points.size() / 3);
#pragma omp master
        values.resize(points.size() / 3);
#pragma omp barrier
      }
#pragma omp single
      {
        for(_it = _leaves.begin(); _it != _leaves.end(); ++_it) {
          if(_it->first.have(entity)) {
            break;
          }
        }
      }
      if(_it != _leaves.end()) {
        _it->second->evaluate(values, points, gaussPoints, elementType, entity);
      }
#pragma omp barrier
    }

    virtual MultaryNode< T_Scalar, T_Degree > *copy() const override
    {
      return new MultaryNode(*this);
    }

    void addFunction(const function::ExecutionTreeWithDegree< T_Scalar, T_Degree > *tree, const domain::Domain &domain) const
    {
      for(auto it = _leaves.begin(); it != _leaves.end(); ++it) {
        if((domain & it->first) != domain::Domain()) {
          msg::warning << "The new domain overlaps an existing domain inside a piecewise function: ignoring" << msg::endl;
          return;
        }
      }
      _leaves.insert(std::make_pair(domain, tree));
    }

    virtual void exportTree(const domain::Domain &domain, common::HTMLNode *node) const override
    {
      for(auto it = _leaves.begin(); it != _leaves.end(); ++it) {
        if((domain & it->first) == domain) {
          node->name = "piecewise function";
          node->scalar = scalar::Name< T_Scalar >::name;
          node->degree = T_Degree;
          node->nodeType = this->nodeType();
          node->constant = true;
          for(auto it = domain.cbegin(); it != domain.cend(); ++it) {
            node->constant &= this->isConstant(*it);
          }
          node->leaves.resize(1);
          it->second->exportTree(domain, &(node->leaves[0]));
        }
      }
    }

    virtual common::Memory peakMemoryByGaussPoint(const domain::Domain &domain) const override
    {
      common::Memory peakMemory;
      for(auto it = _leaves.begin(); it != _leaves.end(); ++it) {
        if((domain & it->first) != domain::Domain()) {
          common::Memory mem = it->second->peakMemoryByGaussPoint(domain);
          if(mem > peakMemory) {
            peakMemory = mem;
          }
        }
      }
      return peakMemory;
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


#endif // H_GMSHFEM_MULTARYNODE
