// GmshFEM - Copyright (C) 2019-2022, A. Royer, E. Béchet, C. Geuzaine, Université de Liège
//
// See the LICENSE.txt file for license information. Please report all
// issues on https://gitlab.onelab.info/gmsh/fem/issues

#ifndef H_GMSHFEM_TERNARYNODE
#define H_GMSHFEM_TERNARYNODE

#include "ExecutionTree.h"
#include "Message.h"
#include "OmpInterface.h"
#include "OperationsInterface.h"
#include "templateUtils.h"

namespace gmshfem::function
{


  template< class T_Scalar, Degree T_Degree, Degree T_DegreeA, Degree T_DegreeB, Degree T_DegreeC >
  class TernaryOperation : public OperationsInterface< T_Scalar, T_Degree >
  {
   public:
    typedef T_Scalar D_Scalar;
    constexpr static Degree D_Degree = T_Degree;
    constexpr static Degree D_DegreeA = T_DegreeA;
    constexpr static Degree D_DegreeB = T_DegreeB;
    constexpr static Degree D_DegreeC = T_DegreeC;

    TernaryOperation() :
      OperationsInterface< T_Scalar, T_Degree >()
    {
    }

    TernaryOperation(const TernaryOperation &other) :
      OperationsInterface< T_Scalar, T_Degree >(other)
    {
    }

    virtual bool operator==(const TernaryOperation &other) const
    {
      return true;
    }
  };

  template< class T_Op >
  class TernaryNode final : public ExecutionTreeWithDegreeAndLeaves< typename T_Op::D_Scalar, T_Op::D_Degree, T_Op::D_DegreeA, T_Op::D_DegreeB, T_Op::D_DegreeC >
  {
   private:
    const std::vector< scalar::Precision< typename T_Op::D_Scalar >, numa::allocator< scalar::Precision< typename T_Op::D_Scalar > > > _zeros{0., 0., 0.};
    T_Op _op;
    std::vector< std::vector< bool > > _isEqual;
    mutable std::vector< bool > _isConstant;
    mutable OutputVector< typename T_Op::D_Scalar, T_Op::D_DegreeA > _tmpA;
    mutable OutputVector< typename T_Op::D_Scalar, T_Op::D_DegreeB > _tmpB;
    mutable OutputVector< typename T_Op::D_Scalar, T_Op::D_DegreeC > _tmpC;

    TernaryNode(const TernaryNode &other) :
      ExecutionTreeWithDegreeAndLeaves< typename T_Op::D_Scalar, T_Op::D_Degree, T_Op::D_DegreeA, T_Op::D_DegreeB, T_Op::D_DegreeC >(other),
      _op(other._op), _isEqual(other._isEqual), _isConstant(3)
    {
      this->_leaves.set(0, other._leaves.get(0)->copy());
      this->_leaves.set(1, other._leaves.get(1)->copy());
      this->_leaves.set(2, other._leaves.get(2)->copy());
    }

   public:
    TernaryNode(const function::ExecutionTreeWithDegree< typename T_Op::D_Scalar, T_Op::D_DegreeA > *const leaf0,
                const function::ExecutionTreeWithDegree< typename T_Op::D_Scalar, T_Op::D_DegreeB > *const leaf1,
                const function::ExecutionTreeWithDegree< typename T_Op::D_Scalar, T_Op::D_DegreeC > *const leaf2) :
      ExecutionTreeWithDegreeAndLeaves< typename T_Op::D_Scalar, T_Op::D_Degree, T_Op::D_DegreeA, T_Op::D_DegreeB, T_Op::D_DegreeC >(),
      _op(), _isEqual(3, std::vector< bool >(3)), _isConstant(3)
    {
      this->_leaves.set(0, leaf0);
      this->_leaves.set(1, leaf1);
      this->_leaves.set(2, leaf2);

      for(int i = 0; i < 3; ++i) {
        for(int j = 0; j < i; ++j) {
          _isEqual[i][j] = _isEqual[j][i] = *(this->_leaves.get(i)) == *(this->_leaves.get(j));
        }
      }
    }

    TernaryNode(const T_Op &other,
                const function::ExecutionTreeWithDegree< typename T_Op::D_Scalar, T_Op::D_DegreeA > *const leaf0,
                const function::ExecutionTreeWithDegree< typename T_Op::D_Scalar, T_Op::D_DegreeB > *const leaf1,
                const function::ExecutionTreeWithDegree< typename T_Op::D_Scalar, T_Op::D_DegreeC > *const leaf2) :
      ExecutionTreeWithDegreeAndLeaves< typename T_Op::D_Scalar, T_Op::D_Degree, T_Op::D_DegreeA, T_Op::D_DegreeB, T_Op::D_DegreeC >(),
      _op(other), _isEqual(3, std::vector< bool >(3)), _isConstant(3)
    {
      this->_leaves.set(0, leaf0);
      this->_leaves.set(1, leaf1);
      this->_leaves.set(2, leaf2);

      for(int i = 0; i < 3; ++i) {
        for(int j = 0; j < i; ++j) {
          _isEqual[i][j] = _isEqual[j][i] = *(this->_leaves.get(i)) == *(this->_leaves.get(j));
        }
      }
    }

    ~TernaryNode()
    {
    }

    virtual NodeType nodeType() const override
    {
      return NodeType::Ternary;
    }

    // EfficientEvaluate

    struct EfficientEvaluate {
      void operator()(const TernaryNode< T_Op > &enclose, OutputVector< typename T_Op::D_Scalar, T_Op::D_Degree > &values, const std::vector< scalar::Precision< typename T_Op::D_Scalar >, numa::allocator< scalar::Precision< typename T_Op::D_Scalar > > > &points, const std::vector< scalar::Precision< typename T_Op::D_Scalar > > &gaussPoints, const int elementType, const std::pair< int, int > &entity) const
      {
        bool isAllConstant = true;
        for(int i = 0; i < 3; ++i) {
#pragma omp critical
          enclose._isConstant[i] = enclose._leaves.get(i)->isConstant(entity);
          if(!enclose._isConstant[i]) {
            isAllConstant = false;
          }
        }

        OutputVector< typename T_Op::D_Scalar, T_Op::D_DegreeA > *A = nullptr;
        OutputVector< typename T_Op::D_Scalar, T_Op::D_DegreeB > *B = nullptr;
        OutputVector< typename T_Op::D_Scalar, T_Op::D_DegreeC > *C = nullptr;

        constexpr std::array< Degree, 3 > degrees{T_Op::D_DegreeA, T_Op::D_DegreeB, T_Op::D_DegreeC};
        auto alias = [&](auto i) constexpr
        {
          if constexpr(i == 0)
            return &A;
          else if constexpr(i == 1)
            return &B;
          else if constexpr(i == 2)
            return &C;
        };
        auto aliasTmp = [&](auto i) constexpr
        {
          if constexpr(i == 0)
            return &enclose._tmpA;
          else if constexpr(i == 1)
            return &enclose._tmpB;
          else if constexpr(i == 2)
            return &enclose._tmpC;
        };

        if(enclose._op.canUseSameVectorsForOutputAndInputs()) {
          bool outputInputSwapFound = false;
          common::constexpr_for< 0, 3, 1 >(
            [&](auto i) constexpr {
              if(outputInputSwapFound) {
                return;
              }
              if constexpr(T_Op::D_Degree == degrees[i]) {
                if((enclose._isConstant[i] && isAllConstant) || !enclose._isConstant[i]) {
                  outputInputSwapFound = true;
                  *alias(i) = &values;
                }
              }
              return;
            });
        }

        common::constexpr_for< 0, 3, 1 >(
          [&](auto i) constexpr {
            bool equalFound = false;
            common::constexpr_for< 0, (int)i, 1 >(
              [&](auto j) constexpr {
                if(enclose._isEqual[i][j] && *alias(i) == nullptr) {
                  if constexpr(degrees[i] == degrees[j]) {
                    *alias(i) = *alias(j);
                    equalFound = true;
                    return;
                  }
                }
                return;
              });

            if(*alias(i) == nullptr) {
              *alias(i) = aliasTmp(i);
            }

            if(!equalFound) {
              enclose._leaves.get(i)->evaluate(**alias(i), enclose._isConstant[i] ? enclose._zeros : points, gaussPoints, elementType, entity);
            }
            return;
          });
#pragma omp barrier
        common::turn_into_tuple< 0, 3 >(
          [&](auto t) {
            enclose._op(values,
                        InputVector< typename T_Op::D_Scalar, degrees[0], typename std::tuple_element< 0, decltype(t) >::type >(*A),
                        InputVector< typename T_Op::D_Scalar, degrees[1], typename std::tuple_element< 1, decltype(t) >::type >(*B),
                        InputVector< typename T_Op::D_Scalar, degrees[2], typename std::tuple_element< 2, decltype(t) >::type >(*C),
                        points, gaussPoints, elementType, entity);
          },
          enclose._isConstant);
      }

      common::Memory localMemory(const TernaryNode< T_Op > &enclose, const domain::Domain &domain)
      {
        bool isAllConstant = true;
        for(auto it = domain.cbegin(); it != domain.cend(); ++it) {
          for(int i = 0; i < 3; ++i) {
            enclose._isConstant[i] = enclose._leaves.get(i)->isConstant(*it);
            if(!enclose._isConstant[i]) {
              isAllConstant = false;
            }
          }
        }

        constexpr std::array< Degree, 3 > degrees{T_Op::D_DegreeA, T_Op::D_DegreeB, T_Op::D_DegreeC};
        std::array< bool, 3 > isNeeded{true, true, true};

        if(enclose._op.canUseSameVectorsForOutputAndInputs()) {
          common::constexpr_for< 0, 3, 1 >(
            [&](auto i) constexpr {
              if constexpr(T_Op::D_Degree == degrees[i]) {
                if((enclose._isConstant[i] && isAllConstant) || !enclose._isConstant[i]) {
                  isNeeded[i] = false;
                }
              }
              return;
            });
        }

        common::constexpr_for< 0, 3, 1 >(
          [&](auto i) constexpr {
            common::constexpr_for< 0, (int)i, 1 >(
              [&](auto j) constexpr {
                if(enclose._isEqual[i][j] && isNeeded[i]) {
                  if constexpr(degrees[i] == degrees[j]) {
                    isNeeded[i] = false;
                    return;
                  }
                }
                return;
              });
            return;
          });

        common::Memory memory;
        common::constexpr_for< 0, 3, 1 >(
          [&](auto i) constexpr {
            if(isNeeded[i]) {
              memory += sizeof(typename MathObject< typename T_Op::D_Scalar, degrees[i] >::Object);
            }
          });

        return memory;
      }

      common::Memory peakMemory(const TernaryNode< T_Op > &enclose, const domain::Domain &domain)
      {
        bool isAllConstant = true;
        for(auto it = domain.cbegin(); it != domain.cend(); ++it) {
          for(int i = 0; i < 3; ++i) {
            enclose._isConstant[i] = enclose._leaves.get(i)->isConstant(*it);
            if(!enclose._isConstant[i]) {
              isAllConstant = false;
            }
          }
        }

        constexpr std::array< Degree, 3 > degrees{T_Op::D_DegreeA, T_Op::D_DegreeB, T_Op::D_DegreeC};
        std::array< bool, 3 > isNeeded{true, true, true};

        if(enclose._op.canUseSameVectorsForOutputAndInputs()) {
          common::constexpr_for< 0, 3, 1 >(
            [&](auto i) constexpr {
              if constexpr(T_Op::D_Degree == degrees[i]) {
                if((enclose._isConstant[i] && isAllConstant) || !enclose._isConstant[i]) {
                  isNeeded[i] = false;
                }
              }
              return;
            });
        }

        common::constexpr_for< 0, 3, 1 >(
          [&](auto i) constexpr {
            common::constexpr_for< 0, (int)i, 1 >(
              [&](auto j) constexpr {
                if(enclose._isEqual[i][j] && isNeeded[i]) {
                  if constexpr(degrees[i] == degrees[j]) {
                    isNeeded[i] = false;
                    return;
                  }
                }
                return;
              });
            return;
          });

        std::array< common::Memory, 3 > memory;
        common::Memory local;
        common::constexpr_for< 0, 3, 1 >(
          [&](auto i) constexpr {
            if(isNeeded[i]) {
              local = sizeof(typename MathObject< typename T_Op::D_Scalar, degrees[i] >::Object);
            }
            memory[i] = local + (i != 0 ? memory[i - 1] : 0);
          });

        local = 0;
        for(int i = 0; i < 3; ++i) {
          if(local < memory[i]) {
            local = memory[i];
          }
        }

        return local;
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
#pragma omp single nowait
      {
        _tmpA.clear();
        _tmpA.shrink_to_fit();
      }
#pragma omp single nowait
      {
        _tmpB.clear();
        _tmpB.shrink_to_fit();
      }
#pragma omp single
      {
        _tmpC.clear();
        _tmpC.shrink_to_fit();
      }
    }

    virtual TernaryNode< T_Op > *copy() const override
    {
      return new TernaryNode< T_Op >(*this);
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

      if(auto cast = dynamic_cast< const TernaryNode< T_Op > * >(&other)) {
        if(_op == cast->_op) {
          for(int i = 0; i < 3; ++i) {
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


#endif // H_GMSHFEM_TERNARYNODE
