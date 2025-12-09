// GmshFEM - Copyright (C) 2019-2022, A. Royer, E. Béchet, C. Geuzaine, Université de Liège
//
// See the LICENSE.txt file for license information. Please report all
// issues on https://gitlab.onelab.info/gmsh/fem/issues

#ifndef H_GMSHFEM_NOVENARYNODE
#define H_GMSHFEM_NOVENARYNODE

#include "ExecutionTree.h"
#include "OmpInterface.h"
#include "OperationsInterface.h"
#include "templateUtils.h"

namespace gmshfem::function
{


  template< class T_Scalar, Degree T_Degree, Degree T_DegreeA, Degree T_DegreeB, Degree T_DegreeC, Degree T_DegreeD, Degree T_DegreeE, Degree T_DegreeF, Degree T_DegreeG, Degree T_DegreeH, Degree T_DegreeI >
  class NovenaryOperation : public OperationsInterface< T_Scalar, T_Degree >
  {
   public:
    typedef T_Scalar D_Scalar;
    constexpr static Degree D_Degree = T_Degree;

    constexpr static Degree D_DegreeA = T_DegreeA;
    constexpr static Degree D_DegreeB = T_DegreeB;
    constexpr static Degree D_DegreeC = T_DegreeC;

    constexpr static Degree D_DegreeD = T_DegreeD;
    constexpr static Degree D_DegreeE = T_DegreeE;
    constexpr static Degree D_DegreeF = T_DegreeF;

    constexpr static Degree D_DegreeG = T_DegreeG;
    constexpr static Degree D_DegreeH = T_DegreeH;
    constexpr static Degree D_DegreeI = T_DegreeI;

    NovenaryOperation() :
      OperationsInterface< T_Scalar, T_Degree >()
    {
    }

    NovenaryOperation(const NovenaryOperation &other) :
      OperationsInterface< T_Scalar, T_Degree >(other)
    {
    }

    virtual bool operator==(const NovenaryOperation &other) const
    {
      return true;
    }
  };

  template< class T_Op >
  class NovenaryNode final : public ExecutionTreeWithDegreeAndLeaves< typename T_Op::D_Scalar, T_Op::D_Degree, T_Op::D_DegreeA, T_Op::D_DegreeB, T_Op::D_DegreeC, T_Op::D_DegreeD, T_Op::D_DegreeE, T_Op::D_DegreeF, T_Op::D_DegreeG, T_Op::D_DegreeH, T_Op::D_DegreeI >
  {
   private:
    const std::vector< scalar::Precision< typename T_Op::D_Scalar >, numa::allocator< scalar::Precision< typename T_Op::D_Scalar > > > _zeros{0., 0., 0.};
    T_Op _op;
    std::vector< std::vector< bool > > _isEqual;
    mutable std::vector< bool > _isConstant;
    mutable OutputVector< typename T_Op::D_Scalar, T_Op::D_DegreeA > _tmpA;
    mutable OutputVector< typename T_Op::D_Scalar, T_Op::D_DegreeB > _tmpB;
    mutable OutputVector< typename T_Op::D_Scalar, T_Op::D_DegreeC > _tmpC;
    mutable OutputVector< typename T_Op::D_Scalar, T_Op::D_DegreeD > _tmpD;
    mutable OutputVector< typename T_Op::D_Scalar, T_Op::D_DegreeE > _tmpE;
    mutable OutputVector< typename T_Op::D_Scalar, T_Op::D_DegreeF > _tmpF;
    mutable OutputVector< typename T_Op::D_Scalar, T_Op::D_DegreeG > _tmpG;
    mutable OutputVector< typename T_Op::D_Scalar, T_Op::D_DegreeH > _tmpH;
    mutable OutputVector< typename T_Op::D_Scalar, T_Op::D_DegreeI > _tmpI;

    NovenaryNode(const NovenaryNode &other) :
      ExecutionTreeWithDegreeAndLeaves< typename T_Op::D_Scalar, T_Op::D_Degree, T_Op::D_DegreeA, T_Op::D_DegreeB, T_Op::D_DegreeC, T_Op::D_DegreeD, T_Op::D_DegreeE, T_Op::D_DegreeF, T_Op::D_DegreeG, T_Op::D_DegreeH, T_Op::D_DegreeI >(other), _op(other._op), _isEqual(other._isEqual), _isConstant(9)
    {
      this->_leaves.set(0, other._leaves.get(0)->copy());
      this->_leaves.set(1, other._leaves.get(1)->copy());
      this->_leaves.set(2, other._leaves.get(2)->copy());
      this->_leaves.set(3, other._leaves.get(3)->copy());
      this->_leaves.set(4, other._leaves.get(4)->copy());
      this->_leaves.set(5, other._leaves.get(5)->copy());
      this->_leaves.set(6, other._leaves.get(6)->copy());
      this->_leaves.set(7, other._leaves.get(7)->copy());
      this->_leaves.set(8, other._leaves.get(8)->copy());
    }

   public:
    NovenaryNode(const function::ExecutionTreeWithDegree< typename T_Op::D_Scalar, T_Op::D_DegreeA > *const leaf0,
                 const function::ExecutionTreeWithDegree< typename T_Op::D_Scalar, T_Op::D_DegreeB > *const leaf1,
                 const function::ExecutionTreeWithDegree< typename T_Op::D_Scalar, T_Op::D_DegreeC > *const leaf2,
                 const function::ExecutionTreeWithDegree< typename T_Op::D_Scalar, T_Op::D_DegreeD > *const leaf3,
                 const function::ExecutionTreeWithDegree< typename T_Op::D_Scalar, T_Op::D_DegreeE > *const leaf4,
                 const function::ExecutionTreeWithDegree< typename T_Op::D_Scalar, T_Op::D_DegreeF > *const leaf5,
                 const function::ExecutionTreeWithDegree< typename T_Op::D_Scalar, T_Op::D_DegreeG > *const leaf6,
                 const function::ExecutionTreeWithDegree< typename T_Op::D_Scalar, T_Op::D_DegreeH > *const leaf7,
                 const function::ExecutionTreeWithDegree< typename T_Op::D_Scalar, T_Op::D_DegreeI > *const leaf8) :
      ExecutionTreeWithDegreeAndLeaves< typename T_Op::D_Scalar, T_Op::D_Degree, T_Op::D_DegreeA, T_Op::D_DegreeB, T_Op::D_DegreeC, T_Op::D_DegreeD, T_Op::D_DegreeE, T_Op::D_DegreeF, T_Op::D_DegreeG, T_Op::D_DegreeH, T_Op::D_DegreeI >(),
      _op(), _isEqual(9, std::vector< bool >(9)), _isConstant(9)
    {
      this->_leaves.set(0, leaf0);
      this->_leaves.set(1, leaf1);
      this->_leaves.set(2, leaf2);
      this->_leaves.set(3, leaf3);
      this->_leaves.set(4, leaf4);
      this->_leaves.set(5, leaf5);
      this->_leaves.set(6, leaf6);
      this->_leaves.set(7, leaf7);
      this->_leaves.set(8, leaf8);

      for(int i = 0; i < 9; ++i) {
        for(int j = 0; j < i; ++j) {
          _isEqual[i][j] = _isEqual[j][i] = *(this->_leaves.get(i)) == *(this->_leaves.get(j));
        }
      }
    }

    NovenaryNode(const T_Op &other,
                 const function::ExecutionTreeWithDegree< typename T_Op::D_Scalar, T_Op::D_DegreeA > *const leaf0,
                 const function::ExecutionTreeWithDegree< typename T_Op::D_Scalar, T_Op::D_DegreeB > *const leaf1,
                 const function::ExecutionTreeWithDegree< typename T_Op::D_Scalar, T_Op::D_DegreeC > *const leaf2,
                 const function::ExecutionTreeWithDegree< typename T_Op::D_Scalar, T_Op::D_DegreeD > *const leaf3,
                 const function::ExecutionTreeWithDegree< typename T_Op::D_Scalar, T_Op::D_DegreeE > *const leaf4,
                 const function::ExecutionTreeWithDegree< typename T_Op::D_Scalar, T_Op::D_DegreeF > *const leaf5,
                 const function::ExecutionTreeWithDegree< typename T_Op::D_Scalar, T_Op::D_DegreeG > *const leaf6,
                 const function::ExecutionTreeWithDegree< typename T_Op::D_Scalar, T_Op::D_DegreeH > *const leaf7,
                 const function::ExecutionTreeWithDegree< typename T_Op::D_Scalar, T_Op::D_DegreeI > *const leaf8) :
      ExecutionTreeWithDegreeAndLeaves< typename T_Op::D_Scalar, T_Op::D_Degree, T_Op::D_DegreeA, T_Op::D_DegreeB, T_Op::D_DegreeC, T_Op::D_DegreeD, T_Op::D_DegreeE, T_Op::D_DegreeF, T_Op::D_DegreeG, T_Op::D_DegreeH, T_Op::D_DegreeI >(),
      _op(other), _isEqual(9, std::vector< bool >(9)), _isConstant(9)
    {
      this->_leaves.set(0, leaf0);
      this->_leaves.set(1, leaf1);
      this->_leaves.set(2, leaf2);
      this->_leaves.set(3, leaf3);
      this->_leaves.set(4, leaf4);
      this->_leaves.set(5, leaf5);
      this->_leaves.set(6, leaf6);
      this->_leaves.set(7, leaf7);
      this->_leaves.set(8, leaf8);

      for(int i = 0; i < 9; ++i) {
        for(int j = 0; j < i; ++j) {
          _isEqual[i][j] = _isEqual[j][i] = *(this->_leaves.get(i)) == *(this->_leaves.get(j));
        }
      }
    }

    ~NovenaryNode()
    {
    }

    virtual NodeType nodeType() const override
    {
      return NodeType::Novenary;
    }

    // EfficientEvaluate

    struct EfficientEvaluate {
      void operator()(const NovenaryNode< T_Op > &enclose, OutputVector< typename T_Op::D_Scalar, T_Op::D_Degree > &values, const std::vector< scalar::Precision< typename T_Op::D_Scalar >, numa::allocator< scalar::Precision< typename T_Op::D_Scalar > > > &points, const std::vector< scalar::Precision< typename T_Op::D_Scalar > > &gaussPoints, const int elementType, const std::pair< int, int > &entity) const
      {
        bool isAllConstant = true;
        for(int i = 0; i < 9; ++i) {
#pragma omp critical
          enclose._isConstant[i] = enclose._leaves.get(i)->isConstant(entity);
          if(!enclose._isConstant[i]) {
            isAllConstant = false;
          }
        }

        OutputVector< typename T_Op::D_Scalar, T_Op::D_DegreeA > *A = nullptr;
        OutputVector< typename T_Op::D_Scalar, T_Op::D_DegreeB > *B = nullptr;
        OutputVector< typename T_Op::D_Scalar, T_Op::D_DegreeC > *C = nullptr;

        OutputVector< typename T_Op::D_Scalar, T_Op::D_DegreeD > *D = nullptr;
        OutputVector< typename T_Op::D_Scalar, T_Op::D_DegreeE > *E = nullptr;
        OutputVector< typename T_Op::D_Scalar, T_Op::D_DegreeF > *F = nullptr;

        OutputVector< typename T_Op::D_Scalar, T_Op::D_DegreeG > *G = nullptr;
        OutputVector< typename T_Op::D_Scalar, T_Op::D_DegreeH > *H = nullptr;
        OutputVector< typename T_Op::D_Scalar, T_Op::D_DegreeI > *I = nullptr;

        constexpr std::array< Degree, 9 > degrees{T_Op::D_DegreeA, T_Op::D_DegreeB, T_Op::D_DegreeC, T_Op::D_DegreeD, T_Op::D_DegreeE, T_Op::D_DegreeF, T_Op::D_DegreeG, T_Op::D_DegreeH, T_Op::D_DegreeI};
        auto alias = [&](auto i) constexpr
        {
          if constexpr(i == 0)
            return &A;
          else if constexpr(i == 1)
            return &B;
          else if constexpr(i == 2)
            return &C;
          else if constexpr(i == 3)
            return &D;
          else if constexpr(i == 4)
            return &E;
          else if constexpr(i == 5)
            return &F;
          else if constexpr(i == 6)
            return &G;
          else if constexpr(i == 7)
            return &H;
          else if constexpr(i == 8)
            return &I;
        };
        auto aliasTmp = [&](auto i) constexpr
        {
          if constexpr(i == 0)
            return &enclose._tmpA;
          else if constexpr(i == 1)
            return &enclose._tmpB;
          else if constexpr(i == 2)
            return &enclose._tmpC;
          else if constexpr(i == 3)
            return &enclose._tmpD;
          else if constexpr(i == 4)
            return &enclose._tmpE;
          else if constexpr(i == 5)
            return &enclose._tmpF;
          else if constexpr(i == 6)
            return &enclose._tmpG;
          else if constexpr(i == 7)
            return &enclose._tmpH;
          else if constexpr(i == 8)
            return &enclose._tmpI;
        };

        if(enclose._op.canUseSameVectorsForOutputAndInputs()) {
          bool outputInputSwapFound = false;
          common::constexpr_for< 0, 9, 1 >(
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

        common::constexpr_for< 0, 9, 1 >(
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
        common::turn_into_tuple< 0, 9 >(
          [&](auto t) {
            enclose._op(values,
                        InputVector< typename T_Op::D_Scalar, degrees[0], typename std::tuple_element< 0, decltype(t) >::type >(*A),
                        InputVector< typename T_Op::D_Scalar, degrees[1], typename std::tuple_element< 1, decltype(t) >::type >(*B),
                        InputVector< typename T_Op::D_Scalar, degrees[2], typename std::tuple_element< 2, decltype(t) >::type >(*C),
                        InputVector< typename T_Op::D_Scalar, degrees[3], typename std::tuple_element< 3, decltype(t) >::type >(*D),
                        InputVector< typename T_Op::D_Scalar, degrees[4], typename std::tuple_element< 4, decltype(t) >::type >(*E),
                        InputVector< typename T_Op::D_Scalar, degrees[5], typename std::tuple_element< 5, decltype(t) >::type >(*F),
                        InputVector< typename T_Op::D_Scalar, degrees[6], typename std::tuple_element< 6, decltype(t) >::type >(*G),
                        InputVector< typename T_Op::D_Scalar, degrees[7], typename std::tuple_element< 7, decltype(t) >::type >(*H),
                        InputVector< typename T_Op::D_Scalar, degrees[8], typename std::tuple_element< 8, decltype(t) >::type >(*I),
                        points, gaussPoints, elementType, entity);
          },
          enclose._isConstant);
      }

      common::Memory localMemory(const NovenaryNode< T_Op > &enclose, const domain::Domain &domain)
      {
        bool isAllConstant = true;
        for(auto it = domain.cbegin(); it != domain.cend(); ++it) {
          for(int i = 0; i < 9; ++i) {
            enclose._isConstant[i] = enclose._leaves.get(i)->isConstant(*it);
            if(!enclose._isConstant[i]) {
              isAllConstant = false;
            }
          }
        }

        constexpr std::array< Degree, 9 > degrees{T_Op::D_DegreeA, T_Op::D_DegreeB, T_Op::D_DegreeC, T_Op::D_DegreeD, T_Op::D_DegreeE, T_Op::D_DegreeF, T_Op::D_DegreeG, T_Op::D_DegreeH, T_Op::D_DegreeI};
        std::array< bool, 9 > isNeeded{true, true, true, true, true, true, true, true, true};

        if(enclose._op.canUseSameVectorsForOutputAndInputs()) {
          common::constexpr_for< 0, 9, 1 >(
            [&](auto i) constexpr {
              if constexpr(T_Op::D_Degree == degrees[i]) {
                if((enclose._isConstant[i] && isAllConstant) || !enclose._isConstant[i]) {
                  isNeeded[i] = false;
                }
              }
              return;
            });
        }

        common::constexpr_for< 0, 9, 1 >(
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
        common::constexpr_for< 0, 9, 1 >(
          [&](auto i) constexpr {
            if(isNeeded[i]) {
              memory += sizeof(typename MathObject< typename T_Op::D_Scalar, degrees[i] >::Object);
            }
          });

        return memory;
      }

      common::Memory peakMemory(const NovenaryNode< T_Op > &enclose, const domain::Domain &domain)
      {
        bool isAllConstant = true;
        for(auto it = domain.cbegin(); it != domain.cend(); ++it) {
          for(int i = 0; i < 9; ++i) {
            enclose._isConstant[i] = enclose._leaves.get(i)->isConstant(*it);
            if(!enclose._isConstant[i]) {
              isAllConstant = false;
            }
          }
        }

        constexpr std::array< Degree, 9 > degrees{T_Op::D_DegreeA, T_Op::D_DegreeB, T_Op::D_DegreeC, T_Op::D_DegreeD, T_Op::D_DegreeE, T_Op::D_DegreeF, T_Op::D_DegreeG, T_Op::D_DegreeH, T_Op::D_DegreeI};
        std::array< bool, 9 > isNeeded{true, true, true, true, true, true, true, true, true};

        if(enclose._op.canUseSameVectorsForOutputAndInputs()) {
          common::constexpr_for< 0, 9, 1 >(
            [&](auto i) constexpr {
              if constexpr(T_Op::D_Degree == degrees[i]) {
                if((enclose._isConstant[i] && isAllConstant) || !enclose._isConstant[i]) {
                  isNeeded[i] = false;
                }
              }
              return;
            });
        }

        common::constexpr_for< 0, 9, 1 >(
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

        std::array< common::Memory, 9 > memory;
        common::Memory local;
        common::constexpr_for< 0, 9, 1 >(
          [&](auto i) constexpr {
            if(isNeeded[i]) {
              local = sizeof(typename MathObject< typename T_Op::D_Scalar, degrees[i] >::Object);
            }
            memory[i] = local + (i != 0 ? memory[i - 1] : 0);
          });

        local = 0;
        for(int i = 0; i < 9; ++i) {
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
#pragma omp single nowait
      {
        _tmpC.clear();
        _tmpC.shrink_to_fit();
      }
#pragma omp single nowait
      {
        _tmpD.clear();
        _tmpD.shrink_to_fit();
      }
#pragma omp single nowait
      {
        _tmpE.clear();
        _tmpE.shrink_to_fit();
      }
#pragma omp single nowait
      {
        _tmpF.clear();
        _tmpF.shrink_to_fit();
      }
#pragma omp single nowait
      {
        _tmpG.clear();
        _tmpG.shrink_to_fit();
      }
#pragma omp single nowait
      {
        _tmpH.clear();
        _tmpH.shrink_to_fit();
      }
#pragma omp single
      {
        _tmpI.clear();
        _tmpI.shrink_to_fit();
      }
    }

    virtual NovenaryNode< T_Op > *copy() const override
    {
      return new NovenaryNode< T_Op >(*this);
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

      if(auto cast = dynamic_cast< const NovenaryNode< T_Op > * >(&other)) {
        if(_op == cast->_op) {
          for(int i = 0; i < 9; ++i) {
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


#endif // H_GMSHFEM_NOVENARYNODE
