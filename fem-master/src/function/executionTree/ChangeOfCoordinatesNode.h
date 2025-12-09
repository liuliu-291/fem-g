// GmshFEM - Copyright (C) 2019-2022, A. Royer, E. Béchet, C. Geuzaine, Université de Liège
//
// See the LICENSE.txt file for license information. Please report all
// issues on https://gitlab.onelab.info/gmsh/fem/issues

#ifndef H_GMSHFEM_CHANGEOFCOORDINATESNODE
#define H_GMSHFEM_CHANGEOFCOORDINATESNODE

#include "ExecutionTree.h"
#include "FieldNode.h"
#include "OmpInterface.h"
#include "OperationsInterface.h"
#include "fieldOperations.h"

namespace gmshfem::function
{


  template< class T_Scalar, Degree T_Degree >
  class ChangeOfCoordinatesNode final : public ExecutionTreeWithDegreeAndLeaves< T_Scalar, T_Degree, T_Degree >
  {
   private:
    const std::vector< scalar::Precision< T_Scalar >, numa::allocator< scalar::Precision< T_Scalar > > > _zeros{0., 0., 0.};
    const function::ExecutionTreeWithDegree< scalar::Precision< T_Scalar >, Degree::Degree0 > *_coordLeaf[3];
    mutable std::vector< scalar::Precision< T_Scalar >, numa::allocator< scalar::Precision< T_Scalar > > > _tmpPoints;
    mutable OutputVector< scalar::Precision< T_Scalar >, Degree::Degree0 > _tmpCoord;

    ChangeOfCoordinatesNode(const ChangeOfCoordinatesNode &other) :
      ExecutionTreeWithDegreeAndLeaves< T_Scalar, T_Degree, T_Degree >(other)
    {
      init(other._coordLeaf[0]->copy(), other._coordLeaf[1]->copy(), other._coordLeaf[2]->copy(), other._leaves.get(0)->copy());
    }

   public:
    ChangeOfCoordinatesNode(const function::ExecutionTreeWithDegree< scalar::Precision< T_Scalar >, Degree::Degree0 > *const leafX,
                            const function::ExecutionTreeWithDegree< scalar::Precision< T_Scalar >, Degree::Degree0 > *const leafY,
                            const function::ExecutionTreeWithDegree< scalar::Precision< T_Scalar >, Degree::Degree0 > *const leafZ,
                            const function::ExecutionTreeWithDegree< T_Scalar, T_Degree > *const leaf0) :
      ExecutionTreeWithDegreeAndLeaves< T_Scalar, T_Degree, T_Degree >()
    {
      init(leafX, leafY, leafZ, leaf0);
    }

    void init(const function::ExecutionTreeWithDegree< scalar::Precision< T_Scalar >, Degree::Degree0 > *const leafX,
              const function::ExecutionTreeWithDegree< scalar::Precision< T_Scalar >, Degree::Degree0 > *const leafY,
              const function::ExecutionTreeWithDegree< scalar::Precision< T_Scalar >, Degree::Degree0 > *const leafZ,
              const function::ExecutionTree< T_Scalar > *const leaf0)
    {
      _coordLeaf[0] = leafX;
      _coordLeaf[1] = leafY;
      _coordLeaf[2] = leafZ;
      // We need to manually activate the point evaluation algorithm if FieldNodes are found
      for(auto it = leaf0->begin(); it != leaf0->end(); ++it) {
        if((*it)->arity() == function::Arity::Nullary) {
          const function::FieldNode< function::None< T_Scalar, field::Form::Form0 > > *castForm0 = dynamic_cast< const function::FieldNode< function::None< T_Scalar, field::Form::Form0 > > * >(*it);
          if(castForm0 != nullptr) {
            castForm0->operation().setPointEvaluation(true);
            continue;
          }
          const function::FieldNode< function::None< T_Scalar, field::Form::Form1 > > *castForm1 = dynamic_cast< const function::FieldNode< function::None< T_Scalar, field::Form::Form1 > > * >(*it);
          if(castForm1 != nullptr) {
            castForm1->operation().setPointEvaluation(true);
            continue;
          }
          const function::FieldNode< function::None< T_Scalar, field::Form::Form2 > > *castForm2 = dynamic_cast< const function::FieldNode< function::None< T_Scalar, field::Form::Form2 > > * >(*it);
          if(castForm2 != nullptr) {
            castForm2->operation().setPointEvaluation(true);
            continue;
          }
          const function::FieldNode< function::None< T_Scalar, field::Form::Form3 > > *castForm3 = dynamic_cast< const function::FieldNode< function::None< T_Scalar, field::Form::Form3 > > * >(*it);
          if(castForm3 != nullptr) {
            castForm3->operation().setPointEvaluation(true);
            continue;
          }

          const function::FieldNode< function::Derivative< T_Scalar, field::Form::Form0 > > *castdForm0 = dynamic_cast< const function::FieldNode< function::Derivative< T_Scalar, field::Form::Form0 > > * >(*it);
          if(castdForm0 != nullptr) {
            castdForm0->operation().setPointEvaluation(true);
            continue;
          }
          const function::FieldNode< function::Derivative< T_Scalar, field::Form::Form1 > > *castdForm1 = dynamic_cast< const function::FieldNode< function::Derivative< T_Scalar, field::Form::Form1 > > * >(*it);
          if(castdForm1 != nullptr) {
            castdForm1->operation().setPointEvaluation(true);
            continue;
          }
          const function::FieldNode< function::Derivative< T_Scalar, field::Form::Form2 > > *castdForm2 = dynamic_cast< const function::FieldNode< function::Derivative< T_Scalar, field::Form::Form2 > > * >(*it);
          if(castdForm2 != nullptr) {
            castdForm2->operation().setPointEvaluation(true);
            continue;
          }
        }
      }
      this->_leaves.set(0, leaf0);
    }

    virtual NodeType nodeType() const override
    {
      return NodeType::ChangeOfCoordinates;
    }

    ~ChangeOfCoordinatesNode()
    {
      delete _coordLeaf[0];
      delete _coordLeaf[1];
      delete _coordLeaf[2];
    }

    virtual void evaluate(OutputVector< T_Scalar, T_Degree > &values, const std::vector< scalar::Precision< T_Scalar >, numa::allocator< scalar::Precision< T_Scalar > > > &points, const std::vector< scalar::Precision< T_Scalar > > &gaussPoints, const int elementType, const std::pair< int, int > &entity) const override
    {
#pragma omp master
      _tmpPoints.resize(points.size());
#pragma omp barrier

      for(auto i = 0; i < 3; ++i) {
        _coordLeaf[i]->evaluate(_tmpCoord, _coordLeaf[i]->isConstant(entity) ? _zeros : points, gaussPoints, elementType, entity);
#pragma omp barrier
        if(_coordLeaf[i]->isConstant(entity)) {
          InputVector< scalar::Precision< T_Scalar >, Degree::Degree0, std::true_type > inputVector(_tmpCoord);
#pragma omp for
          for(auto j = 0ULL; j < points.size() / 3; ++j) {
            _tmpPoints[3 * j + i] = inputVector[j];
          }
        }
        else {
          InputVector< scalar::Precision< T_Scalar >, Degree::Degree0, std::false_type > inputVector(_tmpCoord);
#pragma omp for
          for(auto j = 0ULL; j < points.size() / 3; ++j) {
            _tmpPoints[3 * j + i] = inputVector[j];
          }
        }
#pragma omp single
        _tmpCoord.clear();
      }
#pragma omp single nowait
      {
        _tmpCoord.clear();
        _tmpCoord.shrink_to_fit();
      }
#pragma omp master
      values.resize(points.size() / 3);
#pragma omp barrier

      this->_leaves.get(0)->evaluate(values, _tmpPoints, gaussPoints, elementType, entity);
    }

    virtual ChangeOfCoordinatesNode *copy() const override
    {
      return new ChangeOfCoordinatesNode(*this);
    }

    virtual void exportTree(const domain::Domain &domain, common::HTMLNode *node) const override
    {
      node->name = "change of coordinates (x->a, y->b, z->c) of d";
      node->scalar = scalar::Name< T_Scalar >::name;
      node->localMemory = common::Memory(3 * sizeof(scalar::Precision< T_Scalar >));
      node->degree = T_Degree;
      node->nodeType = this->nodeType();
      node->constant = true;
      for(auto it = domain.cbegin(); it != domain.cend(); ++it) {
        node->constant &= this->isConstant(*it);
      }
      node->leaves.resize(4);
      for(unsigned int i = 0; i < 3; ++i) {
        this->_coordLeaf[i]->exportTree(domain, &(node->leaves[i]));
      }
      this->_leaves.get(0)->exportTree(domain, &(node->leaves[3]));
    }

    virtual common::Memory peakMemoryByGaussPoint(const domain::Domain &domain) const override
    {
      bool isConst = true;
      for(auto it = domain.cbegin(); it != domain.cend(); ++it) {
        for(auto i = 0; i < 3; ++i) {
          if(!this->_coordLeaf[i]->isConstant(*it)) {
            isConst = false;
            break;
          }
        }
      }
      common::Memory memCoord = (isConst ? common::Memory(0) : sizeof(scalar::Precision< T_Scalar >));
      common::Memory memLeaf = this->_leaves.get(0)->peakMemoryByGaussPoint(domain);

      return (memCoord > memLeaf ? memCoord : memLeaf) + common::Memory(3 * sizeof(scalar::Precision< T_Scalar >));
    }

    virtual bool operator==(const ExecutionTreeInterface &other) const override
    {
      if(this->tag() == other.tag()) {
        return true;
      }

      const std::type_info &ext1 = typeid(*this);
      const std::type_info &ext2 = typeid(other);
      if(ext1.hash_code() == ext2.hash_code()) {
        if(*(this->_coordLeaf[0]) == *(static_cast< const ChangeOfCoordinatesNode< T_Scalar, T_Degree > & >(other)._coordLeaf[0]) && *(this->_coordLeaf[1]) == *(static_cast< const ChangeOfCoordinatesNode< T_Scalar, T_Degree > & >(other)._coordLeaf[1]) && *(this->_coordLeaf[2]) == *(static_cast< const ChangeOfCoordinatesNode< T_Scalar, T_Degree > & >(other)._coordLeaf[2]) && *(this->_leaves.get(0)) == *static_cast< const ChangeOfCoordinatesNode< T_Scalar, T_Degree > & >(other)._leaves.get(0)) {
          return true;
        }
      }

      return false;
    }
  };


} // namespace gmshfem::function


#endif // H_GMSHFEM_CHANGEOFCOORDINATESNODE
