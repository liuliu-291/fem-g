// GmshFEM - Copyright (C) 2019-2022, A. Royer, E. Béchet, C. Geuzaine, Université de Liège
//
// See the LICENSE.txt file for license information. Please report all
// issues on https://gitlab.onelab.info/gmsh/fem/issues

#include "Tensor4Function.h"

#include "ExecutionTree.h"
#include "FieldNode.h"
#include "NullaryNode.h"
#include "binaryOperations.h"
#include "fieldNodeFunctions.h"
#include "fieldOperations.h"
#include "gmshfemDefines.h"
#include "instantiate.h"
#include "nullaryOperations.h"
#include "scalarTypeOperations.h"

namespace gmshfem::function
{


  template< class T_Scalar >
  Function< T_Scalar, Degree::Degree4 >::Function() :
    _tree(nullptr)
  {
  }

  template< class T_Scalar >
  Function< T_Scalar, Degree::Degree4 >::Function(const ExecutionTreeWithDegree< T_Scalar, Degree::Degree4 > *tree) :
    _tree(tree)
  {
  }

  template< class T_Scalar >
  Function< T_Scalar, Degree::Degree4 >::Function(const Function &other) :
    _tree(other._tree ? other._tree->copy() : nullptr)
  {
  }

  template< class T_Scalar >
  Function< T_Scalar, Degree::Degree4 >::Function(Function &&other) :
    _tree(other._tree)
  {
    other._tree = nullptr;
  }

  template< class T_Scalar >
  template< class T_SFINAE, class >
  Function< T_Scalar, Degree::Degree4 >::Function(const Function< scalar::Precision< T_SFINAE >, Degree::Degree4 > &other) :
    _tree(new ScalarTypeNode< RealToComplex< scalar::Precision< T_SFINAE >, Degree::Degree4 > >(other.copy()))
  {
  }

  template< class T_Scalar >
  Function< T_Scalar, Degree::Degree4 >::Function(const typename MathObject< T_Scalar, Degree::Degree4 >::Object &value) :
    _tree(new NullaryNode< function::Constant< T_Scalar, Degree::Degree4 > >(function::Constant< T_Scalar, Degree::Degree4 >(value)))
  {
  }

  template< class T_Scalar >
  Function< T_Scalar, Degree::Degree4 > &Function< T_Scalar, Degree::Degree4 >::operator=(const Function &other)
  {
    const ExecutionTreeWithDegree< T_Scalar, Degree::Degree4 > *tmp = _tree;
    _tree = (other._tree ? other._tree->copy() : nullptr);

    if(tmp) {
      delete tmp;
      tmp = nullptr;
    }

    return *this;
  }

  template< class T_Scalar >
  Function< T_Scalar, Degree::Degree4 > &Function< T_Scalar, Degree::Degree4 >::operator=(Function &&other)
  {
    std::swap(_tree, other._tree);
    if(other._tree) {
      delete other._tree;
      other._tree = nullptr;
    }

    return *this;
  }

  template< class T_Scalar >
  template< class T_SFINAE, class >
  Function< T_Scalar, Degree::Degree4 > &Function< T_Scalar, Degree::Degree4 >::operator=(const Function< scalar::Precision< T_SFINAE >, Degree::Degree4 > &other)
  {
    if(_tree) {
      delete _tree;
    }
    const ExecutionTreeWithDegree< scalar::Precision< T_SFINAE >, Degree::Degree4 > *otherTree = other.copy();
    if(otherTree) {
      _tree = new ScalarTypeNode< RealToComplex< scalar::Precision< T_SFINAE >, Degree::Degree4 > >(otherTree);
    }
    else {
      _tree = nullptr;
    }
    return *this;
  }

  template< class T_Scalar >
  Function< T_Scalar, Degree::Degree4 > &Function< T_Scalar, Degree::Degree4 >::operator=(const typename MathObject< T_Scalar, Degree::Degree4 >::Object &value)
  {
    if(_tree) {
      delete _tree;
    }

    _tree = new function::NullaryNode< function::Constant< T_Scalar, Degree::Degree4 > >(function::Constant< T_Scalar, Degree::Degree4 >(value));
    return *this;
  }

  template< class T_Scalar >
  Function< T_Scalar, Degree::Degree4 >::~Function()
  {
    if(_tree) {
      delete _tree;
    }
  }

  template< class T_Scalar >
  ExecutionTreeIterator Function< T_Scalar, Degree::Degree4 >::begin() const
  {
    return _tree->begin();
  }

  template< class T_Scalar >
  ExecutionTreeIterator Function< T_Scalar, Degree::Degree4 >::end() const
  {
    return _tree->end();
  }

  template< class T_Scalar >
  const ExecutionTreeWithDegree< T_Scalar, Degree::Degree4 > *Function< T_Scalar, Degree::Degree4 >::getTree() const
  {
    return _tree;
  }

  template< class T_Scalar >
  const ExecutionTreeWithDegree< T_Scalar, Degree::Degree4 > *Function< T_Scalar, Degree::Degree4 >::copy() const
  {
    if(_tree) {
      return _tree->copy();
    }
    return nullptr;
  }

  template< class T_Scalar >
  void Function< T_Scalar, Degree::Degree4 >::evaluate(std::vector< typename MathObject< T_Scalar, Degree::Degree4 >::Object, numa::allocator< typename MathObject< T_Scalar, Degree::Degree4 >::Object > > &values, const std::vector< scalar::Precision< T_Scalar >, numa::allocator< scalar::Precision< T_Scalar > > > &points, const std::vector< scalar::Precision< T_Scalar > > &gaussPoints, const int elementType, const std::pair< int, int > &entity) const
  {
    if(points.size() == 0) {
      return;
    }
    if(_tree) {
      if(values.size() < points.size() / 3) {
#pragma omp barrier
#pragma omp single
        values.resize(points.size() / 3);
      }
      FunctionAllocator< typename MathObject< T_Scalar, Degree::Degree4 >::Object > allocator(&values[0]);
      OutputVector< T_Scalar, Degree::Degree4 > valuesOV(points.size() / 3, allocator);
#pragma omp barrier
      _tree->evaluate(valuesOV, points, gaussPoints, elementType, entity);
#pragma omp barrier
#pragma omp single
      MemoryPoolAllocator::instance()->shrinkTo(GMSHFEM_FUNCTION_MEMORY_RESIDUAL);
    }
  }

  template< class T_Scalar >
  void Function< T_Scalar, Degree::Degree4 >::evaluate(typename MathObject< T_Scalar, Degree::Degree4 >::Object &value, const scalar::Precision< T_Scalar > x, const scalar::Precision< T_Scalar > y, const scalar::Precision< T_Scalar > z, const std::pair< int, int > &entity) const
  {
    if(_tree) {
      _tree->evaluate(value, x, y, z, entity);
      MemoryPoolAllocator::instance()->shrinkTo(GMSHFEM_FUNCTION_MEMORY_RESIDUAL);
    }
  }

  template< class T_Scalar >
  common::Timer Function< T_Scalar, Degree::Degree4 >::evaluate(typename std::vector< std::vector< typename MathObject< T_Scalar, Degree::Degree4 >::Object, numa::allocator< typename MathObject< T_Scalar, Degree::Degree4 >::Object > > > &values, const domain::Domain &domain, const std::string &gauss) const
  {
    common::Timer timer;
    if(_tree) {
      values.resize(domain.numberOfEntities());
      unsigned int i = 0;
      for(auto itEntity = domain.cbegin(); itEntity != domain.cend(); ++itEntity) {
        const bool needPoints = !isConstant(*itEntity);
        values[i].clear();
        std::vector< int > elementTypes;
        gmsh::model::mesh::getElementTypes(elementTypes, itEntity->first, itEntity->second);
        for(auto typeIndex = 0ULL; typeIndex < elementTypes.size(); ++typeIndex) {
          std::vector< scalar::Precision< T_Scalar >, numa::allocator< scalar::Precision< T_Scalar > > > jacobians;
          std::vector< double > gmshJacobians;
          std::vector< scalar::Precision< T_Scalar >, numa::allocator< scalar::Precision< T_Scalar > > > determinants;
          std::vector< double > gmshDeterminants;
          std::vector< scalar::Precision< T_Scalar >, numa::allocator< scalar::Precision< T_Scalar > > > points;
          std::vector< double > gmshPoints;

          std::vector< scalar::Precision< T_Scalar > > gaussWeights;
          std::vector< scalar::Precision< T_Scalar > > gaussPoints;
          std::vector< double > gmshGaussPoints;

          unsigned int nbrOfGaussPoints = field::FunctionSpaceInterface< scalar::Precision< T_Scalar > >::GetGaussInfo(gauss, elementTypes[typeIndex], gaussWeights, gaussPoints);
          scalar::copy(gmshGaussPoints, gaussPoints);

          gmsh::model::mesh::preallocateJacobians(elementTypes[typeIndex], nbrOfGaussPoints, false, true, needPoints, gmshJacobians, gmshDeterminants, gmshPoints, itEntity->second);

#pragma omp parallel num_threads(omp::getMaxThreads())
          {
            const unsigned int numThreads = omp::getNumThreads();
            const unsigned int myThreadID = omp::getThreadNum();

            gmsh::model::mesh::getJacobians(elementTypes[typeIndex], gmshGaussPoints, gmshJacobians, gmshDeterminants, gmshPoints, itEntity->second, myThreadID, numThreads);
          }
          numa::copy(jacobians, gmshJacobians);
          numa::copy(determinants, gmshDeterminants);
          numa::copy(points, gmshPoints);

          if(domain.haveJacobiansModificators(*itEntity)) {
            domain.applyJacobiansModificator(points, determinants, jacobians, *itEntity);
          }

          if(values[i].size() < points.size() / 3) {
            values[i].resize(points.size() / 3);
          }
          FunctionAllocator< typename MathObject< T_Scalar, Degree::Degree4 >::Object > allocator(&values[i][0]);
          OutputVector< T_Scalar, Degree::Degree4 > valuesOV(points.size() / 3, allocator);
          timer.tick();
#pragma omp parallel num_threads(omp::getMaxThreads())
          _tree->evaluate(valuesOV, points, gaussPoints, elementTypes[typeIndex], *itEntity);
          timer.tock();
          ++i;
        }
      }
      MemoryPoolAllocator::instance()->shrinkTo(GMSHFEM_FUNCTION_MEMORY_RESIDUAL);
    }
    return timer;
  }

  template< class T_Scalar >
  bool Function< T_Scalar, Degree::Degree4 >::isConstant(const std::pair< int, int > &entity) const
  {
    if(_tree) {
      return _tree->isConstant(entity);
    }
    return true;
  }

  template< class T_Scalar >
  void Function< T_Scalar, Degree::Degree4 >::setPointEvaluation() const
  {
    for(auto it = begin(); it != end(); ++it) {
      function::setPointEvaluation(it);
    }
  }

  template< class T_Scalar >
  void Function< T_Scalar, Degree::Degree4 >::exportExecutionTree(const domain::Domain &domain, const std::string &name, const std::string &path)
  {
    common::HTMLTree html(name, path);
    common::HTMLNode *node = html.getHead();

    _tree->exportTree(domain, node);

    html.write();
  }
  
  template< class T_Scalar >
  const Function< T_Scalar, Degree::Degree4 > &Function< T_Scalar, Degree::Degree4 >::getEvaluableFunction() const
  {
    return *this;
  }

  INSTANTIATE_CLASS_2(Function, 4, 1, TEMPLATE_ARGS(std::complex< double >, std::complex< float >, double, float), TEMPLATE_ARGS(Degree::Degree4))


  template Function< std::complex< double >, Degree::Degree4 >::Function(const Function< scalar::Precision< std::complex< double > >, Degree::Degree4 > &other);
  template Function< std::complex< double >, Degree::Degree4 > &Function< std::complex< double >, Degree::Degree4 >::operator=(const Function< scalar::Precision< std::complex< double > >, Degree::Degree4 > &other);

  template Function< std::complex< float >, Degree::Degree4 >::Function(const Function< scalar::Precision< std::complex< float > >, Degree::Degree4 > &other);
  template Function< std::complex< float >, Degree::Degree4 > &Function< std::complex< float >, Degree::Degree4 >::operator=(const Function< scalar::Precision< std::complex< float > >, Degree::Degree4 > &other);


  // operator+
  template< class T_Scalar >
  Function< T_Scalar, Degree::Degree4 > operator+(const Function< T_Scalar, Degree::Degree4 > &a, const Function< T_Scalar, Degree::Degree4 > &b)
  {
    if(a.isConstant(std::make_pair(-1, -1)) && b.isConstant(std::make_pair(-1, -1))) {
      typename MathObject< T_Scalar, Degree::Degree4 >::Object valueA, valueB;
      a.evaluate(valueA, 0., 0., 0., std::make_pair(-1, -1));
      b.evaluate(valueB, 0., 0., 0., std::make_pair(-1, -1));
      return Function< T_Scalar, Degree::Degree4 >(valueA + valueB);
    }
    return Function< T_Scalar, Degree::Degree4 >(new BinaryNode< Add< T_Scalar, Degree::Degree4 > >(a.copy(), b.copy()));
  }

  INSTANTIATE_OPP(TEMPLATE_RETURN(Function< TEMPLATE_PARAM_1, Degree::Degree4 >), , 0, 4, class, TEMPLATE_ARGS(std::complex< double >, std::complex< float >, double, float), TEMPLATE_PARAMS(const Function< TEMPLATE_PARAM_1, Degree::Degree4 > &, const Function< TEMPLATE_PARAM_1, Degree::Degree4 > &))


  template< class T_Scalar >
  Function< T_Scalar, Degree::Degree4 > operator+(const GeneralEvaluableObject< T_Scalar, Degree::Degree4 > &a, const GeneralEvaluableObject< T_Scalar, Degree::Degree4 > &b)
  {
    return a.getEvaluableFunction() + b.getEvaluableFunction();
  }

  INSTANTIATE_OPP(TEMPLATE_RETURN(Function< TEMPLATE_PARAM_1, Degree::Degree4 >), , 1, 4, class, TEMPLATE_ARGS(std::complex< double >, std::complex< float >, double, float), TEMPLATE_PARAMS(const GeneralEvaluableObject< TEMPLATE_PARAM_1, Degree::Degree4 > &, const GeneralEvaluableObject< TEMPLATE_PARAM_1, Degree::Degree4 > &))


  template< class T_Scalar >
  Function< T_Scalar, Degree::Degree4 > operator+(const typename MathObject< T_Scalar, Degree::Degree4 >::Object &a, const GeneralEvaluableObject< T_Scalar, Degree::Degree4 > &b)
  {
    return Function< T_Scalar, Degree::Degree4 >(a) + b.getEvaluableFunction();
  }

  INSTANTIATE_OPP(TEMPLATE_RETURN(Function< TEMPLATE_PARAM_1, Degree::Degree4 >), , 2, 4, class, TEMPLATE_ARGS(std::complex< double >, std::complex< float >, double, float), TEMPLATE_PARAMS(const typename MathObject< TEMPLATE_PARAM_1, Degree::Degree4 >::Object &, const GeneralEvaluableObject< TEMPLATE_PARAM_1, Degree::Degree4 > &))


  template< class T_Scalar >
  Function< T_Scalar, Degree::Degree4 > operator+(const GeneralEvaluableObject< T_Scalar, Degree::Degree4 > &a, const typename MathObject< T_Scalar, Degree::Degree4 >::Object &b)
  {
    return a.getEvaluableFunction() + Function< T_Scalar, Degree::Degree4 >(b);
  }

  INSTANTIATE_OPP(TEMPLATE_RETURN(Function< TEMPLATE_PARAM_1, Degree::Degree4 >), , 3, 4, class, TEMPLATE_ARGS(std::complex< double >, std::complex< float >, double, float), TEMPLATE_PARAMS(const GeneralEvaluableObject< TEMPLATE_PARAM_1, Degree::Degree4 > &, const typename MathObject< TEMPLATE_PARAM_1, Degree::Degree4 >::Object &))


  // operator-
  template< class T_Scalar >
  Function< T_Scalar, Degree::Degree4 > operator-(const Function< T_Scalar, Degree::Degree4 > &a, const Function< T_Scalar, Degree::Degree4 > &b)
  {
    if(a.isConstant(std::make_pair(-1, -1)) && b.isConstant(std::make_pair(-1, -1))) {
      typename MathObject< T_Scalar, Degree::Degree4 >::Object valueA, valueB;
      a.evaluate(valueA, 0., 0., 0., std::make_pair(-1, -1));
      b.evaluate(valueB, 0., 0., 0., std::make_pair(-1, -1));
      return Function< T_Scalar, Degree::Degree4 >(valueA - valueB);
    }
    return Function< T_Scalar, Degree::Degree4 >(new BinaryNode< Sub< T_Scalar, Degree::Degree4 > >(a.copy(), b.copy()));
  }

  INSTANTIATE_OPM(TEMPLATE_RETURN(Function< TEMPLATE_PARAM_1, Degree::Degree4 >), , 0, 4, class, TEMPLATE_ARGS(std::complex< double >, std::complex< float >, double, float), TEMPLATE_PARAMS(const Function< TEMPLATE_PARAM_1, Degree::Degree4 > &, const Function< TEMPLATE_PARAM_1, Degree::Degree4 > &))


  template< class T_Scalar >
  Function< T_Scalar, Degree::Degree4 > operator-(const GeneralEvaluableObject< T_Scalar, Degree::Degree4 > &a, const GeneralEvaluableObject< T_Scalar, Degree::Degree4 > &b)
  {
    return a.getEvaluableFunction() - b.getEvaluableFunction();
  }

  INSTANTIATE_OPM(TEMPLATE_RETURN(Function< TEMPLATE_PARAM_1, Degree::Degree4 >), , 1, 4, class, TEMPLATE_ARGS(std::complex< double >, std::complex< float >, double, float), TEMPLATE_PARAMS(const GeneralEvaluableObject< TEMPLATE_PARAM_1, Degree::Degree4 > &, const GeneralEvaluableObject< TEMPLATE_PARAM_1, Degree::Degree4 > &))


  template< class T_Scalar >
  Function< T_Scalar, Degree::Degree4 > operator-(const typename MathObject< T_Scalar, Degree::Degree4 >::Object &a, const GeneralEvaluableObject< T_Scalar, Degree::Degree4 > &b)
  {
    return Function< T_Scalar, Degree::Degree4 >(a) - b.getEvaluableFunction();
  }

  INSTANTIATE_OPM(TEMPLATE_RETURN(Function< TEMPLATE_PARAM_1, Degree::Degree4 >), , 2, 4, class, TEMPLATE_ARGS(std::complex< double >, std::complex< float >, double, float), TEMPLATE_PARAMS(const typename MathObject< TEMPLATE_PARAM_1, Degree::Degree4 >::Object &, const GeneralEvaluableObject< TEMPLATE_PARAM_1, Degree::Degree4 > &))


  template< class T_Scalar >
  Function< T_Scalar, Degree::Degree4 > operator-(const GeneralEvaluableObject< T_Scalar, Degree::Degree4 > &a, const typename MathObject< T_Scalar, Degree::Degree4 >::Object &b)
  {
    return a.getEvaluableFunction() - Function< T_Scalar, Degree::Degree4 >(b);
  }

  INSTANTIATE_OPM(TEMPLATE_RETURN(Function< TEMPLATE_PARAM_1, Degree::Degree4 >), , 3, 4, class, TEMPLATE_ARGS(std::complex< double >, std::complex< float >, double, float), TEMPLATE_PARAMS(const GeneralEvaluableObject< TEMPLATE_PARAM_1, Degree::Degree4 > &, const typename MathObject< TEMPLATE_PARAM_1, Degree::Degree4 >::Object &))


  // operator*
  template< class T_Scalar >
  Function< T_Scalar, Degree::Degree4 > operator*(const Function< T_Scalar, Degree::Degree0 > &a, const Function< T_Scalar, Degree::Degree4 > &b)
  {
    if(a.isConstant(std::make_pair(-1, -1)) && b.isConstant(std::make_pair(-1, -1))) {
      typename MathObject< T_Scalar, Degree::Degree0 >::Object valueA;
      typename MathObject< T_Scalar, Degree::Degree4 >::Object valueB;
      a.evaluate(valueA, 0., 0., 0., std::make_pair(-1, -1));
      b.evaluate(valueB, 0., 0., 0., std::make_pair(-1, -1));
      typename MathObject< T_Scalar, Degree::Degree4 >::Object tmp;
      for(unsigned int I = 0; I < 3; ++I) {
        for(unsigned int J = 0; J < 3; ++J) {
          tmp(I, J) = valueA * valueB(I, J);
        }
      }
      return Function< T_Scalar, Degree::Degree4 >(tmp);
    }
    return Function< T_Scalar, Degree::Degree4 >(new BinaryNode< Mul< T_Scalar, Degree::Degree0, Degree::Degree4 > >(a.copy(), b.copy()));
  }

  INSTANTIATE_OPT(TEMPLATE_RETURN(Function< TEMPLATE_PARAM_1, Degree::Degree4 >), , 0, 4, class, TEMPLATE_ARGS(std::complex< double >, std::complex< float >, double, float), TEMPLATE_PARAMS(const Function< TEMPLATE_PARAM_1, Degree::Degree0 > &, const Function< TEMPLATE_PARAM_1, Degree::Degree4 > &))


  template< class T_Scalar >
  Function< T_Scalar, Degree::Degree4 > operator*(const GeneralEvaluableObject< T_Scalar, Degree::Degree0 > &a, const GeneralEvaluableObject< T_Scalar, Degree::Degree4 > &b)
  {
    return a.getEvaluableFunction() * b.getEvaluableFunction();
  }

  INSTANTIATE_OPT(TEMPLATE_RETURN(Function< TEMPLATE_PARAM_1, Degree::Degree4 >), , 1, 4, class, TEMPLATE_ARGS(std::complex< double >, std::complex< float >, double, float), TEMPLATE_PARAMS(const GeneralEvaluableObject< TEMPLATE_PARAM_1, Degree::Degree0 > &, const GeneralEvaluableObject< TEMPLATE_PARAM_1, Degree::Degree4 > &))


  template< class T_Scalar >
  Function< T_Scalar, Degree::Degree4 > operator*(const typename MathObject< T_Scalar, Degree::Degree0 >::Object &a, const GeneralEvaluableObject< T_Scalar, Degree::Degree4 > &b)
  {
    return Function< T_Scalar, Degree::Degree0 >(a) * b.getEvaluableFunction();
  }

  INSTANTIATE_OPT(TEMPLATE_RETURN(Function< TEMPLATE_PARAM_1, Degree::Degree4 >), , 2, 4, class, TEMPLATE_ARGS(std::complex< double >, std::complex< float >, double, float), TEMPLATE_PARAMS(const typename MathObject< TEMPLATE_PARAM_1, Degree::Degree0 >::Object &, const GeneralEvaluableObject< TEMPLATE_PARAM_1, Degree::Degree4 > &))


  template< class T_Scalar >
  Function< T_Scalar, Degree::Degree4 > operator*(const GeneralEvaluableObject< T_Scalar, Degree::Degree0 > &a, const typename MathObject< T_Scalar, Degree::Degree4 >::Object &b)
  {
    return a.getEvaluableFunction() * Function< T_Scalar, Degree::Degree4 >(b);
  }

  INSTANTIATE_OPT(TEMPLATE_RETURN(Function< TEMPLATE_PARAM_1, Degree::Degree4 >), , 3, 4, class, TEMPLATE_ARGS(std::complex< double >, std::complex< float >, double, float), TEMPLATE_PARAMS(const GeneralEvaluableObject< TEMPLATE_PARAM_1, Degree::Degree0 > &, const typename MathObject< TEMPLATE_PARAM_1, Degree::Degree4 >::Object &))


  template< class T_Scalar >
  Function< T_Scalar, Degree::Degree4 > operator*(const Function< T_Scalar, Degree::Degree4 > &a, const Function< T_Scalar, Degree::Degree0 > &b)
  {
    if(a.isConstant(std::make_pair(-1, -1)) && b.isConstant(std::make_pair(-1, -1))) {
      typename MathObject< T_Scalar, Degree::Degree4 >::Object valueA;
      typename MathObject< T_Scalar, Degree::Degree0 >::Object valueB;
      a.evaluate(valueA, 0., 0., 0., std::make_pair(-1, -1));
      b.evaluate(valueB, 0., 0., 0., std::make_pair(-1, -1));
      typename MathObject< T_Scalar, Degree::Degree4 >::Object tmp;
      for(unsigned int I = 0; I < 3; ++I) {
        for(unsigned int J = 0; J < 3; ++J) {
          tmp(I, J) = valueA(I, J) * valueB;
        }
      }
      return Function< T_Scalar, Degree::Degree4 >(tmp);
    }
    return Function< T_Scalar, Degree::Degree4 >(new BinaryNode< Mul< T_Scalar, Degree::Degree4, Degree::Degree0 > >(a.copy(), b.copy()));
  }

  INSTANTIATE_OPT(TEMPLATE_RETURN(Function< TEMPLATE_PARAM_1, Degree::Degree4 >), , 4, 4, class, TEMPLATE_ARGS(std::complex< double >, std::complex< float >, double, float), TEMPLATE_PARAMS(const Function< TEMPLATE_PARAM_1, Degree::Degree4 > &, const Function< TEMPLATE_PARAM_1, Degree::Degree0 > &))


  template< class T_Scalar >
  Function< T_Scalar, Degree::Degree4 > operator*(const GeneralEvaluableObject< T_Scalar, Degree::Degree4 > &a, const GeneralEvaluableObject< T_Scalar, Degree::Degree0 > &b)
  {
    return a.getEvaluableFunction() * b.getEvaluableFunction();
  }

  INSTANTIATE_OPT(TEMPLATE_RETURN(Function< TEMPLATE_PARAM_1, Degree::Degree4 >), , 5, 4, class, TEMPLATE_ARGS(std::complex< double >, std::complex< float >, double, float), TEMPLATE_PARAMS(const GeneralEvaluableObject< TEMPLATE_PARAM_1, Degree::Degree4 > &, const GeneralEvaluableObject< TEMPLATE_PARAM_1, Degree::Degree0 > &))


  template< class T_Scalar >
  Function< T_Scalar, Degree::Degree4 > operator*(const typename MathObject< T_Scalar, Degree::Degree4 >::Object &a, const GeneralEvaluableObject< T_Scalar, Degree::Degree0 > &b)
  {
    return Function< T_Scalar, Degree::Degree4 >(a) * b.getEvaluableFunction();
  }

  INSTANTIATE_OPT(TEMPLATE_RETURN(Function< TEMPLATE_PARAM_1, Degree::Degree4 >), , 6, 4, class, TEMPLATE_ARGS(std::complex< double >, std::complex< float >, double, float), TEMPLATE_PARAMS(const typename MathObject< TEMPLATE_PARAM_1, Degree::Degree4 >::Object &, const GeneralEvaluableObject< TEMPLATE_PARAM_1, Degree::Degree0 > &))


  template< class T_Scalar >
  Function< T_Scalar, Degree::Degree4 > operator*(const GeneralEvaluableObject< T_Scalar, Degree::Degree4 > &a, const typename MathObject< T_Scalar, Degree::Degree0 >::Object &b)
  {
    return a.getEvaluableFunction() * Function< T_Scalar, Degree::Degree0 >(b);
  }

  INSTANTIATE_OPT(TEMPLATE_RETURN(Function< TEMPLATE_PARAM_1, Degree::Degree4 >), , 7, 4, class, TEMPLATE_ARGS(std::complex< double >, std::complex< float >, double, float), TEMPLATE_PARAMS(const GeneralEvaluableObject< TEMPLATE_PARAM_1, Degree::Degree4 > &, const typename MathObject< TEMPLATE_PARAM_1, Degree::Degree0 >::Object &))
  
  
  template< class T_Scalar >
  Function< T_Scalar, Degree::Degree4 > operator*(const Function< T_Scalar, Degree::Degree4 > &a, const Function< T_Scalar, Degree::Degree4 > &b)
  {
    if(a.isConstant(std::make_pair(-1, -1)) && b.isConstant(std::make_pair(-1, -1))) {
      typename MathObject< T_Scalar, Degree::Degree4 >::Object valueA;
      typename MathObject< T_Scalar, Degree::Degree4 >::Object valueB;
      a.evaluate(valueA, 0., 0., 0., std::make_pair(-1, -1));
      b.evaluate(valueB, 0., 0., 0., std::make_pair(-1, -1));
      typename MathObject< T_Scalar, Degree::Degree4 >::Object tmp;
      for(unsigned int I = 0; I < 3; ++I) {
        for(unsigned int J = 0; J < 3; ++J) {
          for(unsigned int M = 0; M < 3; ++M) {
            for(unsigned int N = 0; N < 3; ++N) {
              tmp(I, J)(M, N) = 0.;
              for(unsigned int K = 0; K < 3; ++K) {
                for(unsigned int L = 0; L < 3; ++L) {
                  tmp(I, J)(M, N) += valueA(I, J)(K, L) * valueB(K, L)(M, N);
                }
              }
            }
          }
        }
      }
      return Function< T_Scalar, Degree::Degree4 >(tmp);
    }
    return Function< T_Scalar, Degree::Degree4 >(new BinaryNode< Mul< T_Scalar, Degree::Degree4, Degree::Degree4 > >(a.copy(), b.copy()));
  }

  INSTANTIATE_OPT(TEMPLATE_RETURN(Function< TEMPLATE_PARAM_1, Degree::Degree4 >), , 8, 4, class, TEMPLATE_ARGS(std::complex< double >, std::complex< float >, double, float), TEMPLATE_PARAMS(const Function< TEMPLATE_PARAM_1, Degree::Degree4 > &, const Function< TEMPLATE_PARAM_1, Degree::Degree4 > &))


  template< class T_Scalar >
  Function< T_Scalar, Degree::Degree4 > operator*(const GeneralEvaluableObject< T_Scalar, Degree::Degree4 > &a, const GeneralEvaluableObject< T_Scalar, Degree::Degree4 > &b)
  {
    return a.getEvaluableFunction() * b.getEvaluableFunction();
  }

  INSTANTIATE_OPT(TEMPLATE_RETURN(Function< TEMPLATE_PARAM_1, Degree::Degree4 >), , 9, 4, class, TEMPLATE_ARGS(std::complex< double >, std::complex< float >, double, float), TEMPLATE_PARAMS(const GeneralEvaluableObject< TEMPLATE_PARAM_1, Degree::Degree4 > &, const GeneralEvaluableObject< TEMPLATE_PARAM_1, Degree::Degree4 > &))


  template< class T_Scalar >
  Function< T_Scalar, Degree::Degree4 > operator*(const typename MathObject< T_Scalar, Degree::Degree4 >::Object &a, const GeneralEvaluableObject< T_Scalar, Degree::Degree4 > &b)
  {
    return Function< T_Scalar, Degree::Degree4 >(a) * b.getEvaluableFunction();
  }

  INSTANTIATE_OPT(TEMPLATE_RETURN(Function< TEMPLATE_PARAM_1, Degree::Degree4 >), , 10, 4, class, TEMPLATE_ARGS(std::complex< double >, std::complex< float >, double, float), TEMPLATE_PARAMS(const typename MathObject< TEMPLATE_PARAM_1, Degree::Degree4 >::Object &, const GeneralEvaluableObject< TEMPLATE_PARAM_1, Degree::Degree4 > &))


  template< class T_Scalar >
  Function< T_Scalar, Degree::Degree4 > operator*(const GeneralEvaluableObject< T_Scalar, Degree::Degree4 > &a, const typename MathObject< T_Scalar, Degree::Degree4 >::Object &b)
  {
    return a.getEvaluableFunction() * Function< T_Scalar, Degree::Degree4 >(b);
  }

  INSTANTIATE_OPT(TEMPLATE_RETURN(Function< TEMPLATE_PARAM_1, Degree::Degree4 >), , 11, 4, class, TEMPLATE_ARGS(std::complex< double >, std::complex< float >, double, float), TEMPLATE_PARAMS(const GeneralEvaluableObject< TEMPLATE_PARAM_1, Degree::Degree4 > &, const typename MathObject< TEMPLATE_PARAM_1, Degree::Degree4 >::Object &))


  // operator/
  template< class T_Scalar >
  Function< T_Scalar, Degree::Degree4 > operator/(const Function< T_Scalar, Degree::Degree4 > &a, const Function< T_Scalar, Degree::Degree0 > &b)
  {
    if(a.isConstant(std::make_pair(-1, -1)) && b.isConstant(std::make_pair(-1, -1))) {
      typename MathObject< T_Scalar, Degree::Degree4 >::Object valueA;
      typename MathObject< T_Scalar, Degree::Degree0 >::Object valueB;
      a.evaluate(valueA, 0., 0., 0., std::make_pair(-1, -1));
      b.evaluate(valueB, 0., 0., 0., std::make_pair(-1, -1));
      typename MathObject< T_Scalar, Degree::Degree4 >::Object tmp;
      for(unsigned int I = 0; I < 3; ++I) {
        for(unsigned int J = 0; J < 3; ++J) {
          tmp(I, J) = valueA(I, J) / valueB;
        }
      }
      return Function< T_Scalar, Degree::Degree4 >(tmp);
    }
    return Function< T_Scalar, Degree::Degree4 >(new BinaryNode< Div< T_Scalar, Degree::Degree4 > >(a.copy(), b.copy()));
  }

  INSTANTIATE_OPD(TEMPLATE_RETURN(Function< TEMPLATE_PARAM_1, Degree::Degree4 >), , 0, 4, class, TEMPLATE_ARGS(std::complex< double >, std::complex< float >, double, float), TEMPLATE_PARAMS(const Function< TEMPLATE_PARAM_1, Degree::Degree4 > &, const Function< TEMPLATE_PARAM_1, Degree::Degree0 > &))


  template< class T_Scalar >
  Function< T_Scalar, Degree::Degree4 > operator/(const GeneralEvaluableObject< T_Scalar, Degree::Degree4 > &a, const GeneralEvaluableObject< T_Scalar, Degree::Degree0 > &b)
  {
    return a.getEvaluableFunction() / b.getEvaluableFunction();
  }

  INSTANTIATE_OPD(TEMPLATE_RETURN(Function< TEMPLATE_PARAM_1, Degree::Degree4 >), , 1, 4, class, TEMPLATE_ARGS(std::complex< double >, std::complex< float >, double, float), TEMPLATE_PARAMS(const GeneralEvaluableObject< TEMPLATE_PARAM_1, Degree::Degree4 > &, const GeneralEvaluableObject< TEMPLATE_PARAM_1, Degree::Degree0 > &))


  template< class T_Scalar >
  Function< T_Scalar, Degree::Degree4 > operator/(const typename MathObject< T_Scalar, Degree::Degree4 >::Object &a, const GeneralEvaluableObject< T_Scalar, Degree::Degree0 > &b)
  {
    return Function< T_Scalar, Degree::Degree4 >(a) / b.getEvaluableFunction();
  }

  INSTANTIATE_OPD(TEMPLATE_RETURN(Function< TEMPLATE_PARAM_1, Degree::Degree4 >), , 2, 4, class, TEMPLATE_ARGS(std::complex< double >, std::complex< float >, double, float), TEMPLATE_PARAMS(const typename MathObject< TEMPLATE_PARAM_1, Degree::Degree4 >::Object &, const GeneralEvaluableObject< TEMPLATE_PARAM_1, Degree::Degree0 > &))


  template< class T_Scalar >
  Function< T_Scalar, Degree::Degree4 > operator/(const GeneralEvaluableObject< T_Scalar, Degree::Degree4 > &a, const typename MathObject< T_Scalar, Degree::Degree0 >::Object &b)
  {
    return a.getEvaluableFunction() / Function< T_Scalar, Degree::Degree0 >(b);
  }

  INSTANTIATE_OPD(TEMPLATE_RETURN(Function< TEMPLATE_PARAM_1, Degree::Degree4 >), , 3, 4, class, TEMPLATE_ARGS(std::complex< double >, std::complex< float >, double, float), TEMPLATE_PARAMS(const GeneralEvaluableObject< TEMPLATE_PARAM_1, Degree::Degree4 > &a, const typename MathObject< TEMPLATE_PARAM_1, Degree::Degree0 >::Object &))


} // namespace gmshfem::function
