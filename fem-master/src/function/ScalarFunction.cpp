// GmshFEM - Copyright (C) 2019-2022, A. Royer, E. Béchet, C. Geuzaine, Université de Liège
//
// See the LICENSE.txt file for license information. Please report all
// issues on https://gitlab.onelab.info/gmsh/fem/issues

#include "ScalarFunction.h"

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
#include "unaryOperations.h"

namespace gmshfem::function
{


  template< class T_Scalar >
  Function< T_Scalar, Degree::Degree0 >::Function() :
    _tree(nullptr)
  {
  }

  template< class T_Scalar >
  Function< T_Scalar, Degree::Degree0 >::Function(const ExecutionTreeWithDegree< T_Scalar, Degree::Degree0 > *tree) :
    _tree(tree)
  {
  }

  template< class T_Scalar >
  Function< T_Scalar, Degree::Degree0 >::Function(const Function &other) :
    _tree(other._tree ? other._tree->copy() : nullptr)
  {
  }

  template< class T_Scalar >
  Function< T_Scalar, Degree::Degree0 >::Function(Function &&other) :
    _tree(other._tree)
  {
    other._tree = nullptr;
  }

  template< class T_Scalar >
  template< class T_SFINAE, class >
  Function< T_Scalar, Degree::Degree0 >::Function(const Function< scalar::Precision< T_SFINAE >, Degree::Degree0 > &other) :
    _tree(new ScalarTypeNode< RealToComplex< scalar::Precision< T_SFINAE >, Degree::Degree0 > >(other.copy()))
  {
  }

  template< class T_Scalar >
  Function< T_Scalar, Degree::Degree0 >::Function(const typename MathObject< T_Scalar, Degree::Degree0 >::Object &value) :
    _tree(new function::NullaryNode< function::Constant< T_Scalar, Degree::Degree0 > >(function::Constant< T_Scalar, Degree::Degree0 >(value)))
  {
  }

  template< class T_Scalar >
  template< class T_SFINAE, class >
  Function< T_Scalar, Degree::Degree0 >::Function(const typename MathObject< scalar::Precision< T_SFINAE >, Degree::Degree0 >::Object &value) :
    _tree(new function::NullaryNode< function::Constant< T_SFINAE, Degree::Degree0 > >(function::Constant< T_SFINAE, Degree::Degree0 >(T_SFINAE(value))))
  {
  }

  template< class T_Scalar >
  Function< T_Scalar, Degree::Degree0 >::Function(const field::Field< T_Scalar, field::Form::Form0 > &field) :
    _tree(new FieldNode< None< T_Scalar, field::Form::Form0 > >(&field))
  {
  }

  template< class T_Scalar >
  Function< T_Scalar, Degree::Degree0 >::Function(const field::Field< T_Scalar, field::Form::Form3 > &field) :
    _tree(new FieldNode< None< T_Scalar, field::Form::Form3 > >(&field))
  {
  }

  template< class T_Scalar >
  Function< T_Scalar, Degree::Degree0 > &Function< T_Scalar, Degree::Degree0 >::operator=(const Function &other)
  {
    const ExecutionTreeWithDegree< T_Scalar, Degree::Degree0 > *tmp = _tree;
    _tree = (other._tree ? other._tree->copy() : nullptr);

    if(tmp) {
      delete tmp;
      tmp = nullptr;
    }

    return *this;
  }

  template< class T_Scalar >
  Function< T_Scalar, Degree::Degree0 > &Function< T_Scalar, Degree::Degree0 >::operator=(Function &&other)
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
  Function< T_Scalar, Degree::Degree0 > &Function< T_Scalar, Degree::Degree0 >::operator=(const Function< scalar::Precision< T_SFINAE >, Degree::Degree0 > &other)
  {
    if(_tree) {
      delete _tree;
    }
    const ExecutionTreeWithDegree< scalar::Precision< T_SFINAE >, Degree::Degree0 > *otherTree = other.copy();
    if(otherTree) {
      _tree = new ScalarTypeNode< RealToComplex< scalar::Precision< T_SFINAE >, Degree::Degree0 > >(otherTree);
    }
    else {
      _tree = nullptr;
    }
    return *this;
  }

  template< class T_Scalar >
  Function< T_Scalar, Degree::Degree0 > &Function< T_Scalar, Degree::Degree0 >::operator=(const typename MathObject< T_Scalar, Degree::Degree0 >::Object &value)
  {
    if(_tree) {
      delete _tree;
    }
    _tree = new function::NullaryNode< function::Constant< T_Scalar, Degree::Degree0 > >(function::Constant< T_Scalar, Degree::Degree0 >(value));
    return *this;
  }

  template< class T_Scalar >
  template< class T_SFINAE, class >
  Function< T_Scalar, Degree::Degree0 > &Function< T_Scalar, Degree::Degree0 >::operator=(const typename MathObject< scalar::Precision< T_SFINAE >, Degree::Degree0 >::Object &value)
  {
    if(_tree) {
      delete _tree;
    }
    _tree = new function::NullaryNode< function::Constant< T_SFINAE, Degree::Degree0 > >(function::Constant< T_SFINAE, Degree::Degree0 >(T_Scalar(value)));
    return *this;
  }

  template< class T_Scalar >
  Function< T_Scalar, Degree::Degree0 > &Function< T_Scalar, Degree::Degree0 >::operator=(const field::Field< T_Scalar, field::Form::Form0 > &field)
  {
    if(_tree) {
      delete _tree;
    }

    _tree = new FieldNode< None< T_Scalar, field::Form::Form0 > >(&field);
    return *this;
  }

  template< class T_Scalar >
  Function< T_Scalar, Degree::Degree0 > &Function< T_Scalar, Degree::Degree0 >::operator=(const field::Field< T_Scalar, field::Form::Form3 > &field)
  {
    if(_tree) {
      delete _tree;
    }

    _tree = new FieldNode< None< T_Scalar, field::Form::Form3 > >(&field);
    return *this;
  }

  template< class T_Scalar >
  Function< T_Scalar, Degree::Degree0 >::~Function()
  {
    if(_tree) {
      delete _tree;
    }
  }

  template< class T_Scalar >
  ExecutionTreeIterator Function< T_Scalar, Degree::Degree0 >::begin() const
  {
    return _tree->begin();
  }

  template< class T_Scalar >
  ExecutionTreeIterator Function< T_Scalar, Degree::Degree0 >::end() const
  {
    return _tree->end();
  }

  template< class T_Scalar >
  const ExecutionTreeWithDegree< T_Scalar, Degree::Degree0 > *Function< T_Scalar, Degree::Degree0 >::getTree() const
  {
    return _tree;
  }

  template< class T_Scalar >
  const ExecutionTreeWithDegree< T_Scalar, Degree::Degree0 > *Function< T_Scalar, Degree::Degree0 >::copy() const
  {
    if(_tree) {
      return _tree->copy();
    }
    return nullptr;
  }

  template< class T_Scalar >
  void Function< T_Scalar, Degree::Degree0 >::evaluate(std::vector< typename MathObject< T_Scalar, Degree::Degree0 >::Object, numa::allocator< typename MathObject< T_Scalar, Degree::Degree0 >::Object > > &values, const std::vector< scalar::Precision< T_Scalar >, numa::allocator< scalar::Precision< T_Scalar > > > &points, const std::vector< scalar::Precision< T_Scalar > > &gaussPoints, const int elementType, const std::pair< int, int > &entity) const
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
      FunctionAllocator< typename MathObject< T_Scalar, Degree::Degree0 >::Object > allocator(&values[0]);
      OutputVector< T_Scalar, Degree::Degree0 > valuesOV(points.size() / 3, allocator);
#pragma omp barrier
      _tree->evaluate(valuesOV, points, gaussPoints, elementType, entity);
#pragma omp barrier
#pragma omp single
      MemoryPoolAllocator::instance()->shrinkTo(GMSHFEM_FUNCTION_MEMORY_RESIDUAL);
    }
  }

  template< class T_Scalar >
  void Function< T_Scalar, Degree::Degree0 >::evaluate(typename MathObject< T_Scalar, Degree::Degree0 >::Object &value, const scalar::Precision< T_Scalar > x, const scalar::Precision< T_Scalar > y, const scalar::Precision< T_Scalar > z, const std::pair< int, int > &entity) const
  {
    if(_tree) {
      _tree->evaluate(value, x, y, z, entity);
      MemoryPoolAllocator::instance()->shrinkTo(GMSHFEM_FUNCTION_MEMORY_RESIDUAL);
    }
  }

  template< class T_Scalar >
  common::Timer Function< T_Scalar, Degree::Degree0 >::evaluate(typename std::vector< std::vector< typename MathObject< T_Scalar, Degree::Degree0 >::Object, numa::allocator< typename MathObject< T_Scalar, Degree::Degree0 >::Object > > > &values, const domain::Domain &domain, const std::string &gauss) const
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
          FunctionAllocator< typename MathObject< T_Scalar, Degree::Degree0 >::Object > allocator(&values[i][0]);
          OutputVector< T_Scalar, Degree::Degree0 > valuesOV(points.size() / 3, allocator);
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
  bool Function< T_Scalar, Degree::Degree0 >::isConstant(const std::pair< int, int > &entity) const
  {
    if(_tree) {
      return _tree->isConstant(entity);
    }
    return true;
  }

  template< class T_Scalar >
  void Function< T_Scalar, Degree::Degree0 >::setPointEvaluation() const
  {
    for(auto it = begin(); it != end(); ++it) {
      function::setPointEvaluation(it);
    }
  }

  template< class T_Scalar >
  void Function< T_Scalar, Degree::Degree0 >::exportExecutionTree(const domain::Domain &domain, const std::string &name, const std::string &path)
  {
    common::HTMLTree html(name, path);
    common::HTMLNode *node = html.getHead();

    _tree->exportTree(domain, node);
    common::Memory peakMemory = _tree->peakMemoryByGaussPoint(domain);
    peakMemory += sizeof(typename MathObject< T_Scalar, Degree::Degree0 >::Object);

    html.setPeakMemory(peakMemory);
    html.write();
  }

  template< class T_Scalar >
  const Function< T_Scalar, Degree::Degree0 > &Function< T_Scalar, Degree::Degree0 >::getEvaluableFunction() const
  {
    return *this;
  }

  INSTANTIATE_CLASS_2(Function, 4, 1, TEMPLATE_ARGS(std::complex< double >, std::complex< float >, double, float), TEMPLATE_ARGS(Degree::Degree0))

  template Function< std::complex< double >, Degree::Degree0 >::Function(const Function< scalar::Precision< std::complex< double > >, Degree::Degree0 > &other);
  template Function< std::complex< double >, Degree::Degree0 >::Function(const typename MathObject< scalar::Precision< std::complex< double > >, Degree::Degree0 >::Object &value);
  template Function< std::complex< double >, Degree::Degree0 > &Function< std::complex< double >, Degree::Degree0 >::operator=(const typename MathObject< scalar::Precision< std::complex< double > >, Degree::Degree0 >::Object &value);
  template Function< std::complex< double >, Degree::Degree0 > &Function< std::complex< double >, Degree::Degree0 >::operator=(const Function< scalar::Precision< std::complex< double > >, Degree::Degree0 > &other);

  template Function< std::complex< float >, Degree::Degree0 >::Function(const Function< scalar::Precision< std::complex< float > >, Degree::Degree0 > &other);
  template Function< std::complex< float >, Degree::Degree0 >::Function(const typename MathObject< scalar::Precision< std::complex< float > >, Degree::Degree0 >::Object &value);
  template Function< std::complex< float >, Degree::Degree0 > &Function< std::complex< float >, Degree::Degree0 >::operator=(const typename MathObject< scalar::Precision< std::complex< float > >, Degree::Degree0 >::Object &value);
  template Function< std::complex< float >, Degree::Degree0 > &Function< std::complex< float >, Degree::Degree0 >::operator=(const Function< scalar::Precision< std::complex< float > >, Degree::Degree0 > &other);


  namespace literals
  {

    Function< std::complex< double >, Degree::Degree0 > operator"" _cd_sf(const long double value)
    {
      return Function< std::complex< double >, Degree::Degree0 >(value);
    }

    Function< double, Degree::Degree0 > operator"" _d_sf(const long double value)
    {
      return Function< double, Degree::Degree0 >(value);
    }

    Function< std::complex< float >, Degree::Degree0 > operator"" _cf_sf(const long double value)
    {
      return Function< std::complex< float >, Degree::Degree0 >(value);
    }

    Function< float, Degree::Degree0 > operator"" _f_sf(const long double value)
    {
      return Function< float, Degree::Degree0 >(value);
    }

  } // namespace literals


  // operator+
  template< class T_Scalar >
  Function< T_Scalar, Degree::Degree0 > operator+(const Function< T_Scalar, Degree::Degree0 > &a, const Function< T_Scalar, Degree::Degree0 > &b)
  {
    if(a.isConstant(std::make_pair(-1, -1)) && b.isConstant(std::make_pair(-1, -1))) {
      typename MathObject< T_Scalar, Degree::Degree0 >::Object valueA, valueB;
      a.evaluate(valueA, 0., 0., 0., std::make_pair(-1, -1));
      b.evaluate(valueB, 0., 0., 0., std::make_pair(-1, -1));
      return Function< T_Scalar, Degree::Degree0 >(valueA + valueB);
    }
    return Function< T_Scalar, Degree::Degree0 >(new BinaryNode< Add< T_Scalar, Degree::Degree0 > >(a.copy(), b.copy()));
  }

  INSTANTIATE_OPP(TEMPLATE_RETURN(Function< TEMPLATE_PARAM_1, Degree::Degree0 >), , 0, 4, class, TEMPLATE_ARGS(std::complex< double >, std::complex< float >, double, float), TEMPLATE_PARAMS(const Function< TEMPLATE_PARAM_1, Degree::Degree0 > &, const Function< TEMPLATE_PARAM_1, Degree::Degree0 > &))


  template< class T_Scalar >
  Function< T_Scalar, Degree::Degree0 > operator+(const GeneralEvaluableObject< T_Scalar, Degree::Degree0 > &a, const GeneralEvaluableObject< T_Scalar, Degree::Degree0 > &b)
  {
    return a.getEvaluableFunction() + b.getEvaluableFunction();
  }

  INSTANTIATE_OPP(TEMPLATE_RETURN(Function< TEMPLATE_PARAM_1, Degree::Degree0 >), , 1, 4, class, TEMPLATE_ARGS(std::complex< double >, std::complex< float >, double, float), TEMPLATE_PARAMS(const GeneralEvaluableObject< TEMPLATE_PARAM_1, Degree::Degree0 > &, const GeneralEvaluableObject< TEMPLATE_PARAM_1, Degree::Degree0 > &))


  template< class T_Scalar >
  Function< T_Scalar, Degree::Degree0 > operator+(const typename MathObject< T_Scalar, Degree::Degree0 >::Object &a, const GeneralEvaluableObject< T_Scalar, Degree::Degree0 > &b)
  {
    return Function< T_Scalar, Degree::Degree0 >(a) + b.getEvaluableFunction();
  }

  INSTANTIATE_OPP(TEMPLATE_RETURN(Function< TEMPLATE_PARAM_1, Degree::Degree0 >), , 2, 4, class, TEMPLATE_ARGS(std::complex< double >, std::complex< float >, double, float), TEMPLATE_PARAMS(const typename MathObject< TEMPLATE_PARAM_1, Degree::Degree0 >::Object &, const GeneralEvaluableObject< TEMPLATE_PARAM_1, Degree::Degree0 > &))


  template< class T_Scalar >
  Function< T_Scalar, Degree::Degree0 > operator+(const GeneralEvaluableObject< T_Scalar, Degree::Degree0 > &a, const typename MathObject< T_Scalar, Degree::Degree0 >::Object &b)
  {
    return a.getEvaluableFunction() + Function< T_Scalar, Degree::Degree0 >(b);
  }

  INSTANTIATE_OPP(TEMPLATE_RETURN(Function< TEMPLATE_PARAM_1, Degree::Degree0 >), , 3, 4, class, TEMPLATE_ARGS(std::complex< double >, std::complex< float >, double, float), TEMPLATE_PARAMS(const GeneralEvaluableObject< TEMPLATE_PARAM_1, Degree::Degree0 > &, const typename MathObject< TEMPLATE_PARAM_1, Degree::Degree0 >::Object &))


  // operator-
  template< class T_Scalar >
  Function< T_Scalar, Degree::Degree0 > operator-(const Function< T_Scalar, Degree::Degree0 > &a, const Function< T_Scalar, Degree::Degree0 > &b)
  {
    if(a.isConstant(std::make_pair(-1, -1)) && b.isConstant(std::make_pair(-1, -1))) {
      typename MathObject< T_Scalar, Degree::Degree0 >::Object valueA, valueB;
      a.evaluate(valueA, 0., 0., 0., std::make_pair(-1, -1));
      b.evaluate(valueB, 0., 0., 0., std::make_pair(-1, -1));
      return Function< T_Scalar, Degree::Degree0 >(valueA - valueB);
    }
    return Function< T_Scalar, Degree::Degree0 >(new BinaryNode< Sub< T_Scalar, Degree::Degree0 > >(a.copy(), b.copy()));
  }

  INSTANTIATE_OPM(TEMPLATE_RETURN(Function< TEMPLATE_PARAM_1, Degree::Degree0 >), , 0, 4, class, TEMPLATE_ARGS(std::complex< double >, std::complex< float >, double, float), TEMPLATE_PARAMS(const Function< TEMPLATE_PARAM_1, Degree::Degree0 > &, const Function< TEMPLATE_PARAM_1, Degree::Degree0 > &))


  template< class T_Scalar >
  Function< T_Scalar, Degree::Degree0 > operator-(const GeneralEvaluableObject< T_Scalar, Degree::Degree0 > &a, const GeneralEvaluableObject< T_Scalar, Degree::Degree0 > &b)
  {
    return a.getEvaluableFunction() - b.getEvaluableFunction();
  }

  INSTANTIATE_OPM(TEMPLATE_RETURN(Function< TEMPLATE_PARAM_1, Degree::Degree0 >), , 1, 4, class, TEMPLATE_ARGS(std::complex< double >, std::complex< float >, double, float), TEMPLATE_PARAMS(const GeneralEvaluableObject< TEMPLATE_PARAM_1, Degree::Degree0 > &, const GeneralEvaluableObject< TEMPLATE_PARAM_1, Degree::Degree0 > &))


  template< class T_Scalar >
  Function< T_Scalar, Degree::Degree0 > operator-(const typename MathObject< T_Scalar, Degree::Degree0 >::Object &a, const GeneralEvaluableObject< T_Scalar, Degree::Degree0 > &b)
  {
    return Function< T_Scalar, Degree::Degree0 >(a) - b.getEvaluableFunction();
  }

  INSTANTIATE_OPM(TEMPLATE_RETURN(Function< TEMPLATE_PARAM_1, Degree::Degree0 >), , 2, 4, class, TEMPLATE_ARGS(std::complex< double >, std::complex< float >, double, float), TEMPLATE_PARAMS(const typename MathObject< TEMPLATE_PARAM_1, Degree::Degree0 >::Object &, const GeneralEvaluableObject< TEMPLATE_PARAM_1, Degree::Degree0 > &))


  template< class T_Scalar >
  Function< T_Scalar, Degree::Degree0 > operator-(const GeneralEvaluableObject< T_Scalar, Degree::Degree0 > &a, const typename MathObject< T_Scalar, Degree::Degree0 >::Object &b)
  {
    return a.getEvaluableFunction() - Function< T_Scalar, Degree::Degree0 >(b);
  }

  INSTANTIATE_OPM(TEMPLATE_RETURN(Function< TEMPLATE_PARAM_1, Degree::Degree0 >), , 3, 4, class, TEMPLATE_ARGS(std::complex< double >, std::complex< float >, double, float), TEMPLATE_PARAMS(const GeneralEvaluableObject< TEMPLATE_PARAM_1, Degree::Degree0 > &, const typename MathObject< TEMPLATE_PARAM_1, Degree::Degree0 >::Object &))


  // operator*
  template< class T_Scalar >
  Function< T_Scalar, Degree::Degree0 > operator*(const Function< T_Scalar, Degree::Degree0 > &a, const Function< T_Scalar, Degree::Degree0 > &b)
  {
    if(a.isConstant(std::make_pair(-1, -1)) && b.isConstant(std::make_pair(-1, -1))) {
      typename MathObject< T_Scalar, Degree::Degree0 >::Object valueA, valueB;
      a.evaluate(valueA, 0., 0., 0., std::make_pair(-1, -1));
      b.evaluate(valueB, 0., 0., 0., std::make_pair(-1, -1));
      return Function< T_Scalar, Degree::Degree0 >(valueA * valueB);
    }
    return Function< T_Scalar, Degree::Degree0 >(new BinaryNode< Mul< T_Scalar, Degree::Degree0, Degree::Degree0 > >(a.copy(), b.copy()));
  }

  INSTANTIATE_OPT(TEMPLATE_RETURN(Function< TEMPLATE_PARAM_1, Degree::Degree0 >), , 0, 4, class, TEMPLATE_ARGS(std::complex< double >, std::complex< float >, double, float), TEMPLATE_PARAMS(const Function< TEMPLATE_PARAM_1, Degree::Degree0 > &, const Function< TEMPLATE_PARAM_1, Degree::Degree0 > &))

  template< class T_Scalar >
  Function< T_Scalar, Degree::Degree0 > operator*(const GeneralEvaluableObject< T_Scalar, Degree::Degree0 > &a, const GeneralEvaluableObject< T_Scalar, Degree::Degree0 > &b)
  {
    return a.getEvaluableFunction() * b.getEvaluableFunction();
  }

  INSTANTIATE_OPT(TEMPLATE_RETURN(Function< TEMPLATE_PARAM_1, Degree::Degree0 >), , 1, 4, class, TEMPLATE_ARGS(std::complex< double >, std::complex< float >, double, float), TEMPLATE_PARAMS(const GeneralEvaluableObject< TEMPLATE_PARAM_1, Degree::Degree0 > &, const GeneralEvaluableObject< TEMPLATE_PARAM_1, Degree::Degree0 > &))


  template< class T_Scalar >
  Function< T_Scalar, Degree::Degree0 > operator*(const typename MathObject< T_Scalar, Degree::Degree0 >::Object &a, const GeneralEvaluableObject< T_Scalar, Degree::Degree0 > &b)
  {
    return Function< T_Scalar, Degree::Degree0 >(a) * b.getEvaluableFunction();
  }

  INSTANTIATE_OPT(TEMPLATE_RETURN(Function< TEMPLATE_PARAM_1, Degree::Degree0 >), , 2, 4, class, TEMPLATE_ARGS(std::complex< double >, std::complex< float >, double, float), TEMPLATE_PARAMS(const typename MathObject< TEMPLATE_PARAM_1, Degree::Degree0 >::Object &, const GeneralEvaluableObject< TEMPLATE_PARAM_1, Degree::Degree0 > &))


  template< class T_Scalar >
  Function< T_Scalar, Degree::Degree0 > operator*(const GeneralEvaluableObject< T_Scalar, Degree::Degree0 > &a, const typename MathObject< T_Scalar, Degree::Degree0 >::Object &b)
  {
    return a.getEvaluableFunction() * Function< T_Scalar, Degree::Degree0 >(b);
  }

  INSTANTIATE_OPT(TEMPLATE_RETURN(Function< TEMPLATE_PARAM_1, Degree::Degree0 >), , 3, 4, class, TEMPLATE_ARGS(std::complex< double >, std::complex< float >, double, float), TEMPLATE_PARAMS(const GeneralEvaluableObject< TEMPLATE_PARAM_1, Degree::Degree0 > &, const typename MathObject< TEMPLATE_PARAM_1, Degree::Degree0 >::Object &))


  template< class T_Scalar >
  Function< T_Scalar, Degree::Degree0 > operator*(const Function< T_Scalar, Degree::Degree1 > &a, const Function< T_Scalar, Degree::Degree1 > &b)
  {
    if(a.isConstant(std::make_pair(-1, -1)) && b.isConstant(std::make_pair(-1, -1))) {
      typename MathObject< T_Scalar, Degree::Degree1 >::Object valueA, valueB;
      a.evaluate(valueA, 0., 0., 0., std::make_pair(-1, -1));
      b.evaluate(valueB, 0., 0., 0., std::make_pair(-1, -1));
      return Function< T_Scalar, Degree::Degree0 >((valueA.transpose() * valueB)(0));
    }
    return Function< T_Scalar, Degree::Degree0 >(new BinaryNode< Mul< T_Scalar, Degree::Degree1, Degree::Degree1 > >(a.copy(), b.copy()));
  }

  INSTANTIATE_OPT(TEMPLATE_RETURN(Function< TEMPLATE_PARAM_1, Degree::Degree0 >), , 4, 4, class, TEMPLATE_ARGS(std::complex< double >, std::complex< float >, double, float), TEMPLATE_PARAMS(const Function< TEMPLATE_PARAM_1, Degree::Degree1 > &, const Function< TEMPLATE_PARAM_1, Degree::Degree1 > &))


  template< class T_Scalar >
  Function< T_Scalar, Degree::Degree0 > operator*(const GeneralEvaluableObject< T_Scalar, Degree::Degree1 > &a, const GeneralEvaluableObject< T_Scalar, Degree::Degree1 > &b)
  {
    return a.getEvaluableFunction() * b.getEvaluableFunction();
  }

  INSTANTIATE_OPT(TEMPLATE_RETURN(Function< TEMPLATE_PARAM_1, Degree::Degree0 >), , 5, 4, class, TEMPLATE_ARGS(std::complex< double >, std::complex< float >, double, float), TEMPLATE_PARAMS(const GeneralEvaluableObject< TEMPLATE_PARAM_1, Degree::Degree1 > &, const GeneralEvaluableObject< TEMPLATE_PARAM_1, Degree::Degree1 > &))


  template< class T_Scalar >
  Function< T_Scalar, Degree::Degree0 > operator*(const typename MathObject< T_Scalar, Degree::Degree1 >::Object &a, const GeneralEvaluableObject< T_Scalar, Degree::Degree1 > &b)
  {
    return Function< T_Scalar, Degree::Degree1 >(a) * b.getEvaluableFunction();
  }

  INSTANTIATE_OPT(TEMPLATE_RETURN(Function< TEMPLATE_PARAM_1, Degree::Degree0 >), , 6, 4, class, TEMPLATE_ARGS(std::complex< double >, std::complex< float >, double, float), TEMPLATE_PARAMS(const typename MathObject< TEMPLATE_PARAM_1, Degree::Degree1 >::Object &, const GeneralEvaluableObject< TEMPLATE_PARAM_1, Degree::Degree1 > &))


  template< class T_Scalar >
  Function< T_Scalar, Degree::Degree0 > operator*(const GeneralEvaluableObject< T_Scalar, Degree::Degree1 > &a, const typename MathObject< T_Scalar, Degree::Degree1 >::Object &b)
  {
    return a.getEvaluableFunction() * Function< T_Scalar, Degree::Degree1 >(b);
  }

  INSTANTIATE_OPT(TEMPLATE_RETURN(Function< TEMPLATE_PARAM_1, Degree::Degree0 >), , 7, 4, class, TEMPLATE_ARGS(std::complex< double >, std::complex< float >, double, float), TEMPLATE_PARAMS(const GeneralEvaluableObject< TEMPLATE_PARAM_1, Degree::Degree1 > &, const typename MathObject< TEMPLATE_PARAM_1, Degree::Degree1 >::Object &))


  // operator/
  template< class T_Scalar >
  Function< T_Scalar, Degree::Degree0 > operator/(const Function< T_Scalar, Degree::Degree0 > &a, const Function< T_Scalar, Degree::Degree0 > &b)
  {
    if(a.isConstant(std::make_pair(-1, -1)) && b.isConstant(std::make_pair(-1, -1))) {
      typename MathObject< T_Scalar, Degree::Degree0 >::Object valueA, valueB;
      a.evaluate(valueA, 0., 0., 0., std::make_pair(-1, -1));
      b.evaluate(valueB, 0., 0., 0., std::make_pair(-1, -1));
      return Function< T_Scalar, Degree::Degree0 >(valueA / valueB);
    }
    if(a.isConstant(std::make_pair(-1, -1))) {
      typename MathObject< T_Scalar, Degree::Degree0 >::Object value;
      a.evaluate(value, 0., 0., 0., std::make_pair(-1, -1));
      if(value == T_Scalar(1.)) {
        return Function< T_Scalar, Degree::Degree0 >(new UnaryNode< Inv< T_Scalar, Degree::Degree0 > >(b.copy()));
      }
    }
    return Function< T_Scalar, Degree::Degree0 >(new BinaryNode< Div< T_Scalar, Degree::Degree0 > >(a.copy(), b.copy()));
  }

  INSTANTIATE_OPD(TEMPLATE_RETURN(Function< TEMPLATE_PARAM_1, Degree::Degree0 >), , 0, 4, class, TEMPLATE_ARGS(std::complex< double >, std::complex< float >, double, float), TEMPLATE_PARAMS(const Function< TEMPLATE_PARAM_1, Degree::Degree0 > &, const Function< TEMPLATE_PARAM_1, Degree::Degree0 > &))


  template< class T_Scalar >
  Function< T_Scalar, Degree::Degree0 > operator/(const GeneralEvaluableObject< T_Scalar, Degree::Degree0 > &a, const GeneralEvaluableObject< T_Scalar, Degree::Degree0 > &b)
  {
    return a.getEvaluableFunction() / b.getEvaluableFunction();
  }

  INSTANTIATE_OPD(TEMPLATE_RETURN(Function< TEMPLATE_PARAM_1, Degree::Degree0 >), , 1, 4, class, TEMPLATE_ARGS(std::complex< double >, std::complex< float >, double, float), TEMPLATE_PARAMS(const GeneralEvaluableObject< TEMPLATE_PARAM_1, Degree::Degree0 > &, const GeneralEvaluableObject< TEMPLATE_PARAM_1, Degree::Degree0 > &))


  template< class T_Scalar >
  Function< T_Scalar, Degree::Degree0 > operator/(const typename MathObject< T_Scalar, Degree::Degree0 >::Object &a, const GeneralEvaluableObject< T_Scalar, Degree::Degree0 > &b)
  {
    return Function< T_Scalar, Degree::Degree0 >(a) / b.getEvaluableFunction();
  }

  INSTANTIATE_OPD(TEMPLATE_RETURN(Function< TEMPLATE_PARAM_1, Degree::Degree0 >), , 2, 4, class, TEMPLATE_ARGS(std::complex< double >, std::complex< float >, double, float), TEMPLATE_PARAMS(const typename MathObject< TEMPLATE_PARAM_1, Degree::Degree0 >::Object &, const GeneralEvaluableObject< TEMPLATE_PARAM_1, Degree::Degree0 > &))


  template< class T_Scalar >
  Function< T_Scalar, Degree::Degree0 > operator/(const GeneralEvaluableObject< T_Scalar, Degree::Degree0 > &a, const typename MathObject< T_Scalar, Degree::Degree0 >::Object &b)
  {
    return a.getEvaluableFunction() / Function< T_Scalar, Degree::Degree0 >(b);
  }

  INSTANTIATE_OPD(TEMPLATE_RETURN(Function< TEMPLATE_PARAM_1, Degree::Degree0 >), , 3, 4, class, TEMPLATE_ARGS(std::complex< double >, std::complex< float >, double, float), TEMPLATE_PARAMS(const GeneralEvaluableObject< TEMPLATE_PARAM_1, Degree::Degree0 > &, const typename MathObject< TEMPLATE_PARAM_1, Degree::Degree0 >::Object &))


} // namespace gmshfem::function
