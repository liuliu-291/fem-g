// GmshFEM - Copyright (C) 2019-2022, A. Royer, E. Béchet, C. Geuzaine, Université de Liège
//
// See the LICENSE.txt file for license information. Please report all
// issues on https://gitlab.onelab.info/gmsh/fem/issues

#ifndef H_GMSHFEM_FIELDOPERATIONS
#define H_GMSHFEM_FIELDOPERATIONS

#include "FieldEvaluator.h"
#include "FieldModule.h"
#include "FieldNode.h"
#include "MathObject.h"
#include "OmpInterface.h"

namespace gmshfem::function
{


  // Available field operations:
  //  Derivative
  //  DerivativeCompound
  //  None
  //  NoneCompound

  //
  // Derivative
  //

  template< class T_Scalar, field::Form T_Form >
  class Derivative final : public FieldOperation< T_Scalar, field::DegreeOfForm< field::NextForm< T_Form >::value >::value >
  {
   private:
    mutable std::vector< int > _typeKeys;
    mutable std::vector< unsigned long long > _entityKeys;
    mutable std::vector< T_Scalar > _dofsValues;
    mutable std::vector< scalar::Precision< T_Scalar > > _functions;
    mutable std::vector< int > _fsIndex;
    mutable unsigned int _nbrDofsByElements;
    mutable std::vector< double > _gmshJacobians, _gmshDeterminants, _gmshPoints, _gmshGaussPoints;
    mutable std::vector< scalar::Precision< T_Scalar >, numa::allocator< scalar::Precision< T_Scalar > > > _jacobians, _determinants;
    mutable term::evaluator::FieldEvaluator< T_Scalar, field::NextForm< T_Form >::value > _fieldEvaluator;
    mutable bool _pointEvaluation;
    mutable std::string _currentModel;

   public:
    Derivative() :
      _pointEvaluation(false)
    {
    }

    Derivative(const Derivative &other) :
      FieldOperation< T_Scalar, field::DegreeOfForm< field::NextForm< T_Form >::value >::value >(other), _pointEvaluation(other._pointEvaluation)
    {
    }

    void setPointEvaluation(const bool pointEvaluation) const
    {
      _pointEvaluation = pointEvaluation;
    }

    void operator()(const field::FieldInterface< T_Scalar > *const field, OutputVector< T_Scalar, field::DegreeOfForm< field::NextForm< T_Form >::value >::value > &values, const std::vector< scalar::Precision< T_Scalar >, numa::allocator< scalar::Precision< T_Scalar > > > &points, const std::vector< scalar::Precision< T_Scalar > > &gaussPoints, const int elementType, const std::pair< int, int > &entity) const
    {
      // If the current model is not the same on which the field is defined, we use the point evaluation algorithm
      // otherwise, we can interpolate the field as usual.
#pragma omp single
      {
        gmsh::model::getCurrent(_currentModel);
        if(_currentModel != field->model()) {
          _pointEvaluation = true;
          gmsh::model::setCurrent(field->model());
        }
      }
      if(_pointEvaluation) {
        double verbose = 0.;
#pragma omp master
        {
          gmsh::option::getNumber("General.Verbosity", verbose);
          gmsh::option::setNumber("General.Verbosity", 0);
        }
#pragma omp single
        gmsh::model::mesh::rebuildElementCache(true);
#pragma omp for
        for(auto i = 0ULL; i < points.size() / 3; ++i) {
          std::size_t elementTag = 0;
          int elementType = 0;
          std::vector< std::size_t > nodeTags;
          double uGmsh = 0., vGmsh = 0., wGmsh = 0.;
          scalar::Precision< T_Scalar > u = 0., v = 0., w = 0.;
          try {
            gmsh::model::mesh::getElementByCoordinates(points[3 * i], points[3 * i + 1], points[3 * i + 2], elementTag, elementType, nodeTags, uGmsh, vGmsh, wGmsh, field->domain().maxDim(), true);
            u = uGmsh;
            v = vGmsh;
            w = wGmsh;
          }
          catch(...) {
            MathObject< T_Scalar, field::DegreeOfForm< field::NextForm< T_Form >::value >::value >::zero(values[i]);
            continue;
          }

          std::vector< scalar::Precision< T_Scalar > > coordinates;
          std::vector< int > typeKeys;
          std::vector< unsigned long long > entityKeys;
          field->getFunctionSpace()->getKeysOnElement(false, typeKeys, entityKeys, coordinates, elementTag);

          std::vector< T_Scalar > dofsValues;
          dofsValues.resize(typeKeys.size());
          field->getValues(typeKeys, entityKeys, dofsValues, 0, dofsValues.size());

          std::vector< scalar::Precision< T_Scalar > > functions;
          const unsigned int orientation = field->getFunctionSpace()->getOrientationOfElement(elementTag);
          const unsigned int nbrDofsByElements = field->getFunctionSpace()->getBasisFunction_noexcept(true, functions, {u, v, w}, elementType, orientation);

          // Jacobians
          std::vector< double > gmshJacobians, gmshDeterminants, gmshPoints;
          std::vector< scalar::Precision< T_Scalar > > jacobians, determinants, points;
          if(_fieldEvaluator.needJacobians()) {
            gmsh::model::mesh::getJacobian(elementTag, {uGmsh, vGmsh, wGmsh}, gmshJacobians, gmshDeterminants, gmshPoints);
            scalar::move(jacobians, gmshJacobians);
            scalar::move(determinants, gmshDeterminants);
            scalar::move(points, gmshPoints);
            //              Not correct: need to find the entity tag in the field's mesh
            //              if(domain.haveJacobiansModificators(entity)) {
            //                domain.applyJacobiansModificator(points, determinants, jacobians, entity);
            //              }
          }

          // Evaluation
          Eigen::VectorX< T_Scalar > vecDof = Eigen::Map< const Eigen::VectorX< T_Scalar > >(&(dofsValues[0]), nbrDofsByElements);
          Eigen::MatrixX< scalar::Precision< T_Scalar > > fieldExpression(_fieldEvaluator.size(), nbrDofsByElements);
          _fieldEvaluator(fieldExpression, &functions[0], &jacobians[0], &determinants[0], nbrDofsByElements);
          MathObject< T_Scalar, field::DegreeOfForm< field::NextForm< T_Form >::value >::value >::copy(values[i], fieldExpression * vecDof);
        }
#pragma omp master
        {
          gmsh::option::setNumber("General.Verbosity", verbose);
        }
      }
      else {
        const unsigned int nbrOfElements = points.size() / gaussPoints.size();
        const unsigned int numThreads = omp::getNumThreads();
        const unsigned int myThreadID = omp::getThreadNum();

#pragma omp single
        _nbrDofsByElements = field->getFunctionSpace()->getBasisFunctions(true, _functions, _fsIndex, gaussPoints, elementType, entity);

        std::vector< scalar::Precision< T_Scalar > > coordinates;
#pragma omp single
        field->getFunctionSpace()->getKeys(false, _typeKeys, _entityKeys, coordinates, elementType, entity);
        const unsigned int begin = (myThreadID * _typeKeys.size()) / numThreads;
        const unsigned int end = ((myThreadID + 1) * _typeKeys.size()) / numThreads;

        const bool needOffset = !field->getFunctionSpace()->isConstantByElements();

#pragma omp single
        _dofsValues.resize(_typeKeys.size());

        field->getValues(_typeKeys, _entityKeys, _dofsValues, begin, end);
        const unsigned int nbrGaussPoints = gaussPoints.size() / 3;

#pragma omp single
        scalar::copy(_gmshGaussPoints, gaussPoints);

        if(_fieldEvaluator.needJacobians()) {
          const domain::Domain domain = field->domain();
#pragma omp single
          gmsh::model::mesh::preallocateJacobians(elementType, nbrGaussPoints, true, true, false, _gmshJacobians, _gmshDeterminants, _gmshPoints, entity.second);
          gmsh::model::mesh::getJacobians(elementType, _gmshGaussPoints, _gmshJacobians, _gmshDeterminants, _gmshPoints, entity.second, myThreadID, numThreads);
#pragma omp barrier
#pragma omp single
          numa::copy(_jacobians, _gmshJacobians);
#pragma omp single
          numa::copy(_determinants, _gmshDeterminants);

          if(domain.haveJacobiansModificators(entity)) {
            domain.applyJacobiansModificator(points, _determinants, _jacobians, entity);
#pragma omp barrier
          }
        }

        Eigen::MatrixX< scalar::Precision< T_Scalar > > fieldExpression(_fieldEvaluator.size(), _nbrDofsByElements);
#pragma omp for
        for(auto i = 0U; i < nbrOfElements; ++i) {
          Eigen::VectorX< T_Scalar > vecDof = Eigen::Map< const Eigen::VectorX< T_Scalar > >(&(_dofsValues[i * _nbrDofsByElements]), _nbrDofsByElements);
          for(auto j = 0U; j < nbrGaussPoints; ++j) {
            _fieldEvaluator(fieldExpression, &_functions[((needOffset ? _fsIndex[i] * nbrGaussPoints : 0) + j) * field::NumberOfComponentsOfForm< field::NextForm< T_Form >::value >::value * _nbrDofsByElements], &_jacobians[i * 9 * nbrGaussPoints + j * 9], &_determinants[i * nbrGaussPoints + j], _nbrDofsByElements);
            MathObject< T_Scalar, field::DegreeOfForm< field::NextForm< T_Form >::value >::value >::copy(values[i * nbrGaussPoints + j], fieldExpression * vecDof);
          }
        }
      }
      // Restore the current model
#pragma omp single
      {
        if(_currentModel != field->model()) {
          gmsh::model::setCurrent(_currentModel);
        }
      }
    }

    std::string name() const override
    {
      return "derivative of field";
    }
  };

  //
  // DerivativeCompound
  //

  template< class T_Scalar, field::Form T_Form, unsigned int T_NumFields >
  class DerivativeCompound final : public FieldOperation< T_Scalar, field::DegreeOfCompoundForm< field::NextForm< T_Form >::value >::value >
  {
   private:
    mutable std::vector< int > _typeKeys;
    mutable std::vector< unsigned long long > _entityKeys;
    mutable std::vector< T_Scalar > _dofsValues;
    mutable std::vector< scalar::Precision< T_Scalar > > _functions;
    mutable std::vector< int > _fsIndex;
    mutable unsigned int _nbrDofsByElements;
    mutable std::vector< double > _gmshJacobians, _gmshDeterminants, _gmshPoints, _gmshGaussPoints;
    mutable std::vector< scalar::Precision< T_Scalar >, numa::allocator< scalar::Precision< T_Scalar > > > _jacobians, _determinants;
    mutable term::evaluator::CompoundFieldEvaluator< T_Scalar, field::NextForm< T_Form >::value, T_NumFields > _fieldEvaluator;
    mutable bool _pointEvaluation;
    mutable std::string _currentModel;

   public:
    DerivativeCompound() :
      _pointEvaluation(false)
    {
    }

    DerivativeCompound(const DerivativeCompound &other) :
      FieldOperation< T_Scalar, field::DegreeOfCompoundForm< field::NextForm< T_Form >::value >::value >(other), _pointEvaluation(other._pointEvaluation)
    {
    }

    void setPointEvaluation(const bool pointEvaluation) const
    {
      _pointEvaluation = pointEvaluation;
    }

    void operator()(const field::FieldInterface< T_Scalar > *const field, OutputVector< T_Scalar, field::DegreeOfCompoundForm< field::NextForm< T_Form >::value >::value > &values, const std::vector< scalar::Precision< T_Scalar >, numa::allocator< scalar::Precision< T_Scalar > > > &points, const std::vector< scalar::Precision< T_Scalar > > &gaussPoints, const int elementType, const std::pair< int, int > &entity) const
    {
      // If the current model is not the same on which the field is defined, we use the point evaluation algorithm
      // otherwise, we can interpolate the field as usual.
#pragma omp single
      {
        gmsh::model::getCurrent(_currentModel);
        if(_currentModel != field->model()) {
          _pointEvaluation = true;
          gmsh::model::setCurrent(field->model());
        }
      }
      if(_pointEvaluation) {
        double verbose = 0.;
#pragma omp master
        {
          gmsh::option::getNumber("General.Verbosity", verbose);
          gmsh::option::setNumber("General.Verbosity", 0);
        }
#pragma omp single
        gmsh::model::mesh::rebuildElementCache(true);
#pragma omp for
        for(auto i = 0ULL; i < points.size() / 3; ++i) {
          std::size_t elementTag = 0;
          int elementType = 0;
          std::vector< std::size_t > nodeTags;
          double uGmsh = 0., vGmsh = 0., wGmsh = 0.;
          scalar::Precision< T_Scalar > u = 0., v = 0., w = 0.;
          try {
            gmsh::model::mesh::getElementByCoordinates(points[3 * i], points[3 * i + 1], points[3 * i + 2], elementTag, elementType, nodeTags, uGmsh, vGmsh, wGmsh, field->domain().maxDim(), true);
            u = uGmsh;
            v = vGmsh;
            w = wGmsh;
          }
          catch(...) {
            MathObject< T_Scalar, field::DegreeOfCompoundForm< field::NextForm< T_Form >::value >::value >::zero(values[i]);
            continue;
          }

          std::vector< scalar::Precision< T_Scalar > > coordinates;
          std::vector< int > typeKeys;
          std::vector< unsigned long long > entityKeys;
          field->getFunctionSpace()->getKeysOnElement(false, typeKeys, entityKeys, coordinates, elementTag);

          std::vector< T_Scalar > dofsValues;
          dofsValues.resize(T_NumFields * typeKeys.size());
          field->getValues(typeKeys, entityKeys, dofsValues, 0, typeKeys.size());

          std::vector< scalar::Precision< T_Scalar > > functions;
          const unsigned int orientation = field->getFunctionSpace()->getOrientationOfElement(elementTag);
          const unsigned int nbrDofsByElements = field->getFunctionSpace()->getBasisFunction_noexcept(true, functions, {u, v, w}, elementType, orientation);

          // Jacobians
          std::vector< double > gmshJacobians, gmshDeterminants, gmshPoints;
          std::vector< scalar::Precision< T_Scalar > > jacobians, determinants, points;
          if(_fieldEvaluator.needJacobians()) {
            gmsh::model::mesh::getJacobian(elementTag, {uGmsh, vGmsh, wGmsh}, gmshJacobians, gmshDeterminants, gmshPoints);
            scalar::move(jacobians, gmshJacobians);
            scalar::move(determinants, gmshDeterminants);
            scalar::move(points, gmshPoints);
            //              Not correct: need to find the entity tag in the field's mesh
            //              if(domain.haveJacobiansModificators(entity)) {
            //                domain.applyJacobiansModificator(points, determinants, jacobians, entity);
            //              }
          }

          // Evaluation
          Eigen::VectorX< T_Scalar > vecDof = Eigen::Map< const Eigen::VectorX< T_Scalar > >(&(dofsValues[0]), T_NumFields * nbrDofsByElements);
          Eigen::MatrixX< scalar::Precision< T_Scalar > > fieldExpression(_fieldEvaluator.size(), T_NumFields * nbrDofsByElements);
          fieldExpression.setZero();
          _fieldEvaluator(fieldExpression, &functions[0], &jacobians[0], &determinants[0], nbrDofsByElements);
          Eigen::MatrixX< T_Scalar > tmp = fieldExpression * vecDof;
          MathObject< T_Scalar, field::DegreeOfCompoundForm< field::NextForm< T_Form >::value >::value >::copy(values[i], tmp);
        }
#pragma omp master
        {
          gmsh::option::setNumber("General.Verbosity", verbose);
        }
      }
      else {
        const unsigned int nbrOfElements = points.size() / gaussPoints.size();
        const unsigned int numThreads = omp::getNumThreads();
        const unsigned int myThreadID = omp::getThreadNum();

#pragma omp single
        _nbrDofsByElements = field->getFunctionSpace()->getBasisFunctions(true, _functions, _fsIndex, gaussPoints, elementType, entity);

        std::vector< scalar::Precision< T_Scalar > > coordinates;
#pragma omp single
        field->getFunctionSpace()->getKeys(false, _typeKeys, _entityKeys, coordinates, elementType, entity);
        const unsigned int begin = (myThreadID * _typeKeys.size()) / numThreads;
        const unsigned int end = ((myThreadID + 1) * _typeKeys.size()) / numThreads;

        const bool needOffset = !field->getFunctionSpace()->isConstantByElements();

#pragma omp single
        _dofsValues.resize(T_NumFields * _typeKeys.size());

        field->getValues(_typeKeys, _entityKeys, _dofsValues, begin, end);
        const unsigned int nbrGaussPoints = gaussPoints.size() / 3;

#pragma omp single
        scalar::copy(_gmshGaussPoints, gaussPoints);

        if(_fieldEvaluator.needJacobians()) {
          const domain::Domain domain = field->domain();
#pragma omp single
          gmsh::model::mesh::preallocateJacobians(elementType, nbrGaussPoints, true, true, false, _gmshJacobians, _gmshDeterminants, _gmshPoints, entity.second);
          gmsh::model::mesh::getJacobians(elementType, _gmshGaussPoints, _gmshJacobians, _gmshDeterminants, _gmshPoints, entity.second, myThreadID, numThreads);
#pragma omp barrier
#pragma omp single
          numa::copy(_jacobians, _gmshJacobians);
#pragma omp single
          numa::copy(_determinants, _gmshDeterminants);

          if(domain.haveJacobiansModificators(entity)) {
            domain.applyJacobiansModificator(points, _determinants, _jacobians, entity);
#pragma omp barrier
          }
        }

        Eigen::MatrixX< scalar::Precision< T_Scalar > > fieldExpression(_fieldEvaluator.size(), T_NumFields * _nbrDofsByElements);
        fieldExpression.setZero();
#pragma omp for
        for(auto i = 0U; i < nbrOfElements; ++i) {
          Eigen::VectorX< T_Scalar > vecDof = Eigen::Map< const Eigen::VectorX< T_Scalar > >(&(_dofsValues[i * T_NumFields * _nbrDofsByElements]), T_NumFields * _nbrDofsByElements);
          for(auto j = 0U; j < nbrGaussPoints; ++j) {
            _fieldEvaluator(fieldExpression, &_functions[((needOffset ? _fsIndex[i] * nbrGaussPoints : 0) + j) * field::NumberOfComponentsOfForm< field::NextForm< T_Form >::value >::value * _nbrDofsByElements], &_jacobians[i * 9 * nbrGaussPoints + j * 9], &_determinants[i * nbrGaussPoints + j], _nbrDofsByElements);
            Eigen::MatrixX< T_Scalar > tmp = fieldExpression * vecDof;
            MathObject< T_Scalar, field::DegreeOfCompoundForm< field::NextForm< T_Form >::value >::value >::copy(values[i * nbrGaussPoints + j], tmp);
          }
        }
      }
      // Restore the current model
#pragma omp single
      {
        if(_currentModel != field->model()) {
          gmsh::model::setCurrent(_currentModel);
        }
      }
    }

    std::string name() const override
    {
      return "derivative of compound field";
    }
  };

  //
  // None
  //

  template< class T_Scalar, field::Form T_Form >
  class None final : public FieldOperation< T_Scalar, field::DegreeOfForm< T_Form >::value >
  {
   private:
    mutable std::vector< int > _typeKeys;
    mutable std::vector< unsigned long long > _entityKeys;
    mutable std::vector< T_Scalar > _dofsValues;
    mutable std::vector< scalar::Precision< T_Scalar > > _functions;
    mutable std::vector< int > _fsIndex;
    mutable unsigned int _nbrDofsByElements;
    mutable std::vector< double > _gmshJacobians, _gmshDeterminants, _gmshPoints, _gmshGaussPoints;
    mutable std::vector< scalar::Precision< T_Scalar >, numa::allocator< scalar::Precision< T_Scalar > > > _jacobians, _determinants;
    mutable term::evaluator::FieldEvaluator< T_Scalar, T_Form > _fieldEvaluator;
    mutable bool _pointEvaluation;
    mutable std::string _currentModel;

   public:
    None() :
      _pointEvaluation(false)
    {
    }

    None(const None &other) :
      FieldOperation< T_Scalar, field::DegreeOfForm< T_Form >::value >(other), _pointEvaluation(other._pointEvaluation)
    {
    }

    void setPointEvaluation(const bool pointEvaluation) const
    {
      _pointEvaluation = pointEvaluation;
    }

    void operator()(const field::FieldInterface< T_Scalar > *const field, OutputVector< T_Scalar, field::DegreeOfForm< T_Form >::value > &values, const std::vector< scalar::Precision< T_Scalar >, numa::allocator< scalar::Precision< T_Scalar > > > &points, const std::vector< scalar::Precision< T_Scalar > > &gaussPoints, const int elementType, const std::pair< int, int > &entity) const
    {
      // If the current model is not the same on which the field is defined, we use the point evaluation algorithm
      // otherwise, we can interpolate the field as usual.
#pragma omp single
      {
        gmsh::model::getCurrent(_currentModel);
        if(_currentModel != field->model()) {
          _pointEvaluation = true;
          gmsh::model::setCurrent(field->model());
        }
      }
      if(_pointEvaluation) {
        double verbose = 0.;
#pragma omp master
        {
          gmsh::option::getNumber("General.Verbosity", verbose);
          gmsh::option::setNumber("General.Verbosity", 0);
        }
#pragma omp single
        gmsh::model::mesh::rebuildElementCache(true);
#pragma omp for
        for(auto i = 0ULL; i < points.size() / 3; ++i) {
          std::size_t elementTag = 0;
          int elementType = 0;
          std::vector< std::size_t > nodeTags;
          double uGmsh = 0., vGmsh = 0., wGmsh = 0.;
          scalar::Precision< T_Scalar > u = -1., v = -1., w = -1.;
          try {
            gmsh::model::mesh::getElementByCoordinates(points[3 * i], points[3 * i + 1], points[3 * i + 2], elementTag, elementType, nodeTags, uGmsh, vGmsh, wGmsh, field->domain().maxDim(), true);
            u = uGmsh;
            v = vGmsh;
            w = wGmsh;
          }
          catch(...) {
            MathObject< T_Scalar, field::DegreeOfForm< T_Form >::value >::zero(values[i]);
            continue;
          }

          std::vector< scalar::Precision< T_Scalar > > coordinates;
          std::vector< int > typeKeys;
          std::vector< unsigned long long > entityKeys;
          field->getFunctionSpace()->getKeysOnElement(false, typeKeys, entityKeys, coordinates, elementTag);

          std::vector< T_Scalar > dofsValues;
          dofsValues.resize(typeKeys.size());
          field->getValues(typeKeys, entityKeys, dofsValues, 0, dofsValues.size());

          std::vector< scalar::Precision< T_Scalar > > functions;
          const unsigned int orientation = field->getFunctionSpace()->getOrientationOfElement(elementTag);
          const unsigned int nbrDofsByElements = field->getFunctionSpace()->getBasisFunction_noexcept(false, functions, {u, v, w}, elementType, orientation);

          // Jacobians
          std::vector< double > gmshJacobians, gmshDeterminants, gmshPoints;
          std::vector< scalar::Precision< T_Scalar > > jacobians, determinants, points;
          if(_fieldEvaluator.needJacobians()) {
            gmsh::model::mesh::getJacobian(elementTag, {uGmsh, vGmsh, wGmsh}, gmshJacobians, gmshDeterminants, gmshPoints);
            scalar::move(jacobians, gmshJacobians);
            scalar::move(determinants, gmshDeterminants);
            scalar::move(points, gmshPoints);
            //              Not correct: need to find the entity tag in the field's mesh
            //              if(domain.haveJacobiansModificators(entity)) {
            //                domain.applyJacobiansModificator(points, determinants, jacobians, entity);
            //              }
          }

          // Evaluation
          Eigen::VectorX< T_Scalar > vecDof = Eigen::Map< const Eigen::VectorX< T_Scalar > >(&(dofsValues[0]), nbrDofsByElements);
          Eigen::MatrixX< scalar::Precision< T_Scalar > > fieldExpression(_fieldEvaluator.size(), nbrDofsByElements);
          _fieldEvaluator(fieldExpression, &functions[0], &jacobians[0], &determinants[0], nbrDofsByElements);
          MathObject< T_Scalar, field::DegreeOfForm< T_Form >::value >::copy(values[i], fieldExpression * vecDof);
        }
#pragma omp master
        {
          gmsh::option::setNumber("General.Verbosity", verbose);
        }
      }
      else {
        const unsigned int nbrOfElements = points.size() / gaussPoints.size();
        const unsigned int numThreads = omp::getNumThreads();
        const unsigned int myThreadID = omp::getThreadNum();

#pragma omp single
        _nbrDofsByElements = field->getFunctionSpace()->getBasisFunctions(false, _functions, _fsIndex, gaussPoints, elementType, entity);

        std::vector< scalar::Precision< T_Scalar > > coordinates;
#pragma omp single
        field->getFunctionSpace()->getKeys(false, _typeKeys, _entityKeys, coordinates, elementType, entity);
        const unsigned int begin = (myThreadID * _typeKeys.size()) / numThreads;
        const unsigned int end = ((myThreadID + 1) * _typeKeys.size()) / numThreads;

        const bool needOffset = !field->getFunctionSpace()->isConstantByElements();

#pragma omp single
        _dofsValues.resize(_typeKeys.size());

        field->getValues(_typeKeys, _entityKeys, _dofsValues, begin, end);
        const unsigned int nbrGaussPoints = gaussPoints.size() / 3;

#pragma omp single
        scalar::copy(_gmshGaussPoints, gaussPoints);

        if(_fieldEvaluator.needJacobians()) {
          const domain::Domain domain = field->domain();
#pragma omp single
          gmsh::model::mesh::preallocateJacobians(elementType, nbrGaussPoints, true, true, false, _gmshJacobians, _gmshDeterminants, _gmshPoints, entity.second);
          gmsh::model::mesh::getJacobians(elementType, _gmshGaussPoints, _gmshJacobians, _gmshDeterminants, _gmshPoints, entity.second, myThreadID, numThreads);
#pragma omp barrier
#pragma omp single
          numa::copy(_jacobians, _gmshJacobians);
#pragma omp single
          numa::copy(_determinants, _gmshDeterminants);

          if(domain.haveJacobiansModificators(entity)) {
            domain.applyJacobiansModificator(points, _determinants, _jacobians, entity);
#pragma omp barrier
          }
        }

        Eigen::MatrixX< scalar::Precision< T_Scalar > > fieldExpression(_fieldEvaluator.size(), _nbrDofsByElements);
#pragma omp for
        for(auto i = 0U; i < nbrOfElements; ++i) {
          Eigen::VectorX< T_Scalar > vecDof = Eigen::Map< const Eigen::VectorX< T_Scalar > >(&(_dofsValues[i * _nbrDofsByElements]), _nbrDofsByElements);
          for(auto j = 0U; j < nbrGaussPoints; ++j) {
            _fieldEvaluator(fieldExpression, &_functions[((needOffset ? _fsIndex[i] * nbrGaussPoints : 0) + j) * field::NumberOfComponentsOfForm< T_Form >::value * _nbrDofsByElements], &_jacobians[i * 9 * nbrGaussPoints + j * 9], &_determinants[i * nbrGaussPoints + j], _nbrDofsByElements);
            MathObject< T_Scalar, field::DegreeOfForm< T_Form >::value >::copy(values[i * nbrGaussPoints + j], fieldExpression * vecDof);
          }
        }
      }
      // Restore the current model
#pragma omp single
      {
        if(_currentModel != field->model()) {
          gmsh::model::setCurrent(_currentModel);
        }
      }
    }

    std::string name() const override
    {
      return "evaluation of field";
    }
  };


  //
  // NoneCompound
  //

  template< class T_Scalar, field::Form T_Form, unsigned int T_NumFields >
  class NoneCompound final : public FieldOperation< T_Scalar, field::DegreeOfCompoundForm< T_Form >::value >
  {
   private:
    mutable std::vector< int > _typeKeys;
    mutable std::vector< unsigned long long > _entityKeys;
    mutable std::vector< T_Scalar > _dofsValues;
    mutable std::vector< scalar::Precision< T_Scalar > > _functions;
    mutable std::vector< int > _fsIndex;
    mutable unsigned int _nbrDofsByElements;
    mutable std::vector< double > _gmshJacobians, _gmshDeterminants, _gmshPoints, _gmshGaussPoints;
    mutable std::vector< scalar::Precision< T_Scalar >, numa::allocator< scalar::Precision< T_Scalar > > > _jacobians, _determinants;
    mutable term::evaluator::CompoundFieldEvaluator< T_Scalar, T_Form, T_NumFields > _fieldEvaluator;
    mutable bool _pointEvaluation;
    mutable std::string _currentModel;

   public:
    NoneCompound() :
      _pointEvaluation(false)
    {
    }

    NoneCompound(const NoneCompound &other) :
      FieldOperation< T_Scalar, field::DegreeOfCompoundForm< T_Form >::value >(other), _pointEvaluation(other._pointEvaluation)
    {
    }

    void setPointEvaluation(const bool pointEvaluation) const
    {
      _pointEvaluation = pointEvaluation;
    }

    void operator()(const field::FieldInterface< T_Scalar > *const field, OutputVector< T_Scalar, field::DegreeOfCompoundForm< T_Form >::value > &values, const std::vector< scalar::Precision< T_Scalar >, numa::allocator< scalar::Precision< T_Scalar > > > &points, const std::vector< scalar::Precision< T_Scalar > > &gaussPoints, const int elementType, const std::pair< int, int > &entity) const
    {
      // If the current model is not the same on which the field is defined, we use the point evaluation algorithm
      // otherwise, we can interpolate the field as usual.
#pragma omp single
      {
        gmsh::model::getCurrent(_currentModel);
        if(_currentModel != field->model()) {
          _pointEvaluation = true;
          gmsh::model::setCurrent(field->model());
        }
      }
      if(_pointEvaluation) {
        double verbose = 0.;
#pragma omp master
        {
          gmsh::option::getNumber("General.Verbosity", verbose);
          gmsh::option::setNumber("General.Verbosity", 0);
        }
#pragma omp single
        gmsh::model::mesh::rebuildElementCache(true);
#pragma omp for
        for(auto i = 0ULL; i < points.size() / 3; ++i) {
          std::size_t elementTag = 0;
          int elementType = 0;
          std::vector< std::size_t > nodeTags;
          double uGmsh = 0., vGmsh = 0., wGmsh = 0.;
          scalar::Precision< T_Scalar > u = -1., v = -1., w = -1.;
          try {
            gmsh::model::mesh::getElementByCoordinates(points[3 * i], points[3 * i + 1], points[3 * i + 2], elementTag, elementType, nodeTags, uGmsh, vGmsh, wGmsh, field->domain().maxDim(), true);
            u = uGmsh;
            v = vGmsh;
            w = wGmsh;
          }
          catch(...) {
            MathObject< T_Scalar, field::DegreeOfCompoundForm< T_Form >::value >::zero(values[i]);
            continue;
          }

          std::vector< scalar::Precision< T_Scalar > > coordinates;
          std::vector< int > typeKeys;
          std::vector< unsigned long long > entityKeys;
          field->getFunctionSpace()->getKeysOnElement(false, typeKeys, entityKeys, coordinates, elementTag);

          std::vector< T_Scalar > dofsValues;
          dofsValues.resize(T_NumFields * typeKeys.size());
          field->getValues(typeKeys, entityKeys, dofsValues, 0, typeKeys.size());

          std::vector< scalar::Precision< T_Scalar > > functions;
          const unsigned int orientation = field->getFunctionSpace()->getOrientationOfElement(elementTag);
          const unsigned int nbrDofsByElements = field->getFunctionSpace()->getBasisFunction_noexcept(false, functions, {u, v, w}, elementType, orientation);

          // Jacobians
          std::vector< double > gmshJacobians, gmshDeterminants, gmshPoints;
          std::vector< scalar::Precision< T_Scalar > > jacobians, determinants, points;
          if(_fieldEvaluator.needJacobians()) {
            gmsh::model::mesh::getJacobian(elementTag, {uGmsh, vGmsh, wGmsh}, gmshJacobians, gmshDeterminants, gmshPoints);
            scalar::move(jacobians, gmshJacobians);
            scalar::move(determinants, gmshDeterminants);
            scalar::move(points, gmshPoints);
            //              Not correct: need to find the entity tag in the field's mesh
            //              if(domain.haveJacobiansModificators(entity)) {
            //                domain.applyJacobiansModificator(points, determinants, jacobians, entity);
            //              }
          }

          // Evaluation
          Eigen::VectorX< T_Scalar > vecDof = Eigen::Map< const Eigen::VectorX< T_Scalar > >(&(dofsValues[0]), T_NumFields * nbrDofsByElements);
          Eigen::MatrixX< scalar::Precision< T_Scalar > > fieldExpression(_fieldEvaluator.size(), T_NumFields * nbrDofsByElements);
          fieldExpression.setZero();
          _fieldEvaluator(fieldExpression, &functions[0], &jacobians[0], &determinants[0], nbrDofsByElements);
          Eigen::MatrixX< T_Scalar > tmp = fieldExpression * vecDof;
          MathObject< T_Scalar, field::DegreeOfCompoundForm< T_Form >::value >::copy(values[i], tmp);
        }
#pragma omp master
        {
          gmsh::option::setNumber("General.Verbosity", verbose);
        }
      }
      else {
        const unsigned int nbrOfElements = points.size() / gaussPoints.size();
        const unsigned int numThreads = omp::getNumThreads();
        const unsigned int myThreadID = omp::getThreadNum();

#pragma omp single
        _nbrDofsByElements = field->getFunctionSpace()->getBasisFunctions(false, _functions, _fsIndex, gaussPoints, elementType, entity);

        std::vector< scalar::Precision< T_Scalar > > coordinates;
#pragma omp single
        field->getFunctionSpace()->getKeys(false, _typeKeys, _entityKeys, coordinates, elementType, entity);
        const unsigned int begin = (myThreadID * _typeKeys.size()) / numThreads;
        const unsigned int end = ((myThreadID + 1) * _typeKeys.size()) / numThreads;

        const bool needOffset = !field->getFunctionSpace()->isConstantByElements();

#pragma omp single
        _dofsValues.resize(T_NumFields * _typeKeys.size());

        field->getValues(_typeKeys, _entityKeys, _dofsValues, begin, end);
        const unsigned int nbrGaussPoints = gaussPoints.size() / 3;

#pragma omp single
        scalar::copy(_gmshGaussPoints, gaussPoints);

        if(_fieldEvaluator.needJacobians()) {
          const domain::Domain domain = field->domain();
#pragma omp single
          gmsh::model::mesh::preallocateJacobians(elementType, nbrGaussPoints, true, true, false, _gmshJacobians, _gmshDeterminants, _gmshPoints, entity.second);
          gmsh::model::mesh::getJacobians(elementType, _gmshGaussPoints, _gmshJacobians, _gmshDeterminants, _gmshPoints, entity.second, myThreadID, numThreads);
#pragma omp barrier
#pragma omp single
          numa::copy(_jacobians, _gmshJacobians);
#pragma omp single
          numa::copy(_determinants, _gmshDeterminants);

          if(domain.haveJacobiansModificators(entity)) {
            domain.applyJacobiansModificator(points, _determinants, _jacobians, entity);
#pragma omp barrier
          }
        }

        Eigen::MatrixX< scalar::Precision< T_Scalar > > fieldExpression(_fieldEvaluator.size(), T_NumFields * _nbrDofsByElements);
        fieldExpression.setZero();
#pragma omp for
        for(auto i = 0U; i < nbrOfElements; ++i) {
          Eigen::VectorX< T_Scalar > vecDof = Eigen::Map< const Eigen::VectorX< T_Scalar > >(&(_dofsValues[i * T_NumFields * _nbrDofsByElements]), T_NumFields * _nbrDofsByElements);
          for(auto j = 0U; j < nbrGaussPoints; ++j) {
            _fieldEvaluator(fieldExpression, &_functions[((needOffset ? _fsIndex[i] * nbrGaussPoints : 0) + j) * field::NumberOfComponentsOfForm< T_Form >::value * _nbrDofsByElements], &_jacobians[i * 9 * nbrGaussPoints + j * 9], &_determinants[i * nbrGaussPoints + j], _nbrDofsByElements);
            Eigen::MatrixX< T_Scalar > tmp = fieldExpression * vecDof;
            MathObject< T_Scalar, field::DegreeOfCompoundForm< T_Form >::value >::copy(values[i * nbrGaussPoints + j], tmp);
          }
        }
      }
      // Restore the current model
#pragma omp single
      {
        if(_currentModel != field->model()) {
          gmsh::model::setCurrent(_currentModel);
        }
      }
    }

    std::string name() const override
    {
      return "evaluation of compound field";
    }
  };


  // Available field scalar type operations:
  //  Eigenfunction
  //  EigenfunctionCompound

  //
  // Eigenfunction
  //

  template< class T_Scalar, field::Form T_Form >
  class Eigenfunction final : public FieldScalarTypeOperation< scalar::ComplexPrecision< T_Scalar >, T_Scalar, field::DegreeOfForm< T_Form >::value >
  {
   private:
    mutable std::vector< int > _typeKeys;
    mutable std::vector< unsigned long long > _entityKeys;
    mutable std::vector< scalar::ComplexPrecision< T_Scalar > > _dofsValues;
    mutable std::vector< scalar::Precision< scalar::ComplexPrecision< T_Scalar > > > _functions;
    mutable std::vector< int > _fsIndex;
    mutable unsigned int _nbrDofsByElements;
    mutable std::vector< double > _gmshJacobians, _gmshDeterminants, _gmshPoints, _gmshGaussPoints;
    mutable std::vector< scalar::Precision< T_Scalar >, numa::allocator< scalar::Precision< T_Scalar > > > _jacobians, _determinants;
    mutable term::evaluator::FieldEvaluator< scalar::ComplexPrecision< T_Scalar >, T_Form > _fieldEvaluator;
    mutable bool _pointEvaluation;
    mutable unsigned int _eigenTag;
    mutable std::string _currentModel;

   public:
    Eigenfunction(const unsigned int eigenTag) :
      _pointEvaluation(false), _eigenTag(eigenTag)
    {
    }

    Eigenfunction(const Eigenfunction &other) :
      FieldScalarTypeOperation< scalar::ComplexPrecision< T_Scalar >, T_Scalar, field::DegreeOfForm< T_Form >::value >(other), _pointEvaluation(other._pointEvaluation), _eigenTag(other._eigenTag)
    {
    }

    void setPointEvaluation(const bool pointEvaluation) const
    {
      _pointEvaluation = pointEvaluation;
    }

    void operator()(const field::FieldInterface< T_Scalar > *const field, OutputVector< scalar::ComplexPrecision< T_Scalar >, field::DegreeOfForm< T_Form >::value > &values, const std::vector< scalar::Precision< T_Scalar >, numa::allocator< scalar::Precision< T_Scalar > > > &points, const std::vector< scalar::Precision< T_Scalar > > &gaussPoints, const int elementType, const std::pair< int, int > &entity) const
    {
      // If the current model is not the same on which the field is defined, we use the point evaluation algorithm
      // otherwise, we can interpolate the field as usual.
#pragma omp single
      {
        gmsh::model::getCurrent(_currentModel);
        if(_currentModel != field->model()) {
          _pointEvaluation = true;
          gmsh::model::setCurrent(field->model());
        }
      }
      if(_pointEvaluation) {
        double verbose = 0.;
#pragma omp master
        {
          gmsh::option::getNumber("General.Verbosity", verbose);
          gmsh::option::setNumber("General.Verbosity", 0);
        }
#pragma omp single
        gmsh::model::mesh::rebuildElementCache(true);
#pragma omp for
        for(auto i = 0ULL; i < points.size() / 3; ++i) {
          std::size_t elementTag = 0;
          int elementType = 0;
          std::vector< std::size_t > nodeTags;
          double uGmsh = 0., vGmsh = 0., wGmsh = 0.;
          scalar::Precision< T_Scalar > u = 0., v = 0., w = 0.;
          try {
            gmsh::model::mesh::getElementByCoordinates(points[3 * i], points[3 * i + 1], points[3 * i + 2], elementTag, elementType, nodeTags, uGmsh, vGmsh, wGmsh, field->domain().maxDim(), true);
            u = uGmsh;
            v = vGmsh;
            w = wGmsh;
          }
          catch(...) {
            MathObject< scalar::ComplexPrecision< T_Scalar >, field::DegreeOfForm< T_Form >::value >::zero(values[i]);
            continue;
          }

          std::vector< scalar::Precision< T_Scalar > > coordinates;
          std::vector< int > typeKeys;
          std::vector< unsigned long long > entityKeys;
          field->getFunctionSpace()->getKeysOnElement(false, typeKeys, entityKeys, coordinates, elementTag);

          std::vector< scalar::ComplexPrecision< T_Scalar > > dofsValues;
          dofsValues.resize(typeKeys.size());
          field::EigenpairModule< T_Scalar > *module = static_cast< field::EigenpairModule< T_Scalar > * >(field->getModule("Eigenpair"));
          if(module == nullptr) {
            throw common::Exception("There is no eigenpair module in field '" + field->name());
          }
          module->getValues(_eigenTag, typeKeys, entityKeys, dofsValues, 0, dofsValues.size());

          std::vector< scalar::Precision< T_Scalar > > functions;
          const unsigned int orientation = field->getFunctionSpace()->getOrientationOfElement(elementTag);
          const unsigned int nbrDofsByElements = field->getFunctionSpace()->getBasisFunction_noexcept(false, functions, {u, v, w}, elementType, orientation);

          // Jacobians
          std::vector< double > gmshJacobians, gmshDeterminants, gmshPoints;
          std::vector< scalar::Precision< T_Scalar > > jacobians, determinants, points;
          if(_fieldEvaluator.needJacobians()) {
            gmsh::model::mesh::getJacobian(elementTag, {uGmsh, vGmsh, wGmsh}, gmshJacobians, gmshDeterminants, gmshPoints);
            scalar::move(jacobians, gmshJacobians);
            scalar::move(determinants, gmshDeterminants);
            scalar::move(points, gmshPoints);
            //              Not correct: need to find the entity tag in the field's mesh
            //              if(domain.haveJacobiansModificators(entity)) {
            //                domain.applyJacobiansModificator(points, determinants, jacobians, entity);
            //              }
          }

          // Evaluation
          Eigen::VectorX< scalar::ComplexPrecision< T_Scalar > > vecDof = Eigen::Map< const Eigen::VectorX< scalar::ComplexPrecision< T_Scalar > > >(&(dofsValues[0]), nbrDofsByElements);
          Eigen::MatrixX< scalar::Precision< T_Scalar > > fieldExpression(_fieldEvaluator.size(), nbrDofsByElements);
          _fieldEvaluator(fieldExpression, &functions[0], &jacobians[0], &determinants[0], nbrDofsByElements);
          MathObject< scalar::ComplexPrecision< T_Scalar >, field::DegreeOfForm< T_Form >::value >::copy(values[i], fieldExpression * vecDof);
        }
#pragma omp master
        {
          gmsh::option::setNumber("General.Verbosity", verbose);
        }
      }
      else {
        const unsigned int nbrOfElements = points.size() / gaussPoints.size();
        const unsigned int numThreads = omp::getNumThreads();
        const unsigned int myThreadID = omp::getThreadNum();

#pragma omp single
        _nbrDofsByElements = field->getFunctionSpace()->getBasisFunctions(false, _functions, _fsIndex, gaussPoints, elementType, entity);

        std::vector< scalar::Precision< T_Scalar > > coordinates;
#pragma omp single
        field->getFunctionSpace()->getKeys(false, _typeKeys, _entityKeys, coordinates, elementType, entity);
        const unsigned int begin = (myThreadID * _typeKeys.size()) / numThreads;
        const unsigned int end = ((myThreadID + 1) * _typeKeys.size()) / numThreads;

        const bool needOffset = !field->getFunctionSpace()->isConstantByElements();

#pragma omp single
        _dofsValues.resize(_typeKeys.size());

        field::EigenpairModule< T_Scalar > *module = static_cast< field::EigenpairModule< T_Scalar > * >(field->getModule("Eigenpair"));
        if(module == nullptr) {
          throw common::Exception("There is no eigenpair module in field '" + field->name());
        }
        module->getValues(_eigenTag, _typeKeys, _entityKeys, _dofsValues, begin, end);
        const unsigned int nbrGaussPoints = gaussPoints.size() / 3;

#pragma omp single
        scalar::copy(_gmshGaussPoints, gaussPoints);

        if(_fieldEvaluator.needJacobians()) {
          const domain::Domain domain = field->domain();
#pragma omp single
          gmsh::model::mesh::preallocateJacobians(elementType, nbrGaussPoints, true, true, false, _gmshJacobians, _gmshDeterminants, _gmshPoints, entity.second);
          gmsh::model::mesh::getJacobians(elementType, _gmshGaussPoints, _gmshJacobians, _gmshDeterminants, _gmshPoints, entity.second, myThreadID, numThreads);
#pragma omp barrier
#pragma omp single
          numa::copy(_jacobians, _gmshJacobians);
#pragma omp single
          numa::copy(_determinants, _gmshDeterminants);

          if(domain.haveJacobiansModificators(entity)) {
            domain.applyJacobiansModificator(points, _determinants, _jacobians, entity);
#pragma omp barrier
          }
        }

        Eigen::MatrixX< scalar::Precision< T_Scalar > > fieldExpression(_fieldEvaluator.size(), _nbrDofsByElements);
#pragma omp for
        for(auto i = 0U; i < nbrOfElements; ++i) {
          Eigen::VectorX< scalar::ComplexPrecision< T_Scalar > > vecDof = Eigen::Map< const Eigen::VectorX< scalar::ComplexPrecision< T_Scalar > > >(&(_dofsValues[i * _nbrDofsByElements]), _nbrDofsByElements);
          for(auto j = 0U; j < nbrGaussPoints; ++j) {
            _fieldEvaluator(fieldExpression, &_functions[((needOffset ? _fsIndex[i] * nbrGaussPoints : 0) + j) * field::NumberOfComponentsOfForm< T_Form >::value * _nbrDofsByElements], &_jacobians[i * 9 * nbrGaussPoints + j * 9], &_determinants[i * nbrGaussPoints + j], _nbrDofsByElements);
            MathObject< scalar::ComplexPrecision< T_Scalar >, field::DegreeOfForm< T_Form >::value >::copy(values[i * nbrGaussPoints + j], fieldExpression * vecDof);
          }
        }
      }
      // Restore the current model
#pragma omp single
      {
        if(_currentModel != field->model()) {
          gmsh::model::setCurrent(_currentModel);
        }
      }
    }

    std::string name() const override
    {
      return "eigenfunction of field";
    }
  };

  //
  // EigenfunctionCompound
  //

  template< class T_Scalar, field::Form T_Form, unsigned int T_NumFields >
  class EigenfunctionCompound final : public FieldScalarTypeOperation< scalar::ComplexPrecision< T_Scalar >, T_Scalar, field::DegreeOfCompoundForm< T_Form >::value >
  {
   private:
    mutable std::vector< int > _typeKeys;
    mutable std::vector< unsigned long long > _entityKeys;
    mutable std::vector< scalar::ComplexPrecision< T_Scalar > > _dofsValues;
    mutable std::vector< scalar::Precision< scalar::ComplexPrecision< T_Scalar > > > _functions;
    mutable std::vector< int > _fsIndex;
    mutable unsigned int _nbrDofsByElements;
    mutable std::vector< double > _gmshJacobians, _gmshDeterminants, _gmshPoints, _gmshGaussPoints;
    mutable std::vector< scalar::Precision< T_Scalar >, numa::allocator< scalar::Precision< T_Scalar > > > _jacobians, _determinants;
    mutable term::evaluator::CompoundFieldEvaluator< scalar::ComplexPrecision< T_Scalar >, T_Form, T_NumFields > _fieldEvaluator;
    mutable bool _pointEvaluation;
    mutable unsigned int _eigenTag;
    mutable std::string _currentModel;

   public:
    EigenfunctionCompound(const unsigned int eigenTag) :
      _pointEvaluation(false), _eigenTag(eigenTag)
    {
    }

    EigenfunctionCompound(const EigenfunctionCompound &other) :
      FieldScalarTypeOperation< scalar::ComplexPrecision< T_Scalar >, T_Scalar, field::DegreeOfCompoundForm< T_Form >::value >(other), _pointEvaluation(other._pointEvaluation), _eigenTag(other._eigenTag)
    {
    }

    void setPointEvaluation(const bool pointEvaluation) const
    {
      _pointEvaluation = pointEvaluation;
    }

    void operator()(const field::FieldInterface< T_Scalar > *const field, OutputVector< scalar::ComplexPrecision< T_Scalar >, field::DegreeOfCompoundForm< T_Form >::value > &values, const std::vector< scalar::Precision< T_Scalar >, numa::allocator< scalar::Precision< T_Scalar > > > &points, const std::vector< scalar::Precision< T_Scalar > > &gaussPoints, const int elementType, const std::pair< int, int > &entity) const
    {
      // If the current model is not the same on which the field is defined, we use the point evaluation algorithm
      // otherwise, we can interpolate the field as usual.
#pragma omp single
      {
        gmsh::model::getCurrent(_currentModel);
        if(_currentModel != field->model()) {
          _pointEvaluation = true;
          gmsh::model::setCurrent(field->model());
        }
      }
      if(_pointEvaluation) {
        double verbose = 0.;
#pragma omp master
        {
          gmsh::option::getNumber("General.Verbosity", verbose);
          gmsh::option::setNumber("General.Verbosity", 0);
        }
#pragma omp single
        gmsh::model::mesh::rebuildElementCache(true);
#pragma omp for
        for(auto i = 0ULL; i < points.size() / 3; ++i) {
          std::size_t elementTag = 0;
          int elementType = 0;
          std::vector< std::size_t > nodeTags;
          double uGmsh = 0., vGmsh = 0., wGmsh = 0.;
          scalar::Precision< T_Scalar > u = 0., v = 0., w = 0.;
          try {
            gmsh::model::mesh::getElementByCoordinates(points[3 * i], points[3 * i + 1], points[3 * i + 2], elementTag, elementType, nodeTags, uGmsh, vGmsh, wGmsh, field->domain().maxDim(), true);
            u = uGmsh;
            v = vGmsh;
            w = wGmsh;
          }
          catch(...) {
            MathObject< scalar::ComplexPrecision< T_Scalar >, field::DegreeOfCompoundForm< T_Form >::value >::zero(values[i]);
            continue;
          }

          std::vector< scalar::Precision< T_Scalar > > coordinates;
          std::vector< int > typeKeys;
          std::vector< unsigned long long > entityKeys;
          field->getFunctionSpace()->getKeysOnElement(false, typeKeys, entityKeys, coordinates, elementTag);

          std::vector< scalar::ComplexPrecision< T_Scalar > > dofsValues;
          dofsValues.resize(T_NumFields * typeKeys.size());
          field::EigenpairModule< T_Scalar > *module = static_cast< field::EigenpairModule< T_Scalar > * >(field->getModule("Eigenpair"));
          if(module == nullptr) {
            throw common::Exception("There is no eigenpair module in field '" + field->name());
          }
          module->getValues(_eigenTag, typeKeys, entityKeys, dofsValues, 0, typeKeys.size());

          std::vector< scalar::Precision< T_Scalar > > functions;
          const unsigned int orientation = field->getFunctionSpace()->getOrientationOfElement(elementTag);
          const unsigned int nbrDofsByElements = field->getFunctionSpace()->getBasisFunction_noexcept(false, functions, {u, v, w}, elementType, orientation);

          // Jacobians
          std::vector< double > gmshJacobians, gmshDeterminants, gmshPoints;
          std::vector< scalar::Precision< T_Scalar > > jacobians, determinants, points;
          if(_fieldEvaluator.needJacobians()) {
            gmsh::model::mesh::getJacobian(elementTag, {uGmsh, vGmsh, wGmsh}, gmshJacobians, gmshDeterminants, gmshPoints);
            scalar::move(jacobians, gmshJacobians);
            scalar::move(determinants, gmshDeterminants);
            scalar::move(points, gmshPoints);
            //              Not correct: need to find the entity tag in the field's mesh
            //              if(domain.haveJacobiansModificators(entity)) {
            //                domain.applyJacobiansModificator(points, determinants, jacobians, entity);
            //              }
          }

          // Evaluation
          Eigen::VectorX< scalar::ComplexPrecision< T_Scalar > > vecDof = Eigen::Map< const Eigen::VectorX< scalar::ComplexPrecision< T_Scalar > > >(&(dofsValues[0]), T_NumFields * nbrDofsByElements);
          Eigen::MatrixX< scalar::Precision< T_Scalar > > fieldExpression(_fieldEvaluator.size(), T_NumFields * nbrDofsByElements);
          fieldExpression.setZero();
          _fieldEvaluator(fieldExpression, &functions[0], &jacobians[0], &determinants[0], nbrDofsByElements);
          Eigen::MatrixX< scalar::ComplexPrecision< T_Scalar > > tmp = fieldExpression * vecDof;
          MathObject< scalar::ComplexPrecision< T_Scalar >, field::DegreeOfCompoundForm< T_Form >::value >::copy(values[i], tmp);
        }
#pragma omp master
        {
          gmsh::option::setNumber("General.Verbosity", verbose);
        }
      }
      else {
        const unsigned int nbrOfElements = points.size() / gaussPoints.size();
        const unsigned int numThreads = omp::getNumThreads();
        const unsigned int myThreadID = omp::getThreadNum();

#pragma omp single
        _nbrDofsByElements = field->getFunctionSpace()->getBasisFunctions(false, _functions, _fsIndex, gaussPoints, elementType, entity);

        std::vector< scalar::Precision< T_Scalar > > coordinates;
#pragma omp single
        field->getFunctionSpace()->getKeys(false, _typeKeys, _entityKeys, coordinates, elementType, entity);
        const unsigned int begin = (myThreadID * _typeKeys.size()) / numThreads;
        const unsigned int end = ((myThreadID + 1) * _typeKeys.size()) / numThreads;

        const bool needOffset = !field->getFunctionSpace()->isConstantByElements();

#pragma omp single
        _dofsValues.resize(T_NumFields * _typeKeys.size());

        field::EigenpairModule< T_Scalar > *module = static_cast< field::EigenpairModule< T_Scalar > * >(field->getModule("Eigenpair"));
        if(module == nullptr) {
          throw common::Exception("There is no eigenpair module in field '" + field->name());
        }
        module->getValues(_eigenTag, _typeKeys, _entityKeys, _dofsValues, begin, end);
        const unsigned int nbrGaussPoints = gaussPoints.size() / 3;

#pragma omp single
        scalar::copy(_gmshGaussPoints, gaussPoints);

        if(_fieldEvaluator.needJacobians()) {
          const domain::Domain domain = field->domain();
#pragma omp single
          gmsh::model::mesh::preallocateJacobians(elementType, nbrGaussPoints, true, true, false, _gmshJacobians, _gmshDeterminants, _gmshPoints, entity.second);
          gmsh::model::mesh::getJacobians(elementType, _gmshGaussPoints, _gmshJacobians, _gmshDeterminants, _gmshPoints, entity.second, myThreadID, numThreads);
#pragma omp barrier
#pragma omp single
          numa::copy(_jacobians, _gmshJacobians);
#pragma omp single
          numa::copy(_determinants, _gmshDeterminants);

          if(domain.haveJacobiansModificators(entity)) {
            domain.applyJacobiansModificator(points, _determinants, _jacobians, entity);
#pragma omp barrier
          }
        }

        Eigen::MatrixX< scalar::Precision< T_Scalar > > fieldExpression(_fieldEvaluator.size(), T_NumFields * _nbrDofsByElements);
        fieldExpression.setZero();
#pragma omp for
        for(auto i = 0U; i < nbrOfElements; ++i) {
          Eigen::VectorX< scalar::ComplexPrecision< T_Scalar > > vecDof = Eigen::Map< const Eigen::VectorX< scalar::ComplexPrecision< T_Scalar > > >(&(_dofsValues[i * T_NumFields * _nbrDofsByElements]), T_NumFields * _nbrDofsByElements);
          for(auto j = 0U; j < nbrGaussPoints; ++j) {
            _fieldEvaluator(fieldExpression, &_functions[((needOffset ? _fsIndex[i] * nbrGaussPoints : 0) + j) * field::NumberOfComponentsOfForm< T_Form >::value * _nbrDofsByElements], &_jacobians[i * 9 * nbrGaussPoints + j * 9], &_determinants[i * nbrGaussPoints + j], _nbrDofsByElements);
            Eigen::MatrixX< scalar::ComplexPrecision< T_Scalar > > tmp = fieldExpression * vecDof;
            MathObject< scalar::ComplexPrecision< T_Scalar >, field::DegreeOfCompoundForm< T_Form >::value >::copy(values[i * nbrGaussPoints + j], tmp);
          }
        }
      }
      // Restore the current model
#pragma omp single
      {
        if(_currentModel != field->model()) {
          gmsh::model::setCurrent(_currentModel);
        }
      }
    }

    std::string name() const override
    {
      return "eigenfunction of field";
    }
  };


} // namespace gmshfem::function


#endif // H_GMSHFEM_FIELDOPERATIONS
