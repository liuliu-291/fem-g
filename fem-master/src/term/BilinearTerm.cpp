// GmshFEM - Copyright (C) 2019-2022, A. Royer, E. Béchet, C. Geuzaine, Université de Liège
//
// See the LICENSE.txt file for license information. Please report all
// issues on https://gitlab.onelab.info/gmsh/fem/issues

#include "BilinearTerm.h"

#include "EquationEvaluator.h"
#include "FieldEvaluator.h"
#include "FieldInterface.h"
#include "FunctionSpaceBucket.h"
#include "instantiate.h"

//
//  class BilinearTerm : public BilinearTermInterface
//

namespace gmshfem::term
{


  template< class T_Scalar, field::Form T_FormLhs, field::Form T_FormRhs >
  BilinearTerm< T_Scalar, T_FormLhs, T_FormRhs >::BilinearTerm(const domain::GeometricObject &domain, const std::string &integrationType, const equation::EquationInterface< T_Scalar, Degree::Degree0, T_FormLhs > &equationLhs, const equation::EquationInterface< T_Scalar, Degree::Degree0, T_FormRhs > &equationRhs, const ProductType productType) :
    BilinearTermInterface< T_Scalar >(domain, integrationType, equationLhs.getField(), equationRhs.getField(), productType), _equationLhs0(equationLhs.copy()), _equationLhs1(nullptr), _equationLhs2(nullptr), _fieldEvaluatorLhs(equationLhs.getFieldEvaluator()), _equationEvaluatorLhs(equationLhs.getEquationEvaluator()), _isDerivativeLhs(equationLhs.isDerivative()), _equationRhs0(equationRhs.copy()), _equationRhs1(nullptr), _equationRhs2(nullptr), _fieldEvaluatorRhs(equationRhs.getFieldEvaluator()), _equationEvaluatorRhs(equationRhs.getEquationEvaluator()), _isDerivativeRhs(equationRhs.isDerivative()), _degree(Degree::Degree0), _fsLhs(), _fsIndexLhs(nullptr), _needOffsetLhs(false), _fsRhs(), _fsIndexRhs(nullptr), _needOffsetRhs(false), _gaussWeights(), _fieldExpressionLhs(nullptr), _fieldExpressionRhs(nullptr)
  {
  }

  template< class T_Scalar, field::Form T_FormLhs, field::Form T_FormRhs >
  BilinearTerm< T_Scalar, T_FormLhs, T_FormRhs >::BilinearTerm(const domain::GeometricObject &domain, const std::string &integrationType, const equation::EquationInterface< T_Scalar, Degree::Degree1, T_FormLhs > &equationLhs, const equation::EquationInterface< T_Scalar, Degree::Degree1, T_FormRhs > &equationRhs, const ProductType productType) :
    BilinearTermInterface< T_Scalar >(domain, integrationType, equationLhs.getField(), equationRhs.getField(), productType), _equationLhs0(nullptr), _equationLhs1(equationLhs.copy()), _equationLhs2(nullptr), _fieldEvaluatorLhs(equationLhs.getFieldEvaluator()), _equationEvaluatorLhs(equationLhs.getEquationEvaluator()), _isDerivativeLhs(equationLhs.isDerivative()), _equationRhs0(nullptr), _equationRhs1(equationRhs.copy()), _equationRhs2(nullptr), _fieldEvaluatorRhs(equationRhs.getFieldEvaluator()), _equationEvaluatorRhs(equationRhs.getEquationEvaluator()), _isDerivativeRhs(equationRhs.isDerivative()), _degree(Degree::Degree1), _fsLhs(), _fsIndexLhs(nullptr), _needOffsetLhs(false), _fsRhs(), _fsIndexRhs(nullptr), _needOffsetRhs(false), _gaussWeights(), _fieldExpressionLhs(nullptr), _fieldExpressionRhs(nullptr)
  {
  }

  template< class T_Scalar, field::Form T_FormLhs, field::Form T_FormRhs >
  BilinearTerm< T_Scalar, T_FormLhs, T_FormRhs >::BilinearTerm(const domain::GeometricObject &domain, const std::string &integrationType, const equation::EquationInterface< T_Scalar, Degree::Degree2, T_FormLhs > &equationLhs, const equation::EquationInterface< T_Scalar, Degree::Degree2, T_FormRhs > &equationRhs, const ProductType productType) :
    BilinearTermInterface< T_Scalar >(domain, integrationType, equationLhs.getField(), equationRhs.getField(), productType), _equationLhs0(nullptr), _equationLhs1(nullptr), _equationLhs2(equationLhs.copy()), _fieldEvaluatorLhs(equationLhs.getFieldEvaluator()), _equationEvaluatorLhs(equationLhs.getEquationEvaluator()), _isDerivativeLhs(equationLhs.isDerivative()), _equationRhs0(nullptr), _equationRhs1(nullptr), _equationRhs2(equationRhs.copy()), _fieldEvaluatorRhs(equationRhs.getFieldEvaluator()), _equationEvaluatorRhs(equationRhs.getEquationEvaluator()), _isDerivativeRhs(equationRhs.isDerivative()), _degree(Degree::Degree2), _fsLhs(), _fsIndexLhs(nullptr), _needOffsetLhs(false), _fsRhs(), _fsIndexRhs(nullptr), _needOffsetRhs(false), _gaussWeights(), _fieldExpressionLhs(nullptr), _fieldExpressionRhs(nullptr)
  {
  }

  template< class T_Scalar, field::Form T_FormLhs, field::Form T_FormRhs >
  BilinearTerm< T_Scalar, T_FormLhs, T_FormRhs >::BilinearTerm(const domain::GeometricObject &domain, const std::string &integrationType, const equation::EquationInterface< T_Scalar, Degree::Degree0, T_FormLhs > &equationLhs, equation::EquationInterface< T_Scalar, Degree::Degree0, T_FormRhs > *equationRhs, const ProductType productType) :
    BilinearTermInterface< T_Scalar >(domain, integrationType, equationLhs.getField(), equationRhs->getField(), productType), _equationLhs0(equationLhs.copy()), _equationLhs1(nullptr), _equationLhs2(nullptr), _fieldEvaluatorLhs(equationLhs.getFieldEvaluator()), _equationEvaluatorLhs(equationLhs.getEquationEvaluator()), _isDerivativeLhs(equationLhs.isDerivative()), _equationRhs0(equationRhs), _equationRhs1(nullptr), _equationRhs2(nullptr), _fieldEvaluatorRhs(equationRhs->getFieldEvaluator()), _equationEvaluatorRhs(equationRhs->getEquationEvaluator()), _isDerivativeRhs(equationRhs->isDerivative()), _degree(Degree::Degree0), _fsLhs(), _fsIndexLhs(nullptr), _needOffsetLhs(false), _fsRhs(), _fsIndexRhs(nullptr), _needOffsetRhs(false), _gaussWeights(), _fieldExpressionLhs(nullptr), _fieldExpressionRhs(nullptr)
  {
  }

  template< class T_Scalar, field::Form T_FormLhs, field::Form T_FormRhs >
  BilinearTerm< T_Scalar, T_FormLhs, T_FormRhs >::BilinearTerm(const domain::GeometricObject &domain, const std::string &integrationType, const equation::EquationInterface< T_Scalar, Degree::Degree1, T_FormLhs > &equationLhs, equation::EquationInterface< T_Scalar, Degree::Degree1, T_FormRhs > *equationRhs, const ProductType productType) :
    BilinearTermInterface< T_Scalar >(domain, integrationType, equationLhs.getField(), equationRhs->getField(), productType), _equationLhs0(nullptr), _equationLhs1(equationLhs.copy()), _equationLhs2(nullptr), _fieldEvaluatorLhs(equationLhs.getFieldEvaluator()), _equationEvaluatorLhs(equationLhs.getEquationEvaluator()), _isDerivativeLhs(equationLhs.isDerivative()), _equationRhs0(nullptr), _equationRhs1(equationRhs), _equationRhs2(nullptr), _fieldEvaluatorRhs(equationRhs->getFieldEvaluator()), _equationEvaluatorRhs(equationRhs->getEquationEvaluator()), _isDerivativeRhs(equationRhs->isDerivative()), _degree(Degree::Degree1), _fsLhs(), _fsIndexLhs(nullptr), _needOffsetLhs(false), _fsRhs(), _fsIndexRhs(nullptr), _needOffsetRhs(false), _gaussWeights(), _fieldExpressionLhs(nullptr), _fieldExpressionRhs(nullptr)
  {
  }

  template< class T_Scalar, field::Form T_FormLhs, field::Form T_FormRhs >
  BilinearTerm< T_Scalar, T_FormLhs, T_FormRhs >::BilinearTerm(const domain::GeometricObject &domain, const std::string &integrationType, const equation::EquationInterface< T_Scalar, Degree::Degree2, T_FormLhs > &equationLhs, equation::EquationInterface< T_Scalar, Degree::Degree2, T_FormRhs > *equationRhs, const ProductType productType) :
    BilinearTermInterface< T_Scalar >(domain, integrationType, equationLhs.getField(), equationRhs->getField(), productType), _equationLhs0(nullptr), _equationLhs1(nullptr), _equationLhs2(equationLhs.copy()), _fieldEvaluatorLhs(equationLhs.getFieldEvaluator()), _equationEvaluatorLhs(equationLhs.getEquationEvaluator()), _isDerivativeLhs(equationLhs.isDerivative()), _equationRhs0(nullptr), _equationRhs1(nullptr), _equationRhs2(equationRhs), _fieldEvaluatorRhs(equationRhs->getFieldEvaluator()), _equationEvaluatorRhs(equationRhs->getEquationEvaluator()), _isDerivativeRhs(equationRhs->isDerivative()), _degree(Degree::Degree1), _fsLhs(), _fsIndexLhs(nullptr), _needOffsetLhs(false), _fsRhs(), _fsIndexRhs(nullptr), _needOffsetRhs(false), _gaussWeights(), _fieldExpressionLhs(nullptr), _fieldExpressionRhs(nullptr)
  {
  }

  template< class T_Scalar, field::Form T_FormLhs, field::Form T_FormRhs >
  BilinearTerm< T_Scalar, T_FormLhs, T_FormRhs >::~BilinearTerm()
  {
    _freePrecomputedDate();
    this->_numberOfGaussPoints = 0;
    delete _fieldEvaluatorLhs;
    delete _equationEvaluatorLhs;
    delete _fieldEvaluatorRhs;
    delete _equationEvaluatorRhs;
    if(_degree == Degree::Degree0) {
      delete _equationLhs0;
      delete _equationRhs0;
    }
    else if(_degree == Degree::Degree1) {
      delete _equationLhs1;
      delete _equationRhs1;
    }
    else if(_degree == Degree::Degree2) {
      delete _equationLhs2;
      delete _equationRhs2;
    }
  }

  template< class T_Scalar, field::Form T_FormLhs, field::Form T_FormRhs >
  bool BilinearTerm< T_Scalar, T_FormLhs, T_FormRhs >::needJacobians() const
  {
    return _fieldEvaluatorLhs->needJacobians() || _fieldEvaluatorRhs->needJacobians();
  }

  template< class T_Scalar, field::Form T_FormLhs, field::Form T_FormRhs >
  bool BilinearTerm< T_Scalar, T_FormLhs, T_FormRhs >::needGaussCoordinates(const std::pair< int, int > &entity) const
  {
    if(this->_domain.needGaussCoordinates()) {
      return true;
    }

    if(!_equationEvaluatorLhs->isConstant(entity)) {
      return true;
    }

    if(!_equationEvaluatorRhs->isConstant(entity)) {
      return true;
    }

    return false;
  }

  template< class T_Scalar, field::Form T_FormLhs, field::Form T_FormRhs >
  common::Memory BilinearTerm< T_Scalar, T_FormLhs, T_FormRhs >::memory() const
  {
    return _equationEvaluatorLhs->memory() + _equationEvaluatorRhs->memory() + common::Memory(_gaussWeights.size() * sizeof(scalar::Precision< T_Scalar >));
  }

  template< class T_Scalar, field::Form T_FormLhs, field::Form T_FormRhs >
  equation::UnknownFieldType BilinearTerm< T_Scalar, T_FormLhs, T_FormRhs >::unknownFieldType() const
  {
    if(_degree == Degree::Degree0) {
      return _equationLhs0->getUnknownField()->type();
    }
    else if(_degree == Degree::Degree1) {
      return _equationLhs1->getUnknownField()->type();
    }
    else {
      return _equationLhs2->getUnknownField()->type();
    }
  }

  template< class T_Scalar, field::Form T_FormLhs, field::Form T_FormRhs >
  void BilinearTerm< T_Scalar, T_FormLhs, T_FormRhs >::_freePrecomputedDate()
  {
    _fsRhs.clear();
    _fsIndexRhs = nullptr;
    _fsLhs.clear();
    _fsIndexLhs = nullptr;
    _gaussWeights.clear();
    _equationEvaluatorLhs->clear();
    _equationEvaluatorRhs->clear();
    if(_fieldExpressionLhs == _fieldExpressionRhs) {
      if(_fieldExpressionLhs != nullptr) {
        delete _fieldExpressionLhs;
        _fieldExpressionLhs = nullptr;
      }
    }
    else {
      if(_fieldExpressionLhs != nullptr) {
        delete _fieldExpressionLhs;
        _fieldExpressionLhs = nullptr;
      }
      if(_fieldExpressionRhs != nullptr) {
        delete _fieldExpressionRhs;
        _fieldExpressionRhs = nullptr;
      }
    }
    _expressionLhs.clear();
    _expressionRhs.clear();
  }

  template< class T_Scalar, field::Form T_FormLhs, field::Form T_FormRhs >
  void BilinearTerm< T_Scalar, T_FormLhs, T_FormRhs >::assemblyInitialization(problem::FunctionSpaceBucket< scalar::Precision< T_Scalar > > &functionSpaces, const std::vector< scalar::Precision< T_Scalar >, numa::allocator< scalar::Precision< T_Scalar > > > &points, const int elementType, const std::pair< int, int > &entity)
  {
    _freePrecomputedDate();

    field::FunctionSpaceInterface< scalar::Precision< T_Scalar > > *fs[2] = {this->_field[0]->getFunctionSpace(), this->_field[1]->getFunctionSpace()};

    // **************
    // Function space
    // **************

    // load Gauss information
    std::vector< scalar::Precision< T_Scalar > > gaussPoints;
    this->_numberOfGaussPoints = fs[0]->getGaussInfo(this->_integrationType, elementType, _gaussWeights, gaussPoints);

    std::vector< int > orientation;
    std::vector< scalar::Precision< T_Scalar > > functions;

    // load orientations and basis functions Lhs
    if(!functionSpaces.have(fs[0]->getGmshFemName(_isDerivativeLhs), this->_integrationType)) {
      this->_nbrOfDofsByElement[0] = fs[0]->getBasisFunctions(_isDerivativeLhs, functions, orientation, gaussPoints, elementType, entity);
      functionSpaces.functionSpaceOrientation(fs[0]->getGmshFemOrientationName(), orientation);
      functionSpaces.functionSpace(fs[0]->getGmshFemName(_isDerivativeLhs), this->_integrationType, functions);
      functionSpaces.nbrOfDofsByElement(fs[0]->getGmshFemName(_isDerivativeLhs), this->_integrationType, this->_nbrOfDofsByElement[0]);
    }
    _fsIndexLhs = functionSpaces.functionSpaceOrientation(fs[0]->getGmshFemOrientationName());
    _fsLhs.resize(omp::getMaxThreads());
#pragma omp parallel num_threads(omp::getMaxThreads())
    {
      const unsigned int myThreadID = omp::getThreadNum();
      _fsLhs[myThreadID] = *functionSpaces.functionSpace(fs[0]->getGmshFemName(_isDerivativeLhs), this->_integrationType);
    }
    this->_nbrOfDofsByElement[0] = functionSpaces.nbrOfDofsByElement(fs[0]->getGmshFemName(_isDerivativeLhs), this->_integrationType);

    _needOffsetLhs = !fs[0]->isConstantByElements();

    // load orientations and basis functions Rhs
    if(!functionSpaces.have(fs[1]->getGmshFemName(_isDerivativeRhs), this->_integrationType)) {
      this->_nbrOfDofsByElement[1] = fs[1]->getBasisFunctions(_isDerivativeRhs, functions, orientation, gaussPoints, elementType, entity);
      functionSpaces.functionSpaceOrientation(fs[1]->getGmshFemOrientationName(), orientation);
      functionSpaces.functionSpace(fs[1]->getGmshFemName(_isDerivativeRhs), this->_integrationType, functions);
      functionSpaces.nbrOfDofsByElement(fs[1]->getGmshFemName(_isDerivativeRhs), this->_integrationType, this->_nbrOfDofsByElement[1]);
    }
    _fsIndexRhs = functionSpaces.functionSpaceOrientation(fs[1]->getGmshFemOrientationName());
    _fsRhs.resize(omp::getMaxThreads());
    if(this->_field[0] != this->_field[1] || _isDerivativeLhs != _isDerivativeRhs) {
#pragma omp parallel num_threads(omp::getMaxThreads())
      {
        const unsigned int myThreadID = omp::getThreadNum();
        _fsRhs[myThreadID] = *functionSpaces.functionSpace(fs[1]->getGmshFemName(_isDerivativeRhs), this->_integrationType);
      }
    }
    this->_nbrOfDofsByElement[1] = functionSpaces.nbrOfDofsByElement(fs[1]->getGmshFemName(_isDerivativeRhs), this->_integrationType);

    _needOffsetRhs = !fs[1]->isConstantByElements();

    // **************
    // Lhs part
    // **************
    _equationEvaluatorLhs->initialize(points, gaussPoints, elementType, entity);

    // **************
    // Rhs part
    // **************
    _equationEvaluatorRhs->initialize(points, gaussPoints, elementType, entity);

    // **************
    // Instantiate tmp matrices
    // **************
    if(this->_field[0] == this->_field[1] && _isDerivativeLhs == _isDerivativeRhs) {
      _fieldExpressionLhs = new std::vector< Eigen::MatrixX< scalar::Precision< T_Scalar > >, numa::allocator< Eigen::MatrixX< scalar::Precision< T_Scalar > > > >();
      _fieldExpressionRhs = _fieldExpressionLhs;

      _fieldExpressionLhs->resize(omp::getMaxThreads(), Eigen::MatrixX< scalar::Precision< T_Scalar > >(_fieldEvaluatorLhs->size(), this->_field[0]->multiplicity() * this->_nbrOfDofsByElement[0]));
    }
    else {
      _fieldExpressionLhs = new std::vector< Eigen::MatrixX< scalar::Precision< T_Scalar > >, numa::allocator< Eigen::MatrixX< scalar::Precision< T_Scalar > > > >();
      _fieldExpressionRhs = new std::vector< Eigen::MatrixX< scalar::Precision< T_Scalar > >, numa::allocator< Eigen::MatrixX< scalar::Precision< T_Scalar > > > >();

      _fieldExpressionLhs->resize(omp::getMaxThreads(), Eigen::MatrixX< scalar::Precision< T_Scalar > >(_fieldEvaluatorLhs->size(), this->_field[0]->multiplicity() * this->_nbrOfDofsByElement[0]));
      _fieldExpressionRhs->resize(omp::getMaxThreads(), Eigen::MatrixX< scalar::Precision< T_Scalar > >(_fieldEvaluatorRhs->size(), this->_field[1]->multiplicity() * this->_nbrOfDofsByElement[1]));
    }

#pragma omp parallel for
    for(unsigned int i = 0; i < omp::getMaxThreads(); ++i) {
      (*_fieldExpressionLhs)[i].setZero();
      (*_fieldExpressionRhs)[i].setZero();
    }

    unsigned int size1Lhs = _equationEvaluatorLhs->size1(), size2Lhs = _equationEvaluatorLhs->size2();
    if(size1Lhs == 0) {
      size1Lhs = this->_field[0]->multiplicity() * this->_nbrOfDofsByElement[0];
    }
    else if(size1Lhs == 4) {
      size1Lhs = _fieldEvaluatorLhs->size();
    }
    if(size2Lhs == 0) {
      size2Lhs = this->_field[0]->multiplicity() * this->_nbrOfDofsByElement[0];
    }
    else if(size2Lhs == 4) {
      size2Lhs = _fieldEvaluatorLhs->size();
    }
    _expressionLhs.resize(omp::getMaxThreads(), Eigen::MatrixX< T_Scalar >(size1Lhs, size2Lhs));

    unsigned int size1Rhs = _equationEvaluatorRhs->size1(), size2Rhs = _equationEvaluatorRhs->size2();
    if(size1Rhs == 0) {
      size1Rhs = this->_field[1]->multiplicity() * this->_nbrOfDofsByElement[1];
    }
    else if(size1Rhs == 4) {
      size1Rhs = _fieldEvaluatorRhs->size();
    }
    if(size2Rhs == 0) {
      size2Rhs = this->_field[1]->multiplicity() * this->_nbrOfDofsByElement[1];
    }
    else if(size2Rhs == 4) {
      size2Rhs = _fieldEvaluatorRhs->size();
    }
    _expressionRhs.resize(omp::getMaxThreads(), Eigen::MatrixX< T_Scalar >(size1Rhs, size2Rhs));
  }

  template< class T_Scalar, field::Form T_FormLhs, field::Form T_FormRhs >
  BilinearTerm< T_Scalar, T_FormLhs, T_FormRhs > *BilinearTerm< T_Scalar, T_FormLhs, T_FormRhs >::switchTestFunctionField(const field::FieldInterface< T_Scalar > *primal, field::FieldInterface< T_Scalar > *dual) const
  {
    if(_degree == Degree::Degree0) {
      if(_equationRhs0->getField()->tag() != primal->tag()) {
        return nullptr;
      }

      if(auto ptrDom = dynamic_cast< const domain::Domain * >(&this->_domain)) {
        return new BilinearTerm< T_Scalar, T_FormLhs, T_FormRhs >(ptrDom->getSkinLayer(dual->domain()), this->_integrationType, *_equationLhs0, _equationRhs0->changeUnknownField(dual), this->_productType);
      }

      return new BilinearTerm< T_Scalar, T_FormLhs, T_FormRhs >(this->_domain, this->_integrationType, *_equationLhs0, _equationRhs0->changeUnknownField(dual), this->_productType);
    }
    else if(_degree == Degree::Degree1) {
      if(_equationRhs1->getField()->tag() != primal->tag()) {
        return nullptr;
      }

      if(auto ptrDom = dynamic_cast< const domain::Domain * >(&this->_domain)) {
        return new BilinearTerm< T_Scalar, T_FormLhs, T_FormRhs >(ptrDom->getSkinLayer(dual->domain()), this->_integrationType, *_equationLhs1, _equationRhs1->changeUnknownField(dual), this->_productType);
      }

      return new BilinearTerm< T_Scalar, T_FormLhs, T_FormRhs >(this->_domain, this->_integrationType, *_equationLhs1, _equationRhs1->changeUnknownField(dual), this->_productType);
    }
    return nullptr;
  }

  template< class T_Scalar, field::Form T_FormLhs, field::Form T_FormRhs >
  void BilinearTerm< T_Scalar, T_FormLhs, T_FormRhs >::evaluate(Eigen::MatrixX< T_Scalar > &A_e, const scalar::Precision< T_Scalar > *const determinants, const scalar::Precision< T_Scalar > *const jacobians, const unsigned int elementIndex)
  {
    const unsigned int myThreadID = omp::getThreadNum();
    Eigen::MatrixX< scalar::Precision< T_Scalar > > &fieldExpressionRhs = (*_fieldExpressionRhs)[myThreadID];
    Eigen::MatrixX< scalar::Precision< T_Scalar > > &fieldExpressionLhs = (*_fieldExpressionLhs)[myThreadID];
    Eigen::MatrixX< T_Scalar > &expressionRhs = _expressionRhs[myThreadID];
    Eigen::MatrixX< T_Scalar > &expressionLhs = _expressionLhs[myThreadID];
    const std::vector< scalar::Precision< T_Scalar > > &fsLhs = _fsLhs[myThreadID];
    const std::vector< scalar::Precision< T_Scalar > > &fsRhs = _fsRhs[myThreadID];

    if(this->_productType == term::ProductType::Hermitian) {
      for(auto i = 0U; i < this->_numberOfGaussPoints; ++i) {
        (*_fieldEvaluatorLhs)(fieldExpressionLhs, &fsLhs[((_needOffsetLhs ? (*_fsIndexLhs)[elementIndex] * this->_numberOfGaussPoints : 0) + i) * field::NumberOfComponentsOfForm< T_FormLhs >::value * this->_nbrOfDofsByElement[0]], &jacobians[i * 9], &determinants[i], this->_nbrOfDofsByElement[0]);
        if(_fieldExpressionLhs != _fieldExpressionRhs) {
          (*_fieldEvaluatorRhs)(fieldExpressionRhs, &fsRhs[((_needOffsetRhs ? (*_fsIndexRhs)[elementIndex] * this->_numberOfGaussPoints : 0) + i) * field::NumberOfComponentsOfForm< T_FormRhs >::value * this->_nbrOfDofsByElement[1]], &jacobians[i * 9], &determinants[i], this->_nbrOfDofsByElement[1]);
        }

        (*_equationEvaluatorLhs)(expressionLhs, fieldExpressionLhs, elementIndex * this->_numberOfGaussPoints + i);
        (*_equationEvaluatorRhs)(expressionRhs, fieldExpressionRhs, elementIndex * this->_numberOfGaussPoints + i);

        A_e.noalias() += _gaussWeights[i] * determinants[i] * expressionLhs.transpose() * expressionRhs.conjugate();
      }
    }
    else if(this->_productType == term::ProductType::Scalar) {
      for(auto i = 0U; i < this->_numberOfGaussPoints; ++i) {
        (*_fieldEvaluatorLhs)(fieldExpressionLhs, &fsLhs[((_needOffsetLhs ? (*_fsIndexLhs)[elementIndex] * this->_numberOfGaussPoints : 0) + i) * field::NumberOfComponentsOfForm< T_FormLhs >::value * this->_nbrOfDofsByElement[0]], &jacobians[i * 9], &determinants[i], this->_nbrOfDofsByElement[0]);
        if(_fieldExpressionLhs != _fieldExpressionRhs) {
          (*_fieldEvaluatorRhs)(fieldExpressionRhs, &fsRhs[((_needOffsetRhs ? (*_fsIndexRhs)[elementIndex] * this->_numberOfGaussPoints : 0) + i) * field::NumberOfComponentsOfForm< T_FormRhs >::value * this->_nbrOfDofsByElement[1]], &jacobians[i * 9], &determinants[i], this->_nbrOfDofsByElement[1]);
        }

        (*_equationEvaluatorLhs)(expressionLhs, fieldExpressionLhs, elementIndex * this->_numberOfGaussPoints + i);
        (*_equationEvaluatorRhs)(expressionRhs, fieldExpressionRhs, elementIndex * this->_numberOfGaussPoints + i);

        A_e.noalias() += _gaussWeights[i] * determinants[i] * expressionLhs.transpose() * expressionRhs;
      }
    }
    else {
      throw common::Exception("Unknown product type in bilinear term");
    }
  }

  INSTANTIATE_CLASS_3(BilinearTerm, 4, 4, 4, TEMPLATE_ARGS(std::complex< double >, std::complex< float >, double, float), TEMPLATE_ARGS(field::Form::Form0, field::Form::Form1, field::Form::Form2, field::Form::Form3), TEMPLATE_ARGS(field::Form::Form0, field::Form::Form1, field::Form::Form2, field::Form::Form3))


} // namespace gmshfem::term
