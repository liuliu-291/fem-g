// GmshFEM - Copyright (C) 2019-2022, A. Royer, E. Béchet, C. Geuzaine, Université de Liège
//
// See the LICENSE.txt file for license information. Please report all
// issues on https://gitlab.onelab.info/gmsh/fem/issues

#include "LinearTerm.h"

#include "EquationEvaluator.h"
#include "FieldEvaluator.h"
#include "FieldInterface.h"
#include "FunctionSpaceBucket.h"
#include "instantiate.h"

//
// class LinearTerm : public LinearTermInterface
//

namespace gmshfem::term
{


  template< class T_Scalar, Degree T_DegreeLhs, field::Form T_FormRhs >
  LinearTerm< T_Scalar, T_DegreeLhs, T_FormRhs >::LinearTerm(const domain::GeometricObject &domain, const std::string &integrationType, const function::Function< T_Scalar, T_DegreeLhs > &functionLhs, const equation::EquationInterface< T_Scalar, T_DegreeLhs, T_FormRhs > &equationRhs, const ProductType productType, unsigned rhsIdx ) :
    LinearTermInterface< T_Scalar >(domain, integrationType, equationRhs.getField(), productType, rhsIdx), _equationRhs(equationRhs.copy()), _functionLhs(functionLhs), _fieldEvaluatorRhs(equationRhs.getFieldEvaluator()), _equationEvaluatorRhs(equationRhs.getEquationEvaluator()), _isDerivativeRhs(equationRhs.isDerivative()), _valuesLhs(), _lhsIsConstant(false), _fsRhs(nullptr), _fsIndexRhs(nullptr), _needOffsetRhs(false), _gaussWeights()
  {
  }

  template< class T_Scalar, Degree T_DegreeLhs, field::Form T_FormRhs >
  LinearTerm< T_Scalar, T_DegreeLhs, T_FormRhs >::LinearTerm(const domain::GeometricObject &domain, const std::string &integrationType, const function::Function< T_Scalar, T_DegreeLhs > &functionLhs, equation::EquationInterface< T_Scalar, T_DegreeLhs, T_FormRhs > *equationRhs, const ProductType productType, unsigned rhsIdx) :
    LinearTermInterface< T_Scalar >(domain, integrationType, equationRhs->getField(), productType, rhsIdx), _equationRhs(equationRhs), _functionLhs(functionLhs), _fieldEvaluatorRhs(equationRhs->getFieldEvaluator()), _equationEvaluatorRhs(equationRhs->getEquationEvaluator()), _isDerivativeRhs(equationRhs->isDerivative()), _valuesLhs(), _lhsIsConstant(false), _fsRhs(nullptr), _fsIndexRhs(nullptr), _needOffsetRhs(false), _gaussWeights()
  {
  }

  template< class T_Scalar, Degree T_DegreeLhs, field::Form T_FormRhs >
  LinearTerm< T_Scalar, T_DegreeLhs, T_FormRhs >::~LinearTerm()
  {
    _freePrecomputedDate();
    this->_numberOfGaussPoints = 0;
    delete _equationRhs;
    delete _fieldEvaluatorRhs;
    delete _equationEvaluatorRhs;
  }

  template< class T_Scalar, Degree T_DegreeLhs, field::Form T_FormRhs >
  bool LinearTerm< T_Scalar, T_DegreeLhs, T_FormRhs >::needJacobians() const
  {
    return _fieldEvaluatorRhs->needJacobians();
  }

  template< class T_Scalar, Degree T_DegreeLhs, field::Form T_FormRhs >
  bool LinearTerm< T_Scalar, T_DegreeLhs, T_FormRhs >::needGaussCoordinates(const std::pair< int, int > &entity) const
  {
    if(this->_domain.needGaussCoordinates()) {
      return true;
    }

    if(!_functionLhs.isConstant(entity)) {
      return true;
    }

    if(!_equationEvaluatorRhs->isConstant(entity)) {
      return true;
    }

    return false;
  }

  template< class T_Scalar, Degree T_DegreeLhs, field::Form T_FormRhs >
  common::Memory LinearTerm< T_Scalar, T_DegreeLhs, T_FormRhs >::memory() const
  {
    return common::Memory(_valuesLhs.size() * sizeof(typename MathObject< T_Scalar, T_DegreeLhs >::Object) + _gaussWeights.size() * sizeof(scalar::Precision< T_Scalar >));
  }

  template< class T_Scalar, Degree T_DegreeLhs, field::Form T_FormRhs >
  void LinearTerm< T_Scalar, T_DegreeLhs, T_FormRhs >::_freePrecomputedDate()
  {
    _fsRhs = nullptr;
    _fsIndexRhs = nullptr;
    _gaussWeights.clear();
    _valuesLhs.clear();
    _equationEvaluatorRhs->clear();
    _fieldExpressionRhs.clear();
    _expressionRhs.clear();
  }

  template< class T_Scalar, Degree T_DegreeLhs, field::Form T_FormRhs >
  void LinearTerm< T_Scalar, T_DegreeLhs, T_FormRhs >::assemblyInitialization(problem::FunctionSpaceBucket< scalar::Precision< T_Scalar > > &functionSpaces, const std::vector< scalar::Precision< T_Scalar >, numa::allocator< scalar::Precision< T_Scalar > > > &points, const int elementType, const std::pair< int, int > &entity)
  {
    _freePrecomputedDate();

    // **************
    // Function space
    // **************
    field::FunctionSpaceInterface< scalar::Precision< T_Scalar > > *fs = this->_fieldRhs->getFunctionSpace();
    // load Gauss information
    std::vector< scalar::Precision< T_Scalar > > gaussPoints;
    this->_numberOfGaussPoints = fs->getGaussInfo(this->_integrationType, elementType, _gaussWeights, gaussPoints);

    // load orientations and basis functions
    std::vector< int > orientation;
    std::vector< scalar::Precision< T_Scalar > > functions;
    if(!functionSpaces.have(fs->getGmshFemName(_isDerivativeRhs), this->_integrationType)) {
      this->_nbrOfDofsByElement = fs->getBasisFunctions(_isDerivativeRhs, functions, orientation, gaussPoints, elementType, entity);
      functionSpaces.functionSpaceOrientation(fs->getGmshFemOrientationName(), orientation);
      functionSpaces.functionSpace(fs->getGmshFemName(_isDerivativeRhs), this->_integrationType, functions);
      functionSpaces.nbrOfDofsByElement(fs->getGmshFemName(_isDerivativeRhs), this->_integrationType, this->_nbrOfDofsByElement);
    }
    _fsIndexRhs = functionSpaces.functionSpaceOrientation(fs->getGmshFemOrientationName());
    _fsRhs = functionSpaces.functionSpace(fs->getGmshFemName(_isDerivativeRhs), this->_integrationType);
    this->_nbrOfDofsByElement = functionSpaces.nbrOfDofsByElement(fs->getGmshFemName(_isDerivativeRhs), this->_integrationType);

    _needOffsetRhs = !fs->isConstantByElements();

    // **************
    // Lhs part
    // **************
    if(_functionLhs.isConstant(entity)) {
      _valuesLhs.resize(1);
      _functionLhs.evaluate(_valuesLhs[0], 0., 0., 0., entity);
      _lhsIsConstant = true;
    }
    else {
#pragma omp parallel num_threads(omp::getMaxThreads())
      {
        _functionLhs.evaluate(_valuesLhs, points, gaussPoints, elementType, entity);
      }
      _lhsIsConstant = false;
    }

    // **************
    // Rhs part
    // **************
    _equationEvaluatorRhs->initialize(points, gaussPoints, elementType, entity);

    // **************
    // Instantiate tmp matrices
    // **************
    _fieldExpressionRhs.resize(omp::getMaxThreads(), Eigen::MatrixX< scalar::Precision< T_Scalar > >(_fieldEvaluatorRhs->size(), this->_fieldRhs->multiplicity() * this->_nbrOfDofsByElement));
    
  #pragma omp parallel for
    for(unsigned int i = 0; i < omp::getMaxThreads(); ++i) {
      _fieldExpressionRhs[i].setZero();
    }

    unsigned int size1Rhs = _equationEvaluatorRhs->size1(), size2Rhs = _equationEvaluatorRhs->size2();
    if(size1Rhs == 0) {
      size1Rhs = this->_fieldRhs->multiplicity() * this->_nbrOfDofsByElement;
    }
    else if(size1Rhs == 4) {
      size1Rhs = _fieldEvaluatorRhs->size();
    }
    if(size2Rhs == 0) {
      size2Rhs = this->_fieldRhs->multiplicity() * this->_nbrOfDofsByElement;
    }
    else if(size2Rhs == 4) {
      size2Rhs = _fieldEvaluatorRhs->size();
    }
    _expressionRhs.resize(omp::getMaxThreads(), Eigen::MatrixX< T_Scalar >(size1Rhs, size2Rhs));
  }

  template< class T_Scalar, Degree T_DegreeLhs, field::Form T_FormRhs >
  LinearTerm< T_Scalar, T_DegreeLhs, T_FormRhs > *LinearTerm< T_Scalar, T_DegreeLhs, T_FormRhs >::switchTestFunctionField(const field::FieldInterface< T_Scalar > *primal, field::FieldInterface< T_Scalar > *dual) const
  {
    if(_equationRhs->getField()->tag() != primal->tag()) {
      return nullptr;
    }

    if(auto ptrDom = dynamic_cast< const domain::Domain * >(&this->_domain)) {
      return new LinearTerm< T_Scalar, T_DegreeLhs, T_FormRhs >(ptrDom->getSkinLayer(dual->domain()), this->_integrationType, _functionLhs, _equationRhs->changeUnknownField(dual), this->_productType);
    }

    return new LinearTerm< T_Scalar, T_DegreeLhs, T_FormRhs >(this->_domain, this->_integrationType, _functionLhs, _equationRhs->changeUnknownField(dual), this->_productType);
  }

  template< class T_Scalar, Degree T_DegreeLhs, field::Form T_FormRhs >
  void LinearTerm< T_Scalar, T_DegreeLhs, T_FormRhs >::evaluate(Eigen::VectorX< T_Scalar > &b_e, const scalar::Precision< T_Scalar > *const determinants, const scalar::Precision< T_Scalar > *const jacobians, const unsigned int elementIndex)
  {
    const unsigned int myThreadID = omp::getThreadNum();
    Eigen::MatrixX< scalar::Precision< T_Scalar > > &fieldExpressionRhs = _fieldExpressionRhs[myThreadID];
    Eigen::MatrixX< T_Scalar > &expressionRhs = _expressionRhs[myThreadID];

    if(this->_productType == term::ProductType::Hermitian) {
      for(auto i = 0U; i < this->_numberOfGaussPoints; ++i) {
        (*_fieldEvaluatorRhs)(fieldExpressionRhs, &(*_fsRhs)[((_needOffsetRhs ? (*_fsIndexRhs)[elementIndex] * this->_numberOfGaussPoints : 0) + i) * field::NumberOfComponentsOfForm< T_FormRhs >::value * this->_nbrOfDofsByElement], &jacobians[i * 9], &determinants[i], this->_nbrOfDofsByElement);

        (*_equationEvaluatorRhs)(expressionRhs, fieldExpressionRhs, elementIndex * this->_numberOfGaussPoints + i);

        if constexpr(T_DegreeLhs != Degree::Degree2) {
          b_e.noalias() += _gaussWeights[i] * determinants[i] * expressionRhs.adjoint() * _valuesLhs[_lhsIsConstant ? 0 : elementIndex * this->_numberOfGaussPoints + i];
        }
        else {
          for(unsigned int d = 0; d < this->_fieldRhs->multiplicity() * this->_nbrOfDofsByElement; ++d) {
            b_e(d) += _gaussWeights[i] * determinants[i] * (_valuesLhs[_lhsIsConstant ? 0 : elementIndex * this->_numberOfGaussPoints + i] * Eigen::Map< const Eigen::Matrix3< T_Scalar > >(expressionRhs.col(d).data(), 3, 3).adjoint()).trace();
          }
        }
      }
    }
    else if(this->_productType == term::ProductType::Scalar) {
      for(auto i = 0U; i < this->_numberOfGaussPoints; ++i) {
        (*_fieldEvaluatorRhs)(fieldExpressionRhs, &(*_fsRhs)[((_needOffsetRhs ? (*_fsIndexRhs)[elementIndex] * this->_numberOfGaussPoints : 0) + i) * field::NumberOfComponentsOfForm< T_FormRhs >::value * this->_nbrOfDofsByElement], &jacobians[i * 9], &determinants[i], this->_nbrOfDofsByElement);

        (*_equationEvaluatorRhs)(expressionRhs, fieldExpressionRhs, elementIndex * this->_numberOfGaussPoints + i);

        if constexpr(T_DegreeLhs != Degree::Degree2) {
          b_e.noalias() += _gaussWeights[i] * determinants[i] * expressionRhs.transpose() * _valuesLhs[_lhsIsConstant ? 0 : elementIndex * this->_numberOfGaussPoints + i];
        }
        else {
          for(unsigned int d = 0; d < this->_fieldRhs->multiplicity() * this->_nbrOfDofsByElement; ++d) {
            b_e(d) += _gaussWeights[i] * determinants[i] * (_valuesLhs[_lhsIsConstant ? 0 : elementIndex * this->_numberOfGaussPoints + i] * Eigen::Map< const Eigen::Matrix3< T_Scalar > >(expressionRhs.col(d).data(), 3, 3).transpose()).trace();
          }
        }
      }
    }
    else {
      throw common::Exception("Unknown product type (in linear term)");
    }
  }

  INSTANTIATE_CLASS_3(LinearTerm, 4, 3, 4, TEMPLATE_ARGS(std::complex< double >, std::complex< float >, double, float), TEMPLATE_ARGS(Degree::Degree0, Degree::Degree1, Degree::Degree2), TEMPLATE_ARGS(field::Form::Form0, field::Form::Form1, field::Form::Form2, field::Form::Form3))


} // namespace gmshfem::term
