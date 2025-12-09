// GmshFEM - Copyright (C) 2019-2022, A. Royer, E. Béchet, C. Geuzaine, Université de Liège
//
// See the LICENSE.txt file for license information. Please report all
// issues on https://gitlab.onelab.info/gmsh/fem/issues

#ifndef H_GMSHFEM_LINEARTERM
#define H_GMSHFEM_LINEARTERM

#include "Equation.h"
#include "FieldObject.h"
#include "Function.h"
#include "LinearTermInterface.h"

namespace gmshfem::term::evaluator
{
  template< class T_Scalar, field::Form T_Form >
  class FieldEvaluator;

  template< class T_Scalar, field::Form T_Form >
  class EquationEvaluatorInterface;
} // namespace gmshfem::term::evaluator

namespace gmshfem::term
{


  template< class T_Scalar, Degree T_DegreeLhs, field::Form T_FormRhs >
  class LinearTerm final : public LinearTermInterface< T_Scalar >
  {
   protected:
    const equation::EquationInterface< T_Scalar, T_DegreeLhs, T_FormRhs > *const _equationRhs;
    const function::Function< T_Scalar, T_DegreeLhs > _functionLhs;
    const evaluator::FieldEvaluator< T_Scalar, T_FormRhs > *const _fieldEvaluatorRhs;
    const evaluator::EquationEvaluatorInterface< T_Scalar, T_FormRhs > *const _equationEvaluatorRhs;
    const bool _isDerivativeRhs;
    std::vector< typename MathObject< T_Scalar, T_DegreeLhs >::Object, numa::allocator< typename MathObject< T_Scalar, T_DegreeLhs >::Object > > _valuesLhs;
    bool _lhsIsConstant;

    const std::vector< scalar::Precision< T_Scalar > > *_fsRhs;
    const std::vector< int > *_fsIndexRhs;
    bool _needOffsetRhs;
    std::vector< scalar::Precision< T_Scalar > > _gaussWeights;

    // Tmp assembly matrices
    std::vector< Eigen::MatrixX< scalar::Precision< T_Scalar > >, numa::allocator< Eigen::MatrixX< scalar::Precision< T_Scalar > > > > _fieldExpressionRhs;
    std::vector< Eigen::MatrixX< T_Scalar >, numa::allocator< Eigen::MatrixX< T_Scalar > > > _expressionRhs;

    void _freePrecomputedDate();

   public:
    LinearTerm(const domain::GeometricObject &domain, const std::string &integrationType, const function::Function< T_Scalar, T_DegreeLhs > &functionLhs, const equation::EquationInterface< T_Scalar, T_DegreeLhs, T_FormRhs > &equationRhs, const ProductType productType, unsigned rhsIdx = 0);
    LinearTerm(const domain::GeometricObject &domain, const std::string &integrationType, const function::Function< T_Scalar, T_DegreeLhs > &functionLhs, equation::EquationInterface< T_Scalar, T_DegreeLhs, T_FormRhs > *equationRhs, const ProductType productType, unsigned rhsIdx = 0);
    virtual ~LinearTerm();

    virtual bool needJacobians() const override;
    virtual bool needGaussCoordinates(const std::pair< int, int > &entity) const override;
    virtual common::Memory memory() const override;

    virtual void assemblyInitialization(problem::FunctionSpaceBucket< scalar::Precision< T_Scalar > > &functionSpaces, const std::vector< scalar::Precision< T_Scalar >, numa::allocator< scalar::Precision< T_Scalar > > > &points, const int elementType, const std::pair< int, int > &entity) override;
    virtual LinearTerm< T_Scalar, T_DegreeLhs, T_FormRhs > *switchTestFunctionField(const field::FieldInterface< T_Scalar > *primal, field::FieldInterface< T_Scalar > *dual) const override;

    virtual void evaluate(Eigen::VectorX< T_Scalar > &b_e, const scalar::Precision< T_Scalar > *const determinants, const scalar::Precision< T_Scalar > *const jacobians, const unsigned int elementIndex) override;
  };


} // namespace gmshfem::term

#endif // H_GMSHFEM_LINEARTERM
