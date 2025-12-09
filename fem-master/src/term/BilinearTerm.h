// GmshFEM - Copyright (C) 2019-2022, A. Royer, E. Béchet, C. Geuzaine, Université de Liège
//
// See the LICENSE.txt file for license information. Please report all
// issues on https://gitlab.onelab.info/gmsh/fem/issues

#ifndef H_GMSHFEM_BILINEARTERM
#define H_GMSHFEM_BILINEARTERM

#include "BilinearTermInterface.h"
#include "Equation.h"
#include "FieldObject.h"

namespace gmshfem::term::evaluator
{
  template< class T_Scalar, field::Form T_Form >
  class FieldEvaluator;

  template< class T_Scalar, field::Form T_Form >
  class EquationEvaluatorInterface;
} // namespace gmshfem::term::evaluator

namespace gmshfem::term
{


  template< class T_Scalar, field::Form T_FormLhs, field::Form T_FormRhs >
  class BilinearTerm final : public BilinearTermInterface< T_Scalar >
  {
   private:
    const equation::EquationInterface< T_Scalar, Degree::Degree0, T_FormLhs > *const _equationLhs0;
    const equation::EquationInterface< T_Scalar, Degree::Degree1, T_FormLhs > *const _equationLhs1;
    const equation::EquationInterface< T_Scalar, Degree::Degree2, T_FormLhs > *const _equationLhs2;
    const evaluator::FieldEvaluator< T_Scalar, T_FormLhs > *const _fieldEvaluatorLhs;
    const evaluator::EquationEvaluatorInterface< T_Scalar, T_FormLhs > *const _equationEvaluatorLhs;
    const bool _isDerivativeLhs;
    const equation::EquationInterface< T_Scalar, Degree::Degree0, T_FormRhs > *const _equationRhs0;
    const equation::EquationInterface< T_Scalar, Degree::Degree1, T_FormRhs > *const _equationRhs1;
    const equation::EquationInterface< T_Scalar, Degree::Degree2, T_FormRhs > *const _equationRhs2;
    const evaluator::FieldEvaluator< T_Scalar, T_FormRhs > *const _fieldEvaluatorRhs;
    const evaluator::EquationEvaluatorInterface< T_Scalar, T_FormRhs > *const _equationEvaluatorRhs;
    const bool _isDerivativeRhs;

    const Degree _degree;

    std::vector< std::vector< scalar::Precision< T_Scalar > > > _fsLhs;
    const std::vector< int > *_fsIndexLhs;
    bool _needOffsetLhs;
    std::vector< std::vector< scalar::Precision< T_Scalar > > > _fsRhs;
    const std::vector< int > *_fsIndexRhs;
    bool _needOffsetRhs;
    std::vector< scalar::Precision< T_Scalar > > _gaussWeights;

    // Tmp assembly matrices
    std::vector< Eigen::MatrixX< scalar::Precision< T_Scalar > >, numa::allocator< Eigen::MatrixX< scalar::Precision< T_Scalar > > > > *_fieldExpressionLhs;
    std::vector< Eigen::MatrixX< scalar::Precision< T_Scalar > >, numa::allocator< Eigen::MatrixX< scalar::Precision< T_Scalar > > > > *_fieldExpressionRhs;
    std::vector< Eigen::MatrixX< T_Scalar >, numa::allocator< Eigen::MatrixX< T_Scalar > > > _expressionLhs;
    std::vector< Eigen::MatrixX< T_Scalar >, numa::allocator< Eigen::MatrixX< T_Scalar > > > _expressionRhs;

    void _freePrecomputedDate();

   public:
    BilinearTerm(const domain::GeometricObject &domain, const std::string &integrationType, const equation::EquationInterface< T_Scalar, Degree::Degree0, T_FormLhs > &equationLhs, const equation::EquationInterface< T_Scalar, Degree::Degree0, T_FormRhs > &equationRhs, const ProductType productType);
    BilinearTerm(const domain::GeometricObject &domain, const std::string &integrationType, const equation::EquationInterface< T_Scalar, Degree::Degree1, T_FormLhs > &equationLhs, const equation::EquationInterface< T_Scalar, Degree::Degree1, T_FormRhs > &equationRhs, const ProductType productType);
    BilinearTerm(const domain::GeometricObject &domain, const std::string &integrationType, const equation::EquationInterface< T_Scalar, Degree::Degree2, T_FormLhs > &equationLhs, const equation::EquationInterface< T_Scalar, Degree::Degree2, T_FormRhs > &equationRhs, const ProductType productType);
    BilinearTerm(const domain::GeometricObject &domain, const std::string &integrationType, const equation::EquationInterface< T_Scalar, Degree::Degree0, T_FormLhs > &equationLhs, equation::EquationInterface< T_Scalar, Degree::Degree0, T_FormRhs > *equationRhs, const ProductType productType);
    BilinearTerm(const domain::GeometricObject &domain, const std::string &integrationType, const equation::EquationInterface< T_Scalar, Degree::Degree1, T_FormLhs > &equationLhs, equation::EquationInterface< T_Scalar, Degree::Degree1, T_FormRhs > *equationRhs, const ProductType productType);
    BilinearTerm(const domain::GeometricObject &domain, const std::string &integrationType, const equation::EquationInterface< T_Scalar, Degree::Degree2, T_FormLhs > &equationLhs, equation::EquationInterface< T_Scalar, Degree::Degree2, T_FormRhs > *equationRhs, const ProductType productType);
    virtual ~BilinearTerm();

    virtual bool needJacobians() const override;
    bool needGaussCoordinates(const std::pair< int, int > &entity) const override;
    virtual common::Memory memory() const override;

    equation::UnknownFieldType unknownFieldType() const override;

    virtual void assemblyInitialization(problem::FunctionSpaceBucket< scalar::Precision< T_Scalar > > &functionSpaces, const std::vector< scalar::Precision< T_Scalar >, numa::allocator< scalar::Precision< T_Scalar > > > &points, const int elementType, const std::pair< int, int > &entity) override;
    virtual BilinearTerm< T_Scalar, T_FormLhs, T_FormRhs > *switchTestFunctionField(const field::FieldInterface< T_Scalar > *primal, field::FieldInterface< T_Scalar > *dual) const override;

    virtual void evaluate(Eigen::MatrixX< T_Scalar > &A_e, const scalar::Precision< T_Scalar > *const determinants, const scalar::Precision< T_Scalar > *const jacobians, const unsigned int elementIndex) override;
  };


} // namespace gmshfem::term

#endif // H_GMSHFEM_BILINEARTERM
