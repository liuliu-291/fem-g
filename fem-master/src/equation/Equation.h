// GmshFEM - Copyright (C) 2019-2022, A. Royer, E. Béchet, C. Geuzaine, Université de Liège
//
// See the LICENSE.txt file for license information. Please report all
// issues on https://gitlab.onelab.info/gmsh/fem/issues

#ifndef H_GMSHFEM_EQUATION
#define H_GMSHFEM_EQUATION

#include "EquationEvaluator.h"
#include "FieldEvaluator.h"
#include "FieldObject.h"
#include "Function.h"
#include "Product.h"

namespace gmshfem::equation
{
  template< class T_Scalar >
  class UnknownFieldInterface;

  template< class T_Scalar, field::Form T_Form >
  class UnknownField;

  template< class T_Scalar, field::Form T_Form, unsigned int T_NumFields >
  class UnknownCompoundField;
} // namespace gmshfem::equation

namespace gmshfem::field
{
  template< class T_Scalar >
  class FieldInterface;
}

namespace gmshfem::equation
{


  //  Function  : left    : < T_Scalar, T_DegreeL >
  //  Product   : leftPr  : < T_Scalar, T_DegreeL, T_DegreeF, T_DegreeRL >
  //  Field     : field   : < T_Scalar, T_Form = T_DegreeF >
  //  ----------------------------------------
  //  = < T_Scalar, T_DegreeRL >
  //  ----------------------------------------
  //  Product   : rightPr : < T_Scalar, T_DegreeRL, T_DegreeR, T_DegreeRR >
  //  Function  : right   : < T_Scalar, T_DegreeR>
  //  ----------------------------------------
  //  = < T_Scalar, T_DegreeRR >
  //

  // *************************************************
  // EquationInterface< T_Scalar, T_Degree, T_Form >
  // *************************************************


  template< class T_Scalar, Degree T_Degree, field::Form T_Form >
  class EquationInterface
  {
   public:
    EquationInterface() {}
    virtual ~EquationInterface() {}

    virtual field::Form form() const = 0;
    virtual field::FieldInterface< T_Scalar > *getField() const = 0;
    virtual term::evaluator::EquationEvaluatorInterface< T_Scalar, T_Form > *getEquationEvaluator() const = 0;
    virtual term::evaluator::FieldEvaluator< T_Scalar, T_Form > *getFieldEvaluator() const = 0;
    virtual const UnknownFieldInterface< T_Scalar > *getUnknownField() const = 0;

    virtual bool isDerivative() const = 0;
    virtual EquationInterface< T_Scalar, T_Degree, T_Form > *copy() const = 0;
    virtual EquationInterface< T_Scalar, T_Degree, T_Form > *changeUnknownField(field::FieldInterface< T_Scalar > *dual) const = 0;
  };

  // *************************************************
  // EquationLhs< T_Scalar, T_Degree, T_Form, T_Product >
  // *************************************************

  template< class T_Scalar, Degree T_Degree, field::Form T_Form, Product T_Product >
  class EquationLhs final : public EquationInterface< T_Scalar, DegreeOfProduct< T_Degree, T_Product, field::DegreeOfForm< T_Form >::value >::value, T_Form >
  {
   private:
    const function::Function< T_Scalar, T_Degree > _left;
    const UnknownField< T_Scalar, T_Form > *_field;

   public:
    EquationLhs(const EquationLhs &other);
    EquationLhs(EquationLhs &&other);
    ~EquationLhs();

    EquationLhs(const function::Function< T_Scalar, T_Degree > &left, const UnknownField< T_Scalar, T_Form > &field);

    virtual field::Form form() const override;
    virtual field::FieldInterface< T_Scalar > *getField() const override;
    virtual term::evaluator::EquationEvaluator< T_Scalar, T_Form, T_Degree, T_Product, Degree::Empty, Product::Empty, term::evaluator::EvaluationOrder::Left > *getEquationEvaluator() const override;
    virtual term::evaluator::FieldEvaluator< T_Scalar, T_Form > *getFieldEvaluator() const override;
    virtual const UnknownFieldInterface< T_Scalar > *getUnknownField() const override;

    virtual bool isDerivative() const override;
    virtual EquationLhs< T_Scalar, T_Degree, T_Form, T_Product > *copy() const override;
    virtual EquationLhs< T_Scalar, T_Degree, T_Form, T_Product > *changeUnknownField(field::FieldInterface< T_Scalar > *dual) const override;

    function::Function< T_Scalar, T_Degree > getFunction() const;
  };

  // *************************************************
  // EquationLhsCompound< T_Scalar, T_Degree, T_Form, T_Product >
  // *************************************************

  template< class T_Scalar, Degree T_Degree, field::Form T_Form, Product T_Product, unsigned int T_NumFields >
  class EquationLhsCompound final : public EquationInterface< T_Scalar, DegreeOfProduct< T_Degree, T_Product, field::DegreeOfCompoundForm< T_Form >::value >::value, T_Form >
  {
   private:
    const function::Function< T_Scalar, T_Degree > _left;
    const UnknownCompoundField< T_Scalar, T_Form, T_NumFields > *_field;

   public:
    EquationLhsCompound(const EquationLhsCompound &other);
    EquationLhsCompound(EquationLhsCompound &&other);
    ~EquationLhsCompound();

    EquationLhsCompound(const function::Function< T_Scalar, T_Degree > &left, const UnknownCompoundField< T_Scalar, T_Form, T_NumFields > &field);

    virtual field::Form form() const override;
    virtual field::FieldInterface< T_Scalar > *getField() const override;
    virtual term::evaluator::EquationEvaluatorCompound< T_Scalar, T_Form, T_Degree, T_Product, Degree::Empty, Product::Empty, term::evaluator::EvaluationOrder::Left > *getEquationEvaluator() const override;
    virtual term::evaluator::CompoundFieldEvaluator< T_Scalar, T_Form, T_NumFields > *getFieldEvaluator() const override;
    virtual const UnknownFieldInterface< T_Scalar > *getUnknownField() const override;

    virtual bool isDerivative() const override;
    virtual EquationLhsCompound< T_Scalar, T_Degree, T_Form, T_Product, T_NumFields > *copy() const override;
    virtual EquationLhsCompound< T_Scalar, T_Degree, T_Form, T_Product, T_NumFields > *changeUnknownField(field::FieldInterface< T_Scalar > *dual) const override;

    function::Function< T_Scalar, T_Degree > getFunction() const;
  };

  // *************************************************
  // EquationRhs< T_Scalar, T_Degree, T_Form, T_Product >
  // *************************************************

  template< class T_Scalar, Degree T_Degree, field::Form T_Form, Product T_Product >
  class EquationRhs final : public EquationInterface< T_Scalar, DegreeOfProduct< field::DegreeOfForm< T_Form >::value, T_Product, T_Degree >::value, T_Form >
  {
   private:
    const function::Function< T_Scalar, T_Degree > _right;
    const UnknownField< T_Scalar, T_Form > *_field;

   public:
    EquationRhs(const EquationRhs &other);
    EquationRhs(EquationRhs &&other);
    ~EquationRhs();

    EquationRhs(const function::Function< T_Scalar, T_Degree > &right, const UnknownField< T_Scalar, T_Form > &field);

    virtual field::Form form() const override;
    virtual field::FieldInterface< T_Scalar > *getField() const override;
    virtual term::evaluator::EquationEvaluator< T_Scalar, T_Form, Degree::Empty, Product::Empty, T_Degree, T_Product, term::evaluator::EvaluationOrder::Right > *getEquationEvaluator() const override;
    virtual term::evaluator::FieldEvaluator< T_Scalar, T_Form > *getFieldEvaluator() const override;
    virtual const UnknownFieldInterface< T_Scalar > *getUnknownField() const override;

    virtual bool isDerivative() const override;
    virtual EquationRhs< T_Scalar, T_Degree, T_Form, T_Product > *copy() const override;
    virtual EquationRhs< T_Scalar, T_Degree, T_Form, T_Product > *changeUnknownField(field::FieldInterface< T_Scalar > *dual) const override;

    function::Function< T_Scalar, T_Degree > getFunction() const;
  };

  // *************************************************
  // EquationRhsCompound< T_Scalar, T_Degree, T_Form, T_Product >
  // *************************************************

  template< class T_Scalar, Degree T_Degree, field::Form T_Form, Product T_Product, unsigned int T_NumFields >
  class EquationRhsCompound final : public EquationInterface< T_Scalar, DegreeOfProduct< field::DegreeOfCompoundForm< T_Form >::value, T_Product, T_Degree >::value, T_Form >
  {
   private:
    const function::Function< T_Scalar, T_Degree > _right;
    const UnknownCompoundField< T_Scalar, T_Form, T_NumFields > *_field;

   public:
    EquationRhsCompound(const EquationRhsCompound &other);
    EquationRhsCompound(EquationRhsCompound &&other);
    ~EquationRhsCompound();

    EquationRhsCompound(const function::Function< T_Scalar, T_Degree > &right, const UnknownCompoundField< T_Scalar, T_Form, T_NumFields > &field);

    virtual field::Form form() const override;
    virtual field::FieldInterface< T_Scalar > *getField() const override;
    virtual term::evaluator::EquationEvaluatorCompound< T_Scalar, T_Form, Degree::Empty, Product::Empty, T_Degree, T_Product, term::evaluator::EvaluationOrder::Right > *getEquationEvaluator() const override;
    virtual term::evaluator::CompoundFieldEvaluator< T_Scalar, T_Form, T_NumFields > *getFieldEvaluator() const override;
    virtual const UnknownFieldInterface< T_Scalar > *getUnknownField() const override;

    virtual bool isDerivative() const override;
    virtual EquationRhsCompound< T_Scalar, T_Degree, T_Form, T_Product, T_NumFields > *copy() const override;
    virtual EquationRhsCompound< T_Scalar, T_Degree, T_Form, T_Product, T_NumFields > *changeUnknownField(field::FieldInterface< T_Scalar > *dual) const override;

    function::Function< T_Scalar, T_Degree > getFunction() const;
  };

  // *************************************************
  // EquationRhsLhs< T_Scalar, T_DegreeLhs, T_DegreeRhs, T_Form, T_ProductLhs, T_DegreeRhs, EvaluationOrder T_EvaluationOrder >
  // *************************************************

  template< class T_Scalar, Degree T_DegreeLhs, Degree T_DegreeRhs, field::Form T_Form, Product T_ProductLhs, Product T_ProductRhs, term::evaluator::EvaluationOrder T_EvaluationOrder > // T_EvaluationOrder = Left -> (T_DegreeLhs * dof) * T_DegreeRhs, = Right - > T_DegreeLhs * (dof * T_DegreeRhs) >
  class EquationLhsRhs
  {
  };

  // Left

  template< class T_Scalar, Degree T_DegreeLhs, Degree T_DegreeRhs, field::Form T_Form, Product T_ProductLhs, Product T_ProductRhs >
  class EquationLhsRhs< T_Scalar, T_DegreeLhs, T_DegreeRhs, T_Form, T_ProductLhs, T_ProductRhs, term::evaluator::EvaluationOrder::Left > : public EquationInterface< T_Scalar, DegreeOfProduct< DegreeOfProduct< T_DegreeLhs, T_ProductLhs, field::DegreeOfForm< T_Form >::value >::value, T_ProductRhs, T_DegreeRhs >::value, T_Form >
  {
   private:
    const function::Function< T_Scalar, T_DegreeLhs > _left;
    const UnknownField< T_Scalar, T_Form > *_field;
    const function::Function< T_Scalar, T_DegreeRhs > _right;

   public:
    EquationLhsRhs(const EquationLhsRhs &other);
    EquationLhsRhs(EquationLhsRhs &&other);
    ~EquationLhsRhs();

    EquationLhsRhs(const function::Function< T_Scalar, T_DegreeLhs > &left, const UnknownField< T_Scalar, T_Form > &field, const function::Function< T_Scalar, T_DegreeRhs > &right);

    virtual field::Form form() const override;
    virtual field::FieldInterface< T_Scalar > *getField() const override;
    virtual term::evaluator::EquationEvaluator< T_Scalar, T_Form, T_DegreeLhs, T_ProductLhs, T_DegreeRhs, T_ProductRhs, term::evaluator::EvaluationOrder::Left > *getEquationEvaluator() const override;
    virtual term::evaluator::FieldEvaluator< T_Scalar, T_Form > *getFieldEvaluator() const override;
    virtual const UnknownFieldInterface< T_Scalar > *getUnknownField() const override;

    virtual bool isDerivative() const override;
    virtual EquationLhsRhs< T_Scalar, T_DegreeLhs, T_DegreeRhs, T_Form, T_ProductLhs, T_ProductRhs, term::evaluator::EvaluationOrder::Left > *copy() const override;
    virtual EquationLhsRhs< T_Scalar, T_DegreeLhs, T_DegreeRhs, T_Form, T_ProductLhs, T_ProductRhs, term::evaluator::EvaluationOrder::Left > *changeUnknownField(field::FieldInterface< T_Scalar > *dual) const override;
  };

  // Right

  template< class T_Scalar, Degree T_DegreeLhs, Degree T_DegreeRhs, field::Form T_Form, Product T_ProductLhs, Product T_ProductRhs >
  class EquationLhsRhs< T_Scalar, T_DegreeLhs, T_DegreeRhs, T_Form, T_ProductLhs, T_ProductRhs, term::evaluator::EvaluationOrder::Right > : public EquationInterface< T_Scalar, DegreeOfProduct< T_DegreeRhs, T_ProductRhs, DegreeOfProduct< field::DegreeOfForm< T_Form >::value, T_ProductRhs, T_DegreeRhs >::value >::value, T_Form >
  {
   private:
    const function::Function< T_Scalar, T_DegreeLhs > _left;
    const UnknownField< T_Scalar, T_Form > *_field;
    const function::Function< T_Scalar, T_DegreeRhs > _right;

   public:
    EquationLhsRhs(const EquationLhsRhs &other);
    EquationLhsRhs(EquationLhsRhs &&other);
    ~EquationLhsRhs();

    EquationLhsRhs(const function::Function< T_Scalar, T_DegreeLhs > &left, const UnknownField< T_Scalar, T_Form > &field, const function::Function< T_Scalar, T_DegreeRhs > &right);

    virtual field::Form form() const override;
    virtual field::FieldInterface< T_Scalar > *getField() const override;
    virtual term::evaluator::EquationEvaluator< T_Scalar, T_Form, T_DegreeLhs, T_ProductLhs, T_DegreeRhs, T_ProductRhs, term::evaluator::EvaluationOrder::Right > *getEquationEvaluator() const override;
    virtual term::evaluator::FieldEvaluator< T_Scalar, T_Form > *getFieldEvaluator() const override;
    virtual const UnknownFieldInterface< T_Scalar > *getUnknownField() const override;

    virtual bool isDerivative() const override;
    virtual EquationLhsRhs< T_Scalar, T_DegreeLhs, T_DegreeRhs, T_Form, T_ProductLhs, T_ProductRhs, term::evaluator::EvaluationOrder::Right > *copy() const override;
    virtual EquationLhsRhs< T_Scalar, T_DegreeLhs, T_DegreeRhs, T_Form, T_ProductLhs, T_ProductRhs, term::evaluator::EvaluationOrder::Right > *changeUnknownField(field::FieldInterface< T_Scalar > *dual) const override;
  };


  template< class T_Scalar, field::Form T_Form >
  EquationLhs< T_Scalar, Degree::Degree0, T_Form, Product::ScalarProduct > operator+(const UnknownField< T_Scalar, T_Form > &field);
  template< class T_Scalar, field::Form T_Form >
  EquationLhs< T_Scalar, Degree::Degree0, T_Form, Product::ScalarProduct > operator-(const UnknownField< T_Scalar, T_Form > &field);

  template< class T_Scalar, field::Form T_Form, unsigned int T_NumFields >
  EquationLhsCompound< T_Scalar, Degree::Degree0, T_Form, Product::ScalarProduct, T_NumFields > operator+(const UnknownCompoundField< T_Scalar, T_Form, T_NumFields > &field);
  template< class T_Scalar, field::Form T_Form, unsigned int T_NumFields >
  EquationLhsCompound< T_Scalar, Degree::Degree0, T_Form, Product::ScalarProduct, T_NumFields > operator-(const UnknownCompoundField< T_Scalar, T_Form, T_NumFields > &field);

  // Lhs
  template< class T_Scalar, Degree T_Degree, field::Form T_Form >
  EquationLhs< T_Scalar, T_Degree, T_Form, Product::ScalarProduct > operator*(const function::Function< T_Scalar, T_Degree > &function, const UnknownField< T_Scalar, T_Form > &field);
  template< class T_Scalar, field::Form T_Form >
  EquationLhs< T_Scalar, Degree::Degree0, T_Form, Product::ScalarProduct > operator*(const typename MathObject< T_Scalar, Degree::Degree0 >::Object &value, const UnknownField< T_Scalar, T_Form > &field);
  template< class T_Scalar, field::Form T_Form >
  EquationLhs< T_Scalar, Degree::Degree1, T_Form, Product::ScalarProduct > operator*(const typename MathObject< T_Scalar, Degree::Degree1 >::Object &value, const UnknownField< T_Scalar, T_Form > &field);
  template< class T_Scalar, field::Form T_Form >
  EquationLhs< T_Scalar, Degree::Degree2, T_Form, Product::ScalarProduct > operator*(const typename MathObject< T_Scalar, Degree::Degree2 >::Object &value, const UnknownField< T_Scalar, T_Form > &field);

  template< class T_Scalar, Degree T_Degree, field::Form T_Form, unsigned int T_NumFields >
  EquationLhsCompound< T_Scalar, T_Degree, T_Form, Product::ScalarProduct, T_NumFields > operator*(const function::Function< T_Scalar, T_Degree > &function, const UnknownCompoundField< T_Scalar, T_Form, T_NumFields > &field);
  template< class T_Scalar, field::Form T_Form, unsigned int T_NumFields >
  EquationLhsCompound< T_Scalar, Degree::Degree0, T_Form, Product::ScalarProduct, T_NumFields > operator*(const typename MathObject< T_Scalar, Degree::Degree0 >::Object &value, const UnknownCompoundField< T_Scalar, T_Form, T_NumFields > &field);
  template< class T_Scalar, field::Form T_Form, unsigned int T_NumFields >
  EquationLhsCompound< T_Scalar, Degree::Degree1, T_Form, Product::ScalarProduct, T_NumFields > operator*(const typename MathObject< T_Scalar, Degree::Degree1 >::Object &value, const UnknownCompoundField< T_Scalar, T_Form, T_NumFields > &field);
  template< class T_Scalar, field::Form T_Form, unsigned int T_NumFields >
  EquationLhsCompound< T_Scalar, Degree::Degree2, T_Form, Product::ScalarProduct, T_NumFields > operator*(const typename MathObject< T_Scalar, Degree::Degree2 >::Object &value, const UnknownCompoundField< T_Scalar, T_Form, T_NumFields > &field);
  template< class T_Scalar, field::Form T_Form, unsigned int T_NumFields >
  EquationLhsCompound< T_Scalar, Degree::Degree4, T_Form, Product::ScalarProduct, T_NumFields > operator*(const typename MathObject< T_Scalar, Degree::Degree4 >::Object &value, const UnknownCompoundField< T_Scalar, T_Form, T_NumFields > &field);

  template< class T_Scalar, field::Form T_Form >
  EquationLhs< T_Scalar, Degree::Degree1, T_Form, Product::VectorProduct > operator%(const function::Function< T_Scalar, Degree::Degree1 > &function, const UnknownField< T_Scalar, T_Form > &field);
  template< class T_Scalar, field::Form T_Form >
  EquationLhs< T_Scalar, Degree::Degree1, T_Form, Product::VectorProduct > operator%(const typename MathObject< T_Scalar, Degree::Degree1 >::Object &value, const UnknownField< T_Scalar, T_Form > &field);

  // Rhs
  template< class T_Scalar, Degree T_Degree, field::Form T_Form >
  EquationRhs< T_Scalar, T_Degree, T_Form, Product::ScalarProduct > operator*(const UnknownField< T_Scalar, T_Form > &field, const function::Function< T_Scalar, T_Degree > &function);
  template< class T_Scalar, field::Form T_Form >
  EquationRhs< T_Scalar, Degree::Degree0, T_Form, Product::ScalarProduct > operator*(const UnknownField< T_Scalar, T_Form > &field, const typename MathObject< T_Scalar, Degree::Degree0 >::Object &value);
  template< class T_Scalar, field::Form T_Form >
  EquationRhs< T_Scalar, Degree::Degree1, T_Form, Product::ScalarProduct > operator*(const UnknownField< T_Scalar, T_Form > &field, const typename MathObject< T_Scalar, Degree::Degree1 >::Object &value);
  template< class T_Scalar, field::Form T_Form >
  EquationRhs< T_Scalar, Degree::Degree2, T_Form, Product::ScalarProduct > operator*(const UnknownField< T_Scalar, T_Form > &field, const typename MathObject< T_Scalar, Degree::Degree2 >::Object &value);

  template< class T_Scalar, Degree T_Degree, field::Form T_Form, unsigned int T_NumFields >
  EquationRhsCompound< T_Scalar, T_Degree, T_Form, Product::ScalarProduct, T_NumFields > operator*(const UnknownCompoundField< T_Scalar, T_Form, T_NumFields > &field, const function::Function< T_Scalar, T_Degree > &function);
  template< class T_Scalar, field::Form T_Form, unsigned int T_NumFields >
  EquationRhsCompound< T_Scalar, Degree::Degree0, T_Form, Product::ScalarProduct, T_NumFields > operator*(const UnknownCompoundField< T_Scalar, T_Form, T_NumFields > &field, const typename MathObject< T_Scalar, Degree::Degree0 >::Object &value);
  template< class T_Scalar, field::Form T_Form, unsigned int T_NumFields >
  EquationRhsCompound< T_Scalar, Degree::Degree1, T_Form, Product::ScalarProduct, T_NumFields > operator*(const UnknownCompoundField< T_Scalar, T_Form, T_NumFields > &field, const typename MathObject< T_Scalar, Degree::Degree1 >::Object &value);
  template< class T_Scalar, field::Form T_Form, unsigned int T_NumFields >
  EquationRhsCompound< T_Scalar, Degree::Degree2, T_Form, Product::ScalarProduct, T_NumFields > operator*(const UnknownCompoundField< T_Scalar, T_Form, T_NumFields > &field, const typename MathObject< T_Scalar, Degree::Degree2 >::Object &value);

  template< class T_Scalar, field::Form T_Form >
  EquationRhs< T_Scalar, Degree::Degree1, T_Form, Product::VectorProduct > operator%(const UnknownField< T_Scalar, T_Form > &field, const function::Function< T_Scalar, Degree::Degree1 > &function);
  template< class T_Scalar, field::Form T_Form >
  EquationRhs< T_Scalar, Degree::Degree1, T_Form, Product::VectorProduct > operator%(const UnknownField< T_Scalar, T_Form > &field, const typename MathObject< T_Scalar, Degree::Degree1 >::Object &value);

  // LhsRhs
  template< class T_Scalar, field::Form T_Form >
  EquationLhsRhs< T_Scalar, Degree::Degree1, Degree::Degree1, T_Form, Product::VectorProduct, Product::VectorProduct, term::evaluator::EvaluationOrder::Right > operator%(const function::Function< T_Scalar, Degree::Degree1 > &function, const EquationRhs< T_Scalar, Degree::Degree1, T_Form, Product::VectorProduct > &equation);


} // namespace gmshfem::equation

#endif // H_GMSHFEM_EQUATION
