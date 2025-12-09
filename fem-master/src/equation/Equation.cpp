// GmshFEM - Copyright (C) 2019-2022, A. Royer, E. Béchet, C. Geuzaine, Université de Liège
//
// See the LICENSE.txt file for license information. Please report all
// issues on https://gitlab.onelab.info/gmsh/fem/issues

#include "Equation.h"

#include "UnknownField.h"
#include "instantiate.h"

namespace gmshfem::equation
{


  // *************************************************
  // EquationLhs< T_Scalar, T_Degree, T_Form, T_Product >
  // *************************************************

  template< class T_Scalar, Degree T_Degree, field::Form T_Form, Product T_Product >
  EquationLhs< T_Scalar, T_Degree, T_Form, T_Product >::EquationLhs(const EquationLhs &other) :
    _left(other._left), _field(other._field->copy())
  {
  }

  template< class T_Scalar, Degree T_Degree, field::Form T_Form, Product T_Product >
  EquationLhs< T_Scalar, T_Degree, T_Form, T_Product >::EquationLhs(EquationLhs &&other) :
    _left(std::move(other._left)), _field(other._field->copy())
  {
  }

  template< class T_Scalar, Degree T_Degree, field::Form T_Form, Product T_Product >
  EquationLhs< T_Scalar, T_Degree, T_Form, T_Product >::~EquationLhs()
  {
    delete _field;
  }

  template< class T_Scalar, Degree T_Degree, field::Form T_Form, Product T_Product >
  EquationLhs< T_Scalar, T_Degree, T_Form, T_Product >::EquationLhs(const function::Function< T_Scalar, T_Degree > &left, const UnknownField< T_Scalar, T_Form > &field) :
    _left(left), _field(&field)
  {
  }

  template< class T_Scalar, Degree T_Degree, field::Form T_Form, Product T_Product >
  field::Form EquationLhs< T_Scalar, T_Degree, T_Form, T_Product >::form() const
  {
    return T_Form;
  }

  template< class T_Scalar, Degree T_Degree, field::Form T_Form, Product T_Product >
  field::FieldInterface< T_Scalar > *EquationLhs< T_Scalar, T_Degree, T_Form, T_Product >::getField() const
  {
    return _field->getField();
  }

  template< class T_Scalar, Degree T_Degree, field::Form T_Form, Product T_Product >
  term::evaluator::EquationEvaluator< T_Scalar, T_Form, T_Degree, T_Product, Degree::Empty, Product::Empty, term::evaluator::EvaluationOrder::Left > *EquationLhs< T_Scalar, T_Degree, T_Form, T_Product >::getEquationEvaluator() const
  {
    return new term::evaluator::EquationEvaluator< T_Scalar, T_Form, T_Degree, T_Product, Degree::Empty, Product::Empty, term::evaluator::EvaluationOrder::Left >(_left);
  }

  template< class T_Scalar, Degree T_Degree, field::Form T_Form, Product T_Product >
  term::evaluator::FieldEvaluator< T_Scalar, T_Form > *EquationLhs< T_Scalar, T_Degree, T_Form, T_Product >::getFieldEvaluator() const
  {
    return _field->getFieldEvaluator();
  }

  template< class T_Scalar, Degree T_Degree, field::Form T_Form, Product T_Product >
  const UnknownFieldInterface< T_Scalar > *EquationLhs< T_Scalar, T_Degree, T_Form, T_Product >::getUnknownField() const
  {
    return _field;
  }

  template< class T_Scalar, Degree T_Degree, field::Form T_Form, Product T_Product >
  bool EquationLhs< T_Scalar, T_Degree, T_Form, T_Product >::isDerivative() const
  {
    return _field->isDerivative();
  }

  template< class T_Scalar, Degree T_Degree, field::Form T_Form, Product T_Product >
  EquationLhs< T_Scalar, T_Degree, T_Form, T_Product > *EquationLhs< T_Scalar, T_Degree, T_Form, T_Product >::copy() const
  {
    return new EquationLhs< T_Scalar, T_Degree, T_Form, T_Product >(*this);
  }

  template< class T_Scalar, Degree T_Degree, field::Form T_Form, Product T_Product >
  EquationLhs< T_Scalar, T_Degree, T_Form, T_Product > *EquationLhs< T_Scalar, T_Degree, T_Form, T_Product >::changeUnknownField(field::FieldInterface< T_Scalar > *dual) const
  {
    EquationLhs< T_Scalar, T_Degree, T_Form, T_Product > *other = copy();
    UnknownField< T_Scalar, T_Form > *newUnknownField = other->_field->copy();
    newUnknownField->changeField(dual);
    delete other->_field;
    other->_field = newUnknownField;
    return other;
  }

  template< class T_Scalar, Degree T_Degree, field::Form T_Form, Product T_Product >
  function::Function< T_Scalar, T_Degree > EquationLhs< T_Scalar, T_Degree, T_Form, T_Product >::getFunction() const
  {
    return _left;
  }

  INSTANTIATE_CLASS_4(EquationLhs, 4, 3, 4, 1, TEMPLATE_ARGS(std::complex< double >, std::complex< float >, double, float), TEMPLATE_ARGS(Degree::Degree0, Degree::Degree1, Degree::Degree2), TEMPLATE_ARGS(field::Form::Form0, field::Form::Form1, field::Form::Form2, field::Form::Form3), TEMPLATE_ARGS(Product::ScalarProduct))
  INSTANTIATE_CLASS_4(EquationLhs, 4, 1, 2, 1, TEMPLATE_ARGS(std::complex< double >, std::complex< float >, double, float), TEMPLATE_ARGS(Degree::Degree1), TEMPLATE_ARGS(field::Form::Form1, field::Form::Form2), TEMPLATE_ARGS(Product::VectorProduct))

  // *************************************************
  // EquationLhsCompound< T_Scalar, T_Degree, T_Form, T_Product, T_NumFields >
  // *************************************************

  template< class T_Scalar, Degree T_Degree, field::Form T_Form, Product T_Product, unsigned int T_NumFields >
  EquationLhsCompound< T_Scalar, T_Degree, T_Form, T_Product, T_NumFields >::EquationLhsCompound(const EquationLhsCompound &other) :
    _left(other._left), _field(other._field->copy())
  {
  }

  template< class T_Scalar, Degree T_Degree, field::Form T_Form, Product T_Product, unsigned int T_NumFields >
  EquationLhsCompound< T_Scalar, T_Degree, T_Form, T_Product, T_NumFields >::EquationLhsCompound(EquationLhsCompound &&other) :
    _left(std::move(other._left)), _field(other._field->copy())
  {
  }

  template< class T_Scalar, Degree T_Degree, field::Form T_Form, Product T_Product, unsigned int T_NumFields >
  EquationLhsCompound< T_Scalar, T_Degree, T_Form, T_Product, T_NumFields >::~EquationLhsCompound()
  {
    delete _field;
  }

  template< class T_Scalar, Degree T_Degree, field::Form T_Form, Product T_Product, unsigned int T_NumFields >
  EquationLhsCompound< T_Scalar, T_Degree, T_Form, T_Product, T_NumFields >::EquationLhsCompound(const function::Function< T_Scalar, T_Degree > &left, const UnknownCompoundField< T_Scalar, T_Form, T_NumFields > &field) :
    _left(left), _field(&field)
  {
  }

  template< class T_Scalar, Degree T_Degree, field::Form T_Form, Product T_Product, unsigned int T_NumFields >
  field::Form EquationLhsCompound< T_Scalar, T_Degree, T_Form, T_Product, T_NumFields >::form() const
  {
    return T_Form;
  }

  template< class T_Scalar, Degree T_Degree, field::Form T_Form, Product T_Product, unsigned int T_NumFields >
  field::FieldInterface< T_Scalar > *EquationLhsCompound< T_Scalar, T_Degree, T_Form, T_Product, T_NumFields >::getField() const
  {
    return _field->getField();
  }

  template< class T_Scalar, Degree T_Degree, field::Form T_Form, Product T_Product, unsigned int T_NumFields >
  term::evaluator::EquationEvaluatorCompound< T_Scalar, T_Form, T_Degree, T_Product, Degree::Empty, Product::Empty, term::evaluator::EvaluationOrder::Left > *EquationLhsCompound< T_Scalar, T_Degree, T_Form, T_Product, T_NumFields >::getEquationEvaluator() const
  {
    return new term::evaluator::EquationEvaluatorCompound< T_Scalar, T_Form, T_Degree, T_Product, Degree::Empty, Product::Empty, term::evaluator::EvaluationOrder::Left >(_left);
  }

  template< class T_Scalar, Degree T_Degree, field::Form T_Form, Product T_Product, unsigned int T_NumFields >
  term::evaluator::CompoundFieldEvaluator< T_Scalar, T_Form, T_NumFields > *EquationLhsCompound< T_Scalar, T_Degree, T_Form, T_Product, T_NumFields >::getFieldEvaluator() const
  {
    return _field->getFieldEvaluator();
  }

  template< class T_Scalar, Degree T_Degree, field::Form T_Form, Product T_Product, unsigned int T_NumFields >
  const UnknownFieldInterface< T_Scalar > *EquationLhsCompound< T_Scalar, T_Degree, T_Form, T_Product, T_NumFields >::getUnknownField() const
  {
    return _field;
  }

  template< class T_Scalar, Degree T_Degree, field::Form T_Form, Product T_Product, unsigned int T_NumFields >
  bool EquationLhsCompound< T_Scalar, T_Degree, T_Form, T_Product, T_NumFields >::isDerivative() const
  {
    return _field->isDerivative();
  }

  template< class T_Scalar, Degree T_Degree, field::Form T_Form, Product T_Product, unsigned int T_NumFields >
  EquationLhsCompound< T_Scalar, T_Degree, T_Form, T_Product, T_NumFields > *EquationLhsCompound< T_Scalar, T_Degree, T_Form, T_Product, T_NumFields >::copy() const
  {
    return new EquationLhsCompound< T_Scalar, T_Degree, T_Form, T_Product, T_NumFields >(*this);
  }

  template< class T_Scalar, Degree T_Degree, field::Form T_Form, Product T_Product, unsigned int T_NumFields >
  EquationLhsCompound< T_Scalar, T_Degree, T_Form, T_Product, T_NumFields > *EquationLhsCompound< T_Scalar, T_Degree, T_Form, T_Product, T_NumFields >::changeUnknownField(field::FieldInterface< T_Scalar > *dual) const
  {
    EquationLhsCompound< T_Scalar, T_Degree, T_Form, T_Product, T_NumFields > *other = copy();
    UnknownCompoundField< T_Scalar, T_Form, T_NumFields > *newUnknownField = other->_field->copy();
    newUnknownField->changeField(dual);
    delete other->_field;
    other->_field = newUnknownField;
    return other;
  }

  template< class T_Scalar, Degree T_Degree, field::Form T_Form, Product T_Product, unsigned int T_NumFields >
  function::Function< T_Scalar, T_Degree > EquationLhsCompound< T_Scalar, T_Degree, T_Form, T_Product, T_NumFields >::getFunction() const
  {
    return _left;
  }

  INSTANTIATE_CLASS_5(EquationLhsCompound, 4, 3, 4, 1, 2, TEMPLATE_ARGS(std::complex< double >, std::complex< float >, double, float), TEMPLATE_ARGS(Degree::Degree0, Degree::Degree1, Degree::Degree2), TEMPLATE_ARGS(field::Form::Form0, field::Form::Form1, field::Form::Form2, field::Form::Form3), TEMPLATE_ARGS(Product::ScalarProduct), TEMPLATE_ARGS(2, 3))
  INSTANTIATE_CLASS_5(EquationLhsCompound, 4, 1, 2, 1, 2, TEMPLATE_ARGS(std::complex< double >, std::complex< float >, double, float), TEMPLATE_ARGS(Degree::Degree4), TEMPLATE_ARGS(field::Form::Form1, field::Form::Form2), TEMPLATE_ARGS(Product::ScalarProduct), TEMPLATE_ARGS(2, 3))


  // *************************************************
  // EquationRhs< T_Scalar, T_Form, T_Product, T_Degree >
  // *************************************************

  template< class T_Scalar, Degree T_Degree, field::Form T_Form, Product T_Product >
  EquationRhs< T_Scalar, T_Degree, T_Form, T_Product >::EquationRhs(const EquationRhs &other) :
    _right(other._right), _field(other._field->copy())
  {
  }

  template< class T_Scalar, Degree T_Degree, field::Form T_Form, Product T_Product >
  EquationRhs< T_Scalar, T_Degree, T_Form, T_Product >::EquationRhs(EquationRhs &&other) :
    _right(std::move(other._right)), _field(other._field->copy())
  {
  }

  template< class T_Scalar, Degree T_Degree, field::Form T_Form, Product T_Product >
  EquationRhs< T_Scalar, T_Degree, T_Form, T_Product >::~EquationRhs()
  {
    delete _field;
  }

  template< class T_Scalar, Degree T_Degree, field::Form T_Form, Product T_Product >
  EquationRhs< T_Scalar, T_Degree, T_Form, T_Product >::EquationRhs(const function::Function< T_Scalar, T_Degree > &right, const UnknownField< T_Scalar, T_Form > &field) :
    _right(right), _field(&field)
  {
  }

  template< class T_Scalar, Degree T_Degree, field::Form T_Form, Product T_Product >
  field::Form EquationRhs< T_Scalar, T_Degree, T_Form, T_Product >::form() const
  {
    return T_Form;
  }

  template< class T_Scalar, Degree T_Degree, field::Form T_Form, Product T_Product >
  field::FieldInterface< T_Scalar > *EquationRhs< T_Scalar, T_Degree, T_Form, T_Product >::getField() const
  {
    return _field->getField();
  }

  template< class T_Scalar, Degree T_Degree, field::Form T_Form, Product T_Product >
  term::evaluator::EquationEvaluator< T_Scalar, T_Form, Degree::Empty, Product::Empty, T_Degree, T_Product, term::evaluator::EvaluationOrder::Right > *EquationRhs< T_Scalar, T_Degree, T_Form, T_Product >::getEquationEvaluator() const
  {
    return new term::evaluator::EquationEvaluator< T_Scalar, T_Form, Degree::Empty, Product::Empty, T_Degree, T_Product, term::evaluator::EvaluationOrder::Right >(_right);
  }

  template< class T_Scalar, Degree T_Degree, field::Form T_Form, Product T_Product >
  term::evaluator::FieldEvaluator< T_Scalar, T_Form > *EquationRhs< T_Scalar, T_Degree, T_Form, T_Product >::getFieldEvaluator() const
  {
    return _field->getFieldEvaluator();
  }

  template< class T_Scalar, Degree T_Degree, field::Form T_Form, Product T_Product >
  const UnknownFieldInterface< T_Scalar > *EquationRhs< T_Scalar, T_Degree, T_Form, T_Product >::getUnknownField() const
  {
    return _field;
  }

  template< class T_Scalar, Degree T_Degree, field::Form T_Form, Product T_Product >
  bool EquationRhs< T_Scalar, T_Degree, T_Form, T_Product >::isDerivative() const
  {
    return _field->isDerivative();
  }

  template< class T_Scalar, Degree T_Degree, field::Form T_Form, Product T_Product >
  EquationRhs< T_Scalar, T_Degree, T_Form, T_Product > *EquationRhs< T_Scalar, T_Degree, T_Form, T_Product >::copy() const
  {
    return new EquationRhs< T_Scalar, T_Degree, T_Form, T_Product >(*this);
  }

  template< class T_Scalar, Degree T_Degree, field::Form T_Form, Product T_Product >
  EquationRhs< T_Scalar, T_Degree, T_Form, T_Product > *EquationRhs< T_Scalar, T_Degree, T_Form, T_Product >::changeUnknownField(field::FieldInterface< T_Scalar > *dual) const
  {
    EquationRhs< T_Scalar, T_Degree, T_Form, T_Product > *other = copy();
    UnknownField< T_Scalar, T_Form > *newUnknownField = other->_field->copy();
    newUnknownField->changeField(dual);
    delete other->_field;
    other->_field = newUnknownField;
    return other;
  }

  template< class T_Scalar, Degree T_Degree, field::Form T_Form, Product T_Product >
  function::Function< T_Scalar, T_Degree > EquationRhs< T_Scalar, T_Degree, T_Form, T_Product >::getFunction() const
  {
    return _right;
  }

  INSTANTIATE_CLASS_4(EquationRhs, 4, 3, 4, 1, TEMPLATE_ARGS(std::complex< double >, std::complex< float >, double, float), TEMPLATE_ARGS(Degree::Degree0, Degree::Degree1, Degree::Degree2), TEMPLATE_ARGS(field::Form::Form0, field::Form::Form1, field::Form::Form2, field::Form::Form3), TEMPLATE_ARGS(Product::ScalarProduct))
  INSTANTIATE_CLASS_4(EquationRhs, 4, 1, 2, 1, TEMPLATE_ARGS(std::complex< double >, std::complex< float >, double, float), TEMPLATE_ARGS(Degree::Degree1), TEMPLATE_ARGS(field::Form::Form1, field::Form::Form2), TEMPLATE_ARGS(Product::VectorProduct))


  // *************************************************
  // EquationRhsCompound< T_Scalar, T_Degree, T_Form, T_Product, T_NumFields >
  // *************************************************

  template< class T_Scalar, Degree T_Degree, field::Form T_Form, Product T_Product, unsigned int T_NumFields >
  EquationRhsCompound< T_Scalar, T_Degree, T_Form, T_Product, T_NumFields >::EquationRhsCompound(const EquationRhsCompound &other) :
    _right(other._right), _field(other._field->copy())
  {
  }

  template< class T_Scalar, Degree T_Degree, field::Form T_Form, Product T_Product, unsigned int T_NumFields >
  EquationRhsCompound< T_Scalar, T_Degree, T_Form, T_Product, T_NumFields >::EquationRhsCompound(EquationRhsCompound &&other) :
    _right(std::move(other._right)), _field(other._field->copy())
  {
  }

  template< class T_Scalar, Degree T_Degree, field::Form T_Form, Product T_Product, unsigned int T_NumFields >
  EquationRhsCompound< T_Scalar, T_Degree, T_Form, T_Product, T_NumFields >::~EquationRhsCompound()
  {
    delete _field;
  }

  template< class T_Scalar, Degree T_Degree, field::Form T_Form, Product T_Product, unsigned int T_NumFields >
  EquationRhsCompound< T_Scalar, T_Degree, T_Form, T_Product, T_NumFields >::EquationRhsCompound(const function::Function< T_Scalar, T_Degree > &left, const UnknownCompoundField< T_Scalar, T_Form, T_NumFields > &field) :
    _right(left), _field(&field)
  {
  }

  template< class T_Scalar, Degree T_Degree, field::Form T_Form, Product T_Product, unsigned int T_NumFields >
  field::Form EquationRhsCompound< T_Scalar, T_Degree, T_Form, T_Product, T_NumFields >::form() const
  {
    return T_Form;
  }

  template< class T_Scalar, Degree T_Degree, field::Form T_Form, Product T_Product, unsigned int T_NumFields >
  field::FieldInterface< T_Scalar > *EquationRhsCompound< T_Scalar, T_Degree, T_Form, T_Product, T_NumFields >::getField() const
  {
    return _field->getField();
  }

  template< class T_Scalar, Degree T_Degree, field::Form T_Form, Product T_Product, unsigned int T_NumFields >
  term::evaluator::EquationEvaluatorCompound< T_Scalar, T_Form, Degree::Empty, Product::Empty, T_Degree, T_Product, term::evaluator::EvaluationOrder::Right > *EquationRhsCompound< T_Scalar, T_Degree, T_Form, T_Product, T_NumFields >::getEquationEvaluator() const
  {
    return new term::evaluator::EquationEvaluatorCompound< T_Scalar, T_Form, Degree::Empty, Product::Empty, T_Degree, T_Product, term::evaluator::EvaluationOrder::Right >(_right);
  }

  template< class T_Scalar, Degree T_Degree, field::Form T_Form, Product T_Product, unsigned int T_NumFields >
  term::evaluator::CompoundFieldEvaluator< T_Scalar, T_Form, T_NumFields > *EquationRhsCompound< T_Scalar, T_Degree, T_Form, T_Product, T_NumFields >::getFieldEvaluator() const
  {
    return _field->getFieldEvaluator();
  }

  template< class T_Scalar, Degree T_Degree, field::Form T_Form, Product T_Product, unsigned int T_NumFields >
  const UnknownFieldInterface< T_Scalar > *EquationRhsCompound< T_Scalar, T_Degree, T_Form, T_Product, T_NumFields >::getUnknownField() const
  {
    return _field;
  }

  template< class T_Scalar, Degree T_Degree, field::Form T_Form, Product T_Product, unsigned int T_NumFields >
  bool EquationRhsCompound< T_Scalar, T_Degree, T_Form, T_Product, T_NumFields >::isDerivative() const
  {
    return _field->isDerivative();
  }

  template< class T_Scalar, Degree T_Degree, field::Form T_Form, Product T_Product, unsigned int T_NumFields >
  EquationRhsCompound< T_Scalar, T_Degree, T_Form, T_Product, T_NumFields > *EquationRhsCompound< T_Scalar, T_Degree, T_Form, T_Product, T_NumFields >::copy() const
  {
    return new EquationRhsCompound< T_Scalar, T_Degree, T_Form, T_Product, T_NumFields >(*this);
  }

  template< class T_Scalar, Degree T_Degree, field::Form T_Form, Product T_Product, unsigned int T_NumFields >
  EquationRhsCompound< T_Scalar, T_Degree, T_Form, T_Product, T_NumFields > *EquationRhsCompound< T_Scalar, T_Degree, T_Form, T_Product, T_NumFields >::changeUnknownField(field::FieldInterface< T_Scalar > *dual) const
  {
    EquationRhsCompound< T_Scalar, T_Degree, T_Form, T_Product, T_NumFields > *other = copy();
    UnknownCompoundField< T_Scalar, T_Form, T_NumFields > *newUnknownField = other->_field->copy();
    newUnknownField->changeField(dual);
    delete other->_field;
    other->_field = newUnknownField;
    return other;
  }

  template< class T_Scalar, Degree T_Degree, field::Form T_Form, Product T_Product, unsigned int T_NumFields >
  function::Function< T_Scalar, T_Degree > EquationRhsCompound< T_Scalar, T_Degree, T_Form, T_Product, T_NumFields >::getFunction() const
  {
    return _right;
  }

  INSTANTIATE_CLASS_5(EquationRhsCompound, 4, 3, 4, 1, 2, TEMPLATE_ARGS(std::complex< double >, std::complex< float >, double, float), TEMPLATE_ARGS(Degree::Degree0, Degree::Degree1, Degree::Degree2), TEMPLATE_ARGS(field::Form::Form0, field::Form::Form1, field::Form::Form2, field::Form::Form3), TEMPLATE_ARGS(Product::ScalarProduct), TEMPLATE_ARGS(2, 3))


  // *************************************************
  // EquationRhsLhs< T_Scalar, T_DegreeLhs, T_DegreeRhs, T_Form, T_ProductLhs, T_DegreeRhs, EvaluationOrder T_EvaluationOrder >
  // *************************************************

  // Left

  template< class T_Scalar, Degree T_DegreeLhs, Degree T_DegreeRhs, field::Form T_Form, Product T_ProductLhs, Product T_ProductRhs >
  EquationLhsRhs< T_Scalar, T_DegreeLhs, T_DegreeRhs, T_Form, T_ProductLhs, T_ProductRhs, term::evaluator::EvaluationOrder::Left >::EquationLhsRhs(const EquationLhsRhs &other) :
    _left(other._left), _field(other._field->copy()), _right(other._right)
  {
  }

  template< class T_Scalar, Degree T_DegreeLhs, Degree T_DegreeRhs, field::Form T_Form, Product T_ProductLhs, Product T_ProductRhs >
  EquationLhsRhs< T_Scalar, T_DegreeLhs, T_DegreeRhs, T_Form, T_ProductLhs, T_ProductRhs, term::evaluator::EvaluationOrder::Left >::EquationLhsRhs(EquationLhsRhs &&other) :
    _left(std::move(other._left)), _field(other._field->copy()), _right(std::move(other._right))
  {
  }

  template< class T_Scalar, Degree T_DegreeLhs, Degree T_DegreeRhs, field::Form T_Form, Product T_ProductLhs, Product T_ProductRhs >
  EquationLhsRhs< T_Scalar, T_DegreeLhs, T_DegreeRhs, T_Form, T_ProductLhs, T_ProductRhs, term::evaluator::EvaluationOrder::Left >::~EquationLhsRhs()
  {
    delete _field;
  }

  template< class T_Scalar, Degree T_DegreeLhs, Degree T_DegreeRhs, field::Form T_Form, Product T_ProductLhs, Product T_ProductRhs >
  EquationLhsRhs< T_Scalar, T_DegreeLhs, T_DegreeRhs, T_Form, T_ProductLhs, T_ProductRhs, term::evaluator::EvaluationOrder::Left >::EquationLhsRhs(const function::Function< T_Scalar, T_DegreeLhs > &left, const UnknownField< T_Scalar, T_Form > &field, const function::Function< T_Scalar, T_DegreeRhs > &right) :
    _left(left), _field(&field), _right(right)
  {
  }

  template< class T_Scalar, Degree T_DegreeLhs, Degree T_DegreeRhs, field::Form T_Form, Product T_ProductLhs, Product T_ProductRhs >
  field::Form EquationLhsRhs< T_Scalar, T_DegreeLhs, T_DegreeRhs, T_Form, T_ProductLhs, T_ProductRhs, term::evaluator::EvaluationOrder::Left >::form() const
  {
    return T_Form;
  }

  template< class T_Scalar, Degree T_DegreeLhs, Degree T_DegreeRhs, field::Form T_Form, Product T_ProductLhs, Product T_ProductRhs >
  field::FieldInterface< T_Scalar > *EquationLhsRhs< T_Scalar, T_DegreeLhs, T_DegreeRhs, T_Form, T_ProductLhs, T_ProductRhs, term::evaluator::EvaluationOrder::Left >::getField() const
  {
    return _field->getField();
  }

  template< class T_Scalar, Degree T_DegreeLhs, Degree T_DegreeRhs, field::Form T_Form, Product T_ProductLhs, Product T_ProductRhs >
  term::evaluator::EquationEvaluator< T_Scalar, T_Form, T_DegreeLhs, T_ProductLhs, T_DegreeRhs, T_ProductRhs, term::evaluator::EvaluationOrder::Left > *EquationLhsRhs< T_Scalar, T_DegreeLhs, T_DegreeRhs, T_Form, T_ProductLhs, T_ProductRhs, term::evaluator::EvaluationOrder::Left >::getEquationEvaluator() const
  {
    return new term::evaluator::EquationEvaluator< T_Scalar, T_Form, T_DegreeLhs, T_ProductLhs, T_DegreeRhs, T_ProductRhs, term::evaluator::EvaluationOrder::Left >(_left, _right);
  }

  template< class T_Scalar, Degree T_DegreeLhs, Degree T_DegreeRhs, field::Form T_Form, Product T_ProductLhs, Product T_ProductRhs >
  term::evaluator::FieldEvaluator< T_Scalar, T_Form > *EquationLhsRhs< T_Scalar, T_DegreeLhs, T_DegreeRhs, T_Form, T_ProductLhs, T_ProductRhs, term::evaluator::EvaluationOrder::Left >::getFieldEvaluator() const
  {
    return _field->getFieldEvaluator();
  }

  template< class T_Scalar, Degree T_DegreeLhs, Degree T_DegreeRhs, field::Form T_Form, Product T_ProductLhs, Product T_ProductRhs >
  const UnknownFieldInterface< T_Scalar > *EquationLhsRhs< T_Scalar, T_DegreeLhs, T_DegreeRhs, T_Form, T_ProductLhs, T_ProductRhs, term::evaluator::EvaluationOrder::Left >::getUnknownField() const
  {
    return _field;
  }

  template< class T_Scalar, Degree T_DegreeLhs, Degree T_DegreeRhs, field::Form T_Form, Product T_ProductLhs, Product T_ProductRhs >
  bool EquationLhsRhs< T_Scalar, T_DegreeLhs, T_DegreeRhs, T_Form, T_ProductLhs, T_ProductRhs, term::evaluator::EvaluationOrder::Left >::isDerivative() const
  {
    return _field->isDerivative();
  }

  template< class T_Scalar, Degree T_DegreeLhs, Degree T_DegreeRhs, field::Form T_Form, Product T_ProductLhs, Product T_ProductRhs >
  EquationLhsRhs< T_Scalar, T_DegreeLhs, T_DegreeRhs, T_Form, T_ProductLhs, T_ProductRhs, term::evaluator::EvaluationOrder::Left > *EquationLhsRhs< T_Scalar, T_DegreeLhs, T_DegreeRhs, T_Form, T_ProductLhs, T_ProductRhs, term::evaluator::EvaluationOrder::Left >::copy() const
  {
    return new EquationLhsRhs< T_Scalar, T_DegreeLhs, T_DegreeRhs, T_Form, T_ProductLhs, T_ProductRhs, term::evaluator::EvaluationOrder::Left >(*this);
  }

  template< class T_Scalar, Degree T_DegreeLhs, Degree T_DegreeRhs, field::Form T_Form, Product T_ProductLhs, Product T_ProductRhs >
  EquationLhsRhs< T_Scalar, T_DegreeLhs, T_DegreeRhs, T_Form, T_ProductLhs, T_ProductRhs, term::evaluator::EvaluationOrder::Left > *EquationLhsRhs< T_Scalar, T_DegreeLhs, T_DegreeRhs, T_Form, T_ProductLhs, T_ProductRhs, term::evaluator::EvaluationOrder::Left >::changeUnknownField(field::FieldInterface< T_Scalar > *dual) const
  {
    EquationLhsRhs< T_Scalar, T_DegreeLhs, T_DegreeRhs, T_Form, T_ProductLhs, T_ProductRhs, term::evaluator::EvaluationOrder::Left > other = copy();
    UnknownField< T_Scalar, T_Form > *newUnknownField = other->_field->copy();
    newUnknownField->changeField(dual);
    delete other->_field;
    other->_field = newUnknownField;
    return other;
  }

  // Right

  template< class T_Scalar, Degree T_DegreeLhs, Degree T_DegreeRhs, field::Form T_Form, Product T_ProductLhs, Product T_ProductRhs >
  EquationLhsRhs< T_Scalar, T_DegreeLhs, T_DegreeRhs, T_Form, T_ProductLhs, T_ProductRhs, term::evaluator::EvaluationOrder::Right >::EquationLhsRhs(const EquationLhsRhs &other) :
    _left(other._left), _field(other._field->copy()), _right(other._right)
  {
  }

  template< class T_Scalar, Degree T_DegreeLhs, Degree T_DegreeRhs, field::Form T_Form, Product T_ProductLhs, Product T_ProductRhs >
  EquationLhsRhs< T_Scalar, T_DegreeLhs, T_DegreeRhs, T_Form, T_ProductLhs, T_ProductRhs, term::evaluator::EvaluationOrder::Right >::EquationLhsRhs(EquationLhsRhs &&other) :
    _left(std::move(other._left)), _field(other._field->copy()), _right(std::move(other._right))
  {
  }

  template< class T_Scalar, Degree T_DegreeLhs, Degree T_DegreeRhs, field::Form T_Form, Product T_ProductLhs, Product T_ProductRhs >
  EquationLhsRhs< T_Scalar, T_DegreeLhs, T_DegreeRhs, T_Form, T_ProductLhs, T_ProductRhs, term::evaluator::EvaluationOrder::Right >::~EquationLhsRhs()
  {
    delete _field;
  }

  template< class T_Scalar, Degree T_DegreeLhs, Degree T_DegreeRhs, field::Form T_Form, Product T_ProductLhs, Product T_ProductRhs >
  EquationLhsRhs< T_Scalar, T_DegreeLhs, T_DegreeRhs, T_Form, T_ProductLhs, T_ProductRhs, term::evaluator::EvaluationOrder::Right >::EquationLhsRhs(const function::Function< T_Scalar, T_DegreeLhs > &left, const UnknownField< T_Scalar, T_Form > &field, const function::Function< T_Scalar, T_DegreeRhs > &right) :
    _left(left), _field(&field), _right(right)
  {
  }

  template< class T_Scalar, Degree T_DegreeLhs, Degree T_DegreeRhs, field::Form T_Form, Product T_ProductLhs, Product T_ProductRhs >
  field::Form EquationLhsRhs< T_Scalar, T_DegreeLhs, T_DegreeRhs, T_Form, T_ProductLhs, T_ProductRhs, term::evaluator::EvaluationOrder::Right >::form() const
  {
    return T_Form;
  }

  template< class T_Scalar, Degree T_DegreeLhs, Degree T_DegreeRhs, field::Form T_Form, Product T_ProductLhs, Product T_ProductRhs >
  field::FieldInterface< T_Scalar > *EquationLhsRhs< T_Scalar, T_DegreeLhs, T_DegreeRhs, T_Form, T_ProductLhs, T_ProductRhs, term::evaluator::EvaluationOrder::Right >::getField() const
  {
    return _field->getField();
  }

  template< class T_Scalar, Degree T_DegreeLhs, Degree T_DegreeRhs, field::Form T_Form, Product T_ProductLhs, Product T_ProductRhs >
  term::evaluator::EquationEvaluator< T_Scalar, T_Form, T_DegreeLhs, T_ProductLhs, T_DegreeRhs, T_ProductRhs, term::evaluator::EvaluationOrder::Right > *EquationLhsRhs< T_Scalar, T_DegreeLhs, T_DegreeRhs, T_Form, T_ProductLhs, T_ProductRhs, term::evaluator::EvaluationOrder::Right >::getEquationEvaluator() const
  {
    return new term::evaluator::EquationEvaluator< T_Scalar, T_Form, T_DegreeLhs, T_ProductLhs, T_DegreeRhs, T_ProductRhs, term::evaluator::EvaluationOrder::Right >(_left, _right);
  }

  template< class T_Scalar, Degree T_DegreeLhs, Degree T_DegreeRhs, field::Form T_Form, Product T_ProductLhs, Product T_ProductRhs >
  term::evaluator::FieldEvaluator< T_Scalar, T_Form > *EquationLhsRhs< T_Scalar, T_DegreeLhs, T_DegreeRhs, T_Form, T_ProductLhs, T_ProductRhs, term::evaluator::EvaluationOrder::Right >::getFieldEvaluator() const
  {
    return _field->getFieldEvaluator();
  }

  template< class T_Scalar, Degree T_DegreeLhs, Degree T_DegreeRhs, field::Form T_Form, Product T_ProductLhs, Product T_ProductRhs >
  const UnknownFieldInterface< T_Scalar > *EquationLhsRhs< T_Scalar, T_DegreeLhs, T_DegreeRhs, T_Form, T_ProductLhs, T_ProductRhs, term::evaluator::EvaluationOrder::Right >::getUnknownField() const
  {
    return _field;
  }

  template< class T_Scalar, Degree T_DegreeLhs, Degree T_DegreeRhs, field::Form T_Form, Product T_ProductLhs, Product T_ProductRhs >
  bool EquationLhsRhs< T_Scalar, T_DegreeLhs, T_DegreeRhs, T_Form, T_ProductLhs, T_ProductRhs, term::evaluator::EvaluationOrder::Right >::isDerivative() const
  {
    return _field->isDerivative();
  }

  template< class T_Scalar, Degree T_DegreeLhs, Degree T_DegreeRhs, field::Form T_Form, Product T_ProductLhs, Product T_ProductRhs >
  EquationLhsRhs< T_Scalar, T_DegreeLhs, T_DegreeRhs, T_Form, T_ProductLhs, T_ProductRhs, term::evaluator::EvaluationOrder::Right > *EquationLhsRhs< T_Scalar, T_DegreeLhs, T_DegreeRhs, T_Form, T_ProductLhs, T_ProductRhs, term::evaluator::EvaluationOrder::Right >::copy() const
  {
    return new EquationLhsRhs< T_Scalar, T_DegreeLhs, T_DegreeRhs, T_Form, T_ProductLhs, T_ProductRhs, term::evaluator::EvaluationOrder::Right >(*this);
  }

  template< class T_Scalar, Degree T_DegreeLhs, Degree T_DegreeRhs, field::Form T_Form, Product T_ProductLhs, Product T_ProductRhs >
  EquationLhsRhs< T_Scalar, T_DegreeLhs, T_DegreeRhs, T_Form, T_ProductLhs, T_ProductRhs, term::evaluator::EvaluationOrder::Right > *EquationLhsRhs< T_Scalar, T_DegreeLhs, T_DegreeRhs, T_Form, T_ProductLhs, T_ProductRhs, term::evaluator::EvaluationOrder::Right >::changeUnknownField(field::FieldInterface< T_Scalar > *dual) const
  {
    EquationLhsRhs< T_Scalar, T_DegreeLhs, T_DegreeRhs, T_Form, T_ProductLhs, T_ProductRhs, term::evaluator::EvaluationOrder::Right > *other = copy();
    UnknownField< T_Scalar, T_Form > *newUnknownField = other->_field->copy();
    newUnknownField->changeField(dual);
    delete other->_field;
    other->_field = newUnknownField;
    return other;
  }

  INSTANTIATE_CLASS_7(EquationLhsRhs, 4, 1, 1, 2, 1, 1, 1, TEMPLATE_ARGS(std::complex< double >, std::complex< float >, double, float), TEMPLATE_ARGS(Degree::Degree1), TEMPLATE_ARGS(Degree::Degree1), TEMPLATE_ARGS(field::Form::Form1, field::Form::Form2), TEMPLATE_ARGS(Product::VectorProduct), TEMPLATE_ARGS(Product::VectorProduct), TEMPLATE_ARGS(term::evaluator::EvaluationOrder::Right))


  template< class T_Scalar, field::Form T_Form >
  EquationLhs< T_Scalar, Degree::Degree0, T_Form, Product::ScalarProduct > operator+(const UnknownField< T_Scalar, T_Form > &field)
  {
    return EquationLhs< T_Scalar, Degree::Degree0, T_Form, Product::ScalarProduct >(field);
  }

  INSTANTIATE_OPP_2(TEMPLATE_RETURN(EquationLhs< TEMPLATE_PARAM_1, Degree::Degree0, TEMPLATE_PARAM_2, Product::ScalarProduct >), , 0, 4, class, TEMPLATE_ARGS(std::complex< double >, std::complex< float >, double, float), 4, field::Form, TEMPLATE_ARGS(field::Form::Form0, field::Form::Form1, field::Form::Form2, field::Form::Form3), TEMPLATE_PARAMS(const UnknownField< TEMPLATE_PARAM_1, TEMPLATE_PARAM_2 > &))


  template< class T_Scalar, field::Form T_Form >
  EquationLhs< T_Scalar, Degree::Degree0, T_Form, Product::ScalarProduct > operator-(const UnknownField< T_Scalar, T_Form > &field)
  {
    return EquationLhs< T_Scalar, Degree::Degree0, T_Form, Product::ScalarProduct >(function::ScalarFunction< T_Scalar >(-1.), field);
  }

  INSTANTIATE_OPM_2(TEMPLATE_RETURN(EquationLhs< TEMPLATE_PARAM_1, Degree::Degree0, TEMPLATE_PARAM_2, Product::ScalarProduct >), , 0, 4, class, TEMPLATE_ARGS(std::complex< double >, std::complex< float >, double, float), 4, field::Form, TEMPLATE_ARGS(field::Form::Form0, field::Form::Form1, field::Form::Form2, field::Form::Form3), TEMPLATE_PARAMS(const UnknownField< TEMPLATE_PARAM_1, TEMPLATE_PARAM_2 > &))


  template< class T_Scalar, field::Form T_Form, unsigned int T_NumFields >
  EquationLhsCompound< T_Scalar, Degree::Degree0, T_Form, Product::ScalarProduct, T_NumFields > operator+(const UnknownCompoundField< T_Scalar, T_Form, T_NumFields > &field)
  {
    return EquationLhsCompound< T_Scalar, Degree::Degree0, T_Form, Product::ScalarProduct, T_NumFields >(field);
  }

  INSTANTIATE_OPP_3(TEMPLATE_RETURN(EquationLhsCompound< TEMPLATE_PARAM_1, Degree::Degree0, TEMPLATE_PARAM_2, Product::ScalarProduct, TEMPLATE_PARAM_3 >), , 1, 4, class, TEMPLATE_ARGS(std::complex< double >, std::complex< float >, double, float), 4, field::Form, TEMPLATE_ARGS(field::Form::Form0, field::Form::Form1, field::Form::Form2, field::Form::Form3), 2, unsigned int, TEMPLATE_ARGS(2, 3), TEMPLATE_PARAMS(const UnknownCompoundField< TEMPLATE_PARAM_1, TEMPLATE_PARAM_2, TEMPLATE_PARAM_3 > &))


  template< class T_Scalar, field::Form T_Form, unsigned int T_NumFields >
  EquationLhsCompound< T_Scalar, Degree::Degree0, T_Form, Product::ScalarProduct, T_NumFields > operator-(const UnknownCompoundField< T_Scalar, T_Form, T_NumFields > &field)
  {
    return EquationLhsCompound< T_Scalar, Degree::Degree0, T_Form, Product::ScalarProduct, T_NumFields >(function::ScalarFunction< T_Scalar >(-1.), field);
  }

  INSTANTIATE_OPM_3(TEMPLATE_RETURN(EquationLhsCompound< TEMPLATE_PARAM_1, Degree::Degree0, TEMPLATE_PARAM_2, Product::ScalarProduct, TEMPLATE_PARAM_3 >), , 1, 4, class, TEMPLATE_ARGS(std::complex< double >, std::complex< float >, double, float), 4, field::Form, TEMPLATE_ARGS(field::Form::Form0, field::Form::Form1, field::Form::Form2, field::Form::Form3), 2, unsigned int, TEMPLATE_ARGS(2, 3), TEMPLATE_PARAMS(const UnknownCompoundField< TEMPLATE_PARAM_1, TEMPLATE_PARAM_2, TEMPLATE_PARAM_3 > &))


  // Lhs

  template< class T_Scalar, Degree T_Degree, field::Form T_Form >
  EquationLhs< T_Scalar, T_Degree, T_Form, Product::ScalarProduct > operator*(const function::Function< T_Scalar, T_Degree > &function, const UnknownField< T_Scalar, T_Form > &field)
  {
    return EquationLhs< T_Scalar, T_Degree, T_Form, Product::ScalarProduct >(function, field);
  }

  INSTANTIATE_OPT_3(TEMPLATE_RETURN(EquationLhs< TEMPLATE_PARAM_1, TEMPLATE_PARAM_2, TEMPLATE_PARAM_3, Product::ScalarProduct >), , 0, 4, class, TEMPLATE_ARGS(std::complex< double >, std::complex< float >, double, float), 3, Degree, TEMPLATE_ARGS(Degree::Degree0, Degree::Degree1, Degree::Degree2), 4, field::Form, TEMPLATE_ARGS(field::Form::Form0, field::Form::Form1, field::Form::Form2, field::Form::Form3), TEMPLATE_PARAMS(const function::Function< TEMPLATE_PARAM_1, TEMPLATE_PARAM_2 > &, const UnknownField< TEMPLATE_PARAM_1, TEMPLATE_PARAM_3 > &))


  template< class T_Scalar, field::Form T_Form >
  EquationLhs< T_Scalar, Degree::Degree0, T_Form, Product::ScalarProduct > operator*(const typename MathObject< T_Scalar, Degree::Degree0 >::Object &value, const UnknownField< T_Scalar, T_Form > &field)
  {
    return EquationLhs< T_Scalar, Degree::Degree0, T_Form, Product::ScalarProduct >(function::ScalarFunction< T_Scalar >(value), field);
  }

  INSTANTIATE_OPT_2(TEMPLATE_RETURN(EquationLhs< TEMPLATE_PARAM_1, Degree::Degree0, TEMPLATE_PARAM_2, Product::ScalarProduct >), , 1, 4, class, TEMPLATE_ARGS(std::complex< double >, std::complex< float >, double, float), 4, field::Form, TEMPLATE_ARGS(field::Form::Form0, field::Form::Form1, field::Form::Form2, field::Form::Form3), TEMPLATE_PARAMS(const typename MathObject< TEMPLATE_PARAM_1, Degree::Degree0 >::Object &, const UnknownField< TEMPLATE_PARAM_1, TEMPLATE_PARAM_2 > &))


  template< class T_Scalar, field::Form T_Form >
  EquationLhs< T_Scalar, Degree::Degree1, T_Form, Product::ScalarProduct > operator*(const typename MathObject< T_Scalar, Degree::Degree1 >::Object &value, const UnknownField< T_Scalar, T_Form > &field)
  {
    return EquationLhs< T_Scalar, Degree::Degree1, T_Form, Product::ScalarProduct >(function::VectorFunction< T_Scalar >(value), field);
  }

  INSTANTIATE_OPT_2(TEMPLATE_RETURN(EquationLhs< TEMPLATE_PARAM_1, Degree::Degree1, TEMPLATE_PARAM_2, Product::ScalarProduct >), , 2, 4, class, TEMPLATE_ARGS(std::complex< double >, std::complex< float >, double, float), 4, field::Form, TEMPLATE_ARGS(field::Form::Form0, field::Form::Form1, field::Form::Form2, field::Form::Form3), TEMPLATE_PARAMS(const typename MathObject< TEMPLATE_PARAM_1, Degree::Degree1 >::Object &, const UnknownField< TEMPLATE_PARAM_1, TEMPLATE_PARAM_2 > &))


  template< class T_Scalar, field::Form T_Form >
  EquationLhs< T_Scalar, Degree::Degree2, T_Form, Product::ScalarProduct > operator*(const typename MathObject< T_Scalar, Degree::Degree2 >::Object &value, const UnknownField< T_Scalar, T_Form > &field)
  {
    return EquationLhs< T_Scalar, Degree::Degree2, T_Form, Product::ScalarProduct >(function::TensorFunction< T_Scalar >(value), field);
  }

  INSTANTIATE_OPT_2(TEMPLATE_RETURN(EquationLhs< TEMPLATE_PARAM_1, Degree::Degree2, TEMPLATE_PARAM_2, Product::ScalarProduct >), , 3, 4, class, TEMPLATE_ARGS(std::complex< double >, std::complex< float >, double, float), 4, field::Form, TEMPLATE_ARGS(field::Form::Form0, field::Form::Form1, field::Form::Form2, field::Form::Form3), TEMPLATE_PARAMS(const typename MathObject< TEMPLATE_PARAM_1, Degree::Degree2 >::Object &, const UnknownField< TEMPLATE_PARAM_1, TEMPLATE_PARAM_2 > &))


  template< class T_Scalar, Degree T_Degree, field::Form T_Form, unsigned int T_NumFields >
  EquationLhsCompound< T_Scalar, T_Degree, T_Form, Product::ScalarProduct, T_NumFields > operator*(const function::Function< T_Scalar, T_Degree > &function, const UnknownCompoundField< T_Scalar, T_Form, T_NumFields > &field)
  {
    return EquationLhsCompound< T_Scalar, T_Degree, T_Form, Product::ScalarProduct, T_NumFields >(function, field);
  }

  INSTANTIATE_OPT_4(TEMPLATE_RETURN(EquationLhsCompound< TEMPLATE_PARAM_1, TEMPLATE_PARAM_2, TEMPLATE_PARAM_3, Product::ScalarProduct, TEMPLATE_PARAM_4 >), , 4, 4, class, TEMPLATE_ARGS(std::complex< double >, std::complex< float >, double, float), 3, Degree, TEMPLATE_ARGS(Degree::Degree0, Degree::Degree1, Degree::Degree2), 4, field::Form, TEMPLATE_ARGS(field::Form::Form0, field::Form::Form1, field::Form::Form2, field::Form::Form3), 2, unsigned int, TEMPLATE_ARGS(2, 3), TEMPLATE_PARAMS(const function::Function< TEMPLATE_PARAM_1, TEMPLATE_PARAM_2 > &, const UnknownCompoundField< TEMPLATE_PARAM_1, TEMPLATE_PARAM_3, TEMPLATE_PARAM_4 > &))
  INSTANTIATE_OPT_4(TEMPLATE_RETURN(EquationLhsCompound< TEMPLATE_PARAM_1, TEMPLATE_PARAM_2, TEMPLATE_PARAM_3, Product::ScalarProduct, TEMPLATE_PARAM_4 >), , 5, 4, class, TEMPLATE_ARGS(std::complex< double >, std::complex< float >, double, float), 1, Degree, TEMPLATE_ARGS(Degree::Degree4), 2, field::Form, TEMPLATE_ARGS(field::Form::Form1, field::Form::Form2), 2, unsigned int, TEMPLATE_ARGS(2, 3), TEMPLATE_PARAMS(const function::Function< TEMPLATE_PARAM_1, TEMPLATE_PARAM_2 > &, const UnknownCompoundField< TEMPLATE_PARAM_1, TEMPLATE_PARAM_3, TEMPLATE_PARAM_4 > &))


  template< class T_Scalar, field::Form T_Form, unsigned int T_NumFields >
  EquationLhsCompound< T_Scalar, Degree::Degree0, T_Form, Product::ScalarProduct, T_NumFields > operator*(const typename MathObject< T_Scalar, Degree::Degree0 >::Object &value, const UnknownCompoundField< T_Scalar, T_Form, T_NumFields > &field)
  {
    return EquationLhsCompound< T_Scalar, Degree::Degree0, T_Form, Product::ScalarProduct, T_NumFields >(function::ScalarFunction< T_Scalar >(value), field);
  }

  INSTANTIATE_OPT_3(TEMPLATE_RETURN(EquationLhsCompound< TEMPLATE_PARAM_1, Degree::Degree0, TEMPLATE_PARAM_2, Product::ScalarProduct, TEMPLATE_PARAM_3 >), , 6, 4, class, TEMPLATE_ARGS(std::complex< double >, std::complex< float >, double, float), 4, field::Form, TEMPLATE_ARGS(field::Form::Form0, field::Form::Form1, field::Form::Form2, field::Form::Form3), 2, unsigned int, TEMPLATE_ARGS(2, 3), TEMPLATE_PARAMS(const typename MathObject< TEMPLATE_PARAM_1, Degree::Degree0 >::Object &, const UnknownCompoundField< TEMPLATE_PARAM_1, TEMPLATE_PARAM_2, TEMPLATE_PARAM_3 > &))


  template< class T_Scalar, field::Form T_Form, unsigned int T_NumFields >
  EquationLhsCompound< T_Scalar, Degree::Degree1, T_Form, Product::ScalarProduct, T_NumFields > operator*(const typename MathObject< T_Scalar, Degree::Degree1 >::Object &value, const UnknownCompoundField< T_Scalar, T_Form, T_NumFields > &field)
  {
    return EquationLhsCompound< T_Scalar, Degree::Degree1, T_Form, Product::ScalarProduct, T_NumFields >(function::VectorFunction< T_Scalar >(value), field);
  }

  INSTANTIATE_OPT_3(TEMPLATE_RETURN(EquationLhsCompound< TEMPLATE_PARAM_1, Degree::Degree1, TEMPLATE_PARAM_2, Product::ScalarProduct, TEMPLATE_PARAM_3 >), , 7, 4, class, TEMPLATE_ARGS(std::complex< double >, std::complex< float >, double, float), 4, field::Form, TEMPLATE_ARGS(field::Form::Form0, field::Form::Form1, field::Form::Form2, field::Form::Form3), 2, unsigned int, TEMPLATE_ARGS(2, 3), TEMPLATE_PARAMS(const typename MathObject< TEMPLATE_PARAM_1, Degree::Degree1 >::Object &, const UnknownCompoundField< TEMPLATE_PARAM_1, TEMPLATE_PARAM_2, TEMPLATE_PARAM_3 > &))


  template< class T_Scalar, field::Form T_Form, unsigned int T_NumFields >
  EquationLhsCompound< T_Scalar, Degree::Degree2, T_Form, Product::ScalarProduct, T_NumFields > operator*(const typename MathObject< T_Scalar, Degree::Degree2 >::Object &value, const UnknownCompoundField< T_Scalar, T_Form, T_NumFields > &field)
  {
    return EquationLhsCompound< T_Scalar, Degree::Degree2, T_Form, Product::ScalarProduct, T_NumFields >(function::TensorFunction< T_Scalar >(value), field);
  }

  INSTANTIATE_OPT_3(TEMPLATE_RETURN(EquationLhsCompound< TEMPLATE_PARAM_1, Degree::Degree2, TEMPLATE_PARAM_2, Product::ScalarProduct, TEMPLATE_PARAM_3 >), , 8, 4, class, TEMPLATE_ARGS(std::complex< double >, std::complex< float >, double, float), 4, field::Form, TEMPLATE_ARGS(field::Form::Form0, field::Form::Form1, field::Form::Form2, field::Form::Form3), 2, unsigned int, TEMPLATE_ARGS(2, 3), TEMPLATE_PARAMS(const typename MathObject< TEMPLATE_PARAM_1, Degree::Degree2 >::Object &, const UnknownCompoundField< TEMPLATE_PARAM_1, TEMPLATE_PARAM_2, TEMPLATE_PARAM_3 > &))


  template< class T_Scalar, field::Form T_Form, unsigned int T_NumFields >
  EquationLhsCompound< T_Scalar, Degree::Degree4, T_Form, Product::ScalarProduct, T_NumFields > operator*(const typename MathObject< T_Scalar, Degree::Degree4 >::Object &value, const UnknownCompoundField< T_Scalar, T_Form, T_NumFields > &field)
  {
    return EquationLhsCompound< T_Scalar, Degree::Degree4, T_Form, Product::ScalarProduct, T_NumFields >(function::TensorFunction< T_Scalar, 4 >(value), field);
  }

  INSTANTIATE_OPT_3(TEMPLATE_RETURN(EquationLhsCompound< TEMPLATE_PARAM_1, Degree::Degree4, TEMPLATE_PARAM_2, Product::ScalarProduct, TEMPLATE_PARAM_3 >), , 9, 4, class, TEMPLATE_ARGS(std::complex< double >, std::complex< float >, double, float), 2, field::Form, TEMPLATE_ARGS(field::Form::Form1, field::Form::Form2), 2, unsigned int, TEMPLATE_ARGS(2, 3), TEMPLATE_PARAMS(const typename MathObject< TEMPLATE_PARAM_1, Degree::Degree4 >::Object &, const UnknownCompoundField< TEMPLATE_PARAM_1, TEMPLATE_PARAM_2, TEMPLATE_PARAM_3 > &))


  template< class T_Scalar, field::Form T_Form >
  EquationLhs< T_Scalar, Degree::Degree1, T_Form, Product::VectorProduct > operator%(const function::Function< T_Scalar, Degree::Degree1 > &function, const UnknownField< T_Scalar, T_Form > &field)
  {
    return EquationLhs< T_Scalar, Degree::Degree1, T_Form, Product::VectorProduct >(function, field);
  }

  INSTANTIATE_OPV_2(TEMPLATE_RETURN(EquationLhs< TEMPLATE_PARAM_1, Degree::Degree1, TEMPLATE_PARAM_2, Product::VectorProduct >), , 0, 4, class, TEMPLATE_ARGS(std::complex< double >, std::complex< float >, double, float), 2, field::Form, TEMPLATE_ARGS(field::Form::Form1, field::Form::Form2), TEMPLATE_PARAMS(const function::Function< TEMPLATE_PARAM_1, Degree::Degree1 > &, const UnknownField< TEMPLATE_PARAM_1, TEMPLATE_PARAM_2 > &))


  template< class T_Scalar, field::Form T_Form >
  EquationLhs< T_Scalar, Degree::Degree1, T_Form, Product::VectorProduct > operator%(const typename MathObject< T_Scalar, Degree::Degree1 >::Object &value, const UnknownField< T_Scalar, T_Form > &field)
  {
    return EquationLhs< T_Scalar, Degree::Degree1, T_Form, Product::VectorProduct >(function::VectorFunction< T_Scalar >(value), field);
  }

  INSTANTIATE_OPV_2(TEMPLATE_RETURN(EquationLhs< TEMPLATE_PARAM_1, Degree::Degree1, TEMPLATE_PARAM_2, Product::VectorProduct >), , 1, 4, class, TEMPLATE_ARGS(std::complex< double >, std::complex< float >, double, float), 2, field::Form, TEMPLATE_ARGS(field::Form::Form1, field::Form::Form2), TEMPLATE_PARAMS(const typename MathObject< TEMPLATE_PARAM_1, Degree::Degree1 >::Object &, const UnknownField< TEMPLATE_PARAM_1, TEMPLATE_PARAM_2 > &))


  // Rhs

  template< class T_Scalar, Degree T_Degree, field::Form T_Form >
  EquationRhs< T_Scalar, T_Degree, T_Form, Product::ScalarProduct > operator*(const UnknownField< T_Scalar, T_Form > &field, const function::Function< T_Scalar, T_Degree > &function)
  {
    return EquationRhs< T_Scalar, T_Degree, T_Form, Product::ScalarProduct >(function, field);
  }

  INSTANTIATE_OPT_3(TEMPLATE_RETURN(EquationRhs< TEMPLATE_PARAM_1, TEMPLATE_PARAM_2, TEMPLATE_PARAM_3, Product::ScalarProduct >), , 10, 4, class, TEMPLATE_ARGS(std::complex< double >, std::complex< float >, double, float), 3, Degree, TEMPLATE_ARGS(Degree::Degree0, Degree::Degree1, Degree::Degree2), 4, field::Form, TEMPLATE_ARGS(field::Form::Form0, field::Form::Form1, field::Form::Form2, field::Form::Form3), TEMPLATE_PARAMS(const UnknownField< TEMPLATE_PARAM_1, TEMPLATE_PARAM_3 > &, const function::Function< TEMPLATE_PARAM_1, TEMPLATE_PARAM_2 > &))


  template< class T_Scalar, field::Form T_Form >
  EquationRhs< T_Scalar, Degree::Degree0, T_Form, Product::ScalarProduct > operator*(const UnknownField< T_Scalar, T_Form > &field, const typename MathObject< T_Scalar, Degree::Degree0 >::Object &value)
  {
    return EquationRhs< T_Scalar, Degree::Degree0, T_Form, Product::ScalarProduct >(function::ScalarFunction< T_Scalar >(value), field);
  }

  INSTANTIATE_OPT_2(TEMPLATE_RETURN(EquationRhs< TEMPLATE_PARAM_1, Degree::Degree0, TEMPLATE_PARAM_2, Product::ScalarProduct >), , 11, 4, class, TEMPLATE_ARGS(std::complex< double >, std::complex< float >, double, float), 4, field::Form, TEMPLATE_ARGS(field::Form::Form0, field::Form::Form1, field::Form::Form2, field::Form::Form3), TEMPLATE_PARAMS(const UnknownField< TEMPLATE_PARAM_1, TEMPLATE_PARAM_2 > &field, const typename MathObject< TEMPLATE_PARAM_1, Degree::Degree0 >::Object &))


  template< class T_Scalar, field::Form T_Form >
  EquationRhs< T_Scalar, Degree::Degree1, T_Form, Product::ScalarProduct > operator*(const UnknownField< T_Scalar, T_Form > &field, const typename MathObject< T_Scalar, Degree::Degree1 >::Object &value)
  {
    return EquationRhs< T_Scalar, Degree::Degree1, T_Form, Product::ScalarProduct >(function::VectorFunction< T_Scalar >(value), field);
  }

  INSTANTIATE_OPT_2(TEMPLATE_RETURN(EquationRhs< TEMPLATE_PARAM_1, Degree::Degree1, TEMPLATE_PARAM_2, Product::ScalarProduct >), , 12, 4, class, TEMPLATE_ARGS(std::complex< double >, std::complex< float >, double, float), 4, field::Form, TEMPLATE_ARGS(field::Form::Form0, field::Form::Form1, field::Form::Form2, field::Form::Form3), TEMPLATE_PARAMS(const UnknownField< TEMPLATE_PARAM_1, TEMPLATE_PARAM_2 > &field, const typename MathObject< TEMPLATE_PARAM_1, Degree::Degree1 >::Object &))


  template< class T_Scalar, field::Form T_Form >
  EquationRhs< T_Scalar, Degree::Degree2, T_Form, Product::ScalarProduct > operator*(const UnknownField< T_Scalar, T_Form > &field, const typename MathObject< T_Scalar, Degree::Degree2 >::Object &value)
  {
    return EquationRhs< T_Scalar, Degree::Degree2, T_Form, Product::ScalarProduct >(function::TensorFunction< T_Scalar >(value), field);
  }

  INSTANTIATE_OPT_2(TEMPLATE_RETURN(EquationRhs< TEMPLATE_PARAM_1, Degree::Degree2, TEMPLATE_PARAM_2, Product::ScalarProduct >), , 13, 4, class, TEMPLATE_ARGS(std::complex< double >, std::complex< float >, double, float), 4, field::Form, TEMPLATE_ARGS(field::Form::Form0, field::Form::Form1, field::Form::Form2, field::Form::Form3), TEMPLATE_PARAMS(const UnknownField< TEMPLATE_PARAM_1, TEMPLATE_PARAM_2 > &field, const typename MathObject< TEMPLATE_PARAM_1, Degree::Degree2 >::Object &))


  template< class T_Scalar, Degree T_Degree, field::Form T_Form, unsigned int T_NumFields >
  EquationRhsCompound< T_Scalar, T_Degree, T_Form, Product::ScalarProduct, T_NumFields > operator*(const UnknownCompoundField< T_Scalar, T_Form, T_NumFields > &field, const function::Function< T_Scalar, T_Degree > &function)
  {
    return EquationRhsCompound< T_Scalar, T_Degree, T_Form, Product::ScalarProduct, T_NumFields >(function, field);
  }

  INSTANTIATE_OPT_4(TEMPLATE_RETURN(EquationRhsCompound< TEMPLATE_PARAM_1, TEMPLATE_PARAM_2, TEMPLATE_PARAM_3, Product::ScalarProduct, TEMPLATE_PARAM_4 >), , 14, 4, class, TEMPLATE_ARGS(std::complex< double >, std::complex< float >, double, float), 3, Degree, TEMPLATE_ARGS(Degree::Degree0, Degree::Degree1, Degree::Degree2), 4, field::Form, TEMPLATE_ARGS(field::Form::Form0, field::Form::Form1, field::Form::Form2, field::Form::Form3), 2, unsigned int, TEMPLATE_ARGS(2, 3), TEMPLATE_PARAMS(const UnknownCompoundField< TEMPLATE_PARAM_1, TEMPLATE_PARAM_3, TEMPLATE_PARAM_4 > &, const function::Function< TEMPLATE_PARAM_1, TEMPLATE_PARAM_2 > &))


  template< class T_Scalar, field::Form T_Form, unsigned int T_NumFields >
  EquationRhsCompound< T_Scalar, Degree::Degree0, T_Form, Product::ScalarProduct, T_NumFields > operator*(const UnknownCompoundField< T_Scalar, T_Form, T_NumFields > &field, const typename MathObject< T_Scalar, Degree::Degree0 >::Object &value)
  {
    return EquationRhsCompound< T_Scalar, Degree::Degree0, T_Form, Product::ScalarProduct, T_NumFields >(function::ScalarFunction< T_Scalar >(value), field);
  }

  INSTANTIATE_OPT_3(TEMPLATE_RETURN(EquationRhsCompound< TEMPLATE_PARAM_1, Degree::Degree0, TEMPLATE_PARAM_2, Product::ScalarProduct, TEMPLATE_PARAM_3 >), , 15, 4, class, TEMPLATE_ARGS(std::complex< double >, std::complex< float >, double, float), 4, field::Form, TEMPLATE_ARGS(field::Form::Form0, field::Form::Form1, field::Form::Form2, field::Form::Form3), 2, unsigned int, TEMPLATE_ARGS(2, 3), TEMPLATE_PARAMS(const UnknownCompoundField< TEMPLATE_PARAM_1, TEMPLATE_PARAM_2, TEMPLATE_PARAM_3 > &, const typename MathObject< TEMPLATE_PARAM_1, Degree::Degree0 >::Object &))


  template< class T_Scalar, field::Form T_Form, unsigned int T_NumFields >
  EquationRhsCompound< T_Scalar, Degree::Degree1, T_Form, Product::ScalarProduct, T_NumFields > operator*(const UnknownCompoundField< T_Scalar, T_Form, T_NumFields > &field, const typename MathObject< T_Scalar, Degree::Degree1 >::Object &value)
  {
    return EquationRhsCompound< T_Scalar, Degree::Degree1, T_Form, Product::ScalarProduct, T_NumFields >(function::VectorFunction< T_Scalar >(value), field);
  }

  INSTANTIATE_OPT_3(TEMPLATE_RETURN(EquationRhsCompound< TEMPLATE_PARAM_1, Degree::Degree1, TEMPLATE_PARAM_2, Product::ScalarProduct, TEMPLATE_PARAM_3 >), , 16, 4, class, TEMPLATE_ARGS(std::complex< double >, std::complex< float >, double, float), 4, field::Form, TEMPLATE_ARGS(field::Form::Form0, field::Form::Form1, field::Form::Form2, field::Form::Form3), 2, unsigned int, TEMPLATE_ARGS(2, 3), TEMPLATE_PARAMS(const UnknownCompoundField< TEMPLATE_PARAM_1, TEMPLATE_PARAM_2, TEMPLATE_PARAM_3 > &, const typename MathObject< TEMPLATE_PARAM_1, Degree::Degree1 >::Object &))


  template< class T_Scalar, field::Form T_Form, unsigned int T_NumFields >
  EquationRhsCompound< T_Scalar, Degree::Degree2, T_Form, Product::ScalarProduct, T_NumFields > operator*(const UnknownCompoundField< T_Scalar, T_Form, T_NumFields > &field, const typename MathObject< T_Scalar, Degree::Degree2 >::Object &value)
  {
    return EquationRhsCompound< T_Scalar, Degree::Degree2, T_Form, Product::ScalarProduct, T_NumFields >(function::TensorFunction< T_Scalar >(value), field);
  }


  INSTANTIATE_OPT_3(TEMPLATE_RETURN(EquationRhsCompound< TEMPLATE_PARAM_1, Degree::Degree2, TEMPLATE_PARAM_2, Product::ScalarProduct, TEMPLATE_PARAM_3 >), , 17, 4, class, TEMPLATE_ARGS(std::complex< double >, std::complex< float >, double, float), 4, field::Form, TEMPLATE_ARGS(field::Form::Form0, field::Form::Form1, field::Form::Form2, field::Form::Form3), 2, unsigned int, TEMPLATE_ARGS(2, 3), TEMPLATE_PARAMS(const UnknownCompoundField< TEMPLATE_PARAM_1, TEMPLATE_PARAM_2, TEMPLATE_PARAM_3 > &, const typename MathObject< TEMPLATE_PARAM_1, Degree::Degree2 >::Object &))


  template< class T_Scalar, field::Form T_Form >
  EquationRhs< T_Scalar, Degree::Degree1, T_Form, Product::VectorProduct > operator%(const UnknownField< T_Scalar, T_Form > &field, const function::Function< T_Scalar, Degree::Degree1 > &function)
  {
    return EquationRhs< T_Scalar, Degree::Degree1, T_Form, Product::VectorProduct >(function, field);
  }

  INSTANTIATE_OPV_2(TEMPLATE_RETURN(EquationRhs< TEMPLATE_PARAM_1, Degree::Degree1, TEMPLATE_PARAM_2, Product::VectorProduct >), , 2, 4, class, TEMPLATE_ARGS(std::complex< double >, std::complex< float >, double, float), 2, field::Form, TEMPLATE_ARGS(field::Form::Form1, field::Form::Form2), TEMPLATE_PARAMS(const UnknownField< TEMPLATE_PARAM_1, TEMPLATE_PARAM_2 > &, const function::Function< TEMPLATE_PARAM_1, Degree::Degree1 > &))


  template< class T_Scalar, field::Form T_Form >
  EquationRhs< T_Scalar, Degree::Degree1, T_Form, Product::VectorProduct > operator%(const UnknownField< T_Scalar, T_Form > &field, const typename MathObject< T_Scalar, Degree::Degree1 >::Object &value)
  {
    return EquationRhs< T_Scalar, Degree::Degree1, T_Form, Product::VectorProduct >(function::VectorFunction< T_Scalar >(value), field);
  }

  INSTANTIATE_OPV_2(TEMPLATE_RETURN(EquationRhs< TEMPLATE_PARAM_1, Degree::Degree1, TEMPLATE_PARAM_2, Product::VectorProduct >), , 3, 4, class, TEMPLATE_ARGS(std::complex< double >, std::complex< float >, double, float), 2, field::Form, TEMPLATE_ARGS(field::Form::Form1, field::Form::Form2), TEMPLATE_PARAMS(const UnknownField< TEMPLATE_PARAM_1, TEMPLATE_PARAM_2 > &, const typename MathObject< TEMPLATE_PARAM_1, Degree::Degree1 >::Object &))


  // LhsRhs

  template< class T_Scalar, field::Form T_Form >
  EquationLhsRhs< T_Scalar, Degree::Degree1, Degree::Degree1, T_Form, Product::VectorProduct, Product::VectorProduct, term::evaluator::EvaluationOrder::Right > operator%(const function::Function< T_Scalar, Degree::Degree1 > &function, const EquationRhs< T_Scalar, Degree::Degree1, T_Form, Product::VectorProduct > &equation)
  {
    return EquationLhsRhs< T_Scalar, Degree::Degree1, Degree::Degree1, T_Form, Product::VectorProduct, Product::VectorProduct, term::evaluator::EvaluationOrder::Right >(function, *static_cast< equation::UnknownField< T_Scalar, T_Form > * >(equation.getUnknownField()->copy()), equation.getFunction());
  }

  INSTANTIATE_OPV_2(TEMPLATE_RETURN(EquationLhsRhs< TEMPLATE_PARAM_1, Degree::Degree1, Degree::Degree1, TEMPLATE_PARAM_2, Product::VectorProduct, Product::VectorProduct, term::evaluator::EvaluationOrder::Right >), , 4, 4, class, TEMPLATE_ARGS(std::complex< double >, std::complex< float >, double, float), 2, field::Form, TEMPLATE_ARGS(field::Form::Form1, field::Form::Form2), TEMPLATE_PARAMS(const function::Function< TEMPLATE_PARAM_1, Degree::Degree1 > &, const EquationRhs< TEMPLATE_PARAM_1, Degree::Degree1, TEMPLATE_PARAM_2, Product::VectorProduct > &))


} // namespace gmshfem::equation
