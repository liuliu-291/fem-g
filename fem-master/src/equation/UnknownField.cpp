// GmshFEM - Copyright (C) 2019-2022, A. Royer, E. Béchet, C. Geuzaine, Université de Liège
//
// See the LICENSE.txt file for license information. Please report all
// issues on https://gitlab.onelab.info/gmsh/fem/issues

#include "UnknownField.h"

#include "Equation.h"
#include "Exception.h"
#include "FieldEvaluator.h"
#include "Function.h"
#include "instantiate.h"

namespace gmshfem::equation
{


  // ****************************
  // class UnknownFieldInterface
  // ****************************

  template< class T_Scalar >
  UnknownFieldInterface< T_Scalar >::UnknownFieldInterface(const bool isDerivative, field::FieldInterface< T_Scalar > *field, const UnknownFieldType &type) :
    _isDerivative(isDerivative), _type(type), _field(field)
  {
  }

  template< class T_Scalar >
  UnknownFieldInterface< T_Scalar >::~UnknownFieldInterface()
  {
  }

  template< class T_Scalar >
  bool UnknownFieldInterface< T_Scalar >::isDerivative() const
  {
    return _isDerivative;
  }

  template< class T_Scalar >
  UnknownFieldType UnknownFieldInterface< T_Scalar >::type() const
  {
    return _type;
  }

  template< class T_Scalar >
  field::FieldInterface< T_Scalar > *UnknownFieldInterface< T_Scalar >::getField() const
  {
    return _field;
  }

  template< class T_Scalar >
  void UnknownFieldInterface< T_Scalar >::changeField(field::FieldInterface< T_Scalar > *other)
  {
    if(other->form() != _field->form()) {
      throw common::Exception("You switch unknown field '" + _field->name() + "' by a field '" + other->name() + "' that do not have the same form");
    }
    _field = other;
  }

  INSTANTIATE_CLASS(UnknownFieldInterface, 4, TEMPLATE_ARGS(std::complex< double >, std::complex< float >, double, float))


  // ******************
  // class UnknownField
  // ******************

  template< class T_Scalar, field::Form T_Form >
  UnknownField< T_Scalar, T_Form >::UnknownField(const bool isDerivative, field::FieldInterface< T_Scalar > *field, const UnknownFieldType &type) :
    UnknownFieldInterface< T_Scalar >(isDerivative, field, type)
  {
  }

  template< class T_Scalar, field::Form T_Form >
  UnknownField< T_Scalar, T_Form >::~UnknownField()
  {
  }

  template< class T_Scalar, field::Form T_Form >
  UnknownField< T_Scalar, T_Form > *UnknownField< T_Scalar, T_Form >::copy() const
  {
    return new UnknownField< T_Scalar, T_Form >(this->_isDerivative, this->_field, this->_type);
  }

  template< class T_Scalar, field::Form T_Form >
  UnknownField< T_Scalar, field::NextForm< T_Form >::value > *UnknownField< T_Scalar, T_Form >::exteriorDerivative() const
  {
    return new UnknownField< T_Scalar, field::NextForm< T_Form >::value >(true, this->_field, this->_type);
  }

  template< class T_Scalar, field::Form T_Form >
  UnknownField< T_Scalar, T_Form >::operator EquationLhs< T_Scalar, Degree::Degree0, T_Form, Product::ScalarProduct >() const
  {
    return EquationLhs< T_Scalar, Degree::Degree0, T_Form, Product::ScalarProduct >(function::ScalarFunction< T_Scalar >(1.), *this);
  }

  template< class T_Scalar, field::Form T_Form >
  term::evaluator::FieldEvaluator< T_Scalar, T_Form > *UnknownField< T_Scalar, T_Form >::getFieldEvaluator() const
  {
    return new term::evaluator::FieldEvaluator< T_Scalar, T_Form >();
  }

  template< class T_Scalar, field::Form T_Form >
  function::Function< T_Scalar, field::DegreeOfForm< T_Form >::value > UnknownField< T_Scalar, T_Form >::function(const field::FieldInterface< T_Scalar > *const field) const
  {
    if(this->_isDerivative) {
      return function::d(*static_cast< const field::Field< T_Scalar, field::PastForm< T_Form >::value > * >(field ? field : this->_field));
    }

    return function::Function< T_Scalar, field::DegreeOfForm< T_Form >::value >(*static_cast< const field::Field< T_Scalar, T_Form > * >(field ? field : this->_field));
  }

  INSTANTIATE_CLASS_2(UnknownField, 4, 4, TEMPLATE_ARGS(std::complex< double >, std::complex< float >, double, float), TEMPLATE_ARGS(field::Form::Form0, field::Form::Form1, field::Form::Form2, field::Form::Form3))


  // ******************
  // class UnknownCompoundField
  // ******************

  template< class T_Scalar, field::Form T_Form, unsigned int T_NumFields >
  UnknownCompoundField< T_Scalar, T_Form, T_NumFields >::UnknownCompoundField(const bool isDerivative, field::FieldInterface< T_Scalar > *field, const UnknownFieldType &type) :
    UnknownFieldInterface< T_Scalar >(isDerivative, field, type)
  {
  }

  template< class T_Scalar, field::Form T_Form, unsigned int T_NumFields >
  UnknownCompoundField< T_Scalar, T_Form, T_NumFields >::~UnknownCompoundField()
  {
  }

  template< class T_Scalar, field::Form T_Form, unsigned int T_NumFields >
  UnknownCompoundField< T_Scalar, T_Form, T_NumFields > *UnknownCompoundField< T_Scalar, T_Form, T_NumFields >::copy() const
  {
    return new UnknownCompoundField< T_Scalar, T_Form, T_NumFields >(this->_isDerivative, this->_field, this->_type);
  }

  template< class T_Scalar, field::Form T_Form, unsigned int T_NumFields >
  UnknownCompoundField< T_Scalar, field::NextForm< T_Form >::value, T_NumFields > *UnknownCompoundField< T_Scalar, T_Form, T_NumFields >::exteriorDerivative() const
  {
    return new UnknownCompoundField< T_Scalar, field::NextForm< T_Form >::value, T_NumFields >(true, this->_field, this->_type);
  }

  template< class T_Scalar, field::Form T_Form, unsigned int T_NumFields >
  UnknownCompoundField< T_Scalar, T_Form, T_NumFields >::operator EquationLhsCompound< T_Scalar, Degree::Degree0, T_Form, Product::ScalarProduct, T_NumFields >() const
  {
    return EquationLhsCompound< T_Scalar, Degree::Degree0, T_Form, Product::ScalarProduct, T_NumFields >(function::ScalarFunction< T_Scalar >(1.), *this);
  }

  template< class T_Scalar, field::Form T_Form, unsigned int T_NumFields >
  term::evaluator::CompoundFieldEvaluator< T_Scalar, T_Form, T_NumFields > *UnknownCompoundField< T_Scalar, T_Form, T_NumFields >::getFieldEvaluator() const
  {
    return new term::evaluator::CompoundFieldEvaluator< T_Scalar, T_Form, T_NumFields >();
  }

  template< class T_Scalar, field::Form T_Form, unsigned int T_NumFields >
  function::Function< T_Scalar, field::DegreeOfCompoundForm< T_Form >::value > UnknownCompoundField< T_Scalar, T_Form, T_NumFields >::function(const field::FieldInterface< T_Scalar > *const field) const
  {
    if(this->_isDerivative) {
      return function::d(*static_cast< const field::CompoundField< T_Scalar, field::PastForm< T_Form >::value, T_NumFields > * >(field ? field : this->_field));
    }
    return function::Function< T_Scalar, field::DegreeOfCompoundForm< T_Form >::value >(*static_cast< const field::CompoundField< T_Scalar, T_Form, T_NumFields > * >(field ? field : this->_field));
  }

  INSTANTIATE_CLASS_3(UnknownCompoundField, 4, 4, 2, TEMPLATE_ARGS(std::complex< double >, std::complex< float >, double, float), TEMPLATE_ARGS(field::Form::Form0, field::Form::Form1, field::Form::Form2, field::Form::Form3), TEMPLATE_ARGS(2, 3))


  // function

  template< class T_Scalar, field::Form T_Form >
  UnknownField< T_Scalar, T_Form > &dof(field::Field< T_Scalar, T_Form > &field)
  {
    return *(new UnknownField< T_Scalar, T_Form >(false, &field, UnknownFieldType::NotDt));
  }

  INSTANTIATE_FCT_2(TEMPLATE_RETURN(UnknownField< TEMPLATE_PARAM_1, TEMPLATE_PARAM_2 > &), , dof, 0, 4, class, TEMPLATE_ARGS(std::complex< double >, std::complex< float >, double, float), 4, field::Form, TEMPLATE_ARGS(field::Form::Form0, field::Form::Form1, field::Form::Form2, field::Form::Form3), TEMPLATE_PARAMS(field::Field< TEMPLATE_PARAM_1, TEMPLATE_PARAM_2 > &))


  template< class T_Scalar, field::Form T_Form >
  UnknownField< T_Scalar, T_Form > &dt_dof(field::Field< T_Scalar, T_Form > &field)
  {
    return *(new UnknownField< T_Scalar, T_Form >(false, &field, UnknownFieldType::Dt));
  }

  INSTANTIATE_FCT_2(TEMPLATE_RETURN(UnknownField< TEMPLATE_PARAM_1, TEMPLATE_PARAM_2 > &), , dt_dof, 0, 4, class, TEMPLATE_ARGS(std::complex< double >, std::complex< float >, double, float), 4, field::Form, TEMPLATE_ARGS(field::Form::Form0, field::Form::Form1, field::Form::Form2, field::Form::Form3), TEMPLATE_PARAMS(field::Field< TEMPLATE_PARAM_1, TEMPLATE_PARAM_2 > &))


  template< class T_Scalar, field::Form T_Form >
  UnknownField< T_Scalar, T_Form > &dt2_dof(field::Field< T_Scalar, T_Form > &field)
  {
    return *(new UnknownField< T_Scalar, T_Form >(false, &field, UnknownFieldType::DtDt));
  }

  INSTANTIATE_FCT_2(TEMPLATE_RETURN(UnknownField< TEMPLATE_PARAM_1, TEMPLATE_PARAM_2 > &), , dt2_dof, 0, 4, class, TEMPLATE_ARGS(std::complex< double >, std::complex< float >, double, float), 4, field::Form, TEMPLATE_ARGS(field::Form::Form0, field::Form::Form1, field::Form::Form2, field::Form::Form3), TEMPLATE_PARAMS(field::Field< TEMPLATE_PARAM_1, TEMPLATE_PARAM_2 > &))


  template< class T_Scalar, field::Form T_Form >
  UnknownField< T_Scalar, T_Form > &tf(field::Field< T_Scalar, T_Form > &field)
  {
    return *(new UnknownField< T_Scalar, T_Form >(false, &field, UnknownFieldType::Tf));
  }

  INSTANTIATE_FCT_2(TEMPLATE_RETURN(UnknownField< TEMPLATE_PARAM_1, TEMPLATE_PARAM_2 > &), , tf, 0, 4, class, TEMPLATE_ARGS(std::complex< double >, std::complex< float >, double, float), 4, field::Form, TEMPLATE_ARGS(field::Form::Form0, field::Form::Form1, field::Form::Form2, field::Form::Form3), TEMPLATE_PARAMS(field::Field< TEMPLATE_PARAM_1, TEMPLATE_PARAM_2 > &))


  template< class T_Scalar >
  UnknownField< T_Scalar, field::Form::Form1 > &grad(UnknownField< T_Scalar, field::Form::Form0 > &unknownField)
  {
    UnknownField< T_Scalar, field::Form::Form1 > *ret = unknownField.exteriorDerivative();
    delete &unknownField;
    return *ret;
  }

  INSTANTIATE_FCT(TEMPLATE_RETURN(UnknownField< TEMPLATE_PARAM_1, field::Form::Form1 > &), , grad, 0, 4, class, TEMPLATE_ARGS(std::complex< double >, std::complex< float >, double, float), TEMPLATE_PARAMS(UnknownField< TEMPLATE_PARAM_1, field::Form::Form0 > &))


  template< class T_Scalar >
  UnknownField< T_Scalar, field::Form::Form2 > &curl(UnknownField< T_Scalar, field::Form::Form1 > &unknownField)
  {
    if(unknownField.isDerivative()) {
      throw common::Exception("The curl of the gradient of a " + std::string(field::NameOfForm< field::Form::Form0 >::value) + " field is zero");
    }
    UnknownField< T_Scalar, field::Form::Form2 > *ret = unknownField.exteriorDerivative();
    delete &unknownField;
    return *ret;
  }

  INSTANTIATE_FCT(TEMPLATE_RETURN(UnknownField< TEMPLATE_PARAM_1, field::Form::Form2 > &), , curl, 0, 4, class, TEMPLATE_ARGS(std::complex< double >, std::complex< float >, double, float), TEMPLATE_PARAMS(UnknownField< TEMPLATE_PARAM_1, field::Form::Form1 > &))


  template< class T_Scalar >
  UnknownField< T_Scalar, field::Form::Form3 > &div(UnknownField< T_Scalar, field::Form::Form2 > &unknownField)
  {
    if(unknownField.isDerivative()) {
      throw common::Exception("The divergence of the curl of a " + std::string(field::NameOfForm< field::Form::Form1 >::value) + " field is zero");
    }
    UnknownField< T_Scalar, field::Form::Form3 > *ret = unknownField.exteriorDerivative();
    delete &unknownField;
    return *ret;
  }

  INSTANTIATE_FCT(TEMPLATE_RETURN(UnknownField< TEMPLATE_PARAM_1, field::Form::Form3 > &), , div, 0, 4, class, TEMPLATE_ARGS(std::complex< double >, std::complex< float >, double, float), TEMPLATE_PARAMS(UnknownField< TEMPLATE_PARAM_1, field::Form::Form2 > &))


  template< class T_Scalar, field::Form T_Form >
  UnknownField< T_Scalar, field::NextForm< T_Form >::value > &d(UnknownField< T_Scalar, T_Form > &unknownField)
  {
    if(unknownField.isDerivative()) {
      throw common::Exception("An exact differential is closed");
    }

    UnknownField< T_Scalar, field::NextForm< T_Form >::value > *ret = unknownField.exteriorDerivative();
    delete &unknownField;
    return *ret;
  }

  INSTANTIATE_FCT_2(TEMPLATE_RETURN(UnknownField< TEMPLATE_PARAM_1, field::NextForm< TEMPLATE_PARAM_2 >::value > &), , d, 0, 4, class, TEMPLATE_ARGS(std::complex< double >, std::complex< float >, double, float), 4, field::Form, TEMPLATE_ARGS(field::Form::Form0, field::Form::Form1, field::Form::Form2, field::Form::Form3), TEMPLATE_PARAMS(UnknownField< TEMPLATE_PARAM_1, TEMPLATE_PARAM_2 > &))


  template< class T_Scalar, field::Form T_Form, unsigned int T_NumFields >
  UnknownCompoundField< T_Scalar, T_Form, T_NumFields > &dof(field::CompoundField< T_Scalar, T_Form, T_NumFields > &field)
  {
    return *(new UnknownCompoundField< T_Scalar, T_Form, T_NumFields >(false, &field, UnknownFieldType::NotDt));
  }

  INSTANTIATE_FCT_3(TEMPLATE_RETURN(UnknownCompoundField< TEMPLATE_PARAM_1, TEMPLATE_PARAM_2, TEMPLATE_PARAM_3 > &), , dof, 1, 4, class, TEMPLATE_ARGS(std::complex< double >, std::complex< float >, double, float), 4, field::Form, TEMPLATE_ARGS(field::Form::Form0, field::Form::Form1, field::Form::Form2, field::Form::Form3), 2, unsigned int, TEMPLATE_ARGS(2, 3), TEMPLATE_PARAMS(field::CompoundField< TEMPLATE_PARAM_1, TEMPLATE_PARAM_2, TEMPLATE_PARAM_3 > &))


  template< class T_Scalar, field::Form T_Form, unsigned int T_NumFields >
  UnknownCompoundField< T_Scalar, T_Form, T_NumFields > &dt_dof(field::CompoundField< T_Scalar, T_Form, T_NumFields > &field)
  {
    return *(new UnknownCompoundField< T_Scalar, T_Form, T_NumFields >(false, &field, UnknownFieldType::Dt));
  }

  INSTANTIATE_FCT_3(TEMPLATE_RETURN(UnknownCompoundField< TEMPLATE_PARAM_1, TEMPLATE_PARAM_2, TEMPLATE_PARAM_3 > &), , dt_dof, 1, 4, class, TEMPLATE_ARGS(std::complex< double >, std::complex< float >, double, float), 4, field::Form, TEMPLATE_ARGS(field::Form::Form0, field::Form::Form1, field::Form::Form2, field::Form::Form3), 2, unsigned int, TEMPLATE_ARGS(2, 3), TEMPLATE_PARAMS(field::CompoundField< TEMPLATE_PARAM_1, TEMPLATE_PARAM_2, TEMPLATE_PARAM_3 > &))


  template< class T_Scalar, field::Form T_Form, unsigned int T_NumFields >
  UnknownCompoundField< T_Scalar, T_Form, T_NumFields > &dt2_dof(field::CompoundField< T_Scalar, T_Form, T_NumFields > &field)
  {
    return *(new UnknownCompoundField< T_Scalar, T_Form, T_NumFields >(false, &field, UnknownFieldType::DtDt));
  }

  INSTANTIATE_FCT_3(TEMPLATE_RETURN(UnknownCompoundField< TEMPLATE_PARAM_1, TEMPLATE_PARAM_2, TEMPLATE_PARAM_3 > &), , dt2_dof, 1, 4, class, TEMPLATE_ARGS(std::complex< double >, std::complex< float >, double, float), 4, field::Form, TEMPLATE_ARGS(field::Form::Form0, field::Form::Form1, field::Form::Form2, field::Form::Form3), 2, unsigned int, TEMPLATE_ARGS(2, 3), TEMPLATE_PARAMS(field::CompoundField< TEMPLATE_PARAM_1, TEMPLATE_PARAM_2, TEMPLATE_PARAM_3 > &))


  template< class T_Scalar, field::Form T_Form, unsigned int T_NumFields >
  UnknownCompoundField< T_Scalar, T_Form, T_NumFields > &tf(field::CompoundField< T_Scalar, T_Form, T_NumFields > &field)
  {
    return *(new UnknownCompoundField< T_Scalar, T_Form, T_NumFields >(false, &field, UnknownFieldType::Tf));
  }

  INSTANTIATE_FCT_3(TEMPLATE_RETURN(UnknownCompoundField< TEMPLATE_PARAM_1, TEMPLATE_PARAM_2, TEMPLATE_PARAM_3 > &), , tf, 1, 4, class, TEMPLATE_ARGS(std::complex< double >, std::complex< float >, double, float), 4, field::Form, TEMPLATE_ARGS(field::Form::Form0, field::Form::Form1, field::Form::Form2, field::Form::Form3), 2, unsigned int, TEMPLATE_ARGS(2, 3), TEMPLATE_PARAMS(field::CompoundField< TEMPLATE_PARAM_1, TEMPLATE_PARAM_2, TEMPLATE_PARAM_3 > &))


  template< class T_Scalar, unsigned int T_NumFields >
  UnknownCompoundField< T_Scalar, field::Form::Form1, T_NumFields > &grad(UnknownCompoundField< T_Scalar, field::Form::Form0, T_NumFields > &unknownField)
  {
    UnknownCompoundField< T_Scalar, field::Form::Form1, T_NumFields > *ret = unknownField.exteriorDerivative();
    delete &unknownField;
    return *ret;
  }

  INSTANTIATE_FCT_2(TEMPLATE_RETURN(UnknownCompoundField< TEMPLATE_PARAM_1, field::Form::Form1, TEMPLATE_PARAM_2 > &), , grad, 1, 4, class, TEMPLATE_ARGS(std::complex< double >, std::complex< float >, double, float), 2, unsigned int, TEMPLATE_ARGS(2, 3), TEMPLATE_PARAMS(UnknownCompoundField< TEMPLATE_PARAM_1, field::Form::Form0, TEMPLATE_PARAM_2 > &))


  template< class T_Scalar, unsigned int T_NumFields >
  UnknownCompoundField< T_Scalar, field::Form::Form2, T_NumFields > &curl(UnknownCompoundField< T_Scalar, field::Form::Form1, T_NumFields > &unknownField)
  {
    if(unknownField.isDerivative()) {
      throw common::Exception("The curl of the gradient of a " + std::string(field::NameOfForm< field::Form::Form0 >::value) + " field is zero");
    }
    UnknownCompoundField< T_Scalar, field::Form::Form2, T_NumFields > *ret = unknownField.exteriorDerivative();
    delete &unknownField;
    return *ret;
  }

  INSTANTIATE_FCT_2(TEMPLATE_RETURN(UnknownCompoundField< TEMPLATE_PARAM_1, field::Form::Form2, TEMPLATE_PARAM_2 > &), , curl, 1, 4, class, TEMPLATE_ARGS(std::complex< double >, std::complex< float >, double, float), 2, unsigned int, TEMPLATE_ARGS(2, 3), TEMPLATE_PARAMS(UnknownCompoundField< TEMPLATE_PARAM_1, field::Form::Form1, TEMPLATE_PARAM_2 > &))


  template< class T_Scalar, unsigned int T_NumFields >
  UnknownCompoundField< T_Scalar, field::Form::Form3, T_NumFields > &div(UnknownCompoundField< T_Scalar, field::Form::Form2, T_NumFields > &unknownField)
  {
    if(unknownField.isDerivative()) {
      throw common::Exception("The divergence of the curl of a " + std::string(field::NameOfForm< field::Form::Form1 >::value) + " field is zero");
    }
    UnknownCompoundField< T_Scalar, field::Form::Form3, T_NumFields > *ret = unknownField.exteriorDerivative();
    delete &unknownField;
    return *ret;
  }

  INSTANTIATE_FCT_2(TEMPLATE_RETURN(UnknownCompoundField< TEMPLATE_PARAM_1, field::Form::Form3, TEMPLATE_PARAM_2 > &), , div, 1, 4, class, TEMPLATE_ARGS(std::complex< double >, std::complex< float >, double, float), 2, unsigned int, TEMPLATE_ARGS(2, 3), TEMPLATE_PARAMS(UnknownCompoundField< TEMPLATE_PARAM_1, field::Form::Form2, TEMPLATE_PARAM_2 > &))


  template< class T_Scalar, field::Form T_Form, unsigned int T_NumFields >
  UnknownCompoundField< T_Scalar, field::NextForm< T_Form >::value, T_NumFields > &d(UnknownCompoundField< T_Scalar, T_Form, T_NumFields > &unknownField)
  {
    if(unknownField.isDerivative()) {
      throw common::Exception("An exact differential is closed");
    }

    UnknownCompoundField< T_Scalar, field::NextForm< T_Form >::value, T_NumFields > *ret = unknownField.exteriorDerivative();
    delete &unknownField;
    return *ret;
  }

  INSTANTIATE_FCT_3(TEMPLATE_RETURN(UnknownCompoundField< TEMPLATE_PARAM_1, field::NextForm< TEMPLATE_PARAM_2 >::value, TEMPLATE_PARAM_3 > &), , d, 1, 4, class, TEMPLATE_ARGS(std::complex< double >, std::complex< float >, double, float), 4, field::Form, TEMPLATE_ARGS(field::Form::Form0, field::Form::Form1, field::Form::Form2, field::Form::Form3), 2, unsigned int, TEMPLATE_ARGS(2, 3), TEMPLATE_PARAMS(UnknownCompoundField< TEMPLATE_PARAM_1, TEMPLATE_PARAM_2, TEMPLATE_PARAM_3 > &))


} // namespace gmshfem::equation
