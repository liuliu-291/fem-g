// GmshFEM - Copyright (C) 2019-2022, A. Royer, E. Béchet, C. Geuzaine, Université de Liège
//
// See the LICENSE.txt file for license information. Please report all
// issues on https://gitlab.onelab.info/gmsh/fem/issues

#ifndef H_GMSHFEM_UNKNOWNFIELD
#define H_GMSHFEM_UNKNOWNFIELD

#include "FieldInterface.h"
#include "FieldObject.h"
#include "Product.h"

namespace gmshfem::equation
{
  template< class T_Scalar, Degree T_Degree, field::Form T_Form, Product T_Product >
  class EquationLhs;

  template< class T_Scalar, Degree T_Degree, field::Form T_Form, Product T_Product, unsigned int T_NumFields >
  class EquationLhsCompound;
} // namespace gmshfem::equation

namespace gmshfem::term::evaluator
{
  template< class T_Scalar, field::Form T_Form >
  class FieldEvaluator;

  template< class T_Scalar, field::Form T_Form, unsigned int T_NumFields >
  class CompoundFieldEvaluator;
} // namespace gmshfem::term::evaluator

namespace gmshfem::equation
{


  enum class UnknownFieldType {
    Tf,
    NotDt,
    Dt,
    DtDt
  };

  // class UnknownFieldInterface

  template< class T_Scalar >
  class UnknownFieldInterface
  {
   protected:
    const bool _isDerivative;
    const UnknownFieldType _type;
    field::FieldInterface< T_Scalar > *_field;

   public:
    UnknownFieldInterface(const bool isDerivative, field::FieldInterface< T_Scalar > *field, const UnknownFieldType &type);
    virtual ~UnknownFieldInterface();

    virtual bool isDerivative() const;
    virtual UnknownFieldType type() const;

    field::FieldInterface< T_Scalar > *getField() const;
    virtual UnknownFieldInterface< T_Scalar > *copy() const = 0;
    void changeField(field::FieldInterface< T_Scalar > *other);
  };

  // class UnknownField

  template< class T_Scalar, field::Form T_Form >
  class UnknownField final : public UnknownFieldInterface< T_Scalar >
  {
   public:
    UnknownField(const bool isDerivative, field::FieldInterface< T_Scalar > *field, const UnknownFieldType &type);
    virtual ~UnknownField();

    virtual UnknownField< T_Scalar, T_Form > *copy() const override;
    virtual UnknownField< T_Scalar, field::NextForm< T_Form >::value > *exteriorDerivative() const;

    virtual operator EquationLhs< T_Scalar, Degree::Degree0, T_Form, Product::ScalarProduct >() const;

    virtual term::evaluator::FieldEvaluator< T_Scalar, T_Form > *getFieldEvaluator() const;
    virtual function::Function< T_Scalar, field::DegreeOfForm< T_Form >::value > function(const field::FieldInterface< T_Scalar > *const field) const;
  };

  // class UnknownCompoundField

  template< class T_Scalar, field::Form T_Form, unsigned int T_NumFields >
  class UnknownCompoundField final : public UnknownFieldInterface< T_Scalar >
  {
   public:
    UnknownCompoundField(const bool isDerivative, field::FieldInterface< T_Scalar > *field, const UnknownFieldType &type);
    virtual ~UnknownCompoundField();

    virtual UnknownCompoundField< T_Scalar, T_Form, T_NumFields > *copy() const override;
    virtual UnknownCompoundField< T_Scalar, field::NextForm< T_Form >::value, T_NumFields > *exteriorDerivative() const;

    virtual operator EquationLhsCompound< T_Scalar, Degree::Degree0, T_Form, Product::ScalarProduct, T_NumFields >() const;

    virtual term::evaluator::CompoundFieldEvaluator< T_Scalar, T_Form, T_NumFields > *getFieldEvaluator() const;
    virtual function::Function< T_Scalar, field::DegreeOfCompoundForm< T_Form >::value > function(const field::FieldInterface< T_Scalar > *const field) const;
  };

  // function

  template< class T_Scalar, field::Form T_Form >
  UnknownField< T_Scalar, T_Form > &dof(field::Field< T_Scalar, T_Form > &field);
  template< class T_Scalar, field::Form T_Form >
  UnknownField< T_Scalar, T_Form > &dt_dof(field::Field< T_Scalar, T_Form > &field);
  template< class T_Scalar, field::Form T_Form >
  UnknownField< T_Scalar, T_Form > &dt2_dof(field::Field< T_Scalar, T_Form > &field);
  template< class T_Scalar, field::Form T_Form >
  UnknownField< T_Scalar, T_Form > &tf(field::Field< T_Scalar, T_Form > &field);

  template< class T_Scalar >
  UnknownField< T_Scalar, field::Form::Form1 > &grad(UnknownField< T_Scalar, field::Form::Form0 > &unknownField);
  template< class T_Scalar >
  UnknownField< T_Scalar, field::Form::Form2 > &curl(UnknownField< T_Scalar, field::Form::Form1 > &unknownField);
  template< class T_Scalar >
  UnknownField< T_Scalar, field::Form::Form3 > &div(UnknownField< T_Scalar, field::Form::Form2 > &unknownField);

  template< class T_Scalar, field::Form T_Form >
  UnknownField< T_Scalar, field::NextForm< T_Form >::value > &d(UnknownField< T_Scalar, T_Form > &unknownField);


  template< class T_Scalar, field::Form T_Form, unsigned int T_NumFields >
  UnknownCompoundField< T_Scalar, T_Form, T_NumFields > &dof(field::CompoundField< T_Scalar, T_Form, T_NumFields > &field);
  template< class T_Scalar, field::Form T_Form, unsigned int T_NumFields >
  UnknownCompoundField< T_Scalar, T_Form, T_NumFields > &dt_dof(field::CompoundField< T_Scalar, T_Form, T_NumFields > &field);
  template< class T_Scalar, field::Form T_Form, unsigned int T_NumFields >
  UnknownCompoundField< T_Scalar, T_Form, T_NumFields > &dt2_dof(field::CompoundField< T_Scalar, T_Form, T_NumFields > &field);
  template< class T_Scalar, field::Form T_Form, unsigned int T_NumFields >
  UnknownCompoundField< T_Scalar, T_Form, T_NumFields > &tf(field::CompoundField< T_Scalar, T_Form, T_NumFields > &field);

  template< class T_Scalar, unsigned int T_NumFields >
  UnknownCompoundField< T_Scalar, field::Form::Form1, T_NumFields > &grad(UnknownCompoundField< T_Scalar, field::Form::Form0, T_NumFields > &unknownField);
  template< class T_Scalar, unsigned int T_NumFields >
  UnknownCompoundField< T_Scalar, field::Form::Form2, T_NumFields > &curl(UnknownCompoundField< T_Scalar, field::Form::Form1, T_NumFields > &unknownField);
  template< class T_Scalar, unsigned int T_NumFields >
  UnknownCompoundField< T_Scalar, field::Form::Form3, T_NumFields > &div(UnknownCompoundField< T_Scalar, field::Form::Form2, T_NumFields > &unknownField);

  template< class T_Scalar, field::Form T_Form, unsigned int T_NumFields >
  UnknownCompoundField< T_Scalar, field::NextForm< T_Form >::value, T_NumFields > &d(UnknownCompoundField< T_Scalar, T_Form, T_NumFields > &unknownField);


} // namespace gmshfem::equation


#endif // H_GMSHFEM_UNKNOWNFIELD
