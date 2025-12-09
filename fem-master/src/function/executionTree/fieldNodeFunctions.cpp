// GmshFEM - Copyright (C) 2019-2022, A. Royer, E. Béchet, C. Geuzaine, Université de Liège
//
// See the LICENSE.txt file for license information. Please report all
// issues on https://gitlab.onelab.info/gmsh/fem/issues

#include "fieldNodeFunctions.h"

#include "fieldOperations.h"

namespace gmshfem::function
{


  template< class T_Scalar >
  static void s_setPointEvaluation(ExecutionTreeIterator &it)
  {
    if((*it)->nodeType() == NodeType::Field) {
      if(auto cast = dynamic_cast< const function::FieldNode< function::None< T_Scalar, field::Form::Form0 > > * >(*it)) {
        cast->operation().setPointEvaluation(true);
        return;
      }
      if(auto cast = dynamic_cast< const function::FieldNode< function::None< T_Scalar, field::Form::Form1 > > * >(*it)) {
        cast->operation().setPointEvaluation(true);
        return;
      }
      if(auto cast = dynamic_cast< const function::FieldNode< function::None< T_Scalar, field::Form::Form2 > > * >(*it)) {
        cast->operation().setPointEvaluation(true);
        return;
      }
      if(auto cast = dynamic_cast< const function::FieldNode< function::None< T_Scalar, field::Form::Form3 > > * >(*it)) {
        cast->operation().setPointEvaluation(true);
        return;
      }

      if(auto cast = dynamic_cast< const function::FieldNode< function::Derivative< T_Scalar, field::Form::Form0 > > * >(*it)) {
        cast->operation().setPointEvaluation(true);
        return;
      }
      if(auto cast = dynamic_cast< const function::FieldNode< function::Derivative< T_Scalar, field::Form::Form1 > > * >(*it)) {
        cast->operation().setPointEvaluation(true);
        return;
      }
      if(auto cast = dynamic_cast< const function::FieldNode< function::Derivative< T_Scalar, field::Form::Form2 > > * >(*it)) {
        cast->operation().setPointEvaluation(true);
        return;
      }
      if(auto cast = dynamic_cast< const function::FieldNode< function::Derivative< T_Scalar, field::Form::Form3 > > * >(*it)) {
        cast->operation().setPointEvaluation(true);
        return;
      }

      if(auto cast = dynamic_cast< const function::FieldNode< function::NoneCompound< T_Scalar, field::Form::Form0, 2 > > * >(*it)) {
        cast->operation().setPointEvaluation(true);
        return;
      }
      if(auto cast = dynamic_cast< const function::FieldNode< function::NoneCompound< T_Scalar, field::Form::Form1, 2 > > * >(*it)) {
        cast->operation().setPointEvaluation(true);
        return;
      }
      if(auto cast = dynamic_cast< const function::FieldNode< function::NoneCompound< T_Scalar, field::Form::Form2, 2 > > * >(*it)) {
        cast->operation().setPointEvaluation(true);
        return;
      }
      if(auto cast = dynamic_cast< const function::FieldNode< function::NoneCompound< T_Scalar, field::Form::Form3, 2 > > * >(*it)) {
        cast->operation().setPointEvaluation(true);
        return;
      }

      if(auto cast = dynamic_cast< const function::FieldNode< function::NoneCompound< T_Scalar, field::Form::Form0, 3 > > * >(*it)) {
        cast->operation().setPointEvaluation(true);
        return;
      }
      if(auto cast = dynamic_cast< const function::FieldNode< function::NoneCompound< T_Scalar, field::Form::Form1, 3 > > * >(*it)) {
        cast->operation().setPointEvaluation(true);
        return;
      }
      if(auto cast = dynamic_cast< const function::FieldNode< function::NoneCompound< T_Scalar, field::Form::Form2, 3 > > * >(*it)) {
        cast->operation().setPointEvaluation(true);
        return;
      }
      if(auto cast = dynamic_cast< const function::FieldNode< function::NoneCompound< T_Scalar, field::Form::Form3, 3 > > * >(*it)) {
        cast->operation().setPointEvaluation(true);
        return;
      }

      if(auto cast = dynamic_cast< const function::FieldNode< function::DerivativeCompound< T_Scalar, field::Form::Form0, 2 > > * >(*it)) {
        cast->operation().setPointEvaluation(true);
        return;
      }
      if(auto cast = dynamic_cast< const function::FieldNode< function::DerivativeCompound< T_Scalar, field::Form::Form1, 2 > > * >(*it)) {
        cast->operation().setPointEvaluation(true);
        return;
      }
      if(auto cast = dynamic_cast< const function::FieldNode< function::DerivativeCompound< T_Scalar, field::Form::Form2, 2 > > * >(*it)) {
        cast->operation().setPointEvaluation(true);
        return;
      }
      if(auto cast = dynamic_cast< const function::FieldNode< function::DerivativeCompound< T_Scalar, field::Form::Form3, 2 > > * >(*it)) {
        cast->operation().setPointEvaluation(true);
        return;
      }

      if(auto cast = dynamic_cast< const function::FieldNode< function::DerivativeCompound< T_Scalar, field::Form::Form0, 3 > > * >(*it)) {
        cast->operation().setPointEvaluation(true);
        return;
      }
      if(auto cast = dynamic_cast< const function::FieldNode< function::DerivativeCompound< T_Scalar, field::Form::Form1, 3 > > * >(*it)) {
        cast->operation().setPointEvaluation(true);
        return;
      }
      if(auto cast = dynamic_cast< const function::FieldNode< function::DerivativeCompound< T_Scalar, field::Form::Form2, 3 > > * >(*it)) {
        cast->operation().setPointEvaluation(true);
        return;
      }
      if(auto cast = dynamic_cast< const function::FieldNode< function::DerivativeCompound< T_Scalar, field::Form::Form3, 3 > > * >(*it)) {
        cast->operation().setPointEvaluation(true);
        return;
      }
    }
    else if((*it)->nodeType() == NodeType::FieldScalarType) {
      if(auto cast = dynamic_cast< const function::FieldScalarTypeNode< function::Eigenfunction< T_Scalar, field::Form::Form0 > > * >(*it)) {
        cast->operation().setPointEvaluation(true);
        return;
      }
      if(auto cast = dynamic_cast< const function::FieldScalarTypeNode< function::Eigenfunction< T_Scalar, field::Form::Form1 > > * >(*it)) {
        cast->operation().setPointEvaluation(true);
        return;
      }
      if(auto cast = dynamic_cast< const function::FieldScalarTypeNode< function::Eigenfunction< T_Scalar, field::Form::Form2 > > * >(*it)) {
        cast->operation().setPointEvaluation(true);
        return;
      }
      if(auto cast = dynamic_cast< const function::FieldScalarTypeNode< function::Eigenfunction< T_Scalar, field::Form::Form3 > > * >(*it)) {
        cast->operation().setPointEvaluation(true);
        return;
      }

      if(auto cast = dynamic_cast< const function::FieldScalarTypeNode< function::EigenfunctionCompound< T_Scalar, field::Form::Form0, 2 > > * >(*it)) {
        cast->operation().setPointEvaluation(true);
        return;
      }
      if(auto cast = dynamic_cast< const function::FieldScalarTypeNode< function::EigenfunctionCompound< T_Scalar, field::Form::Form1, 2 > > * >(*it)) {
        cast->operation().setPointEvaluation(true);
        return;
      }
      if(auto cast = dynamic_cast< const function::FieldScalarTypeNode< function::EigenfunctionCompound< T_Scalar, field::Form::Form2, 2 > > * >(*it)) {
        cast->operation().setPointEvaluation(true);
        return;
      }
      if(auto cast = dynamic_cast< const function::FieldScalarTypeNode< function::EigenfunctionCompound< T_Scalar, field::Form::Form3, 2 > > * >(*it)) {
        cast->operation().setPointEvaluation(true);
        return;
      }

      if(auto cast = dynamic_cast< const function::FieldScalarTypeNode< function::EigenfunctionCompound< T_Scalar, field::Form::Form0, 3 > > * >(*it)) {
        cast->operation().setPointEvaluation(true);
        return;
      }
      if(auto cast = dynamic_cast< const function::FieldScalarTypeNode< function::EigenfunctionCompound< T_Scalar, field::Form::Form1, 3 > > * >(*it)) {
        cast->operation().setPointEvaluation(true);
        return;
      }
      if(auto cast = dynamic_cast< const function::FieldScalarTypeNode< function::EigenfunctionCompound< T_Scalar, field::Form::Form2, 3 > > * >(*it)) {
        cast->operation().setPointEvaluation(true);
        return;
      }
      if(auto cast = dynamic_cast< const function::FieldScalarTypeNode< function::EigenfunctionCompound< T_Scalar, field::Form::Form3, 3 > > * >(*it)) {
        cast->operation().setPointEvaluation(true);
        return;
      }
    }
  }

  void setPointEvaluation(ExecutionTreeIterator &it)
  {
    s_setPointEvaluation< std::complex< double > >(it);
    s_setPointEvaluation< std::complex< float > >(it);
    s_setPointEvaluation< double >(it);
    s_setPointEvaluation< float >(it);
  }


} // namespace gmshfem::function
