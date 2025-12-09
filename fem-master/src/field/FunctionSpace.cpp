// GmshFEM - Copyright (C) 2019-2022, A. Royer, E. Béchet, C. Geuzaine, Université de Liège
//
// See the LICENSE.txt file for license information. Please report all
// issues on https://gitlab.onelab.info/gmsh/fem/issues

#include "FunctionSpace.h"

#include "instantiate.h"
#include "FieldObject.h"

namespace gmshfem::field
{


  template< class T_PScalar, field::Form T_Form >
  FunctionSpace< T_PScalar, T_Form >::FunctionSpace() :
    FunctionSpaceInterface< T_PScalar >()
  {
  }

  template< class T_PScalar, field::Form T_Form >
  FunctionSpace< T_PScalar, T_Form >::~FunctionSpace()
  {
  }

  template< class T_PScalar, field::Form T_Form >
  field::Form FunctionSpace< T_PScalar, T_Form >::form() const
  {
    return T_Form;
  }

  template< class T_PScalar, field::Form T_Form >
  unsigned int FunctionSpace< T_PScalar, T_Form >::numberOfComponents(const bool derivative) const
  {
    if constexpr (T_Form == field::Form::Form0) {
      return (derivative ? 3 : 1);
    }
    else if  constexpr (T_Form == field::Form::Form1) {
      return (derivative ? 3 : 3);
    }
    else if  constexpr (T_Form == field::Form::Form2) {
      return (derivative ? 1 : 3);
    }
    else if  constexpr (T_Form == field::Form::Form3) {
      return (derivative ? 0 : 3);
    }
  }

  INSTANTIATE_CLASS_2(FunctionSpace, 2, 4, TEMPLATE_ARGS(double, float), TEMPLATE_ARGS(field::Form::Form0, field::Form::Form1, field::Form::Form2, field::Form::Form3))


} // namespace gmshfem::field
