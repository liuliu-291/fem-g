// GmshFEM - Copyright (C) 2019-2022, A. Royer, E. Béchet, C. Geuzaine, Université de Liège
//
// See the LICENSE.txt file for license information. Please report all
// issues on https://gitlab.onelab.info/gmsh/fem/issues

#include "scalarFieldOperators.h"

#include "FieldNode.h"
#include "ScalarTypeNode.h"
#include "fieldOperations.h"
#include "instantiate.h"
#include "scalarTypeOperations.h"

namespace gmshfem::function
{


  // ############################
  // Field
  // ############################

  template< class T_Scalar >
  Function< T_Scalar, Degree::Degree0 > div(const field::Field< T_Scalar, field::Form::Form2 > &field)
  {
    return Function< T_Scalar, Degree::Degree0 >(new FieldNode< Derivative< T_Scalar, field::Form::Form2 > >(&field));
  }

  INSTANTIATE_FCT(TEMPLATE_RETURN(Function< TEMPLATE_PARAM_1, Degree::Degree0 >), , div, 0, 4, class, TEMPLATE_ARGS(std::complex< double >, std::complex< float >, double, float), TEMPLATE_PARAMS(const field::Field< TEMPLATE_PARAM_1, field::Form::Form2 > &))


  template< class T_Scalar >
  Function< T_Scalar, Degree::Degree0 > d(const field::Field< T_Scalar, field::Form::Form2 > &field)
  {
    return Function< T_Scalar, Degree::Degree0 >(new FieldNode< Derivative< T_Scalar, field::Form::Form2 > >(&field));
  }

  INSTANTIATE_FCT(TEMPLATE_RETURN(Function< TEMPLATE_PARAM_1, Degree::Degree0 >), , d, 0, 4, class, TEMPLATE_ARGS(std::complex< double >, std::complex< float >, double, float), TEMPLATE_PARAMS(const field::Field< TEMPLATE_PARAM_1, field::Form::Form2 > &))


  template< class T_Scalar >
  Function< scalar::ComplexPrecision< T_Scalar >, Degree::Degree0 > eigenfunction(const field::Field< T_Scalar, field::Form::Form0 > &field, const unsigned int eigenTag)
  {
    return Function< scalar::ComplexPrecision< T_Scalar >, Degree::Degree0 >(new FieldScalarTypeNode< Eigenfunction< T_Scalar, field::Form::Form0 > >(Eigenfunction< T_Scalar, field::Form::Form0 >(eigenTag), &field));
  }

  INSTANTIATE_FCT(TEMPLATE_RETURN(Function< scalar::ComplexPrecision< TEMPLATE_PARAM_1 >, Degree::Degree0 >), , eigenfunction, 0, 4, class, TEMPLATE_ARGS(std::complex< double >, std::complex< float >, double, float), TEMPLATE_PARAMS(const field::Field< TEMPLATE_PARAM_1, field::Form::Form0 > &, const unsigned int))


  template< class T_Scalar >
  Function< scalar::ComplexPrecision< T_Scalar >, Degree::Degree0 > eigenfunction(const field::Field< T_Scalar, field::Form::Form3 > &field, const unsigned int eigenTag)
  {
    return Function< scalar::ComplexPrecision< T_Scalar >, Degree::Degree0 >(new FieldScalarTypeNode< Eigenfunction< T_Scalar, field::Form::Form3 > >(Eigenfunction< T_Scalar, field::Form::Form3 >(eigenTag), &field));
  }

  INSTANTIATE_FCT(TEMPLATE_RETURN(Function< scalar::ComplexPrecision< TEMPLATE_PARAM_1 >, Degree::Degree0 >), , eigenfunction, 1, 4, class, TEMPLATE_ARGS(std::complex< double >, std::complex< float >, double, float), TEMPLATE_PARAMS(const field::Field< TEMPLATE_PARAM_1, field::Form::Form3 > &, const unsigned int))


} // namespace gmshfem::function
