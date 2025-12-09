// GmshFEM - Copyright (C) 2019-2022, A. Royer, E. Béchet, C. Geuzaine, Université de Liège
//
// See the LICENSE.txt file for license information. Please report all
// issues on https://gitlab.onelab.info/gmsh/fem/issues

#include "vectorFieldOperators.h"

#include "FieldNode.h"
#include "fieldOperations.h"
#include "instantiate.h"

namespace gmshfem::function
{


  // ############################
  // Field
  // ############################

  template< class T_Scalar >
  Function< T_Scalar, Degree::Degree1 > grad(const field::Field< T_Scalar, field::Form::Form0 > &field)
  {
    return Function< T_Scalar, Degree::Degree1 >(new FieldNode< Derivative< T_Scalar, field::Form::Form0 > >(&field));
  }

  INSTANTIATE_FCT(TEMPLATE_RETURN(Function< TEMPLATE_PARAM_1, Degree::Degree1 >), , grad, 0, 4, class, TEMPLATE_ARGS(std::complex< double >, std::complex< float >, double, float), TEMPLATE_PARAMS(const field::Field< TEMPLATE_PARAM_1, field::Form::Form0 > &))


  template< class T_Scalar >
  Function< T_Scalar, Degree::Degree1 > d(const field::Field< T_Scalar, field::Form::Form0 > &field)
  {
    return Function< T_Scalar, Degree::Degree1 >(new FieldNode< Derivative< T_Scalar, field::Form::Form0 > >(&field));
  }

  INSTANTIATE_FCT(TEMPLATE_RETURN(Function< TEMPLATE_PARAM_1, Degree::Degree1 >), , d, 0, 4, class, TEMPLATE_ARGS(std::complex< double >, std::complex< float >, double, float), TEMPLATE_PARAMS(const field::Field< TEMPLATE_PARAM_1, field::Form::Form0 > &))


  template< class T_Scalar >
  Function< T_Scalar, Degree::Degree1 > curl(const field::Field< T_Scalar, field::Form::Form1 > &field)
  {
    return Function< T_Scalar, Degree::Degree1 >(new FieldNode< Derivative< T_Scalar, field::Form::Form1 > >(&field));
  }

  INSTANTIATE_FCT(TEMPLATE_RETURN(Function< TEMPLATE_PARAM_1, Degree::Degree1 >), , curl, 0, 4, class, TEMPLATE_ARGS(std::complex< double >, std::complex< float >, double, float), TEMPLATE_PARAMS(const field::Field< TEMPLATE_PARAM_1, field::Form::Form1 > &))


  template< class T_Scalar >
  Function< T_Scalar, Degree::Degree1 > d(const field::Field< T_Scalar, field::Form::Form1 > &field)
  {
    return Function< T_Scalar, Degree::Degree1 >(new FieldNode< Derivative< T_Scalar, field::Form::Form1 > >(&field));
  }

  INSTANTIATE_FCT(TEMPLATE_RETURN(Function< TEMPLATE_PARAM_1, Degree::Degree1 >), , d, 1, 4, class, TEMPLATE_ARGS(std::complex< double >, std::complex< float >, double, float), TEMPLATE_PARAMS(const field::Field< TEMPLATE_PARAM_1, field::Form::Form1 > &))


  template< class T_Scalar >
  Function< scalar::ComplexPrecision< T_Scalar >, Degree::Degree1 > eigenfunction(const field::Field< T_Scalar, field::Form::Form1 > &field, const unsigned int eigenTag)
  {
    return Function< scalar::ComplexPrecision< T_Scalar >, Degree::Degree1 >(new FieldScalarTypeNode< Eigenfunction< T_Scalar, field::Form::Form1 > >(Eigenfunction< T_Scalar, field::Form::Form1 >(eigenTag), &field));
  }

  INSTANTIATE_FCT(TEMPLATE_RETURN(Function< scalar::ComplexPrecision< TEMPLATE_PARAM_1 >, Degree::Degree1 >), , eigenfunction, 0, 4, class, TEMPLATE_ARGS(std::complex< double >, std::complex< float >, double, float), TEMPLATE_PARAMS(const field::Field< TEMPLATE_PARAM_1, field::Form::Form1 > &, const unsigned int))


  template< class T_Scalar >
  Function< scalar::ComplexPrecision< T_Scalar >, Degree::Degree1 > eigenfunction(const field::Field< T_Scalar, field::Form::Form2 > &field, const unsigned int eigenTag)
  {
    return Function< scalar::ComplexPrecision< T_Scalar >, Degree::Degree1 >(new FieldScalarTypeNode< Eigenfunction< T_Scalar, field::Form::Form2 > >(Eigenfunction< T_Scalar, field::Form::Form2 >(eigenTag), &field));
  }

  INSTANTIATE_FCT(TEMPLATE_RETURN(Function< scalar::ComplexPrecision< TEMPLATE_PARAM_1 >, Degree::Degree1 >), , eigenfunction, 1, 4, class, TEMPLATE_ARGS(std::complex< double >, std::complex< float >, double, float), TEMPLATE_PARAMS(const field::Field< TEMPLATE_PARAM_1, field::Form::Form2 > &, const unsigned int))


  // ############################
  // CompoundField
  // ############################

  template< class T_Scalar, unsigned int T_NumFields >
  Function< T_Scalar, Degree::Degree1 > div(const field::CompoundField< T_Scalar, field::Form::Form2, T_NumFields > &field)
  {
    return Function< T_Scalar, Degree::Degree1 >(new FieldNode< DerivativeCompound< T_Scalar, field::Form::Form2, T_NumFields > >(&field));
  }

  INSTANTIATE_FCT_2(TEMPLATE_RETURN(Function< TEMPLATE_PARAM_1, Degree::Degree1 >), , div, 0, 4, class, TEMPLATE_ARGS(std::complex< double >, std::complex< float >, double, float), 2, unsigned int, TEMPLATE_ARGS(2, 3), TEMPLATE_PARAMS(const field::CompoundField< TEMPLATE_PARAM_1, field::Form::Form2, TEMPLATE_PARAM_2 > &))


  template< class T_Scalar, unsigned int T_NumFields >
  Function< T_Scalar, Degree::Degree1 > d(const field::CompoundField< T_Scalar, field::Form::Form2, T_NumFields > &field)
  {
    return Function< T_Scalar, Degree::Degree1 >(new FieldNode< DerivativeCompound< T_Scalar, field::Form::Form2, T_NumFields > >(&field));
  }

  INSTANTIATE_FCT_2(TEMPLATE_RETURN(Function< TEMPLATE_PARAM_1, Degree::Degree1 >), , d, 2, 4, class, TEMPLATE_ARGS(std::complex< double >, std::complex< float >, double, float), 2, unsigned int, TEMPLATE_ARGS(2, 3), TEMPLATE_PARAMS(const field::CompoundField< TEMPLATE_PARAM_1, field::Form::Form2, TEMPLATE_PARAM_2 > &))


  template< class T_Scalar, unsigned int T_NumFields >
  Function< scalar::ComplexPrecision< T_Scalar >, Degree::Degree1 > eigenfunction(const field::CompoundField< T_Scalar, field::Form::Form0, T_NumFields > &field, const unsigned int eigenTag)
  {
    return Function< scalar::ComplexPrecision< T_Scalar >, Degree::Degree1 >(new FieldScalarTypeNode< EigenfunctionCompound< T_Scalar, field::Form::Form0, T_NumFields > >(EigenfunctionCompound< T_Scalar, field::Form::Form0, T_NumFields >(eigenTag), &field));
  }

  INSTANTIATE_FCT_2(TEMPLATE_RETURN(Function< scalar::ComplexPrecision< TEMPLATE_PARAM_1 >, Degree::Degree1 >), , eigenfunction, 2, 4, class, TEMPLATE_ARGS(std::complex< double >, std::complex< float >, double, float), 2, unsigned int, TEMPLATE_ARGS(2, 3), TEMPLATE_PARAMS(const field::CompoundField< TEMPLATE_PARAM_1, field::Form::Form0, TEMPLATE_PARAM_2 > &, const unsigned int))


  template< class T_Scalar, unsigned int T_NumFields >
  Function< scalar::ComplexPrecision< T_Scalar >, Degree::Degree1 > eigenfunction(const field::CompoundField< T_Scalar, field::Form::Form3, T_NumFields > &field, const unsigned int eigenTag)
  {
    return Function< scalar::ComplexPrecision< T_Scalar >, Degree::Degree1 >(new FieldScalarTypeNode< EigenfunctionCompound< T_Scalar, field::Form::Form3, T_NumFields > >(EigenfunctionCompound< T_Scalar, field::Form::Form3, T_NumFields >(eigenTag), &field));
  }

  INSTANTIATE_FCT_2(TEMPLATE_RETURN(Function< scalar::ComplexPrecision< TEMPLATE_PARAM_1 >, Degree::Degree1 >), , eigenfunction, 3, 4, class, TEMPLATE_ARGS(std::complex< double >, std::complex< float >, double, float), 2, unsigned int, TEMPLATE_ARGS(2, 3), TEMPLATE_PARAMS(const field::CompoundField< TEMPLATE_PARAM_1, field::Form::Form3, TEMPLATE_PARAM_2 > &, const unsigned int))


} // namespace gmshfem::function
