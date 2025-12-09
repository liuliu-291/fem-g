// GmshFEM - Copyright (C) 2019-2022, A. Royer, E. Béchet, C. Geuzaine, Université de Liège
//
// See the LICENSE.txt file for license information. Please report all
// issues on https://gitlab.onelab.info/gmsh/fem/issues

#include "tensorFieldOperators.h"

#include "FieldNode.h"
#include "fieldOperations.h"
#include "instantiate.h"

namespace gmshfem::function
{


  // ############################
  // CompoundField
  // ############################

  template< class T_Scalar, unsigned int T_NumFields >
  Function< T_Scalar, Degree::Degree2 > grad(const field::CompoundField< T_Scalar, field::Form::Form0, T_NumFields > &field)
  {
    return Function< T_Scalar, Degree::Degree2 >(new FieldNode< DerivativeCompound< T_Scalar, field::Form::Form0, T_NumFields > >(&field));
  }

  INSTANTIATE_FCT_2(TEMPLATE_RETURN(Function< TEMPLATE_PARAM_1, Degree::Degree2 >), , grad, 0, 4, class, TEMPLATE_ARGS(std::complex< double >, std::complex< float >, double, float), 2, unsigned int, TEMPLATE_ARGS(2, 3), TEMPLATE_PARAMS(const field::CompoundField< TEMPLATE_PARAM_1, field::Form::Form0, TEMPLATE_PARAM_2 > &))


  template< class T_Scalar, unsigned int T_NumFields >
  Function< T_Scalar, Degree::Degree2 > d(const field::CompoundField< T_Scalar, field::Form::Form0, T_NumFields > &field)
  {
    return Function< T_Scalar, Degree::Degree2 >(new FieldNode< DerivativeCompound< T_Scalar, field::Form::Form0, T_NumFields > >(&field));
  }

  INSTANTIATE_FCT_2(TEMPLATE_RETURN(Function< TEMPLATE_PARAM_1, Degree::Degree2 >), , d, 0, 4, class, TEMPLATE_ARGS(std::complex< double >, std::complex< float >, double, float), 2, unsigned int, TEMPLATE_ARGS(2, 3), TEMPLATE_PARAMS(const field::CompoundField< TEMPLATE_PARAM_1, field::Form::Form0, TEMPLATE_PARAM_2 > &))

  template< class T_Scalar, unsigned int T_NumFields >
  Function< T_Scalar, Degree::Degree2 > curl(const field::CompoundField< T_Scalar, field::Form::Form1, T_NumFields > &field)
  {
    return Function< T_Scalar, Degree::Degree2 >(new FieldNode< DerivativeCompound< T_Scalar, field::Form::Form1, T_NumFields > >(&field));
  }

  INSTANTIATE_FCT_2(TEMPLATE_RETURN(Function< TEMPLATE_PARAM_1, Degree::Degree2 >), , curl, 0, 4, class, TEMPLATE_ARGS(std::complex< double >, std::complex< float >, double, float), 2, unsigned int, TEMPLATE_ARGS(2, 3), TEMPLATE_PARAMS(const field::CompoundField< TEMPLATE_PARAM_1, field::Form::Form1, TEMPLATE_PARAM_2 > &))


  template< class T_Scalar, unsigned int T_NumFields >
  Function< T_Scalar, Degree::Degree2 > d(const field::CompoundField< T_Scalar, field::Form::Form1, T_NumFields > &field)
  {
    return Function< T_Scalar, Degree::Degree2 >(new FieldNode< DerivativeCompound< T_Scalar, field::Form::Form1, T_NumFields > >(&field));
  }

  INSTANTIATE_FCT_2(TEMPLATE_RETURN(Function< TEMPLATE_PARAM_1, Degree::Degree2 >), , d, 1, 4, class, TEMPLATE_ARGS(std::complex< double >, std::complex< float >, double, float), 2, unsigned int, TEMPLATE_ARGS(2, 3), TEMPLATE_PARAMS(const field::CompoundField< TEMPLATE_PARAM_1, field::Form::Form1, TEMPLATE_PARAM_2 > &))

  template< class T_Scalar, unsigned int T_NumFields >
  Function< scalar::ComplexPrecision< T_Scalar >, Degree::Degree2 > eigenfunction(const field::CompoundField< T_Scalar, field::Form::Form1, T_NumFields > &field, const unsigned int eigenTag)
  {
    return Function< scalar::ComplexPrecision< T_Scalar >, Degree::Degree2 >(new FieldScalarTypeNode< EigenfunctionCompound< T_Scalar, field::Form::Form1, T_NumFields > >(EigenfunctionCompound< T_Scalar, field::Form::Form1, T_NumFields >(eigenTag), &field));
  }

  INSTANTIATE_FCT_2(TEMPLATE_RETURN(Function< scalar::ComplexPrecision< TEMPLATE_PARAM_1 >, Degree::Degree2 >), , eigenfunction, 0, 4, class, TEMPLATE_ARGS(std::complex< double >, std::complex< float >, double, float), 2, unsigned int, TEMPLATE_ARGS(2, 3), TEMPLATE_PARAMS(const field::CompoundField< TEMPLATE_PARAM_1, field::Form::Form1, TEMPLATE_PARAM_2 > &, const unsigned int))


  template< class T_Scalar, unsigned int T_NumFields >
  Function< scalar::ComplexPrecision< T_Scalar >, Degree::Degree2 > eigenfunction(const field::CompoundField< T_Scalar, field::Form::Form2, T_NumFields > &field, const unsigned int eigenTag)
  {
    return Function< scalar::ComplexPrecision< T_Scalar >, Degree::Degree2 >(new FieldScalarTypeNode< EigenfunctionCompound< T_Scalar, field::Form::Form2, T_NumFields > >(EigenfunctionCompound< T_Scalar, field::Form::Form2, T_NumFields >(eigenTag), &field));
  }

  INSTANTIATE_FCT_2(TEMPLATE_RETURN(Function< scalar::ComplexPrecision< TEMPLATE_PARAM_1 >, Degree::Degree2 >), , eigenfunction, 1, 4, class, TEMPLATE_ARGS(std::complex< double >, std::complex< float >, double, float), 2, unsigned int, TEMPLATE_ARGS(2, 3), TEMPLATE_PARAMS(const field::CompoundField< TEMPLATE_PARAM_1, field::Form::Form2, TEMPLATE_PARAM_2 > &, const unsigned int))


} // namespace gmshfem::function
