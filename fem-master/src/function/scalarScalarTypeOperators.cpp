// GmshFEM - Copyright (C) 2019-2022, A. Royer, E. Béchet, C. Geuzaine, Université de Liège
//
// See the LICENSE.txt file for license information. Please report all
// issues on https://gitlab.onelab.info/gmsh/fem/issues

#include "scalarScalarTypeOperators.h"

#include "ScalarTypeNode.h"
#include "instantiate.h"
#include "scalarTypeOperations.h"

namespace gmshfem::function
{


  template< class T_Scalar >
  Function< scalar::ComplexPrecision< T_Scalar >, Degree::Degree0 > complex(const GeneralEvaluableObject< T_Scalar, Degree::Degree0 > &a)
  {
    return Function< scalar::ComplexPrecision< T_Scalar >, Degree::Degree0 >(new ScalarTypeNode< RealToComplex< T_Scalar, Degree::Degree0 > >(a.getEvaluableFunction().copy()));
  }

  INSTANTIATE_FCT(TEMPLATE_RETURN(Function< scalar::ComplexPrecision< TEMPLATE_PARAM_1 >, Degree::Degree0 >), , complex, 0, 4, class, TEMPLATE_ARGS(std::complex< double >, std::complex< float >, double, float), TEMPLATE_PARAMS(const GeneralEvaluableObject< TEMPLATE_PARAM_1, Degree::Degree0 > &))


  template< class T_Scalar >
  Function< scalar::Precision< T_Scalar >, Degree::Degree0 > realPart(const GeneralEvaluableObject< T_Scalar, Degree::Degree0 > &a)
  {
    return Function< scalar::Precision< T_Scalar >, Degree::Degree0 >(new ScalarTypeNode< ComplexRealPart< T_Scalar, Degree::Degree0 > >(a.getEvaluableFunction().copy()));
  }

  INSTANTIATE_FCT(TEMPLATE_RETURN(Function< scalar::Precision< TEMPLATE_PARAM_1 >, Degree::Degree0 >), , realPart, 0, 4, class, TEMPLATE_ARGS(std::complex< double >, std::complex< float >, double, float), TEMPLATE_PARAMS(const GeneralEvaluableObject< TEMPLATE_PARAM_1, Degree::Degree0 > &))


  template< class T_Scalar >
  Function< scalar::Precision< T_Scalar >, Degree::Degree0 > imagPart(const GeneralEvaluableObject< T_Scalar, Degree::Degree0 > &a)
  {
    return Function< scalar::Precision< T_Scalar >, Degree::Degree0 >(new ScalarTypeNode< ComplexImaginaryPart< T_Scalar, Degree::Degree0 > >(a.getEvaluableFunction().copy()));
  }

  INSTANTIATE_FCT(TEMPLATE_RETURN(Function< scalar::Precision< TEMPLATE_PARAM_1 >, Degree::Degree0 >), , imagPart, 0, 4, class, TEMPLATE_ARGS(std::complex< double >, std::complex< float >, double, float), TEMPLATE_PARAMS(const GeneralEvaluableObject< TEMPLATE_PARAM_1, Degree::Degree0 > &))


} // namespace gmshfem::function
