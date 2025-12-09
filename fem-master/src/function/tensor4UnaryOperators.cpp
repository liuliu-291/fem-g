// GmshFEM - Copyright (C) 2019-2022, A. Royer, E. Béchet, C. Geuzaine, Université de Liège
//
// See the LICENSE.txt file for license information. Please report all
// issues on https://gitlab.onelab.info/gmsh/fem/issues

#include "tensor4UnaryOperators.h"

#include "ChangeOfCoordinatesNode.h"
#include "UnaryNode.h"
#include "instantiate.h"
#include "unaryOperations.h"

namespace gmshfem::function
{


  template< class T_Scalar >
  Function< T_Scalar, Degree::Degree4 > operator-(const GeneralEvaluableObject< T_Scalar, Degree::Degree4 > &a)
  {
    return Function< T_Scalar, Degree::Degree4 >(new UnaryNode< Minus< T_Scalar, Degree::Degree4 > >(a.getEvaluableFunction().copy()));
  }

  INSTANTIATE_OPMM(TEMPLATE_RETURN(Function< TEMPLATE_PARAM_1, Degree::Degree4 >), , 0, 4, class, TEMPLATE_ARGS(std::complex< double >, std::complex< float >, double, float), TEMPLATE_PARAMS(const GeneralEvaluableObject< TEMPLATE_PARAM_1, Degree::Degree4 > &))


  template< class T_Scalar >
  Function< T_Scalar, Degree::Degree4 > operator+(const GeneralEvaluableObject< T_Scalar, Degree::Degree4 > &a)
  {
    return Function< T_Scalar, Degree::Degree4 >(new UnaryNode< Plus< T_Scalar, Degree::Degree4 > >(a.getEvaluableFunction().copy()));
  }

  INSTANTIATE_OPPP(TEMPLATE_RETURN(Function< TEMPLATE_PARAM_1, Degree::Degree4 >), , 0, 4, class, TEMPLATE_ARGS(std::complex< double >, std::complex< float >, double, float), TEMPLATE_PARAMS(const GeneralEvaluableObject< TEMPLATE_PARAM_1, Degree::Degree4 > &))
  
  
  template< class T_Scalar >
  Function< T_Scalar, Degree::Degree4 > changeOfCoordinates(const GeneralEvaluableObject< T_Scalar, Degree::Degree4 > &a, const GeneralEvaluableObject< scalar::Precision< T_Scalar >, Degree::Degree0 > &x, const GeneralEvaluableObject< scalar::Precision< T_Scalar >, Degree::Degree0 > &y, const GeneralEvaluableObject< scalar::Precision< T_Scalar >, Degree::Degree0 > &z)
  {
    return Function< T_Scalar, Degree::Degree4 >(new ChangeOfCoordinatesNode< T_Scalar, Degree::Degree4 >(x.getEvaluableFunction().copy(), y.getEvaluableFunction().copy(), z.getEvaluableFunction().copy(), a.getEvaluableFunction().copy()));
  }

  INSTANTIATE_FCT(TEMPLATE_RETURN(Function< TEMPLATE_PARAM_1, Degree::Degree4 >), , changeOfCoordinates, 0, 4, class, TEMPLATE_ARGS(std::complex< double >, std::complex< float >, double, float), TEMPLATE_PARAMS(const GeneralEvaluableObject< TEMPLATE_PARAM_1, Degree::Degree4 > &, const GeneralEvaluableObject< scalar::Precision< TEMPLATE_PARAM_1 >, Degree::Degree0 > &, const GeneralEvaluableObject< scalar::Precision< TEMPLATE_PARAM_1 >, Degree::Degree0 > &, const GeneralEvaluableObject< scalar::Precision< TEMPLATE_PARAM_1 >, Degree::Degree0 > &))


} // namespace gmshfem::function
