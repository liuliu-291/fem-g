// GmshFEM - Copyright (C) 2019-2022, A. Royer, E. Béchet, C. Geuzaine, Université de Liège
//
// See the LICENSE.txt file for license information. Please report all
// issues on https://gitlab.onelab.info/gmsh/fem/issues

#include "tensorBinaryOperators.h"

#include "BinaryNode.h"
#include "binaryOperations.h"
#include "instantiate.h"

namespace gmshfem::function
{


  template< class T_Scalar >
  Function< T_Scalar, Degree::Degree2 > dyadic(const GeneralEvaluableObject< T_Scalar, Degree::Degree1 > &a, const GeneralEvaluableObject< T_Scalar, Degree::Degree1 > &b)
  {
    return Function< T_Scalar, Degree::Degree2 >(new BinaryNode< Dyadic< T_Scalar, Degree::Degree1 > >(a.getEvaluableFunction().copy(), b.getEvaluableFunction().copy()));
  }

  INSTANTIATE_FCT(TEMPLATE_RETURN(Function< TEMPLATE_PARAM_1, Degree::Degree2 >), , dyadic, 0, 4, class, TEMPLATE_ARGS(std::complex< double >, std::complex< float >, double, float), TEMPLATE_PARAMS(const GeneralEvaluableObject< TEMPLATE_PARAM_1, Degree::Degree1 > &, const GeneralEvaluableObject< TEMPLATE_PARAM_1, Degree::Degree1 > &))
  
  
  template< class T_Scalar >
  Function< T_Scalar, Degree::Degree2 > hadamard(const GeneralEvaluableObject< T_Scalar, Degree::Degree2 > &a, const GeneralEvaluableObject< T_Scalar, Degree::Degree2 > &b)
  {
    return Function< T_Scalar, Degree::Degree2 >(new BinaryNode< Hadamard< T_Scalar, Degree::Degree2 > >(a.getEvaluableFunction().copy(), b.getEvaluableFunction().copy()));
  }

  INSTANTIATE_FCT(TEMPLATE_RETURN(Function< TEMPLATE_PARAM_1, Degree::Degree2 >), , hadamard, 0, 4, class, TEMPLATE_ARGS(std::complex<double>, std::complex<float>, double, float), TEMPLATE_PARAMS(const GeneralEvaluableObject< TEMPLATE_PARAM_1, Degree::Degree2 > &, const GeneralEvaluableObject< TEMPLATE_PARAM_1, Degree::Degree2 > &))


} // namespace gmshfem::function
