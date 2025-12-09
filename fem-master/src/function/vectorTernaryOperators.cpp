// GmshFEM - Copyright (C) 2019-2022, A. Royer, E. Béchet, C. Geuzaine, Université de Liège
//
// See the LICENSE.txt file for license information. Please report all
// issues on https://gitlab.onelab.info/gmsh/fem/issues

#include "vectorTernaryOperators.h"

#include "NullaryNode.h"
#include "TernaryNode.h"
#include "instantiate.h"
#include "nullaryOperations.h"
#include "ternaryOperations.h"

namespace gmshfem::function
{


  template< class T_Scalar >
  Function< T_Scalar, Degree::Degree1 > vector(const GeneralEvaluableObject< T_Scalar, Degree::Degree0 > &a, const GeneralEvaluableObject< T_Scalar, Degree::Degree0 > &b, const GeneralEvaluableObject< T_Scalar, Degree::Degree0 > &c)
  {
    if(a.getEvaluableFunction().isTrivial() && b.getEvaluableFunction().isTrivial() && c.getEvaluableFunction().isTrivial()) {
      typename MathObject< T_Scalar, Degree::Degree1 >::Object vec;
      vec << a.value(), b.value(), c.value();
      return Function< T_Scalar, Degree::Degree1 >(new NullaryNode< Constant< T_Scalar, Degree::Degree1 > >(Constant< T_Scalar, Degree::Degree1 >(vec)));
    }
    else {
      return Function< T_Scalar, Degree::Degree1 >(new TernaryNode< Vector< T_Scalar > >(a.getEvaluableFunction().copy(), b.getEvaluableFunction().copy(), c.getEvaluableFunction().copy()));
    }
  }

  INSTANTIATE_FCT(TEMPLATE_RETURN(Function< TEMPLATE_PARAM_1, Degree::Degree1 >), , vector, 0, 4, class, TEMPLATE_ARGS(std::complex< double >, std::complex< float >, double, float), TEMPLATE_PARAMS(const GeneralEvaluableObject< TEMPLATE_PARAM_1, Degree::Degree0 > &, const GeneralEvaluableObject< TEMPLATE_PARAM_1, Degree::Degree0 > &, const GeneralEvaluableObject< TEMPLATE_PARAM_1, Degree::Degree0 > &))


} // namespace gmshfem::function
