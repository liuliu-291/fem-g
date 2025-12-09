// GmshFEM - Copyright (C) 2019-2022, A. Royer, E. Béchet, C. Geuzaine, Université de Liège
//
// See the LICENSE.txt file for license information. Please report all
// issues on https://gitlab.onelab.info/gmsh/fem/issues

#include "tensorNovenaryOperators.h"

#include "NovenaryNode.h"
#include "NullaryNode.h"
#include "instantiate.h"
#include "novenaryOperations.h"
#include "nullaryOperations.h"

namespace gmshfem::function
{


  template< class T_Scalar >
  Function< T_Scalar, Degree::Degree2 > tensor(const GeneralEvaluableObject< T_Scalar, Degree::Degree0 > &a,
                                               const GeneralEvaluableObject< T_Scalar, Degree::Degree0 > &b,
                                               const GeneralEvaluableObject< T_Scalar, Degree::Degree0 > &c,
                                               const GeneralEvaluableObject< T_Scalar, Degree::Degree0 > &d,
                                               const GeneralEvaluableObject< T_Scalar, Degree::Degree0 > &e,
                                               const GeneralEvaluableObject< T_Scalar, Degree::Degree0 > &f,
                                               const GeneralEvaluableObject< T_Scalar, Degree::Degree0 > &g,
                                               const GeneralEvaluableObject< T_Scalar, Degree::Degree0 > &h,
                                               const GeneralEvaluableObject< T_Scalar, Degree::Degree0 > &i)
  {
    if(a.getEvaluableFunction().isTrivial() && b.getEvaluableFunction().isTrivial() && c.getEvaluableFunction().isTrivial() && d.getEvaluableFunction().isTrivial() && e.getEvaluableFunction().isTrivial() && f.getEvaluableFunction().isTrivial() && g.getEvaluableFunction().isTrivial() && h.getEvaluableFunction().isTrivial() && i.getEvaluableFunction().isTrivial()) {
      typename MathObject< T_Scalar, Degree::Degree2 >::Object mat;
      mat << a.value(), b.value(), c.value(),
        d.value(), e.value(), f.value(),
        g.value(), h.value(), i.value();
      return Function< T_Scalar, Degree::Degree2 >(new NullaryNode< Constant< T_Scalar, Degree::Degree2 > >(Constant< T_Scalar, Degree::Degree2 >(mat)));
    }
    else {
      return Function< T_Scalar, Degree::Degree2 >(new NovenaryNode< Tensor< T_Scalar > >(a.getEvaluableFunction().copy(), b.getEvaluableFunction().copy(), c.getEvaluableFunction().copy(),
                                                                                          d.getEvaluableFunction().copy(), e.getEvaluableFunction().copy(), f.getEvaluableFunction().copy(),
                                                                                          g.getEvaluableFunction().copy(), h.getEvaluableFunction().copy(), i.getEvaluableFunction().copy()));
    }
  }

  INSTANTIATE_FCT(TEMPLATE_RETURN(Function< TEMPLATE_PARAM_1, Degree::Degree2 >), , tensor, 0, 4, class, TEMPLATE_ARGS(std::complex< double >, std::complex< float >, double, float), TEMPLATE_PARAMS(const GeneralEvaluableObject< TEMPLATE_PARAM_1, Degree::Degree0 > &, const GeneralEvaluableObject< TEMPLATE_PARAM_1, Degree::Degree0 > &, const GeneralEvaluableObject< TEMPLATE_PARAM_1, Degree::Degree0 > &, const GeneralEvaluableObject< TEMPLATE_PARAM_1, Degree::Degree0 > &, const GeneralEvaluableObject< TEMPLATE_PARAM_1, Degree::Degree0 > &, const GeneralEvaluableObject< TEMPLATE_PARAM_1, Degree::Degree0 > &, const GeneralEvaluableObject< TEMPLATE_PARAM_1, Degree::Degree0 > &, const GeneralEvaluableObject< TEMPLATE_PARAM_1, Degree::Degree0 > &, const GeneralEvaluableObject< TEMPLATE_PARAM_1, Degree::Degree0 > &))


} // namespace gmshfem::function
