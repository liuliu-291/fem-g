// GmshFEM - Copyright (C) 2019-2022, A. Royer, E. Béchet, C. Geuzaine, Université de Liège
//
// See the LICENSE.txt file for license information. Please report all
// issues on https://gitlab.onelab.info/gmsh/fem/issues

#include "tensor4NovenaryOperators.h"

#include "NovenaryNode.h"
#include "instantiate.h"
#include "novenaryOperations.h"

namespace gmshfem::function
{


  template< class T_Scalar, unsigned int T_Rank >
  Function< T_Scalar, Degree::Degree4 > tensor(const Function< T_Scalar, Degree::Degree2 > &a,
                                               const Function< T_Scalar, Degree::Degree2 > &b,
                                               const Function< T_Scalar, Degree::Degree2 > &c,
                                               const Function< T_Scalar, Degree::Degree2 > &d,
                                               const Function< T_Scalar, Degree::Degree2 > &e,
                                               const Function< T_Scalar, Degree::Degree2 > &f,
                                               const Function< T_Scalar, Degree::Degree2 > &g,
                                               const Function< T_Scalar, Degree::Degree2 > &h,
                                               const Function< T_Scalar, Degree::Degree2 > &i)
  {
    return Function< T_Scalar, Degree(TensorRank< T_Rank >::value) >(new NovenaryNode< Tensor4< T_Scalar > >(a.copy(), b.copy(), c.copy(),
                                                                                                             d.copy(), e.copy(), f.copy(),
                                                                                                             g.copy(), h.copy(), i.copy()));
  }

  INSTANTIATE_FCT_2(TEMPLATE_RETURN(Function< TEMPLATE_PARAM_1, Degree::Degree4 >), , tensor, 0, 4, class, TEMPLATE_ARGS(std::complex< double >, std::complex< float >, double, float), 1, unsigned int, TEMPLATE_ARGS(4), TEMPLATE_PARAMS(const Function< TEMPLATE_PARAM_1, Degree::Degree2 > &, const Function< TEMPLATE_PARAM_1, Degree::Degree2 > &, const Function< TEMPLATE_PARAM_1, Degree::Degree2 > &, const Function< TEMPLATE_PARAM_1, Degree::Degree2 > &, const Function< TEMPLATE_PARAM_1, Degree::Degree2 > &, const Function< TEMPLATE_PARAM_1, Degree::Degree2 > &, const Function< TEMPLATE_PARAM_1, Degree::Degree2 > &, const Function< TEMPLATE_PARAM_1, Degree::Degree2 > &, const Function< TEMPLATE_PARAM_1, Degree::Degree2 > &))


} // namespace gmshfem::function
