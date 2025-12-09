// GmshFEM - Copyright (C) 2019-2022, A. Royer, E. Béchet, C. Geuzaine, Université de Liège
//
// See the LICENSE.txt file for license information. Please report all
// issues on https://gitlab.onelab.info/gmsh/fem/issues

#include "PointEvaluation.h"

#include "Message.h"
#include "instantiate.h"

#include <gmsh.h>

namespace gmshfem::post
{


  template< class T_Scalar, Degree T_Degree >
  PointEvaluation< T_Scalar, T_Degree >::PointEvaluation(typename MathObject< T_Scalar, T_Degree >::Object *const value, const function::Function< T_Scalar, T_Degree > &function, const scalar::Precision< T_Scalar > &x, const scalar::Precision< T_Scalar > &y, const scalar::Precision< T_Scalar > &z) :
    PostInterface(), _value(value), _function(function), _x(x), _y(y), _z(z)
  {
  }

  template< class T_Scalar, Degree T_Degree >
  PointEvaluation< T_Scalar, T_Degree >::~PointEvaluation()
  {
  }

  template< class T_Scalar, Degree T_Degree >
  int PointEvaluation< T_Scalar, T_Degree >::run()
  {
    // We need to manually activate the point evaluation algorithm because we want to use it even if we evaluate the field on its definition model.
    _function.setPointEvaluation();
    typename MathObject< T_Scalar, T_Degree >::Object value;
    _function.evaluate(value, _x, _y, _z, std::pair< int, int >(0, 0));
    *_value = std::move(value);
    return 0;
  }

  INSTANTIATE_CLASS_2(PointEvaluation, 4, 4, TEMPLATE_ARGS(std::complex< double >, std::complex< float >, double, float), TEMPLATE_ARGS(Degree::Degree0, Degree::Degree1, Degree::Degree2, Degree::Degree4))


} // namespace gmshfem::post
