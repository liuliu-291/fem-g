// GmshFEM - Copyright (C) 2019-2022, A. Royer, E. Béchet, C. Geuzaine, Université de Liège
//
// See the LICENSE.txt file for license information. Please report all
// issues on https://gitlab.onelab.info/gmsh/fem/issues

#include "PiecewiseFunction.h"

#include "Message.h"
#include "MultaryNode.h"
#include "ScalarFunction.h"
#include "TensorFunction.h"
#include "VectorFunction.h"
#include "instantiate.h"

namespace gmshfem::function
{


  template< class T_Scalar, Degree T_Degree >
  PiecewiseFunction< T_Scalar, T_Degree >::PiecewiseFunction() :
    function::Function< T_Scalar, T_Degree >(new function::MultaryNode< T_Scalar, T_Degree >())
  {
  }

  template< class T_Scalar, Degree T_Degree >
  PiecewiseFunction< T_Scalar, T_Degree >::~PiecewiseFunction()
  {
  }

  template< class T_Scalar, Degree T_Degree >
  void PiecewiseFunction< T_Scalar, T_Degree >::addFunction(const Function< T_Scalar, T_Degree > &f, const domain::Domain &domain) const
  {
    static_cast< const function::MultaryNode< T_Scalar, T_Degree > * >(this->_tree)->addFunction(f.copy(), domain);
  }

  INSTANTIATE_CLASS_2(PiecewiseFunction, 4, 4, TEMPLATE_ARGS(std::complex< double >, std::complex< float >, double, float), TEMPLATE_ARGS(Degree::Degree0, Degree::Degree1, Degree::Degree2, Degree::Degree4))


} // namespace gmshfem::function
