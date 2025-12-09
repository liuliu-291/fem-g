// GmshFEM - Copyright (C) 2019-2022, A. Royer, E. Béchet, C. Geuzaine, Université de Liège
//
// See the LICENSE.txt file for license information. Please report all
// issues on https://gitlab.onelab.info/gmsh/fem/issues

#include "vectorUnaryOperators.h"

#include "ChangeOfCoordinatesNode.h"
#include "UnaryNode.h"
#include "instantiate.h"
#include "unaryOperations.h"

namespace gmshfem::function
{


  template< class T_Scalar >
  Function< T_Scalar, Degree::Degree1 > operator-(const GeneralEvaluableObject< T_Scalar, Degree::Degree1 > &a)
  {
    return Function< T_Scalar, Degree::Degree1 >(new UnaryNode< Minus< T_Scalar, Degree::Degree1 > >(a.getEvaluableFunction().copy()));
  }

  INSTANTIATE_OPMM(TEMPLATE_RETURN(Function< TEMPLATE_PARAM_1, Degree::Degree1 >), , 0, 4, class, TEMPLATE_ARGS(std::complex< double >, std::complex< float >, double, float), TEMPLATE_PARAMS(const GeneralEvaluableObject< TEMPLATE_PARAM_1, Degree::Degree1 > &))


  template< class T_Scalar >
  Function< T_Scalar, Degree::Degree1 > operator+(const GeneralEvaluableObject< T_Scalar, Degree::Degree1 > &a)
  {
    return Function< T_Scalar, Degree::Degree1 >(new UnaryNode< Plus< T_Scalar, Degree::Degree1 > >(a.getEvaluableFunction().copy()));
  }

  INSTANTIATE_OPPP(TEMPLATE_RETURN(Function< TEMPLATE_PARAM_1, Degree::Degree1 >), , 0, 4, class, TEMPLATE_ARGS(std::complex< double >, std::complex< float >, double, float), TEMPLATE_PARAMS(const GeneralEvaluableObject< TEMPLATE_PARAM_1, Degree::Degree1 > &))


  template< class T_Scalar >
  Function< T_Scalar, Degree::Degree1 > changeOfCoordinates(const GeneralEvaluableObject< T_Scalar, Degree::Degree1 > &a, const GeneralEvaluableObject< scalar::Precision< T_Scalar >, Degree::Degree0 > &x, const GeneralEvaluableObject< scalar::Precision< T_Scalar >, Degree::Degree0 > &y, const GeneralEvaluableObject< scalar::Precision< T_Scalar >, Degree::Degree0 > &z)
  {
    return Function< T_Scalar, Degree::Degree1 >(new ChangeOfCoordinatesNode< T_Scalar, Degree::Degree1 >(x.getEvaluableFunction().copy(), y.getEvaluableFunction().copy(), z.getEvaluableFunction().copy(), a.getEvaluableFunction().copy()));
  }

  INSTANTIATE_FCT(TEMPLATE_RETURN(Function< TEMPLATE_PARAM_1, Degree::Degree1 >), , changeOfCoordinates, 0, 4, class, TEMPLATE_ARGS(std::complex< double >, std::complex< float >, double, float), TEMPLATE_PARAMS(const GeneralEvaluableObject< TEMPLATE_PARAM_1, Degree::Degree1 > &, const GeneralEvaluableObject< scalar::Precision< TEMPLATE_PARAM_1 >, Degree::Degree0 > &, const GeneralEvaluableObject< scalar::Precision< TEMPLATE_PARAM_1 >, Degree::Degree0 > &, const GeneralEvaluableObject< scalar::Precision< TEMPLATE_PARAM_1 >, Degree::Degree0 > &))


  template< class T_Scalar >
  Function< T_Scalar, Degree::Degree1 > commutator(const GeneralEvaluableObject< T_Scalar, Degree::Degree1 > &a, bool *state)
  {
    return Function< T_Scalar, Degree::Degree1 >(new UnaryNode< Commutator< T_Scalar, Degree::Degree1 > >(Commutator< T_Scalar, Degree::Degree1 >(state), a.getEvaluableFunction().copy()));
  }

  INSTANTIATE_FCT(TEMPLATE_RETURN(Function< TEMPLATE_PARAM_1, Degree::Degree1 >), , commutator, 0, 4, class, TEMPLATE_ARGS(std::complex< double >, std::complex< float >, double, float), TEMPLATE_PARAMS(const GeneralEvaluableObject< TEMPLATE_PARAM_1, Degree::Degree1 > &, bool *))


  template< class T_Scalar >
  Function< T_Scalar, Degree::Degree1 > conj(const GeneralEvaluableObject< T_Scalar, Degree::Degree1 > &a)
  {
    return Function< T_Scalar, Degree::Degree1 >(new UnaryNode< Conj< T_Scalar, Degree::Degree1 > >(a.getEvaluableFunction().copy()));
  }

  INSTANTIATE_FCT(TEMPLATE_RETURN(Function< TEMPLATE_PARAM_1, Degree::Degree1 >), , conj, 0, 4, class, TEMPLATE_ARGS(std::complex< double >, std::complex< float >, double, float), TEMPLATE_PARAMS(const GeneralEvaluableObject< TEMPLATE_PARAM_1, Degree::Degree1 > &))


} // namespace gmshfem::function
