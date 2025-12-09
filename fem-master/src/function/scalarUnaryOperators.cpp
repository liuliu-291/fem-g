// GmshFEM - Copyright (C) 2019-2022, A. Royer, E. Béchet, C. Geuzaine, Université de Liège
//
// See the LICENSE.txt file for license information. Please report all
// issues on https://gitlab.onelab.info/gmsh/fem/issues

#include "scalarUnaryOperators.h"

#include "ChangeOfCoordinatesNode.h"
#include "UnaryNode.h"
#include "instantiate.h"
#include "unaryOperations.h"

namespace gmshfem::function
{


  template< class T_Scalar >
  Function< T_Scalar, Degree::Degree0 > operator-(const GeneralEvaluableObject< T_Scalar, Degree::Degree0 > &a)
  {
    return Function< T_Scalar, Degree::Degree0 >(new UnaryNode< Minus< T_Scalar, Degree::Degree0 > >(a.getEvaluableFunction().copy()));
  }

  INSTANTIATE_OPMM(TEMPLATE_RETURN(Function< TEMPLATE_PARAM_1, Degree::Degree0 >), , 0, 4, class, TEMPLATE_ARGS(std::complex< double >, std::complex< float >, double, float), TEMPLATE_PARAMS(const GeneralEvaluableObject< TEMPLATE_PARAM_1, Degree::Degree0 > &))


  template< class T_Scalar >
  Function< T_Scalar, Degree::Degree0 > operator+(const GeneralEvaluableObject< T_Scalar, Degree::Degree0 > &a)
  {
    return Function< T_Scalar, Degree::Degree0 >(new UnaryNode< Plus< T_Scalar, Degree::Degree0 > >(a.getEvaluableFunction().copy()));
  }

  INSTANTIATE_OPPP(TEMPLATE_RETURN(Function< TEMPLATE_PARAM_1, Degree::Degree0 >), , 0, 4, class, TEMPLATE_ARGS(std::complex< double >, std::complex< float >, double, float), TEMPLATE_PARAMS(const GeneralEvaluableObject< TEMPLATE_PARAM_1, Degree::Degree0 > &))


  template< class T_Scalar >
  Function< T_Scalar, Degree::Degree0 > abs(const GeneralEvaluableObject< T_Scalar, Degree::Degree0 > &a)
  {
    return Function< T_Scalar, Degree::Degree0 >(new UnaryNode< Abs< T_Scalar, Degree::Degree0 > >(a.getEvaluableFunction().copy()));
  }

  INSTANTIATE_FCT(TEMPLATE_RETURN(Function< TEMPLATE_PARAM_1, Degree::Degree0 >), , abs, 0, 4, class, TEMPLATE_ARGS(std::complex< double >, std::complex< float >, double, float), TEMPLATE_PARAMS(const GeneralEvaluableObject< TEMPLATE_PARAM_1, Degree::Degree0 > &))


  template< class T_Scalar >
  Function< T_Scalar, Degree::Degree0 > acos(const GeneralEvaluableObject< T_Scalar, Degree::Degree0 > &a)
  {
    return Function< T_Scalar, Degree::Degree0 >(new UnaryNode< Acos< T_Scalar, Degree::Degree0 > >(a.getEvaluableFunction().copy()));
  }

  INSTANTIATE_FCT(TEMPLATE_RETURN(Function< TEMPLATE_PARAM_1, Degree::Degree0 >), , acos, 0, 4, class, TEMPLATE_ARGS(std::complex< double >, std::complex< float >, double, float), TEMPLATE_PARAMS(const GeneralEvaluableObject< TEMPLATE_PARAM_1, Degree::Degree0 > &))


  template< class T_Scalar >
  Function< T_Scalar, Degree::Degree0 > acosh(const GeneralEvaluableObject< T_Scalar, Degree::Degree0 > &a)
  {
    return Function< T_Scalar, Degree::Degree0 >(new UnaryNode< Acosh< T_Scalar, Degree::Degree0 > >(a.getEvaluableFunction().copy()));
  }

  INSTANTIATE_FCT(TEMPLATE_RETURN(Function< TEMPLATE_PARAM_1, Degree::Degree0 >), , acosh, 0, 4, class, TEMPLATE_ARGS(std::complex< double >, std::complex< float >, double, float), TEMPLATE_PARAMS(const GeneralEvaluableObject< TEMPLATE_PARAM_1, Degree::Degree0 > &))


  template< class T_Scalar >
  Function< T_Scalar, Degree::Degree0 > angularComp(const GeneralEvaluableObject< T_Scalar, Degree::Degree1 > &a)
  {
    return Function< T_Scalar, Degree::Degree0 >(new UnaryNode< AngularComp< T_Scalar > >(a.getEvaluableFunction().copy()));
  }

  INSTANTIATE_FCT(TEMPLATE_RETURN(Function< TEMPLATE_PARAM_1, Degree::Degree0 >), , angularComp, 0, 4, class, TEMPLATE_ARGS(std::complex< double >, std::complex< float >, double, float), TEMPLATE_PARAMS(const GeneralEvaluableObject< TEMPLATE_PARAM_1, Degree::Degree1 > &))


  template< class T_Scalar >
  Function< T_Scalar, Degree::Degree0 > asin(const GeneralEvaluableObject< T_Scalar, Degree::Degree0 > &a)
  {
    return Function< T_Scalar, Degree::Degree0 >(new UnaryNode< Asin< T_Scalar, Degree::Degree0 > >(a.getEvaluableFunction().copy()));
  }

  INSTANTIATE_FCT(TEMPLATE_RETURN(Function< TEMPLATE_PARAM_1, Degree::Degree0 >), , asin, 0, 4, class, TEMPLATE_ARGS(std::complex< double >, std::complex< float >, double, float), TEMPLATE_PARAMS(const GeneralEvaluableObject< TEMPLATE_PARAM_1, Degree::Degree0 > &))


  template< class T_Scalar >
  Function< T_Scalar, Degree::Degree0 > asinh(const GeneralEvaluableObject< T_Scalar, Degree::Degree0 > &a)
  {
    return Function< T_Scalar, Degree::Degree0 >(new UnaryNode< Asinh< T_Scalar, Degree::Degree0 > >(a.getEvaluableFunction().copy()));
  }

  INSTANTIATE_FCT(TEMPLATE_RETURN(Function< TEMPLATE_PARAM_1, Degree::Degree0 >), , asinh, 0, 4, class, TEMPLATE_ARGS(std::complex< double >, std::complex< float >, double, float), TEMPLATE_PARAMS(const GeneralEvaluableObject< TEMPLATE_PARAM_1, Degree::Degree0 > &))


  template< class T_Scalar >
  Function< T_Scalar, Degree::Degree0 > atan(const GeneralEvaluableObject< T_Scalar, Degree::Degree0 > &a)
  {
    return Function< T_Scalar, Degree::Degree0 >(new UnaryNode< Atan< T_Scalar, Degree::Degree0 > >(a.getEvaluableFunction().copy()));
  }

  INSTANTIATE_FCT(TEMPLATE_RETURN(Function< TEMPLATE_PARAM_1, Degree::Degree0 >), , atan, 0, 4, class, TEMPLATE_ARGS(std::complex< double >, std::complex< float >, double, float), TEMPLATE_PARAMS(const GeneralEvaluableObject< TEMPLATE_PARAM_1, Degree::Degree0 > &))


  template< class T_Scalar >
  Function< T_Scalar, Degree::Degree0 > atanh(const GeneralEvaluableObject< T_Scalar, Degree::Degree0 > &a)
  {
    return Function< T_Scalar, Degree::Degree0 >(new UnaryNode< Atanh< T_Scalar, Degree::Degree0 > >(a.getEvaluableFunction().copy()));
  }

  INSTANTIATE_FCT(TEMPLATE_RETURN(Function< TEMPLATE_PARAM_1, Degree::Degree0 >), , atanh, 0, 4, class, TEMPLATE_ARGS(std::complex< double >, std::complex< float >, double, float), TEMPLATE_PARAMS(const GeneralEvaluableObject< TEMPLATE_PARAM_1, Degree::Degree0 > &))


  template< class T_Scalar >
  Function< T_Scalar, Degree::Degree0 > cbrt(const GeneralEvaluableObject< T_Scalar, Degree::Degree0 > &a)
  {
    return Function< T_Scalar, Degree::Degree0 >(new UnaryNode< Cbrt< T_Scalar, Degree::Degree0 > >(a.getEvaluableFunction().copy()));
  }

  INSTANTIATE_FCT(TEMPLATE_RETURN(Function< TEMPLATE_PARAM_1, Degree::Degree0 >), , cbrt, 0, 4, class, TEMPLATE_ARGS(std::complex< double >, std::complex< float >, double, float), TEMPLATE_PARAMS(const GeneralEvaluableObject< TEMPLATE_PARAM_1, Degree::Degree0 > &))


  template< class T_Scalar >
  Function< T_Scalar, Degree::Degree0 > changeOfCoordinates(const GeneralEvaluableObject< T_Scalar, Degree::Degree0 > &a, const GeneralEvaluableObject< scalar::Precision< T_Scalar >, Degree::Degree0 > &x, const GeneralEvaluableObject< scalar::Precision< T_Scalar >, Degree::Degree0 > &y, const GeneralEvaluableObject< scalar::Precision< T_Scalar >, Degree::Degree0 > &z)
  {
    return Function< T_Scalar, Degree::Degree0 >(new ChangeOfCoordinatesNode< T_Scalar, Degree::Degree0 >(x.getEvaluableFunction().copy(), y.getEvaluableFunction().copy(), z.getEvaluableFunction().copy(), a.getEvaluableFunction().copy()));
  }

  INSTANTIATE_FCT(TEMPLATE_RETURN(Function< TEMPLATE_PARAM_1, Degree::Degree0 >), , changeOfCoordinates, 0, 4, class, TEMPLATE_ARGS(std::complex< double >, std::complex< float >, double, float), TEMPLATE_PARAMS(const GeneralEvaluableObject< TEMPLATE_PARAM_1, Degree::Degree0 > &, const GeneralEvaluableObject< scalar::Precision< TEMPLATE_PARAM_1 >, Degree::Degree0 > &, const GeneralEvaluableObject< scalar::Precision< TEMPLATE_PARAM_1 >, Degree::Degree0 > &, const GeneralEvaluableObject< scalar::Precision< TEMPLATE_PARAM_1 >, Degree::Degree0 > &))


  template< class T_Scalar >
  Function< T_Scalar, Degree::Degree0 > commutator(const GeneralEvaluableObject< T_Scalar, Degree::Degree0 > &a, bool *state)
  {
    return Function< T_Scalar, Degree::Degree0 >(new UnaryNode< Commutator< T_Scalar, Degree::Degree0 > >(Commutator< T_Scalar, Degree::Degree0 >(state), a.getEvaluableFunction().copy()));
  }

  INSTANTIATE_FCT(TEMPLATE_RETURN(Function< TEMPLATE_PARAM_1, Degree::Degree0 >), , commutator, 0, 4, class, TEMPLATE_ARGS(std::complex< double >, std::complex< float >, double, float), TEMPLATE_PARAMS(const GeneralEvaluableObject< TEMPLATE_PARAM_1, Degree::Degree0 > &, bool *))


  template< class T_Scalar >
  Function< T_Scalar, Degree::Degree0 > conj(const GeneralEvaluableObject< T_Scalar, Degree::Degree0 > &a)
  {
    return Function< T_Scalar, Degree::Degree0 >(new UnaryNode< Conj< T_Scalar, Degree::Degree0 > >(a.getEvaluableFunction().copy()));
  }

  INSTANTIATE_FCT(TEMPLATE_RETURN(Function< TEMPLATE_PARAM_1, Degree::Degree0 >), , conj, 0, 4, class, TEMPLATE_ARGS(std::complex< double >, std::complex< float >, double, float), TEMPLATE_PARAMS(const GeneralEvaluableObject< TEMPLATE_PARAM_1, Degree::Degree0 > &))


  template< class T_Scalar >
  Function< T_Scalar, Degree::Degree0 > cos(const GeneralEvaluableObject< T_Scalar, Degree::Degree0 > &a)
  {
    return Function< T_Scalar, Degree::Degree0 >(new UnaryNode< Cos< T_Scalar, Degree::Degree0 > >(a.getEvaluableFunction().copy()));
  }

  INSTANTIATE_FCT(TEMPLATE_RETURN(Function< TEMPLATE_PARAM_1, Degree::Degree0 >), , cos, 0, 4, class, TEMPLATE_ARGS(std::complex< double >, std::complex< float >, double, float), TEMPLATE_PARAMS(const GeneralEvaluableObject< TEMPLATE_PARAM_1, Degree::Degree0 > &))


  template< class T_Scalar >
  Function< T_Scalar, Degree::Degree0 > cosh(const GeneralEvaluableObject< T_Scalar, Degree::Degree0 > &a)
  {
    return Function< T_Scalar, Degree::Degree0 >(new UnaryNode< Cosh< T_Scalar, Degree::Degree0 > >(a.getEvaluableFunction().copy()));
  }

  INSTANTIATE_FCT(TEMPLATE_RETURN(Function< TEMPLATE_PARAM_1, Degree::Degree0 >), , cosh, 0, 4, class, TEMPLATE_ARGS(std::complex< double >, std::complex< float >, double, float), TEMPLATE_PARAMS(const GeneralEvaluableObject< TEMPLATE_PARAM_1, Degree::Degree0 > &))


  template< class T_Scalar >
  Function< T_Scalar, Degree::Degree0 > cylBesselJ(const scalar::Precision< T_Scalar > v, const GeneralEvaluableObject< T_Scalar, Degree::Degree0 > &a)
  {
    return Function< T_Scalar, Degree::Degree0 >(new UnaryNode< CylBesselJ< T_Scalar > >(CylBesselJ< T_Scalar >(double(v)), a.getEvaluableFunction().copy()));
  }

  INSTANTIATE_FCT(TEMPLATE_RETURN(Function< TEMPLATE_PARAM_1, Degree::Degree0 >), , cylBesselJ, 0, 4, class, TEMPLATE_ARGS(std::complex< double >, std::complex< float >, double, float), TEMPLATE_PARAMS(const scalar::Precision< TEMPLATE_PARAM_1 >, const GeneralEvaluableObject< TEMPLATE_PARAM_1, Degree::Degree0 > &))


  template< class T_Scalar >
  Function< T_Scalar, Degree::Degree0 > cylNeumann(const scalar::Precision< T_Scalar > v, const GeneralEvaluableObject< T_Scalar, Degree::Degree0 > &a)
  {
    return Function< T_Scalar, Degree::Degree0 >(new UnaryNode< CylNeumann< T_Scalar > >(CylNeumann< T_Scalar >(double(v)), a.getEvaluableFunction().copy()));
  }

  INSTANTIATE_FCT(TEMPLATE_RETURN(Function< TEMPLATE_PARAM_1, Degree::Degree0 >), , cylNeumann, 0, 4, class, TEMPLATE_ARGS(std::complex< double >, std::complex< float >, double, float), TEMPLATE_PARAMS(const scalar::Precision< TEMPLATE_PARAM_1 >, const GeneralEvaluableObject< TEMPLATE_PARAM_1, Degree::Degree0 > &))


  template< class T_Scalar >
  Function< T_Scalar, Degree::Degree0 > exp(const GeneralEvaluableObject< T_Scalar, Degree::Degree0 > &a)
  {
    return Function< T_Scalar, Degree::Degree0 >(new UnaryNode< Exp< T_Scalar, Degree::Degree0 > >(a.getEvaluableFunction().copy()));
  }

  INSTANTIATE_FCT(TEMPLATE_RETURN(Function< TEMPLATE_PARAM_1, Degree::Degree0 >), , exp, 0, 4, class, TEMPLATE_ARGS(std::complex< double >, std::complex< float >, double, float), TEMPLATE_PARAMS(const GeneralEvaluableObject< TEMPLATE_PARAM_1, Degree::Degree0 > &))


  template< class T_Scalar >
  Function< T_Scalar, Degree::Degree0 > heaviside(const GeneralEvaluableObject< T_Scalar, Degree::Degree0 > &a)
  {
    return Function< T_Scalar, Degree::Degree0 >(new UnaryNode< Heaviside< T_Scalar, Degree::Degree0 > >(a.getEvaluableFunction().copy()));
  }

  INSTANTIATE_FCT(TEMPLATE_RETURN(Function< TEMPLATE_PARAM_1, Degree::Degree0 >), , heaviside, 0, 4, class, TEMPLATE_ARGS(std::complex< double >, std::complex< float >, double, float), TEMPLATE_PARAMS(const GeneralEvaluableObject< TEMPLATE_PARAM_1, Degree::Degree0 > &))


  template< class T_Scalar >
  Function< T_Scalar, Degree::Degree0 > imag(const GeneralEvaluableObject< T_Scalar, Degree::Degree0 > &a)
  {
    return Function< T_Scalar, Degree::Degree0 >(new UnaryNode< Imag< T_Scalar, Degree::Degree0 > >(a.getEvaluableFunction().copy()));
  }

  INSTANTIATE_FCT(TEMPLATE_RETURN(Function< TEMPLATE_PARAM_1, Degree::Degree0 >), , imag, 0, 4, class, TEMPLATE_ARGS(std::complex< double >, std::complex< float >, double, float), TEMPLATE_PARAMS(const GeneralEvaluableObject< TEMPLATE_PARAM_1, Degree::Degree0 > &))


  template< class T_Scalar >
  Function< T_Scalar, Degree::Degree0 > log(const GeneralEvaluableObject< T_Scalar, Degree::Degree0 > &a)
  {
    return Function< T_Scalar, Degree::Degree0 >(new UnaryNode< Log< T_Scalar, Degree::Degree0 > >(a.getEvaluableFunction().copy()));
  }

  INSTANTIATE_FCT(TEMPLATE_RETURN(Function< TEMPLATE_PARAM_1, Degree::Degree0 >), , log, 0, 4, class, TEMPLATE_ARGS(std::complex< double >, std::complex< float >, double, float), TEMPLATE_PARAMS(const GeneralEvaluableObject< TEMPLATE_PARAM_1, Degree::Degree0 > &))


  template< class T_Scalar >
  Function< T_Scalar, Degree::Degree0 > ln(const GeneralEvaluableObject< T_Scalar, Degree::Degree0 > &a)
  {
    return Function< T_Scalar, Degree::Degree0 >(new UnaryNode< Ln< T_Scalar, Degree::Degree0 > >(a.getEvaluableFunction().copy()));
  }

  INSTANTIATE_FCT(TEMPLATE_RETURN(Function< TEMPLATE_PARAM_1, Degree::Degree0 >), , ln, 0, 4, class, TEMPLATE_ARGS(std::complex< double >, std::complex< float >, double, float), TEMPLATE_PARAMS(const GeneralEvaluableObject< TEMPLATE_PARAM_1, Degree::Degree0 > &))


  template< class T_Scalar >
  Function< T_Scalar, Degree::Degree0 > norm(const GeneralEvaluableObject< T_Scalar, Degree::Degree0 > &a)
  {
    return Function< T_Scalar, Degree::Degree0 >(new UnaryNode< Norm< T_Scalar, Degree::Degree0 > >(Norm< T_Scalar, Degree::Degree0 >(1), a.getEvaluableFunction().copy()));
  }

  INSTANTIATE_FCT(TEMPLATE_RETURN(Function< TEMPLATE_PARAM_1, Degree::Degree0 >), , norm, 0, 4, class, TEMPLATE_ARGS(std::complex< double >, std::complex< float >, double, float), TEMPLATE_PARAMS(const GeneralEvaluableObject< TEMPLATE_PARAM_1, Degree::Degree0 > &))


  template< class T_Scalar >
  Function< T_Scalar, Degree::Degree0 > norm(const GeneralEvaluableObject< T_Scalar, Degree::Degree1 > &a, const int order)
  {
    return Function< T_Scalar, Degree::Degree0 >(new UnaryNode< Norm< T_Scalar, Degree::Degree1 > >(Norm< T_Scalar, Degree::Degree1 >(order), a.getEvaluableFunction().copy()));
  }

  INSTANTIATE_FCT(TEMPLATE_RETURN(Function< TEMPLATE_PARAM_1, Degree::Degree0 >), , norm, 1, 4, class, TEMPLATE_ARGS(std::complex< double >, std::complex< float >, double, float), TEMPLATE_PARAMS(const GeneralEvaluableObject< TEMPLATE_PARAM_1, Degree::Degree1 > &, const int))


  template< class T_Scalar >
  Function< T_Scalar, Degree::Degree0 > norm(const GeneralEvaluableObject< T_Scalar, Degree::Degree2 > &a, const int order)
  {
    return Function< T_Scalar, Degree::Degree0 >(new UnaryNode< Norm< T_Scalar, Degree::Degree2 > >(Norm< T_Scalar, Degree::Degree2 >(order), a.getEvaluableFunction().copy()));
  }

  INSTANTIATE_FCT(TEMPLATE_RETURN(Function< TEMPLATE_PARAM_1, Degree::Degree0 >), , norm, 2, 4, class, TEMPLATE_ARGS(std::complex< double >, std::complex< float >, double, float), TEMPLATE_PARAMS(const GeneralEvaluableObject< TEMPLATE_PARAM_1, Degree::Degree2 > &, const int))


  template< class T_Scalar >
  Function< T_Scalar, Degree::Degree0 > pow(const GeneralEvaluableObject< T_Scalar, Degree::Degree0 > &a, const unsigned int n)
  {
    return Function< T_Scalar, Degree::Degree0 >(new UnaryNode< Pow< T_Scalar, Degree::Degree0 > >(Pow< T_Scalar, Degree::Degree0 >(n), a.getEvaluableFunction().copy()));
  }

  INSTANTIATE_FCT(TEMPLATE_RETURN(Function< TEMPLATE_PARAM_1, Degree::Degree0 >), , pow, 0, 4, class, TEMPLATE_ARGS(std::complex< double >, std::complex< float >, double, float), TEMPLATE_PARAMS(const GeneralEvaluableObject< TEMPLATE_PARAM_1, Degree::Degree0 > &, const unsigned int))


  template< class T_Scalar >
  Function< T_Scalar, Degree::Degree0 > r2dComp(const GeneralEvaluableObject< T_Scalar, Degree::Degree1 > &a)
  {
    return Function< T_Scalar, Degree::Degree0 >(new UnaryNode< R2Dcomp< T_Scalar > >(a.getEvaluableFunction().copy()));
  }

  INSTANTIATE_FCT(TEMPLATE_RETURN(Function< TEMPLATE_PARAM_1, Degree::Degree0 >), , r2dComp, 0, 4, class, TEMPLATE_ARGS(std::complex< double >, std::complex< float >, double, float), TEMPLATE_PARAMS(const GeneralEvaluableObject< TEMPLATE_PARAM_1, Degree::Degree1 > &))


  template< class T_Scalar >
  Function< T_Scalar, Degree::Degree0 > real(const GeneralEvaluableObject< T_Scalar, Degree::Degree0 > &a)
  {
    return Function< T_Scalar, Degree::Degree0 >(new UnaryNode< Real< T_Scalar, Degree::Degree0 > >(a.getEvaluableFunction().copy()));
  }

  INSTANTIATE_FCT(TEMPLATE_RETURN(Function< TEMPLATE_PARAM_1, Degree::Degree0 >), , real, 0, 4, class, TEMPLATE_ARGS(std::complex< double >, std::complex< float >, double, float), TEMPLATE_PARAMS(const GeneralEvaluableObject< TEMPLATE_PARAM_1, Degree::Degree0 > &))


  template< class T_Scalar >
  Function< T_Scalar, Degree::Degree0 > sin(const GeneralEvaluableObject< T_Scalar, Degree::Degree0 > &a)
  {
    return Function< T_Scalar, Degree::Degree0 >(new UnaryNode< Sin< T_Scalar, Degree::Degree0 > >(a.getEvaluableFunction().copy()));
  }

  INSTANTIATE_FCT(TEMPLATE_RETURN(Function< TEMPLATE_PARAM_1, Degree::Degree0 >), , sin, 0, 4, class, TEMPLATE_ARGS(std::complex< double >, std::complex< float >, double, float), TEMPLATE_PARAMS(const GeneralEvaluableObject< TEMPLATE_PARAM_1, Degree::Degree0 > &))


  template< class T_Scalar >
  Function< T_Scalar, Degree::Degree0 > sinh(const GeneralEvaluableObject< T_Scalar, Degree::Degree0 > &a)
  {
    return Function< T_Scalar, Degree::Degree0 >(new UnaryNode< Sinh< T_Scalar, Degree::Degree0 > >(a.getEvaluableFunction().copy()));
  }

  INSTANTIATE_FCT(TEMPLATE_RETURN(Function< TEMPLATE_PARAM_1, Degree::Degree0 >), , sinh, 0, 4, class, TEMPLATE_ARGS(std::complex< double >, std::complex< float >, double, float), TEMPLATE_PARAMS(const GeneralEvaluableObject< TEMPLATE_PARAM_1, Degree::Degree0 > &))


  template< class T_Scalar >
  Function< T_Scalar, Degree::Degree0 > sqrt(const GeneralEvaluableObject< T_Scalar, Degree::Degree0 > &a)
  {
    return Function< T_Scalar, Degree::Degree0 >(new UnaryNode< Sqrt< T_Scalar, Degree::Degree0 > >(a.getEvaluableFunction().copy()));
  }

  INSTANTIATE_FCT(TEMPLATE_RETURN(Function< TEMPLATE_PARAM_1, Degree::Degree0 >), , sqrt, 0, 4, class, TEMPLATE_ARGS(std::complex< double >, std::complex< float >, double, float), TEMPLATE_PARAMS(const GeneralEvaluableObject< TEMPLATE_PARAM_1, Degree::Degree0 > &))


  template< class T_Scalar >
  Function< T_Scalar, Degree::Degree0 > tan(const GeneralEvaluableObject< T_Scalar, Degree::Degree0 > &a)
  {
    return Function< T_Scalar, Degree::Degree0 >(new UnaryNode< Tan< T_Scalar, Degree::Degree0 > >(a.getEvaluableFunction().copy()));
  }

  INSTANTIATE_FCT(TEMPLATE_RETURN(Function< TEMPLATE_PARAM_1, Degree::Degree0 >), , tan, 0, 4, class, TEMPLATE_ARGS(std::complex< double >, std::complex< float >, double, float), TEMPLATE_PARAMS(const GeneralEvaluableObject< TEMPLATE_PARAM_1, Degree::Degree0 > &))


  template< class T_Scalar >
  Function< T_Scalar, Degree::Degree0 > tanh(const GeneralEvaluableObject< T_Scalar, Degree::Degree0 > &a)
  {
    return Function< T_Scalar, Degree::Degree0 >(new UnaryNode< Tanh< T_Scalar, Degree::Degree0 > >(a.getEvaluableFunction().copy()));
  }

  INSTANTIATE_FCT(TEMPLATE_RETURN(Function< TEMPLATE_PARAM_1, Degree::Degree0 >), , tanh, 0, 4, class, TEMPLATE_ARGS(std::complex< double >, std::complex< float >, double, float), TEMPLATE_PARAMS(const GeneralEvaluableObject< TEMPLATE_PARAM_1, Degree::Degree0 > &))
  
  
  template< class T_Scalar >
  Function< T_Scalar, Degree::Degree0 > tr(const GeneralEvaluableObject< T_Scalar, Degree::Degree2 > &a)
  {
    return Function< T_Scalar, Degree::Degree0 >(new UnaryNode< Trace< T_Scalar > >(a.getEvaluableFunction().copy()));
  }

  INSTANTIATE_FCT(TEMPLATE_RETURN(Function< TEMPLATE_PARAM_1, Degree::Degree0 >), , tr, 0, 4, class, TEMPLATE_ARGS(std::complex< double >, std::complex< float >, double, float), TEMPLATE_PARAMS(const GeneralEvaluableObject< TEMPLATE_PARAM_1, Degree::Degree2 > &))


  template< class T_Scalar >
  Function< T_Scalar, Degree::Degree0 > xComp(const GeneralEvaluableObject< T_Scalar, Degree::Degree1 > &a)
  {
    return Function< T_Scalar, Degree::Degree0 >(new UnaryNode< Xcomp< T_Scalar > >(a.getEvaluableFunction().copy()));
  }

  INSTANTIATE_FCT(TEMPLATE_RETURN(Function< TEMPLATE_PARAM_1, Degree::Degree0 >), , xComp, 0, 4, class, TEMPLATE_ARGS(std::complex< double >, std::complex< float >, double, float), TEMPLATE_PARAMS(const GeneralEvaluableObject< TEMPLATE_PARAM_1, Degree::Degree1 > &))


  template< class T_Scalar >
  Function< T_Scalar, Degree::Degree0 > xxComp(const GeneralEvaluableObject< T_Scalar, Degree::Degree2 > &a)
  {
    return Function< T_Scalar, Degree::Degree0 >(new UnaryNode< XXcomp< T_Scalar > >(a.getEvaluableFunction().copy()));
  }

  INSTANTIATE_FCT(TEMPLATE_RETURN(Function< TEMPLATE_PARAM_1, Degree::Degree0 >), , xxComp, 0, 4, class, TEMPLATE_ARGS(std::complex< double >, std::complex< float >, double, float), TEMPLATE_PARAMS(const GeneralEvaluableObject< TEMPLATE_PARAM_1, Degree::Degree2 > &))


  template< class T_Scalar >
  Function< T_Scalar, Degree::Degree0 > xyComp(const GeneralEvaluableObject< T_Scalar, Degree::Degree2 > &a)
  {
    return Function< T_Scalar, Degree::Degree0 >(new UnaryNode< XYcomp< T_Scalar > >(a.getEvaluableFunction().copy()));
  }

  INSTANTIATE_FCT(TEMPLATE_RETURN(Function< TEMPLATE_PARAM_1, Degree::Degree0 >), , xyComp, 0, 4, class, TEMPLATE_ARGS(std::complex< double >, std::complex< float >, double, float), TEMPLATE_PARAMS(const GeneralEvaluableObject< TEMPLATE_PARAM_1, Degree::Degree2 > &))


  template< class T_Scalar >
  Function< T_Scalar, Degree::Degree0 > xzComp(const GeneralEvaluableObject< T_Scalar, Degree::Degree2 > &a)
  {
    return Function< T_Scalar, Degree::Degree0 >(new UnaryNode< XZcomp< T_Scalar > >(a.getEvaluableFunction().copy()));
  }

  INSTANTIATE_FCT(TEMPLATE_RETURN(Function< TEMPLATE_PARAM_1, Degree::Degree0 >), , xzComp, 0, 4, class, TEMPLATE_ARGS(std::complex< double >, std::complex< float >, double, float), TEMPLATE_PARAMS(const GeneralEvaluableObject< TEMPLATE_PARAM_1, Degree::Degree2 > &))


  template< class T_Scalar >
  Function< T_Scalar, Degree::Degree0 > yComp(const GeneralEvaluableObject< T_Scalar, Degree::Degree1 > &a)
  {
    return Function< T_Scalar, Degree::Degree0 >(new UnaryNode< Ycomp< T_Scalar > >(a.getEvaluableFunction().copy()));
  }

  INSTANTIATE_FCT(TEMPLATE_RETURN(Function< TEMPLATE_PARAM_1, Degree::Degree0 >), , yComp, 0, 4, class, TEMPLATE_ARGS(std::complex< double >, std::complex< float >, double, float), TEMPLATE_PARAMS(const GeneralEvaluableObject< TEMPLATE_PARAM_1, Degree::Degree1 > &))


  template< class T_Scalar >
  Function< T_Scalar, Degree::Degree0 > yxComp(const GeneralEvaluableObject< T_Scalar, Degree::Degree2 > &a)
  {
    return Function< T_Scalar, Degree::Degree0 >(new UnaryNode< YXcomp< T_Scalar > >(a.getEvaluableFunction().copy()));
  }

  INSTANTIATE_FCT(TEMPLATE_RETURN(Function< TEMPLATE_PARAM_1, Degree::Degree0 >), , yxComp, 0, 4, class, TEMPLATE_ARGS(std::complex< double >, std::complex< float >, double, float), TEMPLATE_PARAMS(const GeneralEvaluableObject< TEMPLATE_PARAM_1, Degree::Degree2 > &))


  template< class T_Scalar >
  Function< T_Scalar, Degree::Degree0 > yyComp(const GeneralEvaluableObject< T_Scalar, Degree::Degree2 > &a)
  {
    return Function< T_Scalar, Degree::Degree0 >(new UnaryNode< YYcomp< T_Scalar > >(a.getEvaluableFunction().copy()));
  }

  INSTANTIATE_FCT(TEMPLATE_RETURN(Function< TEMPLATE_PARAM_1, Degree::Degree0 >), , yyComp, 0, 4, class, TEMPLATE_ARGS(std::complex< double >, std::complex< float >, double, float), TEMPLATE_PARAMS(const GeneralEvaluableObject< TEMPLATE_PARAM_1, Degree::Degree2 > &))


  template< class T_Scalar >
  Function< T_Scalar, Degree::Degree0 > yzComp(const GeneralEvaluableObject< T_Scalar, Degree::Degree2 > &a)
  {
    return Function< T_Scalar, Degree::Degree0 >(new UnaryNode< YZcomp< T_Scalar > >(a.getEvaluableFunction().copy()));
  }

  INSTANTIATE_FCT(TEMPLATE_RETURN(Function< TEMPLATE_PARAM_1, Degree::Degree0 >), , yzComp, 0, 4, class, TEMPLATE_ARGS(std::complex< double >, std::complex< float >, double, float), TEMPLATE_PARAMS(const GeneralEvaluableObject< TEMPLATE_PARAM_1, Degree::Degree2 > &))


  template< class T_Scalar >
  Function< T_Scalar, Degree::Degree0 > zComp(const GeneralEvaluableObject< T_Scalar, Degree::Degree1 > &a)
  {
    return Function< T_Scalar, Degree::Degree0 >(new UnaryNode< Zcomp< T_Scalar > >(a.getEvaluableFunction().copy()));
  }

  INSTANTIATE_FCT(TEMPLATE_RETURN(Function< TEMPLATE_PARAM_1, Degree::Degree0 >), , zComp, 0, 4, class, TEMPLATE_ARGS(std::complex< double >, std::complex< float >, double, float), TEMPLATE_PARAMS(const GeneralEvaluableObject< TEMPLATE_PARAM_1, Degree::Degree1 > &))


  template< class T_Scalar >
  Function< T_Scalar, Degree::Degree0 > zxComp(const GeneralEvaluableObject< T_Scalar, Degree::Degree2 > &a)
  {
    return Function< T_Scalar, Degree::Degree0 >(new UnaryNode< ZXcomp< T_Scalar > >(a.getEvaluableFunction().copy()));
  }

  INSTANTIATE_FCT(TEMPLATE_RETURN(Function< TEMPLATE_PARAM_1, Degree::Degree0 >), , zxComp, 0, 4, class, TEMPLATE_ARGS(std::complex< double >, std::complex< float >, double, float), TEMPLATE_PARAMS(const GeneralEvaluableObject< TEMPLATE_PARAM_1, Degree::Degree2 > &))


  template< class T_Scalar >
  Function< T_Scalar, Degree::Degree0 > zyComp(const GeneralEvaluableObject< T_Scalar, Degree::Degree2 > &a)
  {
    return Function< T_Scalar, Degree::Degree0 >(new UnaryNode< ZYcomp< T_Scalar > >(a.getEvaluableFunction().copy()));
  }

  INSTANTIATE_FCT(TEMPLATE_RETURN(Function< TEMPLATE_PARAM_1, Degree::Degree0 >), , zyComp, 0, 4, class, TEMPLATE_ARGS(std::complex< double >, std::complex< float >, double, float), TEMPLATE_PARAMS(const GeneralEvaluableObject< TEMPLATE_PARAM_1, Degree::Degree2 > &))


  template< class T_Scalar >
  Function< T_Scalar, Degree::Degree0 > zzComp(const GeneralEvaluableObject< T_Scalar, Degree::Degree2 > &a)
  {
    return Function< T_Scalar, Degree::Degree0 >(new UnaryNode< ZZcomp< T_Scalar > >(a.getEvaluableFunction().copy()));
  }

  INSTANTIATE_FCT(TEMPLATE_RETURN(Function< TEMPLATE_PARAM_1, Degree::Degree0 >), , zzComp, 0, 4, class, TEMPLATE_ARGS(std::complex< double >, std::complex< float >, double, float), TEMPLATE_PARAMS(const GeneralEvaluableObject< TEMPLATE_PARAM_1, Degree::Degree2 > &))


} // namespace gmshfem::function
