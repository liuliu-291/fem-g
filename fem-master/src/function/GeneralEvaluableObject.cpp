// GmshFEM - Copyright (C) 2019-2022, A. Royer, E. Béchet, C. Geuzaine, Université de Liège
//
// See the LICENSE.txt file for license information. Please report all
// issues on https://gitlab.onelab.info/gmsh/fem/issues

#include "GeneralEvaluableObject.h"

#include "Function.h"
#include "instantiate.h"

namespace gmshfem::function
{


  template< class T_Scalar, Degree T_Degree >
  GeneralEvaluableObject< T_Scalar, T_Degree >::GeneralEvaluableObject() :
    _myself(nullptr), _value(), _isTrivial(false)
  {
  }

  template< class T_Scalar, Degree T_Degree >
  GeneralEvaluableObject< T_Scalar, T_Degree >::GeneralEvaluableObject(const typename MathObject< T_Scalar, T_Degree >::Object &value) :
    _myself(new Function< T_Scalar, T_Degree >(value)), _value(value), _isTrivial(true)
  {
  }

  template< class T_Scalar >
  static typename MathObject< T_Scalar, Degree::Degree0 >::Object toObject(const typename MathObject< scalar::Precision< T_Scalar >, Degree::Degree0 >::Object &value)
  {
    return T_Scalar(value);
  }

  template< class T_Scalar >
  static typename MathObject< T_Scalar, Degree::Degree1 >::Object toObject(const typename MathObject< scalar::Precision< T_Scalar >, Degree::Degree1 >::Object &value)
  {
    typename MathObject< T_Scalar, Degree::Degree1 >::Object vec;
    vec << T_Scalar(value(0)), T_Scalar(value(1)), T_Scalar(value(2));
    return vec;
  }

  template< class T_Scalar >
  static typename MathObject< T_Scalar, Degree::Degree2 >::Object toObject(const typename MathObject< scalar::Precision< T_Scalar >, Degree::Degree2 >::Object &value)
  {
    typename MathObject< T_Scalar, Degree::Degree2 >::Object mat;
    mat << T_Scalar(value(0)), T_Scalar(value(1)), T_Scalar(value(2)), T_Scalar(value(3)), T_Scalar(value(4)), T_Scalar(value(5)), T_Scalar(value(6)), T_Scalar(value(7)), T_Scalar(value(8));
    return mat;
  }

  template< class T_Scalar >
  static typename MathObject< T_Scalar, Degree::Degree4 >::Object toObject(const typename MathObject< scalar::Precision< T_Scalar >, Degree::Degree4 >::Object &value)
  {
    typename MathObject< T_Scalar, Degree::Degree4 >::Object mat;
    for(unsigned int i = 0; i < 3; ++i) {
      for(unsigned int j = 0; j < 3; ++j) {
        mat(i, j) << T_Scalar(value(i, j)(0)), T_Scalar(value(i, j)(1)), T_Scalar(value(i, j)(2)), T_Scalar(value(i, j)(3)), T_Scalar(value(i, j)(4)), T_Scalar(value(i, j)(5)), T_Scalar(value(i, j)(6)), T_Scalar(value(i, j)(7)), T_Scalar(value(i, j)(8));
      }
    }
    return mat;
  }

  template< class T_Scalar, Degree T_Degree >
  template< class T_SFINAE, class >
  GeneralEvaluableObject< T_Scalar, T_Degree >::GeneralEvaluableObject(const typename MathObject< scalar::Precision< T_SFINAE >, T_Degree >::Object &value) :
    _myself(new Function< T_Scalar, T_Degree >(toObject< T_SFINAE >(value))), _value(toObject< T_SFINAE >(value)), _isTrivial(true)
  {
  }

  template GeneralEvaluableObject< std::complex< double >, Degree::Degree0 >::GeneralEvaluableObject(const typename MathObject< scalar::Precision< std::complex< double > >, Degree::Degree0 >::Object &value);
  template GeneralEvaluableObject< std::complex< double >, Degree::Degree1 >::GeneralEvaluableObject(const typename MathObject< scalar::Precision< std::complex< double > >, Degree::Degree1 >::Object &value);
  template GeneralEvaluableObject< std::complex< double >, Degree::Degree2 >::GeneralEvaluableObject(const typename MathObject< scalar::Precision< std::complex< double > >, Degree::Degree2 >::Object &value);
  template GeneralEvaluableObject< std::complex< double >, Degree::Degree4 >::GeneralEvaluableObject(const typename MathObject< scalar::Precision< std::complex< double > >, Degree::Degree4 >::Object &value);

  template GeneralEvaluableObject< std::complex< float >, Degree::Degree0 >::GeneralEvaluableObject(const typename MathObject< scalar::Precision< std::complex< float > >, Degree::Degree0 >::Object &value);
  template GeneralEvaluableObject< std::complex< float >, Degree::Degree1 >::GeneralEvaluableObject(const typename MathObject< scalar::Precision< std::complex< float > >, Degree::Degree1 >::Object &value);
  template GeneralEvaluableObject< std::complex< float >, Degree::Degree2 >::GeneralEvaluableObject(const typename MathObject< scalar::Precision< std::complex< float > >, Degree::Degree2 >::Object &value);
  template GeneralEvaluableObject< std::complex< float >, Degree::Degree4 >::GeneralEvaluableObject(const typename MathObject< scalar::Precision< std::complex< float > >, Degree::Degree4 >::Object &value);

  template< class T_Scalar, Degree T_Degree >
  GeneralEvaluableObject< T_Scalar, T_Degree >::~GeneralEvaluableObject()
  {
    if(_myself != nullptr) {
      delete _myself;
    }
  }

  template< class T_Scalar, Degree T_Degree >
  const Function< T_Scalar, T_Degree > &GeneralEvaluableObject< T_Scalar, T_Degree >::getEvaluableFunction() const
  {
    return *_myself;
  }

  template< class T_Scalar, Degree T_Degree >
  bool GeneralEvaluableObject< T_Scalar, T_Degree >::isTrivial() const
  {
    return _isTrivial;
  }

  template< class T_Scalar, Degree T_Degree >
  typename MathObject< T_Scalar, T_Degree >::Object GeneralEvaluableObject< T_Scalar, T_Degree >::value() const
  {
    return _value;
  }

  INSTANTIATE_CLASS_2(GeneralEvaluableObject, 4, 4, TEMPLATE_ARGS(std::complex< double >, std::complex< float >, double, float), TEMPLATE_ARGS(Degree::Degree0, Degree::Degree1, Degree::Degree2, Degree::Degree4))


} // namespace gmshfem::function
