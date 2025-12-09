// GmshFEM - Copyright (C) 2019-2022, A. Royer, E. Béchet, C. Geuzaine, Université de Liège
//
// See the LICENSE.txt file for license information. Please report all
// issues on https://gitlab.onelab.info/gmsh/fem/issues

#ifndef H_GMSHFEM_GENERALEVALUABLEOBJECT
#define H_GMSHFEM_GENERALEVALUABLEOBJECT

#include "MathObject.h"

namespace gmshfem::function
{
  template< class T_Scalar, Degree T_Degree >
  class Function;
}

namespace gmshfem::function
{


  template< class T_Scalar, Degree T_Degree >
  class GeneralEvaluableObject
  {
   private:
    Function< T_Scalar, T_Degree > *_myself;
    typename MathObject< T_Scalar, T_Degree >::Object _value;
    bool _isTrivial;

   public:
    GeneralEvaluableObject();
    GeneralEvaluableObject(const typename MathObject< T_Scalar, T_Degree >::Object &value);
    template< class T_SFINAE = T_Scalar, class = typename std::enable_if< scalar::IsComplex< T_SFINAE >::value >::type >
    GeneralEvaluableObject(const typename MathObject< scalar::Precision< T_SFINAE >, T_Degree >::Object &value);
    virtual ~GeneralEvaluableObject();

    virtual const Function< T_Scalar, T_Degree > &getEvaluableFunction() const;

    bool isTrivial() const;
    typename MathObject< T_Scalar, T_Degree >::Object value() const;
  };


} // namespace gmshfem::function


#endif // H_GMSHFEM_GENERALEVALUABLEOBJECT
