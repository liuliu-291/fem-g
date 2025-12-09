// GmshFEM - Copyright (C) 2019-2022, A. Royer, E. Béchet, C. Geuzaine, Université de Liège
//
// See the LICENSE.txt file for license information. Please report all
// issues on https://gitlab.onelab.info/gmsh/fem/issues

#ifndef H_GMSHFEM_POINTEVALUTATION
#define H_GMSHFEM_POINTEVALUTATION

#include "Function.h"
#include "MathObject.h"
#include "Post.h"
#include "scalar.h"

#include <string>

namespace gmshfem::post
{


  template< class T_Scalar, Degree T_Degree >
  class PointEvaluation final : public PostInterface
  {
   private:
    typename MathObject< T_Scalar, T_Degree >::Object *const _value;
    const function::Function< T_Scalar, T_Degree > _function;
    const scalar::Precision< T_Scalar > _x;
    const scalar::Precision< T_Scalar > _y;
    const scalar::Precision< T_Scalar > _z;

   public:
    PointEvaluation(typename MathObject< T_Scalar, T_Degree >::Object *const value, const function::Function< T_Scalar, T_Degree > &function, const scalar::Precision< T_Scalar > &x, const scalar::Precision< T_Scalar > &y, const scalar::Precision< T_Scalar > &z);
    virtual ~PointEvaluation();

    virtual int run() override;
  };


} // namespace gmshfem::post


#endif // H_GMSHFEM_POINTEVALUTATION
