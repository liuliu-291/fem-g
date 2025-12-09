// GmshFEM - Copyright (C) 2019-2022, A. Royer, E. Béchet, C. Geuzaine, Université de Liège
//
// See the LICENSE.txt file for license information. Please report all
// issues on https://gitlab.onelab.info/gmsh/fem/issues

#ifndef H_GMSHFEM_INTEGRATE
#define H_GMSHFEM_INTEGRATE

#include "Domain.h"
#include "Function.h"
#include "MathObject.h"
#include "Post.h"
#include "scalar.h"

#include <string>

namespace gmshfem::post
{


  template< class T_Scalar, Degree T_Degree >
  class Integrate final : public PostInterface
  {
   private:
    typename MathObject< T_Scalar, T_Degree >::Object *const _value;
    const function::Function< T_Scalar, T_Degree > _function;
    domain::GeometricObject _domain;
    const std::string _integrationType;

   public:
    Integrate(typename MathObject< T_Scalar, T_Degree >::Object *const value, const function::Function< T_Scalar, T_Degree > &function, const domain::GeometricObject &domain, const std::string &integrationType);
    virtual ~Integrate();

    virtual int run() override;
  };


} // namespace gmshfem::post


#endif // H_GMSHFEM_INTEGRATE
