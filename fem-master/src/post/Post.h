// GmshFEM - Copyright (C) 2019-2022, A. Royer, E. Béchet, C. Geuzaine, Université de Liège
//
// See the LICENSE.txt file for license information. Please report all
// issues on https://gitlab.onelab.info/gmsh/fem/issues

#ifndef H_GMSHFEM_POST
#define H_GMSHFEM_POST

#include "FieldInterface.h"
#include "Function.h"
#include "GeneralEvaluableObject.h"
#include "Mesh.h"

namespace gmshfem::post
{


  class PostInterface
  {
   public:
    PostInterface();
    virtual ~PostInterface();

    virtual int run() = 0;
  };

  // Functions

  template< class T_Scalar, field::Form T_Form >
  int save(const field::Field< T_Scalar, T_Form > &field, const std::string &format = "msh", const std::string &path = "", bool memorize = false, const int step = 0, const double time = 0., const int partition = 0);

  template< class T_Scalar, Degree T_Degree >
  int save(const function::GeneralEvaluableObject< T_Scalar, T_Degree > &function, const domain::GeometricObject &domain, const std::string &name, const std::string &format = "msh", const std::string &path = "", bool memorize = false, const int step = 0, const double time = 0., const int partition = 0);
  template< class T_Scalar, Degree T_Degree >
  int save(const function::GeneralEvaluableObject< T_Scalar, T_Degree > &function, const Mesh &mesh, const std::string &name, const std::string &format = "msh", const std::string &path = "", bool memorize = false, const int step = 0, const double time = 0., const int partition = 0);

  template< class T_Scalar, Degree T_Degree >
  typename MathObject< T_Scalar, T_Degree >::Object integrate(const function::GeneralEvaluableObject< T_Scalar, T_Degree > &function, const domain::GeometricObject &domain, const std::string &integrationType);
  template< class T_Scalar, Degree T_Degree >
  typename MathObject< T_Scalar, T_Degree >::Object integrate(const function::GeneralEvaluableObject< T_Scalar, T_Degree > &function, const Mesh &mesh, const std::string &integrationType);

  template< class T_Scalar, Degree T_Degree >
  typename MathObject< T_Scalar, T_Degree >::Object evaluate(const function::GeneralEvaluableObject< T_Scalar, T_Degree > &function, const scalar::Precision< T_Scalar > &x, const scalar::Precision< T_Scalar > &y, const scalar::Precision< T_Scalar > &z);


} // namespace gmshfem::post

#include "PostproMap.h"

#endif // H_GMSHFEM_POST
