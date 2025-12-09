// GmshFEM - Copyright (C) 2019-2022, A. Royer, E. Béchet, C. Geuzaine, Université de Liège
//
// See the LICENSE.txt file for license information. Please report all
// issues on https://gitlab.onelab.info/gmsh/fem/issues

#include "Post.h"

#include "Integrate.h"
#include "PointEvaluation.h"
#include "SaveField.h"
#include "SaveFunction.h"
#include "instantiate.h"

namespace gmshfem::post
{


  PostInterface::PostInterface()
  {
  }

  PostInterface::~PostInterface()
  {
  }

  // **********
  // Functions
  // **********

  template< class T_Scalar, field::Form T_Form >
  int save(const field::Field< T_Scalar, T_Form > &field, const std::string &format, const std::string &path, const bool memorize, const int step, const double time, const int partition)
  {
    post::SaveField< T_Scalar, T_Form > post(field, format, path, memorize, step, time, partition);
    return post.run();
  }

  INSTANTIATE_FCT_2(TEMPLATE_RETURN(int), , save, 0, 4, class, TEMPLATE_ARGS(std::complex< double >, std::complex< float >, double, float), 1, field::Form, TEMPLATE_ARGS(field::Form::Form0), TEMPLATE_PARAMS(const field::Field< TEMPLATE_PARAM_1, TEMPLATE_PARAM_2 > &, const std::string &, const std::string &, const bool, const int, const double, const int))


  template< class T_Scalar, Degree T_Degree >
  int save(const function::GeneralEvaluableObject< T_Scalar, T_Degree > &function, const domain::GeometricObject &domain, const std::string &name, const std::string &format, const std::string &path, const bool memorize, const int step, const double time, const int partition)
  {
    post::SaveFunction< T_Scalar, T_Degree > post(function.getEvaluableFunction(), domain, name, format, path, memorize, step, time, partition);
    return post.run();
  }

  INSTANTIATE_FCT_2(TEMPLATE_RETURN(int), , save, 1, 4, class, TEMPLATE_ARGS(std::complex< double >, std::complex< float >, double, float), 3, Degree, TEMPLATE_ARGS(Degree::Degree0, Degree::Degree1, Degree::Degree2), TEMPLATE_PARAMS(const function::GeneralEvaluableObject< TEMPLATE_PARAM_1, TEMPLATE_PARAM_2 > &, const domain::GeometricObject &, const std::string &, const std::string &, const std::string &, const bool, const int, const double, const int))


  template< class T_Scalar, Degree T_Degree >
  int save(const function::GeneralEvaluableObject< T_Scalar, T_Degree > &function, const Mesh &mesh, const std::string &name, const std::string &format, const std::string &path, const bool memorize, const int step, const double time, const int partition)
  {
    mesh.generateMesh();
    mesh.setCurrent();
    post::SaveFunction< T_Scalar, T_Degree > post(function.getEvaluableFunction(), mesh.getDomain(), name, format, path, memorize, step, time, partition);
    int viewTag = post.run();
    mesh.restoreCurrent();
    return viewTag;
  }

  INSTANTIATE_FCT_2(TEMPLATE_RETURN(int), , save, 2, 4, class, TEMPLATE_ARGS(std::complex< double >, std::complex< float >, double, float), 3, Degree, TEMPLATE_ARGS(Degree::Degree0, Degree::Degree1, Degree::Degree2), TEMPLATE_PARAMS(const function::GeneralEvaluableObject< TEMPLATE_PARAM_1, TEMPLATE_PARAM_2 > &, const Mesh &, const std::string &, const std::string &, const std::string &, const bool, const int, const double, const int))


  template< class T_Scalar, Degree T_Degree >
  typename MathObject< T_Scalar, T_Degree >::Object integrate(const function::GeneralEvaluableObject< T_Scalar, T_Degree > &function, const domain::GeometricObject &domain, const std::string &integrationType)
  {
    typename MathObject< T_Scalar, T_Degree >::Object postQuantity;
    post::Integrate< T_Scalar, T_Degree > post(&postQuantity, function.getEvaluableFunction(), domain, integrationType);
    post.run();
    return postQuantity;
  }

  INSTANTIATE_FCT_2(TEMPLATE_RETURN(typename MathObject< TEMPLATE_PARAM_1, TEMPLATE_PARAM_2 >::Object), , integrate, 0, 4, class, TEMPLATE_ARGS(std::complex< double >, std::complex< float >, double, float), 3, Degree, TEMPLATE_ARGS(Degree::Degree0, Degree::Degree1, Degree::Degree2), TEMPLATE_PARAMS(const function::GeneralEvaluableObject< TEMPLATE_PARAM_1, TEMPLATE_PARAM_2 > &, const domain::GeometricObject &, const std::string &))


  template< class T_Scalar, Degree T_Degree >
  typename MathObject< T_Scalar, T_Degree >::Object integrate(const function::GeneralEvaluableObject< T_Scalar, T_Degree > &function, const Mesh &mesh, const std::string &integrationType)
  {
    mesh.generateMesh();
    mesh.setCurrent();
    typename MathObject< T_Scalar, T_Degree >::Object postQuantity;
    post::Integrate< T_Scalar, T_Degree > post(&postQuantity, function.getEvaluableFunction(), mesh.getDomain(), integrationType);
    post.run();
    mesh.restoreCurrent();
    return postQuantity;
  }

  INSTANTIATE_FCT_2(TEMPLATE_RETURN(typename MathObject< TEMPLATE_PARAM_1, TEMPLATE_PARAM_2 >::Object), , integrate, 1, 4, class, TEMPLATE_ARGS(std::complex< double >, std::complex< float >, double, float), 3, Degree, TEMPLATE_ARGS(Degree::Degree0, Degree::Degree1, Degree::Degree2), TEMPLATE_PARAMS(const function::GeneralEvaluableObject< TEMPLATE_PARAM_1, TEMPLATE_PARAM_2 > &, const Mesh &, const std::string &))


  template< class T_Scalar, Degree T_Degree >
  typename MathObject< T_Scalar, T_Degree >::Object evaluate(const function::GeneralEvaluableObject< T_Scalar, T_Degree > &function, const scalar::Precision< T_Scalar > &x, const scalar::Precision< T_Scalar > &y, const scalar::Precision< T_Scalar > &z)
  {
    typename MathObject< T_Scalar, T_Degree >::Object postQuantity;
    post::PointEvaluation< T_Scalar, T_Degree > post(&postQuantity, function.getEvaluableFunction(), x, y, z);
    post.run();
    return postQuantity;
  }

  INSTANTIATE_FCT_2(TEMPLATE_RETURN(typename MathObject< TEMPLATE_PARAM_1, TEMPLATE_PARAM_2 >::Object), , evaluate, 0, 4, class, TEMPLATE_ARGS(std::complex< double >, std::complex< float >, double, float), 4, Degree, TEMPLATE_ARGS(Degree::Degree0, Degree::Degree1, Degree::Degree2, Degree::Degree4), TEMPLATE_PARAMS(const function::GeneralEvaluableObject< TEMPLATE_PARAM_1, TEMPLATE_PARAM_2 > &, const scalar::Precision< TEMPLATE_PARAM_1 > &, const scalar::Precision< TEMPLATE_PARAM_1 > &, const scalar::Precision< TEMPLATE_PARAM_1 > &))


} // namespace gmshfem::post
