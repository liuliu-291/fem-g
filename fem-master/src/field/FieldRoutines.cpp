// GmshFEM - Copyright (C) 2019-2022, A. Royer, E. Béchet, C. Geuzaine, Université de Liège
//
// See the LICENSE.txt file for license information. Please report all
// issues on https://gitlab.onelab.info/gmsh/fem/issues

#include "FieldRoutines.h"

#include "AlgebraicFunctions.h"
#include "Dof.h"
#include "Exception.h"
#include "Matrix.h"
#include "Message.h"
#include "gmshfemDefines.h"
#include "instantiate.h"

namespace gmshfem::field
{


  template< class T_Scalar >
  Field< T_Scalar, Form::Form0 > multiplyAdd(const T_Scalar &m, const Field< T_Scalar, Form::Form0 > &a, const T_Scalar &n, const Field< T_Scalar, Form::Form0 > &b)
  {
    if(*a.getFunctionSpace() != *b.getFunctionSpace()) {
      throw common::Exception("Unable to apply 'multiplyAdd' to fields having different function spaces");
    }

    Field< T_Scalar, Form::Form0 > output = a;

    for(auto it = a.begin(); it != a.end(); ++it) {
      const dofs::Dof *dof = it->first;
      T_Scalar value = it->second;
      dofs::SearchDof sdof(dof->numType() - GMSHFEM_DOF_FIELD_OFFSET * a.tag() + GMSHFEM_DOF_FIELD_OFFSET * output.tag(), dof->entity());
      auto itFind = output.findDof(&sdof);
      itFind->second = m * value;
    }

    for(auto it = b.begin(); it != b.end(); ++it) {
      const dofs::Dof *dof = it->first;
      T_Scalar value = it->second;
      dofs::SearchDof sdof(dof->numType() - GMSHFEM_DOF_FIELD_OFFSET * b.tag() + GMSHFEM_DOF_FIELD_OFFSET * output.tag(), dof->entity());
      auto itFind = output.findDof(&sdof);
      if(itFind == output.end()) {
        output.setValue(new(output.getNextFixedDofMemoryPlace()) dofs::FixedDof(dof->numType() - GMSHFEM_DOF_FIELD_OFFSET * b.tag() + GMSHFEM_DOF_FIELD_OFFSET * output.tag(), dof->entity()), n * value);
      }
      else {
        itFind->second += n * value;
      }
    }

    output.name("multiplyAdd(" + a.name() + ", " + b.name() + ")");

    return output;
  }

  INSTANTIATE_FCT(TEMPLATE_RETURN(Field< TEMPLATE_PARAM_1, Form::Form0 >), , multiplyAdd, 0, 4, class, TEMPLATE_ARGS(std::complex< double >, std::complex< float >, double, float), TEMPLATE_PARAMS(const TEMPLATE_PARAM_1 &, const Field< TEMPLATE_PARAM_1, Form::Form0 > &, const TEMPLATE_PARAM_1 &, const Field< TEMPLATE_PARAM_1, Form::Form0 > &))


  template< class T_Scalar >
  void multiplyAddAssignement(const T_Scalar &m, Field< T_Scalar, Form::Form0 > &a, const T_Scalar &n, const Field< T_Scalar, Form::Form0 > &b)
  {
    if(*a.getFunctionSpace() != *b.getFunctionSpace()) {
      throw common::Exception("Unable to apply 'multiplyAddAssignement' to fields having different function spaces");
    }

    for(auto it = a.begin(); it != a.end(); ++it) {
      it->second *= m;
    }

    for(auto it = b.begin(); it != b.end(); ++it) {
      const dofs::Dof *dof = it->first;
      const T_Scalar value = it->second;
      dofs::SearchDof sdof(dof->numType() - GMSHFEM_DOF_FIELD_OFFSET * b.tag() + GMSHFEM_DOF_FIELD_OFFSET * a.tag(), dof->entity());
      auto itFind = a.findDof(&sdof);
      if(itFind == a.end()) {
        a.setValue(new(a.getNextFixedDofMemoryPlace()) dofs::FixedDof(dof->numType() - GMSHFEM_DOF_FIELD_OFFSET * b.tag() + GMSHFEM_DOF_FIELD_OFFSET * a.tag(), dof->entity()), n * value);
      }
      else {
        itFind->second += n * value;
      }
    }
  }

  INSTANTIATE_FCT(TEMPLATE_RETURN(void), , multiplyAddAssignement, 0, 4, class, TEMPLATE_ARGS(std::complex< double >, std::complex< float >, double, float), TEMPLATE_PARAMS(const TEMPLATE_PARAM_1 &, Field< TEMPLATE_PARAM_1, Form::Form0 > &, const TEMPLATE_PARAM_1 &, const Field< TEMPLATE_PARAM_1, Form::Form0 > &))


  template< class T_Scalar >
  Field< T_Scalar, Form::Form0 > conjField(const Field< T_Scalar, Form::Form0 > &a)
  {
    return a;
  }

  template< class T_ComplexScalar >
  Field< std::complex< T_ComplexScalar >, Form::Form0 > conjField(const Field< std::complex< T_ComplexScalar >, Form::Form0 > &a)
  {
    Field< std::complex< T_ComplexScalar >, Form::Form0 > output = a;
    for(auto it = output.begin(); it != output.end(); ++it) {
      it->second = std::conj(it->second);
    }
    return output;
  }

  INSTANTIATE_FCT(TEMPLATE_RETURN(Field< TEMPLATE_PARAM_1, Form::Form0 >), , conjField, 0, 4, class, TEMPLATE_ARGS(std::complex< double >, std::complex< float >, double, float), TEMPLATE_PARAMS(const Field< TEMPLATE_PARAM_1, Form::Form0 > &))


  template< class T_Scalar >
  Field< T_Scalar, Form::Form0 > realField(const Field< T_Scalar, Form::Form0 > &a)
  {
    return a;
  }

  template< class T_ComplexScalar >
  Field< std::complex< T_ComplexScalar >, Form::Form0 > realField(const Field< std::complex< T_ComplexScalar >, Form::Form0 > &a)
  {
    Field< std::complex< T_ComplexScalar >, Form::Form0 > output = a;
    for(auto it = output.begin(); it != output.end(); ++it) {
      it->second = ((std::complex< T_ComplexScalar >)std::real(it->second));
    }
    return output;
  }

  INSTANTIATE_FCT(TEMPLATE_RETURN(Field< TEMPLATE_PARAM_1, Form::Form0 >), , realField, 0, 4, class, TEMPLATE_ARGS(std::complex< double >, std::complex< float >, double, float), TEMPLATE_PARAMS(const Field< TEMPLATE_PARAM_1, Form::Form0 > &))


  template< class T_Scalar >
  T_Scalar bilinearProduct(const Field< T_Scalar, Form::Form0 > &a, const Field< T_Scalar, Form::Form0 > &b, const algebra::Matrix< T_Scalar > &M)
  {
    algebra::Vector< T_Scalar > x, y;
    a.getAllVector(x);
    b.getAllVector(y);

    return bilinear(x, M, y);
  }

  INSTANTIATE_FCT(TEMPLATE_RETURN(TEMPLATE_PARAM_1), , bilinearProduct, 0, 4, class, TEMPLATE_ARGS(std::complex< double >, std::complex< float >, double, float), TEMPLATE_PARAMS(const Field< TEMPLATE_PARAM_1, Form::Form0 > &, const Field< TEMPLATE_PARAM_1, Form::Form0 > &, const algebra::Matrix< TEMPLATE_PARAM_1 > &))


} // namespace gmshfem::field
