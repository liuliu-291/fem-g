// GmshFEM - Copyright (C) 2019-2022, A. Royer, E. Béchet, C. Geuzaine, Université de Liège
//
// See the LICENSE.txt file for license information. Please report all
// issues on https://gitlab.onelab.info/gmsh/fem/issues

#include "LinearTerm.h"
#include "instantiate.h"

//
// class LinearTerm : public Term
//

namespace gmshfem::term
{


  template< class T_Scalar >
  LinearTermInterface< T_Scalar >::LinearTermInterface(const domain::GeometricObject &domain, const std::string &integrationType, const field::FieldInterface< T_Scalar > *const fieldRhs, const ProductType productType, unsigned rhsIdx) :
    Term< T_Scalar >(domain, integrationType, productType), _nbrOfDofsByElement(0), _fieldRhs(fieldRhs), _rhsIdx(rhsIdx)
  {
  }

  template< class T_Scalar >
  LinearTermInterface< T_Scalar >::~LinearTermInterface()
  {
  }

  template< class T_Scalar >
  std::string LinearTermInterface< T_Scalar >::termName() const
  {
    return "linear term";
  }

  template< class T_Scalar >
  unsigned int LinearTermInterface< T_Scalar >::nbrDofsByElement(unsigned int field) const
  {
    return this->_fieldRhs->multiplicity() * _nbrOfDofsByElement;
  }

  template< class T_Scalar >
  bool LinearTermInterface< T_Scalar >::isLinear() const
  {
    return true;
  }

  template< class T_Scalar >
  bool LinearTermInterface< T_Scalar >::isBilinear() const
  {
    return false;
  }

  template< class T_Scalar >
  const field::FieldInterface< T_Scalar > *LinearTermInterface< T_Scalar >::field() const
  {
    return _fieldRhs;
  }

  INSTANTIATE_CLASS(LinearTermInterface, 4, TEMPLATE_ARGS(std::complex< double >, std::complex< float >, double, float))


} // namespace gmshfem::term
