// GmshFEM - Copyright (C) 2019-2022, A. Royer, E. Béchet, C. Geuzaine, Université de Liège
//
// See the LICENSE.txt file for license information. Please report all
// issues on https://gitlab.onelab.info/gmsh/fem/issues

#include "BilinearTermInterface.h"

#include "instantiate.h"

//
// class BilinearTermInterface : public Term
//

namespace gmshfem::term
{


  template< class T_Scalar >
  BilinearTermInterface< T_Scalar >::BilinearTermInterface(const domain::GeometricObject &domain, const std::string &integrationType, const field::FieldInterface< T_Scalar > *const fieldLhs, const field::FieldInterface< T_Scalar > *const fieldRhs, const ProductType productType) :
    Term< T_Scalar >(domain, integrationType, productType), _nbrOfDofsByElement{0, 0}, _field{fieldLhs, fieldRhs}
  {
  }

  template< class T_Scalar >
  BilinearTermInterface< T_Scalar >::~BilinearTermInterface()
  {
  }

  template< class T_Scalar >
  std::string BilinearTermInterface< T_Scalar >::termName() const
  {
    return "bilinear term";
  }

  template< class T_Scalar >
  unsigned int BilinearTermInterface< T_Scalar >::nbrDofsByElement(unsigned int field) const
  {
    return this->_field[field]->multiplicity() * _nbrOfDofsByElement[field];
  }

  template< class T_Scalar >
  bool BilinearTermInterface< T_Scalar >::isLinear() const
  {
    return false;
  }

  template< class T_Scalar >
  bool BilinearTermInterface< T_Scalar >::isBilinear() const
  {
    return true;
  }

  template< class T_Scalar >
  const field::FieldInterface< T_Scalar > *BilinearTermInterface< T_Scalar >::field(const unsigned int field) const
  {
    return _field[field];
  }

  INSTANTIATE_CLASS(BilinearTermInterface, 4, TEMPLATE_ARGS(std::complex< double >, std::complex< float >, double, float))


} // namespace gmshfem::term
