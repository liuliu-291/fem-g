// GmshFEM - Copyright (C) 2019-2022, A. Royer, E. Béchet, C. Geuzaine, Université de Liège
//
// See the LICENSE.txt file for license information. Please report all
// issues on https://gitlab.onelab.info/gmsh/fem/issues

#include "Term.h"

#include "ElementBucket.h"
#include "Exception.h"
#include "FieldInterface.h"
#include "IndiceBucket.h"
#include "MathObject.h"
#include "Message.h"
#include "OmpInterface.h"
#include "instantiate.h"

#include <utility>

namespace gmshfem::term
{


  static unsigned int nbrterm = 0;

  template< class T_Scalar >
  Term< T_Scalar >::Term(const domain::GeometricObject &domain, const std::string &integrationType, const ProductType productType) :
    _tag(nbrterm++), _domain(domain), _integrationType(integrationType), _productType(productType), _numberOfGaussPoints(0), _activated(true)
  {
  }

  template< class T_Scalar >
  Term< T_Scalar >::~Term()
  {
  }

  template< class T_Scalar >
  unsigned int Term< T_Scalar >::tag() const
  {
    return _tag;
  }

  template< class T_Scalar >
  std::string Term< T_Scalar >::integrationType() const
  {
    return _integrationType;
  }

  template< class T_Scalar >
  unsigned int Term< T_Scalar >::nbrGaussPoints() const
  {
    return _numberOfGaussPoints;
  }

  template< class T_Scalar >
  bool Term< T_Scalar >::isActivated() const
  {
    return _activated;
  }

  template< class T_Scalar >
  void Term< T_Scalar >::activate()
  {
    _activated = true;
  }

  template< class T_Scalar >
  void Term< T_Scalar >::deactivate()
  {
    _activated = false;
  }

  template< class T_Scalar >
  const domain::GeometricObject &Term< T_Scalar >::domain() const
  {
    return _domain;
  }

  INSTANTIATE_CLASS(Term, 4, TEMPLATE_ARGS(std::complex< double >, std::complex< float >, double, float))


} // namespace gmshfem::term
