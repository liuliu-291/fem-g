// GmshFEM - Copyright (C) 2019-2022, A. Royer, E. Béchet, C. Geuzaine, Université de Liège
//
// See the LICENSE.txt file for license information. Please report all
// issues on https://gitlab.onelab.info/gmsh/fem/issues

#include "GlobalQuantity.h"

#include "FieldInterface.h"
#include "Message.h"
#include "instantiate.h"

#include <gmsh.h>

namespace gmshfem::field
{


  static unsigned int s_nbrGlobalQuantity = 0;
  static std::vector< bool > s_listOfGlobalQuantity;

  template< class T_Scalar >
  GlobalQuantity< T_Scalar >::GlobalQuantity() :
    _name(), _model(), _tag(s_nbrGlobalQuantity++), _domain(), _primal(0.), _dual(0.), _fixedComponent(FixedComponent::None), _associatedPrimalField(nullptr), _associatedDualField(nullptr), _isActivated(false), _associatedPrimalDof(nullptr), _associatedDualDof(nullptr)
  {
    gmsh::model::getCurrent(_model);
    s_listOfGlobalQuantity.resize(s_nbrGlobalQuantity);
    s_listOfGlobalQuantity[_tag] = true;
  }

  template< class T_Scalar >
  GlobalQuantity< T_Scalar >::GlobalQuantity(const std::string &name, const domain::Domain &domain, const std::string &model) :
    _name(name), _model(model), _tag(s_nbrGlobalQuantity++), _domain(domain), _primal(0.), _dual(0.), _fixedComponent(FixedComponent::None), _associatedPrimalField(nullptr), _associatedDualField(nullptr), _isActivated(false), _associatedPrimalDof(nullptr), _associatedDualDof(nullptr)
  {
    if(_model == "") {
      gmsh::model::getCurrent(_model);
    }
    s_listOfGlobalQuantity.resize(s_nbrGlobalQuantity);
    s_listOfGlobalQuantity[_tag] = true;
  }

  template< class T_Scalar >
  GlobalQuantity< T_Scalar >::GlobalQuantity(const GlobalQuantity &other) :
    _name(other._name), _model(other._model), _tag(s_nbrGlobalQuantity++), _domain(other._domain), _primal(other._primal), _dual(other._dual), _fixedComponent(other._fixedComponent), _associatedPrimalField(nullptr), _associatedDualField(nullptr), _isActivated(other._isActivated), _associatedPrimalDof(nullptr), _associatedDualDof(nullptr)
  {
    s_listOfGlobalQuantity.resize(s_nbrGlobalQuantity);
    s_listOfGlobalQuantity[_tag] = true;
  }

  template< class T_Scalar >
  GlobalQuantity< T_Scalar >::~GlobalQuantity()
  {
    s_listOfGlobalQuantity[_tag] = false;
  }

  template< class T_Scalar >
  GlobalQuantity< T_Scalar > &GlobalQuantity< T_Scalar >::operator=(const GlobalQuantity &other)
  {
    _name = other._name;
    _model = other._model;
    _domain = other._domain;
    _primal = other._primal;
    _dual = other._dual;
    _fixedComponent = other._fixedComponent;
    _associatedPrimalField = other._associatedPrimalField;
    _associatedDualField = other._associatedDualField;
    _isActivated = other._isActivated;
    _associatedPrimalDof = other._associatedPrimalDof;
    _associatedDualDof = other._associatedDualDof;

    return *this;
  }

  template< class T_Scalar >
  unsigned int GlobalQuantity< T_Scalar >::tag() const
  {
    return _tag;
  }

  template< class T_Scalar >
  std::string GlobalQuantity< T_Scalar >::name() const
  {
    return _name;
  }

  template< class T_Scalar >
  std::string GlobalQuantity< T_Scalar >::model() const
  {
    return _model;
  }

  template< class T_Scalar >
  void GlobalQuantity< T_Scalar >::name(const std::string &name)
  {
    _name = name;
  }

  template< class T_Scalar >
  domain::Domain GlobalQuantity< T_Scalar >::domain() const
  {
    return _domain;
  }

  template< class T_Scalar >
  void GlobalQuantity< T_Scalar >::domain(const domain::Domain &domain)
  {
    _domain = domain;
  }

  template< class T_Scalar >
  T_Scalar GlobalQuantity< T_Scalar >::getPrimalValue() const
  {
    return _primal;
  }

  template< class T_Scalar >
  T_Scalar GlobalQuantity< T_Scalar >::getDualValue() const
  {
    return _dual;
  }

  template< class T_Scalar >
  void GlobalQuantity< T_Scalar >::_setPrimalValue(const T_Scalar &value)
  {
    _primal = value;
    _isActivated = true;
  }

  template< class T_Scalar >
  void GlobalQuantity< T_Scalar >::_setDualValue(const T_Scalar &value)
  {
    _dual = value;
    _isActivated = true;
  }

  template< class T_Scalar >
  void GlobalQuantity< T_Scalar >::_setFixedComponent(const FixedComponent fixedComponent)
  {
    _fixedComponent = fixedComponent;
  }

  template< class T_Scalar >
  FixedComponent GlobalQuantity< T_Scalar >::fixedComponent() const
  {
    return _fixedComponent;
  }

  template< class T_Scalar >
  bool GlobalQuantity< T_Scalar >::isActivated() const
  {
    return _isActivated;
  }

  template< class T_Scalar >
  void GlobalQuantity< T_Scalar >::setAssociatedPrimalDof(const dofs::Dof *dof)
  {
    _associatedPrimalDof = dof;
  }

  template< class T_Scalar >
  unsigned int GlobalQuantity< T_Scalar >::getNumPrimalDof() const
  {
    return _associatedPrimalDof->numDof() - 1;
  }

  template< class T_Scalar >
  void GlobalQuantity< T_Scalar >::setAssociatedDualDof(const dofs::Dof *dof)
  {
    _associatedDualDof = dof;
  }

  template< class T_Scalar >
  unsigned int GlobalQuantity< T_Scalar >::getNumDualDof() const
  {
    return _associatedDualDof->numDof() - 1;
  }

  template< class T_Scalar >
  void GlobalQuantity< T_Scalar >::setAssociatedDualField(FieldInterface< T_Scalar > *const field)
  {
    _associatedDualField = field;
  }

  template< class T_Scalar >
  FieldInterface< T_Scalar > *GlobalQuantity< T_Scalar >::getAssociatedDualField() const
  {
    return _associatedDualField;
  }

  template< class T_Scalar >
  void GlobalQuantity< T_Scalar >::setAssociatedPrimalField(FieldInterface< T_Scalar > *const field)
  {
    _associatedPrimalField = field;
  }

  template< class T_Scalar >
  FieldInterface< T_Scalar > *GlobalQuantity< T_Scalar >::getAssociatedPrimalField() const
  {
    return _associatedPrimalField;
  }

  INSTANTIATE_CLASS(GlobalQuantity, 4, TEMPLATE_ARGS(std::complex< double >, std::complex< float >, double, float))


  bool globalValueIsStillValid(const unsigned int tag)
  {
    if(tag >= s_listOfGlobalQuantity.size()) {
      return false;
    }

    return s_listOfGlobalQuantity[tag];
  }


} // namespace gmshfem::field
