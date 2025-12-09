// GmshFEM - Copyright (C) 2019-2022, A. Royer, E. Béchet, C. Geuzaine, Université de Liège
//
// See the LICENSE.txt file for license information. Please report all
// issues on https://gitlab.onelab.info/gmsh/fem/issues

#include "FunctionSpaceDivHDiv.h"

#include "Exception.h"
#include "FieldInterface.h"
#include "Message.h"
#include "OmpInterface.h"
#include "instantiate.h"

#include <gmsh.h>

namespace gmshfem::field
{


  template< class T_PScalar >
  FunctionSpaceDivHDiv< T_PScalar >::FunctionSpaceDivHDiv(const FunctionSpace< T_PScalar, Form::Form2 > &fs) :
    FunctionSpace< T_PScalar, Form::Form3 >(), _fs(fs.copy())
  {
  }

  template< class T_PScalar >
  FunctionSpaceDivHDiv< T_PScalar >::~FunctionSpaceDivHDiv()
  {
    delete _fs;
  }

  template< class T_PScalar >
  std::string FunctionSpaceDivHDiv< T_PScalar >::_getGmshName(bool derivative) const
  {
    return _fs->_getGmshName(derivative);
  }

  template< class T_PScalar >
  std::string FunctionSpaceDivHDiv< T_PScalar >::getGmshFemName(bool derivative) const
  {
    return "Div " + _fs->getGmshFemName(derivative);
  }

  template< class T_PScalar >
  std::string FunctionSpaceDivHDiv< T_PScalar >::getGmshFemOrientationName() const
  {
    return _fs->getGmshFemOrientationName();
  }

  template< class T_PScalar >
  FunctionSpaceTypeForm3 FunctionSpaceDivHDiv< T_PScalar >::type() const
  {
    switch(_fs->type()) {
    case FunctionSpaceTypeForm2::P_HierarchicalHCurl:
      return FunctionSpaceTypeForm3::DP_HierarchicalHCurl;
      break;
    default:
      break;
    }
    return FunctionSpaceTypeForm3::DP_HierarchicalHCurl;
  }

  template< class T_PScalar >
  unsigned int FunctionSpaceDivHDiv< T_PScalar >::order() const
  {
    return _fs->order();
  }

  template< class T_PScalar >
  bool FunctionSpaceDivHDiv< T_PScalar >::isConstantByElements() const
  {
    return _fs->isConstantByElements();
  }

  template< class T_PScalar >
  bool FunctionSpaceDivHDiv< T_PScalar >::isFormP() const
  {
    return _fs->isFormP();
  }

  template< class T_PScalar >
  FunctionSpaceDivHDiv< T_PScalar > *FunctionSpaceDivHDiv< T_PScalar >::copy() const
  {
    return new FunctionSpaceDivHDiv< T_PScalar >(*_fs);
  }

  template< class T_PScalar >
  unsigned int FunctionSpaceDivHDiv< T_PScalar >::getBasisFunctions(const bool derivative, std::vector< T_PScalar > &functions, std::vector< int > &orientations, const std::vector< T_PScalar > &integrationPoints, const int elementType, const std::pair< int, int > &entity) const
  {
    if(getDimOfElementType(elementType) < static_cast< unsigned int >(this->form()) - (isFormP() ? 1 : 0)) {
      throw common::Exception("You try to compute basis function on elements where function space is not defined");
    }

    const unsigned int nbrDofsByElement = _fs->getBasisFunctions(true, functions, orientations, integrationPoints, elementType, entity);
    if(derivative) {
      for(auto i = 0U; i < functions.size(); ++i) {
        functions[i] = 0.;
      }
    }

    return nbrDofsByElement;
  }

  template< class T_PScalar >
  unsigned int FunctionSpaceDivHDiv< T_PScalar >::getBasisFunction_noexcept(const bool derivative, std::vector< T_PScalar > &functions, const std::vector< T_PScalar > &integrationPoints, const int elementType, const unsigned int orientation) const noexcept
  {
    try {
      const unsigned int nbrDofsByElement = _fs->getBasisFunction_noexcept(true, functions, integrationPoints, elementType, orientation);
      if(derivative) {
        for(auto i = 0U; i < functions.size(); ++i) {
          functions[i] = 0.;
        }
      }
      return nbrDofsByElement;
    }
    catch(...) {
      msg::error << "Unexpected error in 'FunctionSpaceDivHDiv::getBasisFunction_noexcept'" << msg::endl;
      exit(1);
    }
  }

  INSTANTIATE_CLASS(FunctionSpaceDivHDiv, 2, TEMPLATE_ARGS(double, float))


} // namespace gmshfem::field
