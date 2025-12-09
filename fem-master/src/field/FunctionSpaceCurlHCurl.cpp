// GmshFEM - Copyright (C) 2019-2022, A. Royer, E. Béchet, C. Geuzaine, Université de Liège
//
// See the LICENSE.txt file for license information. Please report all
// issues on https://gitlab.onelab.info/gmsh/fem/issues

#include "FunctionSpaceCurlHCurl.h"

#include "Exception.h"
#include "FieldInterface.h"
#include "Message.h"
#include "OmpInterface.h"
#include "instantiate.h"

#include <gmsh.h>

namespace gmshfem::field
{


  template< class T_PScalar >
  FunctionSpaceCurlHCurl< T_PScalar >::FunctionSpaceCurlHCurl(const FunctionSpace< T_PScalar, Form::Form1 > &fs) :
    FunctionSpace< T_PScalar, Form::Form2 >(), _fs(fs.copy())
  {
  }

  template< class T_PScalar >
  FunctionSpaceCurlHCurl< T_PScalar >::~FunctionSpaceCurlHCurl()
  {
    delete _fs;
  }

  template< class T_PScalar >
  std::string FunctionSpaceCurlHCurl< T_PScalar >::_getGmshName(bool derivative) const
  {
    return _fs->_getGmshName(derivative);
  }

  template< class T_PScalar >
  std::string FunctionSpaceCurlHCurl< T_PScalar >::getGmshFemName(bool derivative) const
  {
    return "Curl " + _fs->getGmshFemName(derivative);
  }

  template< class T_PScalar >
  std::string FunctionSpaceCurlHCurl< T_PScalar >::getGmshFemOrientationName() const
  {
    return _fs->getGmshFemOrientationName();
  }

  template< class T_PScalar >
  FunctionSpaceTypeForm2 FunctionSpaceCurlHCurl< T_PScalar >::type() const
  {
    switch(_fs->type()) {
    case FunctionSpaceTypeForm1::HierarchicalHCurl:
      return FunctionSpaceTypeForm2::D_HierarchicalHCurl;
      break;
    case FunctionSpaceTypeForm1::P_Lagrange:
      return FunctionSpaceTypeForm2::DP_Lagrange;
      break;
    case FunctionSpaceTypeForm1::P_HierarchicalH1:
      return FunctionSpaceTypeForm2::DP_HierarchicalH1;
      break;
    default:
      break;
    }
    return FunctionSpaceTypeForm2::D_HierarchicalHCurl;
  }

  template< class T_PScalar >
  unsigned int FunctionSpaceCurlHCurl< T_PScalar >::order() const
  {
    return _fs->order();
  }

  template< class T_PScalar >
  bool FunctionSpaceCurlHCurl< T_PScalar >::isConstantByElements() const
  {
    return _fs->isConstantByElements();
  }

  template< class T_PScalar >
  bool FunctionSpaceCurlHCurl< T_PScalar >::isFormP() const
  {
    return _fs->isFormP();
  }

  template< class T_PScalar >
  FunctionSpaceCurlHCurl< T_PScalar > *FunctionSpaceCurlHCurl< T_PScalar >::copy() const
  {
    return new FunctionSpaceCurlHCurl< T_PScalar >(*_fs);
  }

  template< class T_PScalar >
  unsigned int FunctionSpaceCurlHCurl< T_PScalar >::getBasisFunctions(const bool derivative, std::vector< T_PScalar > &functions, std::vector< int > &orientations, const std::vector< T_PScalar > &integrationPoints, const int elementType, const std::pair< int, int > &entity) const
  {
    if(getDimOfElementType(elementType) < static_cast< unsigned int >(this->form()) - (isFormP() ? 1 : 0)) {
      throw common::Exception("You try to compute basis function on elements where function space is not defined");
    }

    const unsigned int nbrDofsByElement = _fs->getBasisFunctions(true, functions, orientations, integrationPoints, elementType, entity);
    if(derivative) {
      functions.resize(functions.size() / 3);
      for(auto i = 0U; i < functions.size(); ++i) {
        functions[i] = 0.;
      }
    }

    return nbrDofsByElement;
  }

  template< class T_PScalar >
  unsigned int FunctionSpaceCurlHCurl< T_PScalar >::getBasisFunction_noexcept(const bool derivative, std::vector< T_PScalar > &functions, const std::vector< T_PScalar > &integrationPoints, const int elementType, const unsigned int orientation) const noexcept
  {
    try {
      const unsigned int nbrDofsByElement = _fs->getBasisFunction_noexcept(true, functions, integrationPoints, elementType, orientation);
      if(derivative) {
        functions.resize(functions.size() / 3);
        for(auto i = 0U; i < functions.size(); ++i) {
          functions[i] = 0.;
        }
      }
      return nbrDofsByElement;
    }
    catch(...) {
      msg::error << "Unexpected error in 'FunctionSpaceCurlHCurl::getBasisFunction_noexcept'" << msg::endl;
      exit(1);
    }
  }

  INSTANTIATE_CLASS(FunctionSpaceCurlHCurl, 2, TEMPLATE_ARGS(double, float))


} // namespace gmshfem::field
