// GmshFEM - Copyright (C) 2019-2022, A. Royer, E. Béchet, C. Geuzaine, Université de Liège
//
// See the LICENSE.txt file for license information. Please report all
// issues on https://gitlab.onelab.info/gmsh/fem/issues

#include "FunctionSpaceGradH1.h"

#include "Exception.h"
#include "FieldInterface.h"
#include "Message.h"
#include "OmpInterface.h"
#include "instantiate.h"

#include <gmsh.h>


namespace gmshfem::field
{


  template< class T_PScalar >
  FunctionSpaceGradH1< T_PScalar >::FunctionSpaceGradH1(const FunctionSpace< T_PScalar, Form::Form0 > &fs) :
    FunctionSpace< T_PScalar, Form::Form1 >(), _fs(fs.copy())
  {
  }

  template< class T_PScalar >
  FunctionSpaceGradH1< T_PScalar >::~FunctionSpaceGradH1()
  {
    delete _fs;
  }

  template< class T_PScalar >
  std::string FunctionSpaceGradH1< T_PScalar >::_getGmshName(bool derivative) const
  {
    return _fs->_getGmshName(derivative);
  }

  template< class T_PScalar >
  std::string FunctionSpaceGradH1< T_PScalar >::getGmshFemName(bool derivative) const
  {
    return "Grad " + _fs->getGmshFemName(derivative);
  }

  template< class T_PScalar >
  std::string FunctionSpaceGradH1< T_PScalar >::getGmshFemOrientationName() const
  {
    return _fs->getGmshFemOrientationName();
  }

  template< class T_PScalar >
  FunctionSpaceTypeForm1 FunctionSpaceGradH1< T_PScalar >::type() const
  {
    switch(_fs->type()) {
    case FunctionSpaceTypeForm0::Lagrange:
      return FunctionSpaceTypeForm1::D_Lagrange;
      break;
    case FunctionSpaceTypeForm0::HierarchicalH1:
      return FunctionSpaceTypeForm1::D_HierarchicalH1;
      break;
    default:
      break;
    }
    return FunctionSpaceTypeForm1::D_Lagrange;
  }

  template< class T_PScalar >
  unsigned int FunctionSpaceGradH1< T_PScalar >::order() const
  {
    return _fs->order();
  }

  template< class T_PScalar >
  bool FunctionSpaceGradH1< T_PScalar >::isConstantByElements() const
  {
    return _fs->isConstantByElements();
  }

  template< class T_PScalar >
  bool FunctionSpaceGradH1< T_PScalar >::isFormP() const
  {
    return _fs->isFormP();
  }

  template< class T_PScalar >
  FunctionSpaceGradH1< T_PScalar > *FunctionSpaceGradH1< T_PScalar >::copy() const
  {
    return new FunctionSpaceGradH1< T_PScalar >(*_fs);
  }

  template< class T_PScalar >
  unsigned int FunctionSpaceGradH1< T_PScalar >::getBasisFunctions(const bool derivative, std::vector< T_PScalar > &functions, std::vector< int > &orientations, const std::vector< T_PScalar > &integrationPoints, const int elementType, const std::pair< int, int > &entity) const
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
  unsigned int FunctionSpaceGradH1< T_PScalar >::getBasisFunction_noexcept(const bool derivative, std::vector< T_PScalar > &functions, const std::vector< T_PScalar > &integrationPoints, const int elementType, const unsigned int orientation) const noexcept
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
      msg::error << "Unexpected error in 'FunctionSpaceGradH1::getBasisFunction_noexcept'" << msg::endl;
      exit(1);
    }
  }

  INSTANTIATE_CLASS(FunctionSpaceGradH1, 2, TEMPLATE_ARGS(double, float))


} // namespace gmshfem::field
