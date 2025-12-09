// GmshFEM - Copyright (C) 2019-2022, A. Royer, E. Béchet, C. Geuzaine, Université de Liège
//
// See the LICENSE.txt file for license information. Please report all
// issues on https://gitlab.onelab.info/gmsh/fem/issues

#include "FunctionSpacePerpendicular.h"

#include "Exception.h"
#include "FieldInterface.h"
#include "OmpInterface.h"
#include "instantiate.h"

#include <gmsh.h>


namespace gmshfem::field
{


  template< class T_PScalar >
  FunctionSpacePerpendicular< T_PScalar >::FunctionSpacePerpendicular(const FunctionSpace< T_PScalar, Form::Form0 > &fs) :
    FunctionSpace< T_PScalar, Form::Form1 >(), _fs(fs.copy())
  {
  }

  template< class T_PScalar >
  FunctionSpacePerpendicular< T_PScalar >::~FunctionSpacePerpendicular()
  {
    delete _fs;
  }

  template< class T_PScalar >
  std::string FunctionSpacePerpendicular< T_PScalar >::_getGmshName(bool derivative) const
  {
    return _fs->_getGmshName(derivative);
  }

  template< class T_PScalar >
  std::string FunctionSpacePerpendicular< T_PScalar >::getGmshFemName(bool derivative) const
  {
    return "Perpendicular " + _fs->getGmshFemName(derivative);
  }

  template< class T_PScalar >
  std::string FunctionSpacePerpendicular< T_PScalar >::getGmshFemOrientationName() const
  {
    return _fs->getGmshFemOrientationName();
  }

  template< class T_PScalar >
  FunctionSpaceTypeForm1 FunctionSpacePerpendicular< T_PScalar >::type() const
  {
    switch(_fs->type()) {
    case FunctionSpaceTypeForm0::Lagrange:
      return FunctionSpaceTypeForm1::P_Lagrange;
      break;
    case FunctionSpaceTypeForm0::HierarchicalH1:
      return FunctionSpaceTypeForm1::P_HierarchicalH1;
      break;
    default:
      break;
    }
    return FunctionSpaceTypeForm1::P_Lagrange;
  }

  template< class T_PScalar >
  unsigned int FunctionSpacePerpendicular< T_PScalar >::order() const
  {
    return _fs->order();
  }

  template< class T_PScalar >
  bool FunctionSpacePerpendicular< T_PScalar >::isConstantByElements() const
  {
    return _fs->isConstantByElements();
  }

  template< class T_PScalar >
  bool FunctionSpacePerpendicular< T_PScalar >::isFormP() const
  {
    return true;
  }

  template< class T_PScalar >
  FunctionSpacePerpendicular< T_PScalar > *FunctionSpacePerpendicular< T_PScalar >::copy() const
  {
    return new FunctionSpacePerpendicular< T_PScalar >(*_fs);
  }

  template< class T_PScalar >
  unsigned int FunctionSpacePerpendicular< T_PScalar >::getBasisFunctions(const bool derivative, std::vector< T_PScalar > &functions, std::vector< int > &orientations, const std::vector< T_PScalar > &integrationPoints, const int elementType, const std::pair< int, int > &entity) const
  {
    if(getDimOfElementType(elementType) < static_cast< unsigned int >(this->form()) - (isFormP() ? 1 : 0)) {
      msg::info << isFormP() << msg::endl;
      throw common::Exception("Trying to compute basis functions on elements where no function space is defined");
    }

    unsigned int nbrDofsByElement = _fs->getBasisFunctions(derivative, functions, orientations, integrationPoints, elementType, entity);
    if(derivative) {
      const unsigned int size = functions.size() / 3;
      for(auto i = 0U; i < size; ++i) {
        const T_PScalar tmp = functions[3 * i + 0];
        functions[3 * i + 0] = functions[3 * i + 1];
        functions[3 * i + 1] = -tmp;
      }
    }
    else {
      const unsigned int size = functions.size();
      functions.resize(3 * size, 0.);
      for(auto i = size - 1; i >= 0 && i < size; --i) {
        functions[3 * i + 2] = functions[i];
        functions[i] = 0.;
      }
    }
    return nbrDofsByElement;
  }

  template< class T_PScalar >
  unsigned int FunctionSpacePerpendicular< T_PScalar >::getBasisFunction_noexcept(const bool derivative, std::vector< T_PScalar > &functions, const std::vector< T_PScalar > &integrationPoints, const int elementType, const unsigned int orientation) const noexcept
  {
    try {
      unsigned int nbrDofsByElement = _fs->getBasisFunction_noexcept(derivative, functions, integrationPoints, elementType, orientation);
      if(derivative) {
        const unsigned int size = functions.size() / 3;
        for(auto i = 0U; i < size; ++i) {
          const T_PScalar tmp = functions[3 * i + 0];
          functions[3 * i + 0] = functions[3 * i + 1];
          functions[3 * i + 1] = -tmp;
        }
      }
      else {
        const unsigned int size = functions.size();
        functions.resize(3 * size, 0.);
        for(auto i = size - 1; i >= 0 && i < size; --i) {
          functions[3 * i + 2] = functions[i];
          functions[i] = 0.;
        }
      }
      return nbrDofsByElement;
    }
    catch(...) {
      msg::error << "Unexpected error in 'FunctionSpacePerpendicular::getBasisFunction_noexcept'" << msg::endl;
      exit(1);
    }
  }

  INSTANTIATE_CLASS(FunctionSpacePerpendicular, 2, TEMPLATE_ARGS(double, float))


} // namespace gmshfem::field
