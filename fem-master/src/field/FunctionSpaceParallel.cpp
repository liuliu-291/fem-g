// GmshFEM - Copyright (C) 2019-2022, A. Royer, E. Béchet, C. Geuzaine, Université de Liège
//
// See the LICENSE.txt file for license information. Please report all
// issues on https://gitlab.onelab.info/gmsh/fem/issues

#include "FunctionSpaceParallel.h"

#include "Exception.h"
#include "FieldInterface.h"
#include "OmpInterface.h"
#include "instantiate.h"

#include <gmsh.h>


namespace gmshfem::field
{


  template< class T_PScalar >
  FunctionSpaceParallel< T_PScalar >::FunctionSpaceParallel(const FunctionSpace< T_PScalar, Form::Form1 > &fs) :
    FunctionSpace< T_PScalar, Form::Form2 >(), _fs(fs.copy())
  {
  }

  template< class T_PScalar >
  FunctionSpaceParallel< T_PScalar >::~FunctionSpaceParallel()
  {
    delete _fs;
  }

  template< class T_PScalar >
  std::string FunctionSpaceParallel< T_PScalar >::_getGmshName(bool derivative) const
  {
    return _fs->_getGmshName(derivative);
  }

  template< class T_PScalar >
  std::string FunctionSpaceParallel< T_PScalar >::getGmshFemName(bool derivative) const
  {
    return "Parallel " + _fs->getGmshFemName(derivative);
  }

  template< class T_PScalar >
  std::string FunctionSpaceParallel< T_PScalar >::getGmshFemOrientationName() const
  {
    return _fs->getGmshFemOrientationName();
  }

  template< class T_PScalar >
  FunctionSpaceTypeForm2 FunctionSpaceParallel< T_PScalar >::type() const
  {
    switch(_fs->type()) {
    case FunctionSpaceTypeForm1::HierarchicalHCurl:
      return FunctionSpaceTypeForm2::P_HierarchicalHCurl;
      break;
    default:
      break;
    }
    return FunctionSpaceTypeForm2::P_HierarchicalHCurl;
  }

  template< class T_PScalar >
  unsigned int FunctionSpaceParallel< T_PScalar >::order() const
  {
    return _fs->order();
  }

  template< class T_PScalar >
  bool FunctionSpaceParallel< T_PScalar >::isConstantByElements() const
  {
    return _fs->isConstantByElements();
  }

  template< class T_PScalar >
  bool FunctionSpaceParallel< T_PScalar >::isFormP() const
  {
    return true;
  }

  template< class T_PScalar >
  FunctionSpaceParallel< T_PScalar > *FunctionSpaceParallel< T_PScalar >::copy() const
  {
    return new FunctionSpaceParallel< T_PScalar >(*_fs);
  }

  template< class T_PScalar >
  unsigned int FunctionSpaceParallel< T_PScalar >::getBasisFunctions(const bool derivative, std::vector< T_PScalar > &functions, std::vector< int > &orientations, const std::vector< T_PScalar > &integrationPoints, const int elementType, const std::pair< int, int > &entity) const
  {
    if(getDimOfElementType(elementType) < static_cast< unsigned int >(this->form()) - (isFormP() ? 1 : 0)) {
      throw common::Exception("Trying to compute basis functions on elements where no function space is defined");
    }

    unsigned int nbrDofsByElement = _fs->getBasisFunctions(derivative, functions, orientations, integrationPoints, elementType, entity);
    if(derivative) {
      const unsigned int size = functions.size() / 3;
      for(auto i = 0U; i < size; ++i) {
        functions[i] = -functions[3 * i + 2];
      }
      functions.resize(size);
      functions.shrink_to_fit();
    }
    else {
      const unsigned int size = functions.size() / 3;
      for(auto i = 0U; i < size; ++i) {
        const T_PScalar tmp = functions[3 * i + 0];
        functions[3 * i + 0] = -functions[3 * i + 1];
        functions[3 * i + 1] = tmp;
      }
    }
    return nbrDofsByElement;
  }

  template< class T_PScalar >
  unsigned int FunctionSpaceParallel< T_PScalar >::getBasisFunction_noexcept(const bool derivative, std::vector< T_PScalar > &functions, const std::vector< T_PScalar > &integrationPoints, const int elementType, const unsigned int orientation) const noexcept
  {
    try {
      unsigned int nbrDofsByElement = _fs->getBasisFunction_noexcept(derivative, functions, integrationPoints, elementType, orientation);
      if(derivative) {
        const unsigned int size = functions.size() / 3;
        for(auto i = 0U; i < size; ++i) {
          functions[i] = -functions[3 * i + 2];
        }
        functions.resize(size);
        functions.shrink_to_fit();
      }
      else {
        const unsigned int size = functions.size() / 3;
        for(auto i = 0U; i < size; ++i) {
          const T_PScalar tmp = functions[3 * i + 0];
          functions[3 * i + 0] = -functions[3 * i + 1];
          functions[3 * i + 1] = tmp;
        }
      }
      return nbrDofsByElement;
    }
    catch(...) {
      msg::error << "Unexpected error in 'FunctionSpaceParallel::getBasisFunction_noexcept'" << msg::endl;
      exit(1);
    }
  }

  INSTANTIATE_CLASS(FunctionSpaceParallel, 2, TEMPLATE_ARGS(double, float))


} // namespace gmshfem::field
