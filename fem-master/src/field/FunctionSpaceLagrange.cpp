// GmshFEM - Copyright (C) 2019-2022, A. Royer, E. Béchet, C. Geuzaine, Université de Liège
//
// See the LICENSE.txt file for license information. Please report all
// issues on https://gitlab.onelab.info/gmsh/fem/issues

#include "FunctionSpaceLagrange.h"

#include "Exception.h"
#include "FieldInterface.h"
#include "OmpInterface.h"
#include "instantiate.h"

#include <gmsh.h>
#include <limits>


namespace gmshfem::field
{


  template< class T_PScalar >
  FunctionSpaceLagrange< T_PScalar >::FunctionSpaceLagrange(const unsigned int order) :
    FunctionSpace< T_PScalar, Form::Form0 >(), _order(order)
  {
    if(order != 0) {
      std::vector< int > elementTypes;
      gmsh::model::mesh::getElementTypes(elementTypes);
      std::map< int, std::string > elementsList;
      bool wrongOrder = false;
      for(auto elementType : elementTypes) {
        std::string elementName;
        int dim, elmOrder, numNodes, numPrimaryNodes;
        std::vector< double > localNodeCoord;
        gmsh::model::mesh::getElementProperties(elementType, elementName, dim, elmOrder, numNodes, localNodeCoord, numPrimaryNodes);
        if(dim != 0) {
          elementsList.insert(std::pair(elmOrder, elementName));
          if(static_cast<unsigned int>(elmOrder) != order) {
            wrongOrder = true;
          }
        }
      }
      if(wrongOrder) {
        msg::warning << "Tried to construct a Lagrange function space of order " << order << " on mesh made by elements:" << msg::endl;
        for(auto elm : elementsList) {
          msg::warning << " - " << elm.second << " (order: " << elm.first << ")" << msg::endl;
        }
      }
    }
  }

  template< class T_PScalar >
  FunctionSpaceLagrange< T_PScalar >::~FunctionSpaceLagrange()
  {
  }

  template< class T_PScalar >
  std::string FunctionSpaceLagrange< T_PScalar >::_getGmshName(bool derivative) const
  {
    return (derivative ? "GradLagrange" : "Lagrange");
  }

  template< class T_PScalar >
  std::string FunctionSpaceLagrange< T_PScalar >::getGmshFemName(bool derivative) const
  {
    return std::string(derivative ? "d_" : "") + "Lagrange";
  }

  template< class T_PScalar >
  std::string FunctionSpaceLagrange< T_PScalar >::getGmshFemOrientationName() const
  {
    return "None";
  }

  template< class T_PScalar >
  FunctionSpaceTypeForm0 FunctionSpaceLagrange< T_PScalar >::type() const
  {
    return FunctionSpaceTypeForm0::Lagrange;
  }

  template< class T_PScalar >
  unsigned int FunctionSpaceLagrange< T_PScalar >::order() const
  {
    return 0;
  }

  template< class T_PScalar >
  bool FunctionSpaceLagrange< T_PScalar >::isConstantByElements() const
  {
    return true;
  }

  template< class T_PScalar >
  bool FunctionSpaceLagrange< T_PScalar >::isFormP() const
  {
    return false;
  }

  template< class T_PScalar >
  FunctionSpaceLagrange< T_PScalar > *FunctionSpaceLagrange< T_PScalar >::copy() const
  {
    return new FunctionSpaceLagrange< T_PScalar >(_order);
  }

  template< class T_PScalar >
  unsigned int FunctionSpaceLagrange< T_PScalar >::getBasisFunctions(const bool derivative, std::vector< T_PScalar > &functions, std::vector< int > &orientations, const std::vector< T_PScalar > &integrationPoints, const int elementType, const std::pair< int, int > &entity) const
  {
    if(getDimOfElementType(elementType) < static_cast< unsigned int >(this->form()) - (isFormP() ? 1 : 0)) {
      throw common::Exception("Trying to compute basis functions on elements where no function space is defined");
    }

    return this->getBasisFunction(derivative, functions, integrationPoints, elementType);
  }

  template< class T_PScalar >
  unsigned int FunctionSpaceLagrange< T_PScalar >::getBasisFunction_noexcept(const bool derivative, std::vector< T_PScalar > &functions, const std::vector< T_PScalar > &integrationPoints, const int elementType, const unsigned int orientation) const noexcept
  {
    try {
      functions.clear();
      int fsNumComponents = 0;
      int numOrientation = 0;
      std::vector< double > gmshIntegrationPoints;
      scalar::copy(gmshIntegrationPoints, integrationPoints);
      std::vector< double > gmshFsData;

      gmsh::model::mesh::getBasisFunctions(elementType, gmshIntegrationPoints, _getGmshName(derivative), fsNumComponents, gmshFsData, numOrientation, {static_cast< int >(orientation)});
      const int nbrDofsByElement = gmshFsData.size() / (numOrientation * (integrationPoints.size() / 3) * fsNumComponents);

      scalar::move(functions, gmshFsData);

      this->_errorIfToManyDofsByElements(nbrDofsByElement);
      return nbrDofsByElement;
    }
    catch(...) {
      msg::error << "Unexpected error in 'FunctionSpaceLagrange::getBasisFunction_noexcept'" << msg::endl;
      exit(1);
    }
  }

  INSTANTIATE_CLASS(FunctionSpaceLagrange, 2, TEMPLATE_ARGS(double, float))


} // namespace gmshfem::field
