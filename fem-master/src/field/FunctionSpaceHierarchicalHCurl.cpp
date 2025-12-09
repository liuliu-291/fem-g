// GmshFEM - Copyright (C) 2019-2022, A. Royer, E. Béchet, C. Geuzaine, Université de Liège
//
// See the LICENSE.txt file for license information. Please report all
// issues on https://gitlab.onelab.info/gmsh/fem/issues

#include "FunctionSpaceHierarchicalHCurl.h"

#include "Exception.h"
#include "FieldInterface.h"
#include "Message.h"
#include "OmpInterface.h"
#include "instantiate.h"

#include <gmsh.h>


namespace gmshfem::field
{


  template< class T_PScalar >
  FunctionSpaceHierarchicalHCurl< T_PScalar >::FunctionSpaceHierarchicalHCurl(const unsigned int order) :
    FunctionSpace< T_PScalar, Form::Form1 >(), _order(order)
  {
  }

  template< class T_PScalar >
  FunctionSpaceHierarchicalHCurl< T_PScalar >::~FunctionSpaceHierarchicalHCurl()
  {
  }

  template< class T_PScalar >
  std::string FunctionSpaceHierarchicalHCurl< T_PScalar >::_getGmshName(bool derivative) const
  {
    return (derivative ? "CurlHcurlLegendre" + std::to_string(_order) : "HcurlLegendre" + std::to_string(_order));
  }

  template< class T_PScalar >
  std::string FunctionSpaceHierarchicalHCurl< T_PScalar >::getGmshFemName(bool derivative) const
  {
    return std::string(derivative ? "d_" : "") + "HierarchicalHCurl(" + std::to_string(_order) + ")";
  }

  template< class T_PScalar >
  std::string FunctionSpaceHierarchicalHCurl< T_PScalar >::getGmshFemOrientationName() const
  {
    return "HierarchicalSolin";
  }

  template< class T_PScalar >
  FunctionSpaceTypeForm1 FunctionSpaceHierarchicalHCurl< T_PScalar >::type() const
  {
    return FunctionSpaceTypeForm1::HierarchicalHCurl;
  }

  template< class T_PScalar >
  unsigned int FunctionSpaceHierarchicalHCurl< T_PScalar >::order() const
  {
    return _order;
  }

  template< class T_PScalar >
  bool FunctionSpaceHierarchicalHCurl< T_PScalar >::isConstantByElements() const
  {
    return false;
  }

  template< class T_PScalar >
  bool FunctionSpaceHierarchicalHCurl< T_PScalar >::isFormP() const
  {
    return false;
  }

  template< class T_PScalar >
  FunctionSpaceHierarchicalHCurl< T_PScalar > *FunctionSpaceHierarchicalHCurl< T_PScalar >::copy() const
  {
    return new FunctionSpaceHierarchicalHCurl< T_PScalar >(_order);
  }

  template< class T_PScalar >
  unsigned int FunctionSpaceHierarchicalHCurl< T_PScalar >::getBasisFunctions(const bool derivative, std::vector< T_PScalar > &functions, std::vector< int > &orientations, const std::vector< T_PScalar > &integrationPoints, const int elementType, const std::pair< int, int > &entity) const
  {
    if(getDimOfElementType(elementType) < static_cast< unsigned int >(this->form()) - (isFormP() ? 1 : 0)) {
      throw common::Exception("Trying to compute basis functions on elements where no function space is defined");
    }

    gmsh::model::mesh::preallocateBasisFunctionsOrientation(elementType, orientations, entity.second);
    gmsh::model::mesh::getBasisFunctionsOrientation(elementType, _getGmshName(derivative), orientations, entity.second);

    const unsigned int numberOfOrientations = gmsh::model::mesh::getNumberOfOrientations(elementType, _getGmshName(derivative));
    std::vector< int > uniqueOrientationsNumbering = getDenseOrientations(orientations, numberOfOrientations);

    functions.clear();
    int fsNumComponents = 0;
    int numOrientation = 0;
    std::vector< double > gmshIntegrationPoints;
    scalar::copy(gmshIntegrationPoints, integrationPoints);
    std::vector< double > gmshFsData;

    gmsh::model::mesh::getBasisFunctions(elementType, gmshIntegrationPoints, _getGmshName(derivative), fsNumComponents, gmshFsData, numOrientation, uniqueOrientationsNumbering);
    const int nbrDofsByElement = gmshFsData.size() / (uniqueOrientationsNumbering.size() * (integrationPoints.size() / 3) * fsNumComponents);

    scalar::move(functions, gmshFsData);

    this->_errorIfToManyDofsByElements(nbrDofsByElement);
    return nbrDofsByElement;
  }

  template< class T_PScalar >
  unsigned int FunctionSpaceHierarchicalHCurl< T_PScalar >::getBasisFunction_noexcept(const bool derivative, std::vector< T_PScalar > &functions, const std::vector< T_PScalar > &integrationPoints, const int elementType, const unsigned int orientation) const noexcept
  {
    try {
      functions.clear();
      int fsNumComponents = 0;
      int numOrientation = 0;
      std::vector< double > gmshIntegrationPoints;
      scalar::copy(gmshIntegrationPoints, integrationPoints);
      std::vector< double > gmshFsData;

      gmsh::model::mesh::getBasisFunctions(elementType, gmshIntegrationPoints, _getGmshName(derivative), fsNumComponents, gmshFsData, numOrientation, {static_cast< int >(orientation)});
      const int nbrDofsByElement = gmshFsData.size() / ((integrationPoints.size() / 3) * fsNumComponents);

      scalar::move(functions, gmshFsData);

      this->_errorIfToManyDofsByElements(nbrDofsByElement);
      return nbrDofsByElement;
    }
    catch(...) {
      msg::error << "Unexpected error in 'FunctionSpaceHierarchicalHCurl::getBasisFunction_noexcept'" << msg::endl;
      exit(1);
    }
  }

  INSTANTIATE_CLASS(FunctionSpaceHierarchicalHCurl, 2, TEMPLATE_ARGS(double, float))


} // namespace gmshfem::field
