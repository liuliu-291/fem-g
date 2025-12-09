// GmshFEM - Copyright (C) 2019-2022, A. Royer, E. Béchet, C. Geuzaine, Université de Liège
//
// See the LICENSE.txt file for license information. Please report all
// issues on https://gitlab.onelab.info/gmsh/fem/issues

#include "FunctionSpaceHierarchicalH1.h"

#include "Exception.h"
#include "FieldInterface.h"
#include "OmpInterface.h"
#include "instantiate.h"

#include <gmsh.h>


namespace gmshfem::field
{


  template< class T_PScalar >
  FunctionSpaceHierarchicalH1< T_PScalar >::FunctionSpaceHierarchicalH1(const unsigned int order) :
    FunctionSpace< T_PScalar, Form::Form0 >(), _order(order)
  {
  }

  template< class T_PScalar >
  FunctionSpaceHierarchicalH1< T_PScalar >::~FunctionSpaceHierarchicalH1()
  {
  }

  template< class T_PScalar >
  std::string FunctionSpaceHierarchicalH1< T_PScalar >::_getGmshName(bool derivative) const
  {
    return (derivative ? "GradH1Legendre" + std::to_string(_order) : "H1Legendre" + std::to_string(_order));
  }

  template< class T_PScalar >
  std::string FunctionSpaceHierarchicalH1< T_PScalar >::getGmshFemName(bool derivative) const
  {
    return std::string(derivative ? "d_" : "") + "HierarchicalH1(" + std::to_string(_order) + ")";
  }

  template< class T_PScalar >
  std::string FunctionSpaceHierarchicalH1< T_PScalar >::getGmshFemOrientationName() const
  {
    if(_order == 1) {
      return "None";
    }
    return "HierarchicalSolin";
  }

  template< class T_PScalar >
  FunctionSpaceHierarchicalH1< T_PScalar > *FunctionSpaceHierarchicalH1< T_PScalar >::copy() const
  {
    return new FunctionSpaceHierarchicalH1< T_PScalar >(_order);
  }

  template< class T_PScalar >
  FunctionSpaceTypeForm0 FunctionSpaceHierarchicalH1< T_PScalar >::type() const
  {
    return FunctionSpaceTypeForm0::HierarchicalH1;
  }

  template< class T_PScalar >
  unsigned int FunctionSpaceHierarchicalH1< T_PScalar >::order() const
  {
    return _order;
  }

  template< class T_PScalar >
  bool FunctionSpaceHierarchicalH1< T_PScalar >::isConstantByElements() const
  {
    return false;
  }

  template< class T_PScalar >
  bool FunctionSpaceHierarchicalH1< T_PScalar >::isFormP() const
  {
    return false;
  }

  template< class T_PScalar >
  unsigned int FunctionSpaceHierarchicalH1< T_PScalar >::getBasisFunctions(const bool derivative, std::vector< T_PScalar > &functions, std::vector< int > &orientations, const std::vector< T_PScalar > &integrationPoints, const int elementType, const std::pair< int, int > &entity) const
  {
    if(getDimOfElementType(elementType) < static_cast< unsigned int >(this->form()) - (isFormP() ? 1 : 0)) {
      throw common::Exception("Trying to compute basis function on elements where no function space is defined");
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
  unsigned int FunctionSpaceHierarchicalH1< T_PScalar >::getBasisFunction_noexcept(const bool derivative, std::vector< T_PScalar > &functions, const std::vector< T_PScalar > &integrationPoints, const int elementType, const unsigned int orientation) const noexcept
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
      msg::error << "Unexpected error in 'FunctionSpaceHierarchicalH1::getBasisFunction_noexcept'" << msg::endl;
      exit(1);
    }
  }

  INSTANTIATE_CLASS(FunctionSpaceHierarchicalH1, 2, TEMPLATE_ARGS(double, float))


} // namespace gmshfem::field
