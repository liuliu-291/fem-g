// GmshFEM - Copyright (C) 2019-2022, A. Royer, E. Béchet, C. Geuzaine, Université de Liège
//
// See the LICENSE.txt file for license information. Please report all
// issues on https://gitlab.onelab.info/gmsh/fem/issues

#include "FunctionSpaceConstant.h"

#include "Exception.h"
#include "FieldInterface.h"
#include "OmpInterface.h"
#include "instantiate.h"

#include <gmsh.h>

namespace gmshfem::field
{


  template< class T_PScalar >
  FunctionSpaceConstant< T_PScalar >::FunctionSpaceConstant(const bool pForm) :
    FunctionSpace< T_PScalar, Form::Form3 >(), _pForm(pForm)
  {
  }

  template< class T_PScalar >
  FunctionSpaceConstant< T_PScalar >::~FunctionSpaceConstant()
  {
  }

  template< class T_PScalar >
  std::string FunctionSpaceConstant< T_PScalar >::getGmshFemName(bool derivative) const
  {
    return (_pForm ? std::string("P_") : std::string("")) + "Constant";
  }

  template< class T_PScalar >
  std::string FunctionSpaceConstant< T_PScalar >::getGmshFemOrientationName() const
  {
    return "None";
  }

  template< class T_PScalar >
  FunctionSpaceTypeForm3 FunctionSpaceConstant< T_PScalar >::type() const
  {
    return (_pForm ? FunctionSpaceTypeForm3::P_Constant : FunctionSpaceTypeForm3::Constant);
  }

  template< class T_PScalar >
  unsigned int FunctionSpaceConstant< T_PScalar >::order() const
  {
    return 0;
  }

  template< class T_PScalar >
  bool FunctionSpaceConstant< T_PScalar >::isConstantByElements() const
  {
    return true;
  }

  template< class T_PScalar >
  bool FunctionSpaceConstant< T_PScalar >::isFormP() const
  {
    return _pForm;
  }

  template< class T_PScalar >
  FunctionSpaceConstant< T_PScalar > *FunctionSpaceConstant< T_PScalar >::copy() const
  {
    return new FunctionSpaceConstant< T_PScalar >(_pForm);
  }

  template< class T_PScalar >
  bool FunctionSpaceConstant< T_PScalar >::operator==(const FunctionSpaceInterface< T_PScalar > &other) const
  {
    if(this->form() == other.form()) {
      return (type() == static_cast< const FunctionSpace< T_PScalar, Form::Form3 > & >(other).type());
    }
    return false;
  }

  template< class T_PScalar >
  unsigned int FunctionSpaceConstant< T_PScalar >::getBasisFunctions(const bool derivative, std::vector< T_PScalar > &functions, std::vector< int > &orientations, const std::vector< T_PScalar > &integrationPoints, const int elementType, const std::pair< int, int > &entity) const
  {
    if(getDimOfElementType(elementType) < static_cast< unsigned int >(this->form()) - (isFormP() ? 1 : 0)) {
      throw common::Exception("Trying to compute basis function on elements where no function space is defined");
    }

    return this->getBasisFunction(derivative, functions, integrationPoints, elementType);
  }

  template< class T_PScalar >
  unsigned int FunctionSpaceConstant< T_PScalar >::getBasisFunction_noexcept(const bool derivative, std::vector< T_PScalar > &functions, const std::vector< T_PScalar > &integrationPoints, const int elementType, const unsigned int orientation) const noexcept
  {
    try {
      functions.clear();

      if(derivative) {
        functions.resize(integrationPoints.size(), 0.);
      }
      else {
        functions.resize(integrationPoints.size(), 1.);
      }

      this->_errorIfToManyDofsByElements(1);
      return 1;
    }
    catch(...) {
      msg::error << "Unexpected error in 'FunctionSpaceConstant::getBasisFunction_noexcept'" << msg::endl;
      exit(1);
    }
  }

  template< class T_PScalar >
  void FunctionSpaceConstant< T_PScalar >::getKeys(const bool withCoordinates, std::vector< int > &typeKeys, std::vector< unsigned long long > &entityKeys, std::vector< T_PScalar > &coordinates, const int elementType, const std::pair< int, int > &entity) const
  {
    std::vector< std::size_t > elements;
    std::vector< std::size_t > nodes;
    gmsh::model::mesh::getElementsByType(elementType, elements, nodes, entity.second);

    typeKeys.resize(elements.size());
    entityKeys.resize(elements.size());
    std::vector< double > gmshCoord;
#pragma omp parallel num_threads(omp::getMaxThreads())
    {
#pragma omp for
      for(auto i = 0ULL; i < elements.size(); ++i) {
        typeKeys[i] = 0;
        entityKeys[i] = elements[i];
      }

      if(withCoordinates) {
        const unsigned int numThreads = omp::getNumThreads();
        const unsigned int myThreadID = omp::getThreadNum();

#pragma omp single
        gmsh::model::mesh::preallocateBarycenters(elementType, gmshCoord, entity.second);

        gmsh::model::mesh::getBarycenters(elementType, entity.second, false, true, gmshCoord, myThreadID, numThreads);
#pragma omp barrier
#pragma omp single
        scalar::move(coordinates, gmshCoord);
      }
    }
  }

  template< class T_PScalar >
  void FunctionSpaceConstant< T_PScalar >::getKeysOnElement(const bool withCoordinates, std::vector< int > &typeKeys, std::vector< unsigned long long > &entityKeys, std::vector< T_PScalar > &coordinates, const unsigned int elementTag) const
  {
    typeKeys.clear();
    entityKeys.clear();
    coordinates.clear();

    typeKeys.push_back(0);
    entityKeys.push_back(elementTag);
    if(withCoordinates) {
      std::vector< double > gmshCoord;
      int elementType, dim, tag;
      std::vector< std::size_t > nodeTags;
      gmsh::model::mesh::getElement(elementTag, elementType, nodeTags, dim, tag);

      coordinates.resize(3, 0.);
      for(auto i = 0ULL; i < nodeTags.size(); ++i) {
        std::vector< double > coord;
        std::vector< double > parametricCoord;
        gmsh::model::mesh::getNode(nodeTags[i], coord, parametricCoord, dim, tag);
        coordinates[0] += coord[0];
        coordinates[1] += coord[1];
        coordinates[2] += coord[2];
      }
      coordinates[0] /= nodeTags.size();
      coordinates[1] /= nodeTags.size();
      coordinates[2] /= nodeTags.size();
    }
  }

  template< class T_PScalar >
  bool FunctionSpaceConstant< T_PScalar >::getPeriodicKeys(const bool withCoordinates, std::vector< int > &slaveTypeKeys, std::vector< unsigned long long > &slaveEntityKeys, std::vector< T_PScalar > &slaveCoordinates, std::vector< int > &masterTypeKeys, std::vector< unsigned long long > &masterEntityKeys, std::vector< T_PScalar > &masterCoordinates, const int elementType, const std::pair< int, int > &entity) const
  {
    throw common::Exception("'getPeriodicKeys' is not implemented for " + getGmshFemName(false) + " basis function");
    return true;
  }

  template< class T_PScalar >
  unsigned int FunctionSpaceConstant< T_PScalar >::getNumberOfKeysByElement(const int elementType) const
  {
    if(getDimOfElementType(elementType) < static_cast< unsigned int >(this->form()) - (isFormP() ? 1 : 0)) {
      throw common::Exception("Trying to get the number of keys on elements where no function space is defined");
    }
    return 1;
  }

  template< class T_PScalar >
  unsigned int FunctionSpaceConstant< T_PScalar >::getNumberOfKeysByNode() const
  {
    return 0;
  }

  template< class T_PScalar >
  unsigned int FunctionSpaceConstant< T_PScalar >::getNumberOfKeysByEdge() const
  {
    return 0;
  }

  template< class T_PScalar >
  unsigned int FunctionSpaceConstant< T_PScalar >::getNumberOfKeysByTriangularFace() const
  {
    return _pForm ? 1 : 0;
  }

  template< class T_PScalar >
  unsigned int FunctionSpaceConstant< T_PScalar >::getNumberOfKeysByQuadrangularFace() const
  {
    return _pForm ? 1 : 0;
  }

  template< class T_PScalar >
  unsigned int FunctionSpaceConstant< T_PScalar >::getNumberOfOrientations(const int elementType) const
  {
    if(getDimOfElementType(elementType) < static_cast< unsigned int >(this->form()) - (isFormP() ? 1 : 0)) {
      throw common::Exception("Trying to get the number of orientations on elements where no function space is defined");
    }
    return 1;
  }

  template< class T_PScalar >
  void FunctionSpaceConstant< T_PScalar >::getKeyInformation(const std::vector< int > &typeKeys, const std::vector< unsigned long long > &entityKeys, const int elementType, std::vector< std::pair< int, int > > &infoKey) const
  {
    if(getDimOfElementType(elementType) < static_cast< unsigned int >(this->form()) - (isFormP() ? 1 : 0)) {
      throw common::Exception("Trying to get key information on elements where no function space is defined");
    }
    infoKey.resize(typeKeys.size());
    std::string elementName;
    int dim = 0, order = 0, numNodes = 0, numPrimaryNodes = 0;
    std::vector< double > nodeCoord;
    gmsh::model::mesh::getElementProperties(elementType, elementName, dim, order, numNodes, nodeCoord, numPrimaryNodes);
    for(auto i = 0ULL; i < typeKeys.size(); ++i) {
      infoKey[i] = std::pair< int, int >(dim, 0);
    }
  }

  template< class T_PScalar >
  unsigned int FunctionSpaceConstant< T_PScalar >::getOrientationOfElement(const unsigned int elementTag) const
  {
    return 0;
  }

  INSTANTIATE_CLASS(FunctionSpaceConstant, 2, TEMPLATE_ARGS(double, float))


} // namespace gmshfem::field
