// GmshFEM - Copyright (C) 2019-2022, A. Royer, E. Béchet, C. Geuzaine, Université de Liège
//
// See the LICENSE.txt file for license information. Please report all
// issues on https://gitlab.onelab.info/gmsh/fem/issues

#include "FunctionSpaceInterface.h"

#include "Exception.h"
#include "Message.h"
#include "Options.h"
#include "gmshfemDefines.h"
#include "instantiate.h"

#include <gmsh.h>

namespace gmshfem::field
{


  template< class T_PScalar >
  void FunctionSpaceInterface< T_PScalar >::_errorIfToManyDofsByElements(const unsigned int nbrDofsByElement) const
  {
    if(nbrDofsByElement > GMSHFEM_DOF_FIELD_OFFSET) {
      msg::error << "The number of degrees of freedom (" << nbrDofsByElement << ") are bigger than 'GMSHFEM_DOF_FIELD_OFFSET'. Recompile GmshFEM with higher 'GMSHFEM_DOF_FIELD_OFFSET' constant (at least higher than " << nbrDofsByElement << " )" << msg::endl;
    }
  }

  template< class T_PScalar >
  std::string FunctionSpaceInterface< T_PScalar >::_getGmshName(bool derivative) const
  {
    return "Undefined";
  }

  template< class T_PScalar >
  FunctionSpaceInterface< T_PScalar >::FunctionSpaceInterface()
  {
  }

  template< class T_PScalar >
  FunctionSpaceInterface< T_PScalar >::~FunctionSpaceInterface()
  {
  }

  static void createMesh(const std::string &elementName, const double meshSize)
  {
    gmsh::model::add(elementName);
    gmsh::model::setCurrent(elementName);

    if(elementName == "point") {
      gmsh::model::geo::addPoint(0., 0., 0., meshSize);
    }
    else if(elementName == "line") {
      int p[2];
      p[0] = gmsh::model::geo::addPoint(-1., 0., 0., meshSize);
      p[1] = gmsh::model::geo::addPoint(1., 0., 0., meshSize);

      gmsh::model::geo::addLine(p[0], p[1]);
    }
    else if(elementName == "triangle") {
      int p[3];
      p[0] = gmsh::model::geo::addPoint(0., 0., 0., meshSize);
      p[1] = gmsh::model::geo::addPoint(1., 0., 0., meshSize);
      p[2] = gmsh::model::geo::addPoint(0., 1., 0., meshSize);

      int l[3];
      l[0] = gmsh::model::geo::addLine(p[0], p[1]);
      l[1] = gmsh::model::geo::addLine(p[1], p[2]);
      l[2] = gmsh::model::geo::addLine(p[2], p[0]);

      int ll = gmsh::model::geo::addCurveLoop({l[0], l[1], l[2]});
      gmsh::model::geo::addPlaneSurface({ll});
    }
    else if(elementName == "quadrangle") {
      int p[4];
      p[0] = gmsh::model::geo::addPoint(-1., -1., 0., meshSize);
      p[1] = gmsh::model::geo::addPoint(1., -1., 0., meshSize);
      p[2] = gmsh::model::geo::addPoint(1., 1., 0., meshSize);
      p[3] = gmsh::model::geo::addPoint(-1., 1., 0., meshSize);

      int l[4];
      l[0] = gmsh::model::geo::addLine(p[0], p[1]);
      l[1] = gmsh::model::geo::addLine(p[1], p[2]);
      l[2] = gmsh::model::geo::addLine(p[2], p[3]);
      l[3] = gmsh::model::geo::addLine(p[3], p[0]);

      int ll = gmsh::model::geo::addCurveLoop({l[0], l[1], l[2], l[3]});
      gmsh::model::geo::addPlaneSurface({ll});
    }
    else if(elementName == "tetrahedron") {
      int p[4];
      p[0] = gmsh::model::geo::addPoint(0., 0., 0., meshSize);
      p[1] = gmsh::model::geo::addPoint(1., 0., 0., meshSize);
      p[2] = gmsh::model::geo::addPoint(0., 1., 0., meshSize);
      p[3] = gmsh::model::geo::addPoint(0., 0., 1., meshSize);

      int l[6];
      l[0] = gmsh::model::geo::addLine(p[0], p[1]);
      l[1] = gmsh::model::geo::addLine(p[1], p[2]);
      l[2] = gmsh::model::geo::addLine(p[2], p[0]);
      l[3] = gmsh::model::geo::addLine(p[3], p[0]);
      l[4] = gmsh::model::geo::addLine(p[3], p[2]);
      l[5] = gmsh::model::geo::addLine(p[3], p[1]);

      int ll[4];
      ll[0] = gmsh::model::geo::addCurveLoop({-l[2], -l[1], -l[0]});
      ll[1] = gmsh::model::geo::addCurveLoop({l[0], -l[5], l[3]});
      ll[2] = gmsh::model::geo::addCurveLoop({-l[3], l[4], l[2]});
      ll[3] = gmsh::model::geo::addCurveLoop({l[5], l[1], -l[4]});

      int s[4];
      s[0] = gmsh::model::geo::addPlaneSurface({ll[0]});
      s[1] = gmsh::model::geo::addPlaneSurface({ll[1]});
      s[2] = gmsh::model::geo::addPlaneSurface({ll[2]});
      s[3] = gmsh::model::geo::addPlaneSurface({ll[3]});

      int sl = gmsh::model::geo::addSurfaceLoop({s[0], s[1], s[2], s[3]});
      gmsh::model::geo::addVolume({sl});
    }
    else if(elementName == "pyramid") {
      int p[5];
      p[0] = gmsh::model::geo::addPoint(-1., -1., 0., meshSize);
      p[1] = gmsh::model::geo::addPoint(1., -1., 0., meshSize);
      p[2] = gmsh::model::geo::addPoint(1., 1., 0., meshSize);
      p[3] = gmsh::model::geo::addPoint(-1., 1., 0., meshSize);
      p[4] = gmsh::model::geo::addPoint(0., 0., 1., meshSize);

      int l[8];
      l[0] = gmsh::model::geo::addLine(p[0], p[1]);
      l[1] = gmsh::model::geo::addLine(p[0], p[3]);
      l[2] = gmsh::model::geo::addLine(p[0], p[4]);
      l[3] = gmsh::model::geo::addLine(p[1], p[2]);
      l[4] = gmsh::model::geo::addLine(p[1], p[4]);
      l[5] = gmsh::model::geo::addLine(p[2], p[3]);
      l[6] = gmsh::model::geo::addLine(p[2], p[4]);
      l[7] = gmsh::model::geo::addLine(p[3], p[4]);

      int ll[5];
      ll[0] = gmsh::model::geo::addCurveLoop({l[0], l[4], -l[2]});
      ll[1] = gmsh::model::geo::addCurveLoop({-l[1], l[2], -l[7]});
      ll[2] = gmsh::model::geo::addCurveLoop({l[3], l[6], -l[4]});
      ll[3] = gmsh::model::geo::addCurveLoop({l[5], l[7], -l[6]});
      ll[4] = gmsh::model::geo::addCurveLoop({l[1], -l[5], -l[3], -l[0]});

      int s[5];
      s[0] = gmsh::model::geo::addPlaneSurface({ll[0]});
      s[1] = gmsh::model::geo::addPlaneSurface({ll[1]});
      s[2] = gmsh::model::geo::addPlaneSurface({ll[2]});
      s[3] = gmsh::model::geo::addPlaneSurface({ll[3]});
      s[4] = gmsh::model::geo::addPlaneSurface({ll[4]});

      int sl = gmsh::model::geo::addSurfaceLoop({s[0], s[1], s[2], s[3], s[4]});
      gmsh::model::geo::addVolume({sl});
    }
    else if(elementName == "prism") {
      int p[6];
      p[0] = gmsh::model::geo::addPoint(0., 0., -1., meshSize);
      p[1] = gmsh::model::geo::addPoint(1., 0., -1., meshSize);
      p[2] = gmsh::model::geo::addPoint(0., 1., -1., meshSize);
      p[3] = gmsh::model::geo::addPoint(0., 0., 1., meshSize);
      p[4] = gmsh::model::geo::addPoint(1., 0., 1., meshSize);
      p[5] = gmsh::model::geo::addPoint(0., 1., 1., meshSize);

      int l[9];
      l[0] = gmsh::model::geo::addLine(p[0], p[1]);
      l[1] = gmsh::model::geo::addLine(p[0], p[2]);
      l[2] = gmsh::model::geo::addLine(p[0], p[3]);
      l[3] = gmsh::model::geo::addLine(p[1], p[2]);
      l[4] = gmsh::model::geo::addLine(p[1], p[4]);
      l[5] = gmsh::model::geo::addLine(p[2], p[5]);
      l[6] = gmsh::model::geo::addLine(p[3], p[4]);
      l[7] = gmsh::model::geo::addLine(p[3], p[5]);
      l[8] = gmsh::model::geo::addLine(p[4], p[5]);

      int ll[5];
      ll[0] = gmsh::model::geo::addCurveLoop({l[1], -l[3], -l[0]});
      ll[1] = gmsh::model::geo::addCurveLoop({l[6], l[8], -l[7]});
      ll[2] = gmsh::model::geo::addCurveLoop({l[0], l[4], -l[6], -l[2]});
      ll[3] = gmsh::model::geo::addCurveLoop({l[2], l[7], -l[5], -l[1]});
      ll[4] = gmsh::model::geo::addCurveLoop({l[3], l[5], -l[8], -l[4]});

      int s[5];
      s[0] = gmsh::model::geo::addPlaneSurface({ll[0]});
      s[1] = gmsh::model::geo::addPlaneSurface({ll[1]});
      s[2] = gmsh::model::geo::addPlaneSurface({ll[2]});
      s[3] = gmsh::model::geo::addPlaneSurface({ll[3]});
      s[4] = gmsh::model::geo::addPlaneSurface({ll[4]});

      int sl = gmsh::model::geo::addSurfaceLoop({s[0], s[1], s[2], s[3], s[4]});
      gmsh::model::geo::addVolume({sl});
    }
    else if(elementName == "hexahedron") {
      int p[8];
      p[0] = gmsh::model::geo::addPoint(-1., -1., -1., meshSize);
      p[1] = gmsh::model::geo::addPoint(1., -1., -1., meshSize);
      p[2] = gmsh::model::geo::addPoint(1., 1., -1., meshSize);
      p[3] = gmsh::model::geo::addPoint(-1., 1., -1., meshSize);
      p[4] = gmsh::model::geo::addPoint(-1., -1., 1., meshSize);
      p[5] = gmsh::model::geo::addPoint(1., -1., 1., meshSize);
      p[6] = gmsh::model::geo::addPoint(1., 1., 1., meshSize);
      p[7] = gmsh::model::geo::addPoint(-1., 1., 1., meshSize);

      int l[12];
      l[0] = gmsh::model::geo::addLine(p[0], p[1]);
      l[1] = gmsh::model::geo::addLine(p[0], p[3]);
      l[2] = gmsh::model::geo::addLine(p[0], p[4]);
      l[3] = gmsh::model::geo::addLine(p[1], p[2]);
      l[4] = gmsh::model::geo::addLine(p[1], p[5]);
      l[5] = gmsh::model::geo::addLine(p[2], p[3]);
      l[6] = gmsh::model::geo::addLine(p[2], p[6]);
      l[7] = gmsh::model::geo::addLine(p[3], p[7]);
      l[8] = gmsh::model::geo::addLine(p[4], p[5]);
      l[9] = gmsh::model::geo::addLine(p[4], p[7]);
      l[10] = gmsh::model::geo::addLine(p[5], p[6]);
      l[11] = gmsh::model::geo::addLine(p[6], p[7]);

      int ll[6];
      ll[0] = gmsh::model::geo::addCurveLoop({l[1], -l[5], -l[3], -l[0]});
      ll[1] = gmsh::model::geo::addCurveLoop({l[0], l[4], -l[8], -l[2]});
      ll[2] = gmsh::model::geo::addCurveLoop({l[2], l[9], -l[7], -l[1]});
      ll[3] = gmsh::model::geo::addCurveLoop({l[3], l[6], -l[10], -l[4]});
      ll[4] = gmsh::model::geo::addCurveLoop({l[5], l[7], -l[11], -l[6]});
      ll[5] = gmsh::model::geo::addCurveLoop({l[8], l[10], l[11], -l[9]});

      int s[6];
      s[0] = gmsh::model::geo::addPlaneSurface({ll[0]});
      s[1] = gmsh::model::geo::addPlaneSurface({ll[1]});
      s[2] = gmsh::model::geo::addPlaneSurface({ll[2]});
      s[3] = gmsh::model::geo::addPlaneSurface({ll[3]});
      s[4] = gmsh::model::geo::addPlaneSurface({ll[4]});
      s[5] = gmsh::model::geo::addPlaneSurface({ll[5]});

      int sl = gmsh::model::geo::addSurfaceLoop({s[0], s[1], s[2], s[3], s[4], s[5]});
      gmsh::model::geo::addVolume({sl});
    }
    else {
      msg::error << "Unknown element named: " << elementName << "" << msg::endl;
    }

    gmsh::model::geo::synchronize();
    gmsh::model::mesh::generate();
  }

  template< class T_PScalar >
  void FunctionSpaceInterface< T_PScalar >::draw(const bool derivative, const std::string &elementName, const unsigned int orientation, const double meshSize) const
  {
    std::string currentName;
    gmsh::model::getCurrent(currentName);

    int order = 1;
    if(getGmshFemName(false) == "Lagrange") {
      std::vector< int > elementTypes;
      gmsh::model::mesh::getElementTypes(elementTypes);
      for(auto i = 0ULL; i < elementTypes.size(); ++i) {
        std::string name;
        int dim, gmshOrder, numNodes, numPrimaryNodes;
        std::vector< double > nodeCoord;
        gmsh::model::mesh::getElementProperties(elementTypes[i], name, dim, gmshOrder, numNodes, nodeCoord, numPrimaryNodes);
        if(gmshOrder > order) {
          order = gmshOrder;
        }
      }
    }

    double verbose = 0.;
    gmsh::option::getNumber("General.Verbosity", verbose);
    gmsh::option::setNumber("General.Verbosity", 1);
    createMesh(elementName, meshSize);

    std::vector< std::size_t > nodeTags;
    std::vector< double > gmshCoord;
    std::vector< T_PScalar > coord;
    std::vector< double > parametricCoord;

    gmsh::model::mesh::getNodes(nodeTags, gmshCoord, parametricCoord, -1, -1, true, false);
    scalar::move(coord, gmshCoord);

    std::vector< T_PScalar > functions;
    int elementType = gmsh::model::mesh::getElementType(elementName, order, false);
    const unsigned int nbrDof = getBasisFunction(derivative, functions, coord, elementType, orientation);
    const unsigned int nbrComponents = numberOfComponents(derivative);

    msg::info << "Writing '" << getGmshFemName(derivative) << "_" << elementName << "_[0," << std::to_string(nbrDof - 1) << "].pos'..." << msg::endl;
    for(auto i = 0U; i < nbrDof; ++i) {
      std::vector< std::vector< double > > gmshData(nodeTags.size());
      for(auto j = 0ULL; j < gmshData.size(); ++j) {
        gmshData[j].reserve(nbrComponents);
        for(auto k = 0U; k < nbrComponents; ++k) {
          gmshData[j].push_back(functions[j * nbrDof * nbrComponents + i * nbrComponents + k]);
        }
      }

      const int tag = gmsh::view::add(getGmshFemName(derivative) + "_" + elementName + "_" + std::to_string(i));
      gmsh::view::addModelData(tag, 0, "", "NodeData", nodeTags, gmshData, 0., nbrComponents);
      gmsh::view::write(tag, getGmshFemName(derivative) + "_" + elementName + "_" + std::to_string(i) + ".pos");

      if(!common::Options::instance()->interface) {
        gmsh::view::remove(tag);
      }
    }
    msg::info << "Done writing '" << getGmshFemName(derivative) << "_" << elementName << "_[0," << std::to_string(nbrDof - 1) << "].pos'" << msg::endl;

    gmsh::model::setCurrent(currentName);
    gmsh::option::setNumber("General.Verbosity", verbose);
  }

  template< class T_PScalar >
  unsigned int FunctionSpaceInterface< T_PScalar >::getBasisFunction(const bool derivative, std::vector< T_PScalar > &functions, const std::vector< T_PScalar > &integrationPoints, const int elementType, const unsigned int orientation) const
  {
    if(getDimOfElementType(elementType) < static_cast< unsigned int >(this->form()) - (isFormP() ? 1 : 0)) {
      throw common::Exception("Trying to compute basis functions on elements where no function space is defined");
    }
    if(orientation >= this->getNumberOfOrientations(elementType)) {
      throw common::Exception("Trying to compute orientation '" + std::to_string(orientation) + "' while there are only " + std::to_string(this->getNumberOfOrientations(elementType)));
    }

    return getBasisFunction_noexcept(derivative, functions, integrationPoints, elementType, orientation);
  }

  template< class T_PScalar >
  void FunctionSpaceInterface< T_PScalar >::getKeys(const bool withCoordinates, std::vector< int > &typeKeys, std::vector< unsigned long long > &entityKeys, std::vector< T_PScalar > &coordinates, const int elementType, const std::pair< int, int > &entity) const
  {
    if(entity.first < static_cast< int >(form()) - (isFormP() ? 1 : 0)) {
      throw common::Exception("Trying to get keys on elements where no function space is defined");
    }
    typeKeys.clear();
    entityKeys.clear();
    coordinates.clear();
    std::vector< double > gmshCoord;
    gmsh::model::mesh::getKeys(elementType, _getGmshName(false), typeKeys, *reinterpret_cast< std::vector< std::size_t > * >(&entityKeys), gmshCoord, entity.second, withCoordinates);
    scalar::move(coordinates, gmshCoord);
  }

  template< class T_PScalar >
  void FunctionSpaceInterface< T_PScalar >::getKeysOnElement(const bool withCoordinates, std::vector< int > &typeKeys, std::vector< unsigned long long > &entityKeys, std::vector< T_PScalar > &coordinates, const unsigned int elementTag) const
  {
    typeKeys.clear();
    entityKeys.clear();
    coordinates.clear();
    std::vector< double > gmshCoord;
    gmsh::model::mesh::getKeysForElement(elementTag, _getGmshName(false), typeKeys, *reinterpret_cast< std::vector< std::size_t > * >(&entityKeys), gmshCoord, withCoordinates);
    scalar::move(coordinates, gmshCoord);
  }

  template< class T_PScalar >
  bool FunctionSpaceInterface< T_PScalar >::getPeriodicKeys(const bool withCoordinates, std::vector< int > &slaveTypeKeys, std::vector< unsigned long long > &slaveEntityKeys, std::vector< T_PScalar > &slaveCoordinates, std::vector< int > &masterTypeKeys, std::vector< unsigned long long > &masterEntityKeys, std::vector< T_PScalar > &masterCoordinates, const int elementType, const std::pair< int, int > &entity) const
  {
    if(entity.first < static_cast< int >(form()) - (isFormP() ? 1 : 0)) {
      throw common::Exception("Trying to get periodic keys on elements where no function space is defined");
    }
    slaveTypeKeys.clear();
    slaveEntityKeys.clear();
    slaveCoordinates.clear();
    masterTypeKeys.clear();
    masterEntityKeys.clear();
    masterCoordinates.clear();
    std::vector< double > slaveGmshCoord, masterGmshCoord;
    int masterTag;
    gmsh::model::mesh::getPeriodicKeys(elementType, _getGmshName(false), entity.second, masterTag, slaveTypeKeys, masterTypeKeys, *reinterpret_cast< std::vector< std::size_t > * >(&slaveEntityKeys), *reinterpret_cast< std::vector< std::size_t > * >(&masterEntityKeys), slaveGmshCoord, masterGmshCoord, withCoordinates);
    scalar::move(slaveCoordinates, slaveGmshCoord);
    scalar::move(masterCoordinates, masterGmshCoord);
    if(entity.second == masterTag) {
      return false;
    }
    return true;
  }

  template< class T_PScalar >
  bool FunctionSpaceInterface< T_PScalar >::operator==(const FunctionSpaceInterface< T_PScalar > &other) const
  {
    return getGmshFemName(false) == other.getGmshFemName(false);
  }

  template< class T_PScalar >
  bool FunctionSpaceInterface< T_PScalar >::operator!=(const FunctionSpaceInterface< T_PScalar > &other) const
  {
    return getGmshFemName(false) != other.getGmshFemName(false);
  }

  template< class T_PScalar >
  bool FunctionSpaceInterface< T_PScalar >::operator<(const FunctionSpaceInterface< T_PScalar > &other) const
  {
    return getGmshFemName(false) < other.getGmshFemName(false);
  }

  template< class T_PScalar >
  unsigned int FunctionSpaceInterface< T_PScalar >::getNumberOfKeysByElement(const int elementType) const
  {
    if(getDimOfElementType(elementType) < static_cast< unsigned int >(form()) - (isFormP() ? 1 : 0)) {
      throw common::Exception("Trying to get the number of keys on elements where no function space is defined");
    }
    return gmsh::model::mesh::getNumberOfKeys(elementType, _getGmshName(false));
  }

  template< class T_PScalar >
  unsigned int FunctionSpaceInterface< T_PScalar >::getNumberOfKeysByNode() const
  {
    try {
      return getNumberOfKeysByElement(gmsh::model::mesh::getElementType("Point", 1));
    }
    catch(const common::Exception &exc) {
      return 0;
    }
  }

  template< class T_PScalar >
  unsigned int FunctionSpaceInterface< T_PScalar >::getNumberOfKeysByEdge() const
  {
    try {
      return getNumberOfKeysByElement(gmsh::model::mesh::getElementType("Line", 1)) - 2 * getNumberOfKeysByNode();
    }
    catch(const common::Exception &exc) {
      return 0;
    }
  }

  template< class T_PScalar >
  unsigned int FunctionSpaceInterface< T_PScalar >::getNumberOfKeysByTriangularFace() const
  {
    try {
      return getNumberOfKeysByElement(gmsh::model::mesh::getElementType("Triangle", 1)) - 3 * getNumberOfKeysByEdge() - 3 * getNumberOfKeysByNode();
    }
    catch(const common::Exception &exc) {
      return 0;
    }
  }

  template< class T_PScalar >
  unsigned int FunctionSpaceInterface< T_PScalar >::getNumberOfKeysByQuadrangularFace() const
  {
    try {
      return getNumberOfKeysByElement(gmsh::model::mesh::getElementType("Quadrangle", 1)) - 4 * getNumberOfKeysByEdge() - 4 * getNumberOfKeysByNode();
    }
    catch(const common::Exception &exc) {
      return 0;
    }
  }

  template< class T_PScalar >
  unsigned int FunctionSpaceInterface< T_PScalar >::getNumberOfOrientations(const int elementType) const
  {
    if(getDimOfElementType(elementType) < static_cast< unsigned int >(form()) - (isFormP() ? 1 : 0)) {
      throw common::Exception("Trying to get the number of orientations on elements where no function space is defined");
    }
    return gmsh::model::mesh::getNumberOfOrientations(elementType, _getGmshName(false));
  }

  template< class T_PScalar >
  void FunctionSpaceInterface< T_PScalar >::getKeyInformation(const std::vector< int > &typeKeys, const std::vector< unsigned long long > &entityKeys, const int elementType, std::vector< std::pair< int, int > > &infoKey) const
  {
    if(getDimOfElementType(elementType) < static_cast< unsigned int >(form()) - (isFormP() ? 1 : 0)) {
      throw common::Exception("Trying to get key information on elements where no function space is defined");
    }
    gmsh::model::mesh::getKeysInformation(typeKeys, *reinterpret_cast< const std::vector< std::size_t > * >(&entityKeys), elementType, _getGmshName(false), infoKey);
  }

  template< class T_PScalar >
  unsigned int FunctionSpaceInterface< T_PScalar >::getOrientationOfElement(const unsigned int elementTag) const
  {
    int basisFunctionsOrientation = 0;
    gmsh::model::mesh::getBasisFunctionsOrientationForElement(elementTag, _getGmshName(false), basisFunctionsOrientation);
    return basisFunctionsOrientation;
  }

  template< class T_PScalar >
  unsigned int FunctionSpaceInterface< T_PScalar >::getGaussInfo(const std::string &integrationType, const int elementType, std::vector< T_PScalar > &weights, std::vector< T_PScalar > &points) const
  {
    return GetGaussInfo(integrationType, elementType, weights, points);
  }

  template< class T_PScalar >
  unsigned int FunctionSpaceInterface< T_PScalar >::GetGaussInfo(const std::string &integrationType, const int elementType, std::vector< T_PScalar > &weights, std::vector< T_PScalar > &points)
  {
    std::vector< double > gmshPoints;
    std::vector< double > gmshWeights;
    gmsh::model::mesh::getIntegrationPoints(elementType, integrationType, gmshPoints, gmshWeights);

    scalar::move(weights, gmshWeights);
    scalar::move(points, gmshPoints);

    return weights.size();
  }

  INSTANTIATE_CLASS(FunctionSpaceInterface, 2, TEMPLATE_ARGS(double, float))


  unsigned int getDimOfElementType(const int elementType)
  {
    std::string elementName;
    int dim, order, numNodes, numPrimaryNodes;
    std::vector< double > localNodeCoord;
    gmsh::model::mesh::getElementProperties(elementType, elementName, dim, order, numNodes, localNodeCoord, numPrimaryNodes);
    return dim;
  }

  std::vector< int > getDenseOrientations(std::vector< int > &orientations, const unsigned int numberOfOrientations)
  {
    std::vector< int > uniqueOrientationsNumbering;
    uniqueOrientationsNumbering.reserve(numberOfOrientations / 20);
    std::pair< int, int > previousPair(-1, -1);
    for(auto i = 0ULL; i < orientations.size(); ++i) {
      if(orientations[i] == previousPair.first) {
        orientations[i] = previousPair.second;
      }
      else {
        auto itFind = std::find(uniqueOrientationsNumbering.begin(), uniqueOrientationsNumbering.end(), orientations[i]);
        if(itFind == uniqueOrientationsNumbering.end()) {
          uniqueOrientationsNumbering.push_back(orientations[i]);
          previousPair = std::make_pair(orientations[i], uniqueOrientationsNumbering.size() - 1);
          orientations[i] = previousPair.second;
        }
        else {
          previousPair = std::make_pair(orientations[i], itFind - uniqueOrientationsNumbering.begin());
          orientations[i] = previousPair.second;
        }
      }
    }

    return uniqueOrientationsNumbering;
  }


} // namespace gmshfem::field
