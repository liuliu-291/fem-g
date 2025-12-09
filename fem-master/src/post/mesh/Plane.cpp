// GmshFEM - Copyright (C) 2019-2022, A. Royer, E. Béchet, C. Geuzaine, Université de Liège
//
// See the LICENSE.txt file for license information. Please report all
// issues on https://gitlab.onelab.info/gmsh/fem/issues

#include "Plane.h"

#include "Message.h"

#include <gmsh.h>

namespace gmshfem::post
{


  Plane::Plane(const std::string &name) :
    Mesh(name), _origin{0., 0., 0.}, _vectorX{0., 0., 0.}, _vectorY{0., 0., 0.}, _numPoints{0, 0}
  {
  }

  Plane::Plane(const std::string &name, const double xOrigin, const double yOrigin, const double zOrigin, const double vecXX, const double vecXY, const double vecXZ, const double vecYX, const double vecYY, const double vecYZ, const unsigned int numPointsX, const unsigned int numPointsY) :
    Mesh(name), _origin{xOrigin, yOrigin, zOrigin}, _vectorX{vecXX, vecXY, vecXZ}, _vectorY{vecYX, vecYY, vecYZ}, _numPoints{numPointsX, numPointsY}
  {
  }

  Plane::Plane(Plane &&other) :
    Mesh(std::move(other)), _origin{other._origin[0], other._origin[1], other._origin[2]}, _vectorX{other._vectorX[0], other._vectorX[1], other._vectorX[2]}, _vectorY{other._vectorY[0], other._vectorY[1], other._vectorY[2]}, _numPoints{other._numPoints[0], other._numPoints[1]}
  {
  }

  Plane::~Plane()
  {
  }

  void Plane::generateMesh() const
  {
    if(_haveMesh) {
      return;
    }

    double verbose = 0.;
    gmsh::option::getNumber("General.Verbosity", verbose);
    gmsh::option::setNumber("General.Verbosity", 1);

    gmsh::model::add(_name);
    gmsh::model::setCurrent(_name);

    int tag = gmsh::model::geo::addPoint(_origin[0], _origin[1], _origin[2]);
    std::vector< std::pair< int, int > > outDimTags;
    gmsh::model::geo::extrude({std::make_pair(0, tag)}, _vectorX[0], _vectorX[1], _vectorX[2], outDimTags, {static_cast< int >(_numPoints[0] - 1)}, {}, true);
    gmsh::model::geo::extrude({outDimTags[1]}, _vectorY[0], _vectorY[1], _vectorY[2], outDimTags, {static_cast< int >(_numPoints[1] - 1)}, {}, true);

    gmsh::model::addPhysicalGroup(outDimTags[1].first, {outDimTags[1].second}, 1);
    gmsh::model::setPhysicalName(outDimTags[1].first, 1, "omega_" + _name);

    gmsh::model::geo::synchronize();
    gmsh::model::mesh::generate();

    _domain = domain::Domain(outDimTags[1].first, 1);

    gmsh::option::setNumber("General.Verbosity", verbose);
    _haveMesh = true;
  }


} // namespace gmshfem::post
