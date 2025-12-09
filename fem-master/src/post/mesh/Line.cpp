// GmshFEM - Copyright (C) 2019-2022, A. Royer, E. Béchet, C. Geuzaine, Université de Liège
//
// See the LICENSE.txt file for license information. Please report all
// issues on https://gitlab.onelab.info/gmsh/fem/issues

#include "Line.h"

#include <gmsh.h>

namespace gmshfem::post
{


  Line::Line(const std::string &name) :
    Mesh(name), _min{0., 0., 0.}, _max{0., 0., 0.}, _numPoints(0)
  {
  }

  Line::Line(const std::string &name, const double xMin, const double yMin, const double zMin, const double xMax, const double yMax, const double zMax, const unsigned int numPoints) :
    Mesh(name), _min{xMin, yMin, zMin}, _max{xMax, yMax, zMax}, _numPoints(numPoints)
  {
  }

  Line::Line(Line &&other) :
    Mesh(std::move(other)), _min{other._min[0], other._min[1], other._min[2]}, _max{other._max[0], other._max[1], other._max[2]}, _numPoints(other._numPoints)
  {
  }

  Line::~Line()
  {
  }

  void Line::generateMesh() const
  {
    if(_haveMesh) {
      return;
    }

    double verbose = 0.;
    gmsh::option::getNumber("General.Verbosity", verbose);
    gmsh::option::setNumber("General.Verbosity", 1);

    gmsh::model::add(_name);
    gmsh::model::setCurrent(_name);

    int tag = gmsh::model::geo::addPoint(_min[0], _min[1], _min[2]);
    std::vector< std::pair< int, int > > outDimTags;
    gmsh::model::geo::extrude({std::make_pair(0, tag)}, _max[0] - _min[0], _max[1] - _min[1], _max[2] - _min[2], outDimTags, {static_cast< int >(_numPoints - 1)}, {}, true);

    gmsh::model::addPhysicalGroup(outDimTags[1].first, {outDimTags[1].second}, 1);
    gmsh::model::setPhysicalName(outDimTags[1].first, 1, "omega_" + _name);

    gmsh::model::geo::synchronize();
    gmsh::model::mesh::generate();

    _domain = domain::Domain(outDimTags[1].first, 1);

    gmsh::option::setNumber("General.Verbosity", verbose);
    _haveMesh = true;
  }


} // namespace gmshfem::post
