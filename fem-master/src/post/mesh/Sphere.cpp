// GmshFEM - Copyright (C) 2019-2022, A. Royer, E. Béchet, C. Geuzaine, Université de Liège
//
// See the LICENSE.txt file for license information. Please report all
// issues on https://gitlab.onelab.info/gmsh/fem/issues

#include "Sphere.h"

#include <gmsh.h>
#include <vector>

namespace gmshfem::post
{


  Sphere::Sphere(const std::string &name) :
    Mesh(name), _center{0., 0., 0.}, _numPoints(0)
  {
  }

  Sphere::Sphere(const std::string &name, const double xCenter, const double yCenter, const double zCenter, const double radius, const unsigned int numPoints) :
    Mesh(name), _center{xCenter, yCenter, zCenter}, _radius(radius), _numPoints(numPoints)
  {
  }

  Sphere::Sphere(Sphere &&other) :
    Mesh(std::move(other)), _center{other._center[0], other._center[1], other._center[2]}, _radius(other._radius), _numPoints(other._numPoints)
  {
  }

  Sphere::~Sphere()
  {
  }

  void Sphere::generateMesh() const
  {
    if(_haveMesh) {
      return;
    }

    double verbose = 0.;
    gmsh::option::getNumber("General.Verbosity", verbose);
    gmsh::option::setNumber("General.Verbosity", 1);

    gmsh::model::add(_name);
    gmsh::model::setCurrent(_name);

    int center = gmsh::model::geo::addPoint(_center[0], _center[1], _center[2]);
    std::vector< int > points(6, 0);
    points[0] = gmsh::model::geo::addPoint(_center[0] + _radius, _center[1], _center[2]);
    points[1] = gmsh::model::geo::addPoint(_center[0], _center[1] + _radius, _center[2]);
    points[2] = gmsh::model::geo::addPoint(_center[0] - _radius, _center[1], _center[2]);
    points[3] = gmsh::model::geo::addPoint(_center[0], _center[1] - _radius, _center[2]);
    points[4] = gmsh::model::geo::addPoint(_center[0], _center[1], _center[2] + _radius);
    points[5] = gmsh::model::geo::addPoint(_center[0], _center[1], _center[2] - _radius);

    std::vector< int > lines(12, 0);
    for(auto i = 0; i < 4; ++i) {
      lines[i] = gmsh::model::geo::addCircleArc(points[i], center, points[(i + 1) % 4]);
      lines[i + 4] = gmsh::model::geo::addCircleArc(points[4], center, points[i]);
      lines[i + 8] = gmsh::model::geo::addCircleArc(points[5], center, points[i]);
    }

    std::vector< int > lineLoops(8, 0);
    for(auto i = 0; i < 4; ++i) {
      lineLoops[i] = gmsh::model::geo::addCurveLoop({lines[i], -lines[(i + 1) % 4 + 4], lines[i + 4]});
      lineLoops[i + 4] = gmsh::model::geo::addCurveLoop({-lines[i], -lines[i + 8], lines[(i + 1) % 4 + 8]});
    }

    std::vector< int > surfaces(8, 0);
    for(auto i = 0; i < 8; ++i) {
      surfaces[i] = gmsh::model::geo::addSurfaceFilling({lineLoops[i]}, -1, center);
    }

    gmsh::model::addPhysicalGroup(2, surfaces, 1);
    gmsh::model::setPhysicalName(2, 1, "omega_" + _name);

    gmsh::model::geo::synchronize();

    for(auto i = 0; i < 12; ++i) {
      gmsh::model::mesh::setTransfiniteCurve(lines[i], _numPoints / 4 + 1);
    }

    gmsh::model::mesh::generate();

    _domain = domain::Domain(2, 1);

    gmsh::option::setNumber("General.Verbosity", verbose);
    _haveMesh = true;
  }


} // namespace gmshfem::post
