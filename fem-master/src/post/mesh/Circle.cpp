// GmshFEM - Copyright (C) 2019-2022, A. Royer, E. Béchet, C. Geuzaine, Université de Liège
//
// See the LICENSE.txt file for license information. Please report all
// issues on https://gitlab.onelab.info/gmsh/fem/issues

#include "Circle.h"

#include <gmsh.h>
#include <vector>

namespace gmshfem::post
{


  Circle::Circle(const std::string &name) :
    Mesh(name), _center{0., 0., 0.}, _normal{0., 0., 0.}, _numPoints(0)
  {
  }

  Circle::Circle(const std::string &name, const double xCenter, const double yCenter, const double zCenter, const double radius, const double xNormal, const double yNormal, const double zNormal, const unsigned int numPoints) :
    Mesh(name), _center{xCenter, yCenter, zCenter}, _normal{xNormal, yNormal, zNormal}, _radius(radius), _numPoints(numPoints)
  {
    const double norm = std::sqrt(_normal[0] * _normal[0] + _normal[1] * _normal[1] + _normal[2] * _normal[2]);
    _normal[0] /= norm;
    _normal[1] /= norm;
    _normal[2] /= norm;
  }

  Circle::Circle(Circle &&other) :
    Mesh(std::move(other)), _center{other._center[0], other._center[1], other._center[2]}, _normal{other._normal[0], other._normal[1], other._normal[2]}, _radius(other._radius), _numPoints(other._numPoints)
  {
  }

  Circle::~Circle()
  {
  }

  void Circle::generateMesh() const
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
    std::vector< int > points(4, 0);
    if(_normal[0] == 0. && _normal[1] == 0.) { // normal -> z
      const double sign = _normal[2] < 0. ? -1. : 1.;
      points[0] = gmsh::model::geo::addPoint(_center[0] + sign * _radius, _center[1], _center[2]);
      points[1] = gmsh::model::geo::addPoint(_center[0], _center[1] + sign * _radius, _center[2]);
      points[2] = gmsh::model::geo::addPoint(_center[0] - sign * _radius, _center[1], _center[2]);
      points[3] = gmsh::model::geo::addPoint(_center[0], _center[1] - sign * _radius, _center[2]);
    }
    else if(_normal[1] == 0. && _normal[2] == 0.) { // normal -> x
      const double sign = _normal[0] < 0. ? -1. : 1.;
      points[0] = gmsh::model::geo::addPoint(_center[0], _center[1], _center[2] - sign * _radius);
      points[1] = gmsh::model::geo::addPoint(_center[0], _center[1] + sign * _radius, _center[2]);
      points[2] = gmsh::model::geo::addPoint(_center[0], _center[1], _center[2] + sign * _radius);
      points[3] = gmsh::model::geo::addPoint(_center[0], _center[1] - sign * _radius, _center[2]);
    }
    else if(_normal[0] == 0. && _normal[2] == 0.) { // normal -> y
      const double sign = _normal[0] < 0. ? -1. : 1.;
      points[0] = gmsh::model::geo::addPoint(_center[0] + sign * _radius, _center[1], _center[2]);
      points[1] = gmsh::model::geo::addPoint(_center[0], _center[1], _center[2] - sign * _radius);
      points[2] = gmsh::model::geo::addPoint(_center[0] - sign * _radius, _center[1], _center[2]);
      points[3] = gmsh::model::geo::addPoint(_center[0], _center[1], _center[2] + sign * _radius);
    }
    else { // normal -> general case
      const double A = _normal[2] / std::cos(std::asin(_normal[0]));
      const double B = -_normal[1] / std::cos(std::asin(_normal[0]));
      const double C = std::cos(std::asin(_normal[0]));
      const double D = _normal[0];

      std::vector< std::vector< double > > mat{{C, 0., D},
                                               {B * D, A, -B * C},
                                               {-A * D, B, A * C}};
      std::vector< double > v1{mat[0][0] * _radius, mat[1][0] * _radius, mat[2][0] * _radius};
      std::vector< double > v2{mat[0][1] * _radius, mat[1][1] * _radius, mat[2][1] * _radius};
      points[0] = gmsh::model::geo::addPoint(_center[0] + v1[0], _center[1] + v1[1], _center[2] + v1[2]);
      points[1] = gmsh::model::geo::addPoint(_center[0] + v2[0], _center[1] + v2[1], _center[2] + v2[2]);
      points[2] = gmsh::model::geo::addPoint(_center[0] - v1[0], _center[1] - v1[1], _center[2] - v1[2]);
      points[3] = gmsh::model::geo::addPoint(_center[0] - v2[0], _center[1] - v2[1], _center[2] - v2[2]);
    }

    std::vector< int > lines(4, 0);
    for(auto i = 0; i < 4; ++i) {
      lines[i] = gmsh::model::geo::addCircleArc(points[i], center, points[(i + 1) % 4]);
    }

    gmsh::model::addPhysicalGroup(1, lines, 1);
    gmsh::model::setPhysicalName(1, 1, "omega_" + _name);

    gmsh::model::geo::synchronize();

    for(auto i = 0; i < 4; ++i) {
      gmsh::model::mesh::setTransfiniteCurve(lines[i], _numPoints / 4 + 1);
    }

    gmsh::model::mesh::generate();

    _domain = domain::Domain(1, 1);

    gmsh::option::setNumber("General.Verbosity", verbose);
    _haveMesh = true;
  }


} // namespace gmshfem::post
