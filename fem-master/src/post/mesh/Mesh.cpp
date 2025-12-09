// GmshFEM - Copyright (C) 2019-2022, A. Royer, E. Béchet, C. Geuzaine, Université de Liège
//
// See the LICENSE.txt file for license information. Please report all
// issues on https://gitlab.onelab.info/gmsh/fem/issues

#include "Mesh.h"

#include <gmsh.h>

namespace gmshfem::post
{


  Mesh::Mesh(const std::string &name) :
    _name(name), _currentMesh(), _domain(), _haveMesh(false)
  {
    gmsh::model::getCurrent(_currentMesh);
  }

  Mesh::Mesh(Mesh &&other) :
    _name(std::move(other._name)), _currentMesh(std::move(other._currentMesh)), _domain(std::move(other._domain)), _haveMesh(other._haveMesh)
  {
  }

  Mesh::~Mesh()
  {
    destroyMesh();
  }

  std::string Mesh::name() const
  {
    return _name;
  }

  void Mesh::destroyMesh()
  {
    if(_haveMesh) {
      gmsh::model::setCurrent(_name);
      gmsh::model::remove();
      _haveMesh = false;
    }
  }

  void Mesh::setCurrent() const
  {
    if(_haveMesh) {
      gmsh::model::setCurrent(_name);
    }
  }

  void Mesh::restoreCurrent() const
  {
    gmsh::model::setCurrent(_currentMesh);
  }

  domain::Domain Mesh::getDomain() const
  {
    return _domain;
  }


} // namespace gmshfem::post
