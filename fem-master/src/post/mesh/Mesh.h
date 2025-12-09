// GmshFEM - Copyright (C) 2019-2022, A. Royer, E. Béchet, C. Geuzaine, Université de Liège
//
// See the LICENSE.txt file for license information. Please report all
// issues on https://gitlab.onelab.info/gmsh/fem/issues

#ifndef H_GMSHFEM_MESH
#define H_GMSHFEM_MESH

#include "Domain.h"

#include <string>

namespace gmshfem::post
{


  class Mesh
  {
   protected:
    const std::string _name;
    std::string _currentMesh;
    mutable domain::Domain _domain;
    mutable bool _haveMesh;

   public:
    Mesh(const std::string &name);
    Mesh(Mesh &&other);
    virtual ~Mesh();

    std::string name() const;

    virtual void generateMesh() const = 0;
    void destroyMesh();
    void setCurrent() const;
    void restoreCurrent() const;

    domain::Domain getDomain() const;
  };


} // namespace gmshfem::post

#include "Circle.h"
#include "Disk.h"
#include "Line.h"
#include "Plane.h"
#include "Sphere.h"

#endif // H_GMSHFEM_MESH
