// GmshFEM - Copyright (C) 2019-2022, A. Royer, E. Béchet, C. Geuzaine, Université de Liège
//
// See the LICENSE.txt file for license information. Please report all
// issues on https://gitlab.onelab.info/gmsh/fem/issues

#ifndef H_GMSHFEM_SPHERE
#define H_GMSHFEM_SPHERE

#include "Mesh.h"

namespace gmshfem::post
{


  class Sphere : public Mesh
  {
   private:
    double _center[3];
    double _radius;
    unsigned int _numPoints;

   public:
    Sphere(const std::string &name);
    Sphere(const std::string &name, const double xCenter, const double yCenter, const double zCenter, const double radius, const unsigned int numPoints);
    Sphere(Sphere &&other);
    ~Sphere();

    void generateMesh() const override;
  };


} // namespace gmshfem::post


#endif // H_GMSHFEM_SPHERE
