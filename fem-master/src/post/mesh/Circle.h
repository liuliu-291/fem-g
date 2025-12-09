// GmshFEM - Copyright (C) 2019-2022, A. Royer, E. Béchet, C. Geuzaine, Université de Liège
//
// See the LICENSE.txt file for license information. Please report all
// issues on https://gitlab.onelab.info/gmsh/fem/issues

#ifndef H_GMSHFEM_CIRCLE
#define H_GMSHFEM_CIRCLE

#include "Mesh.h"

namespace gmshfem::post
{


  class Circle : public Mesh
  {
   private:
    double _center[3];
    double _normal[3];
    double _radius;
    unsigned int _numPoints;

   public:
    Circle(const std::string &name);
    Circle(const std::string &name, const double xCenter, const double yCenter, const double zCenter, const double radius, const double xNormal, const double yNormal, const double zNormal, const unsigned int numPoints);
    Circle(Circle &&other);
    ~Circle();

    void generateMesh() const override;
  };


} // namespace gmshfem::post


#endif // H_GMSHFEM_CIRCLE
