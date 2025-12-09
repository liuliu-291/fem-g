// GmshFEM - Copyright (C) 2019-2022, A. Royer, E. Béchet, C. Geuzaine, Université de Liège
//
// See the LICENSE.txt file for license information. Please report all
// issues on https://gitlab.onelab.info/gmsh/fem/issues

#ifndef H_GMSHFEM_PLANE
#define H_GMSHFEM_PLANE

#include "Mesh.h"

namespace gmshfem::post
{


  class Plane : public Mesh
  {
   private:
    double _origin[3];
    double _vectorX[3];
    double _vectorY[3];
    unsigned int _numPoints[2];

   public:
    Plane(const std::string &name);
    Plane(const std::string &name, const double xOrigin, const double yOrigin, const double zOrigin, const double vecXX, const double vecXY, const double vecXZ, const double vecYX, const double vecYY, const double vecYZ, const unsigned int numPointsX, const unsigned int numPointsY);
    Plane(Plane &&other);
    ~Plane();

    void generateMesh() const override;
  };


} // namespace gmshfem::post


#endif // H_GMSHFEM_PLANE
