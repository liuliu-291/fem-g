// GmshFEM - Copyright (C) 2019-2022, A. Royer, E. Béchet, C. Geuzaine, Université de Liège
//
// See the LICENSE.txt file for license information. Please report all
// issues on https://gitlab.onelab.info/gmsh/fem/issues

#ifndef H_GMSHFEM_LINE
#define H_GMSHFEM_LINE

#include "Mesh.h"

namespace gmshfem::post
{


  class Line : public Mesh
  {
   private:
    double _min[3];
    double _max[3];
    unsigned int _numPoints;

   public:
    Line(const std::string &name);
    Line(const std::string &name, const double xMin, const double yMin, const double zMin, const double xMax, const double yMax, const double zMax, const unsigned int numPoints);
    Line(Line &&other);
    ~Line();

    void generateMesh() const override;
  };


} // namespace gmshfem::post


#endif // H_GMSHFEM_LINE
