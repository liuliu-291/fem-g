// GmshFEM - Copyright (C) 2019-2021 Université de Liège
//
// See the LICENSE.txt file for license information. Please report all
// issues on https://gitlab.onelab.info/gmsh/fem/issues

#ifndef H_GMSHFEM_TEST_GEO
#define H_GMSHFEM_TEST_GEO


void removeGeo();

namespace Geo1D {


  void line();
  
  
}

namespace Geo2D {


  void triangle();
  void quadrangle();
  void periodic();
  
  
}

namespace Geo3D {


  void tetrahedra();
  
  
}

#endif // H_GMSHFEM_TEST_GEO
