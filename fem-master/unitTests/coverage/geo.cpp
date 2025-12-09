// GmshFEM - Copyright (C) 2019-2021 Université de Liège
//
// See the LICENSE.txt file for license information. Please report all
// issues on https://gitlab.onelab.info/gmsh/fem/issues

#include <gmsh.h>

#include "geo.h"

void removeGeo()
{
  gmsh::model::remove();
}

namespace Geo1D {


  void line()
  {
    gmsh::model::add("line");
    
    const int p0 = gmsh::model::geo::addPoint(0., 0., 0., 0.05);
    const int p1 = gmsh::model::geo::addPoint(1., 0., 0., 0.05);
    
    const int l = gmsh::model::geo::addLine(p0, p1);
    
    gmsh::model::geo::synchronize();

    gmsh::model::addPhysicalGroup(1, {l}, 1);
    gmsh::model::setPhysicalName(1, 1, "omega");
    gmsh::model::addPhysicalGroup(0, {p0}, 1);
    gmsh::model::setPhysicalName(0, 1, "dirichlet");
    gmsh::model::addPhysicalGroup(0, {p1}, 2);
    gmsh::model::setPhysicalName(0, 2, "neumann");
      
    gmsh::model::geo::synchronize();
    gmsh::model::mesh::generate();
  }
  
  
}

namespace Geo2D {


  void triangle()
  {
    gmsh::model::add("triangle");
    
    const int p0 = gmsh::model::geo::addPoint(0., 0., 0., 0.1);
    const int p1 = gmsh::model::geo::addPoint(1., 0., 0., 0.1);
    const int p2 = gmsh::model::geo::addPoint(1., 1., 0., 0.1);
    const int p3 = gmsh::model::geo::addPoint(0., 1., 0., 0.1);
    
    const int l0 = gmsh::model::geo::addLine(p0, p1);
    const int l1 = gmsh::model::geo::addLine(p1, p2);
    const int l2 = gmsh::model::geo::addLine(p2, p3);
    const int l3 = gmsh::model::geo::addLine(p3, p0);
    
    const int ll = gmsh::model::geo::addCurveLoop({l0, l1, l2, l3});
    const int s = gmsh::model::geo::addPlaneSurface({ll});
    
    gmsh::model::geo::synchronize();

    gmsh::model::addPhysicalGroup(2, {s}, 1);
    gmsh::model::setPhysicalName(2, 1, "omega");
    gmsh::model::addPhysicalGroup(1, {l0}, 1);
    gmsh::model::setPhysicalName(1, 1, "dirichlet");
    gmsh::model::addPhysicalGroup(1, {l2}, 2);
    gmsh::model::setPhysicalName(1, 2, "neumann");
      
    gmsh::model::geo::synchronize();
    gmsh::model::mesh::generate();
  }
  
  void quadrangle()
  {
    gmsh::model::add("triangle");
    
    const int p0 = gmsh::model::geo::addPoint(0., 0., 0., 0.1);
    const int p1 = gmsh::model::geo::addPoint(1., 0., 0., 0.1);
    const int p2 = gmsh::model::geo::addPoint(1., 1., 0., 0.1);
    const int p3 = gmsh::model::geo::addPoint(0., 1., 0., 0.1);
    
    const int l0 = gmsh::model::geo::addLine(p0, p1);
    const int l1 = gmsh::model::geo::addLine(p1, p2);
    const int l2 = gmsh::model::geo::addLine(p2, p3);
    const int l3 = gmsh::model::geo::addLine(p3, p0);
    
    const int ll = gmsh::model::geo::addCurveLoop({l0, l1, l2, l3});
    const int s = gmsh::model::geo::addPlaneSurface({ll});
    
    gmsh::model::geo::synchronize();

    gmsh::model::addPhysicalGroup(2, {s}, 1);
    gmsh::model::setPhysicalName(2, 1, "omega");
    gmsh::model::addPhysicalGroup(1, {l0}, 1);
    gmsh::model::setPhysicalName(1, 1, "dirichlet");
    gmsh::model::addPhysicalGroup(1, {l2}, 2);
    gmsh::model::setPhysicalName(1, 2, "neumann");
      
    gmsh::model::geo::synchronize();
    gmsh::model::mesh::generate();
    gmsh::model::mesh::recombine();
  }
  
  void periodic()
  {
    const double lc = 0.05;
    
    gmsh::model::add("periodic");
    
    int p[9];
    p[0] = gmsh::model::geo::addPoint(0., 0., 0., lc);
    p[1] = gmsh::model::geo::addPoint(0., 1., 0., lc);
    p[2] = gmsh::model::geo::addPoint(0., 2., 0., lc);
    
    p[3] = gmsh::model::geo::addPoint(1., 0., 0., lc);
    p[4] = gmsh::model::geo::addPoint(1., 1., 0., lc);
    p[5] = gmsh::model::geo::addPoint(1., 2., 0., lc);
    
    p[6] = gmsh::model::geo::addPoint(2., 0., 0., lc);
    p[7] = gmsh::model::geo::addPoint(2., 1., 0., lc);
    p[8] = gmsh::model::geo::addPoint(2., 2., 0., lc);
    
    int l[12];
    l[0] = gmsh::model::geo::addLine(p[0], p[1]);
    l[1] = gmsh::model::geo::addLine(p[1], p[2]);
    
    l[2] = gmsh::model::geo::addLine(p[3], p[4]);
    l[3] = gmsh::model::geo::addLine(p[4], p[5]);
    
    l[4] = gmsh::model::geo::addLine(p[6], p[7]);
    l[5] = gmsh::model::geo::addLine(p[7], p[8]);
    
    l[6] = gmsh::model::geo::addLine(p[0], p[3]);
    l[7] = gmsh::model::geo::addLine(p[1], p[4]);
    l[8] = gmsh::model::geo::addLine(p[2], p[5]);
    
    l[9]  = gmsh::model::geo::addLine(p[3], p[6]);
    l[10] = gmsh::model::geo::addLine(p[4], p[7]);
    l[11] = gmsh::model::geo::addLine(p[5], p[8]);
    
    
    int ll[4];
    ll[0] = gmsh::model::geo::addCurveLoop({l[6], l[2], -l[7], -l[0]});
    ll[1] = gmsh::model::geo::addCurveLoop({l[7], l[3], -l[8], -l[1]});
    ll[2] = gmsh::model::geo::addCurveLoop({l[9], l[4], -l[10], -l[2]});
    ll[3] = gmsh::model::geo::addCurveLoop({l[10], l[5], -l[11], -l[3]});
    
    int s[4];
    s[0] = gmsh::model::geo::addPlaneSurface({ll[0]});
    s[1] = gmsh::model::geo::addPlaneSurface({ll[1]});
    s[2] = gmsh::model::geo::addPlaneSurface({ll[2]});
    s[3] = gmsh::model::geo::addPlaneSurface({ll[3]});
    
    gmsh::model::geo::synchronize();
    
    std::vector< double > translation({0, -1, 0, 2, 1, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1});
    gmsh::model::mesh::setPeriodic(1, {l[10]}, {l[2]}, translation);
    
    // physicals
    gmsh::model::addPhysicalGroup(2, {s[0]}, 1);
    gmsh::model::setPhysicalName(2, 1, "omega_1");
    gmsh::model::addPhysicalGroup(2, {s[1]}, 2);
    gmsh::model::setPhysicalName(2, 2, "omega_2");
    gmsh::model::addPhysicalGroup(2, {s[2]}, 3);
    gmsh::model::setPhysicalName(2, 3, "omega_3");
    gmsh::model::addPhysicalGroup(2, {s[3]}, 4);
    gmsh::model::setPhysicalName(2, 4, "omega_4");
    
    gmsh::model::addPhysicalGroup(1, {l[6]}, 1);
    gmsh::model::setPhysicalName(1, 1, "gammaBottom_1");
    gmsh::model::addPhysicalGroup(1, {l[2]}, 2);
    gmsh::model::setPhysicalName(1, 2, "gammaRight_1");
    gmsh::model::addPhysicalGroup(1, {l[7]}, 3);
    gmsh::model::setPhysicalName(1, 3, "gammaTop_1");
    gmsh::model::addPhysicalGroup(1, {l[0]}, 4);
    gmsh::model::setPhysicalName(1, 4, "gammaLeft_1");
    
    gmsh::model::addPhysicalGroup(1, {l[7]}, 5);
    gmsh::model::setPhysicalName(1, 5, "gammaBottom_2");
    gmsh::model::addPhysicalGroup(1, {l[3]}, 6);
    gmsh::model::setPhysicalName(1, 6, "gammaRight_2");
    gmsh::model::addPhysicalGroup(1, {l[8]}, 7);
    gmsh::model::setPhysicalName(1, 7, "gammaTop_2");
    gmsh::model::addPhysicalGroup(1, {l[1]}, 8);
    gmsh::model::setPhysicalName(1, 8, "gammaLeft_2");
    
    gmsh::model::addPhysicalGroup(1, {l[9]}, 9);
    gmsh::model::setPhysicalName(1, 9, "gammaBottom_3");
    gmsh::model::addPhysicalGroup(1, {l[4]}, 10);
    gmsh::model::setPhysicalName(1, 10, "gammaRight_3");
    gmsh::model::addPhysicalGroup(1, {l[10]}, 11);
    gmsh::model::setPhysicalName(1, 11, "gammaTop_3");
    gmsh::model::addPhysicalGroup(1, {l[2]}, 12);
    gmsh::model::setPhysicalName(1, 12, "gammaLeft_3");
    
    gmsh::model::addPhysicalGroup(1, {l[10]}, 13);
    gmsh::model::setPhysicalName(1, 13, "gammaBottom_4");
    gmsh::model::addPhysicalGroup(1, {l[5]}, 14);
    gmsh::model::setPhysicalName(1, 14, "gammaRight_4");
    gmsh::model::addPhysicalGroup(1, {l[11]}, 15);
    gmsh::model::setPhysicalName(1, 15, "gammaTop_4");
    gmsh::model::addPhysicalGroup(1, {l[3]}, 16);
    gmsh::model::setPhysicalName(1, 16, "gammaLeft_4");
    
    gmsh::model::mesh::setTransfiniteSurface(s[0]);
    gmsh::model::mesh::setTransfiniteSurface(s[1]);
    gmsh::model::mesh::setTransfiniteSurface(s[2]);
    gmsh::model::mesh::setTransfiniteSurface(s[3]);
    
    gmsh::model::mesh::generate();
    gmsh::model::mesh::setOrder(2);
  }
  
  
}

namespace Geo3D {


  void tetrahedra()
  {
    gmsh::model::add("tetrahedra");
    
    const int p0 = gmsh::model::geo::addPoint(0., 0., 0., 1./3.);
    const int p1 = gmsh::model::geo::addPoint(1., 0., 0., 1./3.);
    const int p2 = gmsh::model::geo::addPoint(1., 1., 0., 1./3.);
    const int p3 = gmsh::model::geo::addPoint(0., 1., 0., 1./3.);
    
    const int p4 = gmsh::model::geo::addPoint(0., 0., 1., 1./3.);
    const int p5 = gmsh::model::geo::addPoint(1., 0., 1., 1./3.);
    const int p6 = gmsh::model::geo::addPoint(1., 1., 1., 1./3.);
    const int p7 = gmsh::model::geo::addPoint(0., 1., 1., 1./3.);
    
    
    const int l0 = gmsh::model::geo::addLine(p0, p1);
    const int l1 = gmsh::model::geo::addLine(p1, p2);
    const int l2 = gmsh::model::geo::addLine(p2, p3);
    const int l3 = gmsh::model::geo::addLine(p3, p0);
    
    const int l4 = gmsh::model::geo::addLine(p4, p5);
    const int l5 = gmsh::model::geo::addLine(p5, p6);
    const int l6 = gmsh::model::geo::addLine(p6, p7);
    const int l7 = gmsh::model::geo::addLine(p7, p4);
    
    const int l8 = gmsh::model::geo::addLine(p0, p4);
    const int l9 = gmsh::model::geo::addLine(p1, p5);
    const int l10 = gmsh::model::geo::addLine(p2, p6);
    const int l11 = gmsh::model::geo::addLine(p3, p7);
    
    
    const int ll1 = gmsh::model::geo::addCurveLoop({l0, l1, l2, l3}, -1, true);
    const int s1 = gmsh::model::geo::addPlaneSurface({ll1});
    
    const int ll2 = gmsh::model::geo::addCurveLoop({l4, l5, l6, l7}, -1, true);
    const int s2 = gmsh::model::geo::addPlaneSurface({ll2});
    
    const int ll3 = gmsh::model::geo::addCurveLoop({l0, l9, l4, l8}, -1, true);
    const int s3 = gmsh::model::geo::addPlaneSurface({ll3});
    
    const int ll4 = gmsh::model::geo::addCurveLoop({l1, l10, l5, l9}, -1, true);
    const int s4 = gmsh::model::geo::addPlaneSurface({ll4});
    
    const int ll5 = gmsh::model::geo::addCurveLoop({l2, l11, l6, l10}, -1, true);
    const int s5 = gmsh::model::geo::addPlaneSurface({ll5});
    
    const int ll6 = gmsh::model::geo::addCurveLoop({l3, l8, l7, l11}, -1, true);
    const int s6 = gmsh::model::geo::addPlaneSurface({ll6});
    
    
    const int sl = gmsh::model::geo::addSurfaceLoop({s1, s2, s3, s4, s5, s6});
    const int v = gmsh::model::geo::addVolume({sl});
    
    
    gmsh::model::geo::synchronize();

    gmsh::model::addPhysicalGroup(3, {v}, 1);
    gmsh::model::setPhysicalName(3, 1, "omega");
    gmsh::model::addPhysicalGroup(2, {s1}, 1);
    gmsh::model::setPhysicalName(2, 1, "dirichlet");
    gmsh::model::addPhysicalGroup(2, {l2}, 2);
    gmsh::model::setPhysicalName(2, 2, "neumann");
    gmsh::model::addPhysicalGroup(1, {l0}, 1);
    gmsh::model::setPhysicalName(1, 1, "line");
      
    gmsh::model::geo::synchronize();
    gmsh::model::mesh::generate();
  }
  
  
}
