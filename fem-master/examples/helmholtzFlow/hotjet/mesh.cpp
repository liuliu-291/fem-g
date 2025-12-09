#include "mesh.h"

void meshJet(const double sigma, const double SourceMeshRefinement, const int NrOfLayers, const int MeshElemOrder, std::pair<double, double> &xlim, std::pair<double, double> &ylim, double &lc)
{
  // gmsh::initialize();
  std::cout << " - Sigma = " << sigma << std::endl;
  std::cout << " - Number of PML layers = " << NrOfLayers << std::endl;
  std::cout << " - Source refinement factor = " << SourceMeshRefinement << std::endl;
  gmsh::model::add("geometry");
  
  // parameters
  lc = 0.6*sigma;
  double Lx = 7.*sigma;
  double Ly = 6.*sigma;

  double Cy = 3.*sigma;
  double Lx_down = 9.*sigma;

  int Nx = std::ceil((2.0*Lx+Lx_down)/lc);
  std::cout << "Nx = " << Nx << std::endl;
  // int Ny_long = std::ceil(Ly/lc);
  int Ny_short = std::ceil(((Ly-Cy))/lc);
  std::cout << "Ny_short = " << Ny_short << std::endl;
  double Lpml = NrOfLayers*lc;
  double MeshSizeSource = lc/SourceMeshRefinement;

  // geometry definition
  gmsh::model::geo::addPoint(0., 0., 0., MeshSizeSource, 1);
  gmsh::model::geo::addPoint(-4*sigma, 0., 0., MeshSizeSource, 501);
  gmsh::model::geo::addPoint(-3*sigma, 0., 0., MeshSizeSource, 502);
  gmsh::model::geo::addPoint(-2*sigma, 0., 0., MeshSizeSource, 503);
  for (int j=0; j<=15 ; j++) {
    gmsh::model::geo::addPoint((j+2)*sigma, 0., 0., MeshSizeSource, j+504);
  }
  gmsh::model::geo::addPoint(-sigma, 0., 0., MeshSizeSource, 520);
  gmsh::model::geo::addPoint(sigma, 0., 0., MeshSizeSource, 521);
  
  gmsh::model::geo::addPoint(-Lx, -Cy, 0., MeshSizeSource, 523);
  gmsh::model::geo::addPoint(Lx+Lx_down, -Cy, 0., MeshSizeSource, 524);
 
  gmsh::model::geo::addPoint(-Lx, -Ly, 0., lc, 2);
  gmsh::model::geo::addPoint(Lx+Lx_down, -Ly, 0, lc, 3);
  gmsh::model::geo::addPoint(Lx+Lx_down, Ly, 0, lc, 4);
  gmsh::model::geo::addPoint(-Lx, Ly, 0, lc, 5);

  gmsh::model::geo::addPoint(-Lx-Lpml, -Ly, 0, lc, 6);
  gmsh::model::geo::addPoint(-Lx, -Ly-Lpml, 0, lc, 7);
  gmsh::model::geo::addPoint(Lx+Lpml+Lx_down, -Ly, 0, lc, 8);
  gmsh::model::geo::addPoint(Lx+Lx_down, -Ly-Lpml, 0, lc, 9);
  
  gmsh::model::geo::addPoint(Lx+Lpml+Lx_down, Ly, 0, lc, 10);
  gmsh::model::geo::addPoint(Lx+Lx_down, Ly+Lpml, 0, lc, 11);
  gmsh::model::geo::addPoint(-Lx-Lpml, Ly, 0, lc, 12);
  gmsh::model::geo::addPoint(-Lx, Ly+Lpml, 0, lc, 13);

  gmsh::model::geo::addPoint(-Lx-Lpml, -Ly-Lpml, 0, lc, 14);
  gmsh::model::geo::addPoint(Lx+Lpml+Lx_down, -Ly-Lpml, 0, lc, 15);
  gmsh::model::geo::addPoint(Lx+Lpml+Lx_down, Ly+Lpml, 0, lc, 16);
  gmsh::model::geo::addPoint(-Lx-Lpml, Ly+Lpml, 0, lc, 17);

  gmsh::model::geo::addPoint(-Lx-Lpml, Cy, 0, lc, 49);
  gmsh::model::geo::addPoint(-Lx, Cy, 0, lc, 50);
  gmsh::model::geo::addPoint(Lx+Lx_down, Cy, 0, lc, 51);
  gmsh::model::geo::addPoint(Lx+Lpml+Lx_down, Cy, 0, lc, 52);
  
  // adding the lines
  int l50 = gmsh::model::geo::addLine(49,50);
  int l51 = gmsh::model::geo::addLine(50,51);
  int l52 = gmsh::model::geo::addLine(51,52);
  
  int l1 = gmsh::model::geo::addLine(2,3);
  int l2 = gmsh::model::geo::addLine(3,524);
  int l2222 = gmsh::model::geo::addLine(524,51);
  int l222 = gmsh::model::geo::addLine(51,4);
  int l3 = gmsh::model::geo::addLine(4,5);
  
  int l4 = gmsh::model::geo::addLine(50,523);
  int l4444 = gmsh::model::geo::addLine(523,2);
  int l124 = gmsh::model::geo::addLine(5,50);
  int l5 = gmsh::model::geo::addLine(7,9);
  int l6 = gmsh::model::geo::addLine(9,15);
  int l7 = gmsh::model::geo::addLine(15,8);
  int l8 = gmsh::model::geo::addLine(8,52);
  
  int l888 = gmsh::model::geo::addLine(52,10);
  int l9 = gmsh::model::geo::addLine(10,16);
  int l10 = gmsh::model::geo::addLine(16,11);
  int l11 = gmsh::model::geo::addLine(11,13);
  int l12 = gmsh::model::geo::addLine(13,17);
  int l13 = gmsh::model::geo::addLine(17,12);

  int l14 = gmsh::model::geo::addLine(6,49);
  int l114 = gmsh::model::geo::addLine(49,12);
  int l15 = gmsh::model::geo::addLine(6,14);
  int l16 = gmsh::model::geo::addLine(14,7);
  int l17 = gmsh::model::geo::addLine(3,9);
  int l18 = gmsh::model::geo::addLine(3,8);
  int l19 = gmsh::model::geo::addLine(4,10);
  int l20 = gmsh::model::geo::addLine(4,11);

  int l21 = gmsh::model::geo::addLine(5,13);
  int l22 = gmsh::model::geo::addLine(5,12);
  int l23 = gmsh::model::geo::addLine(2,6);
  int l24 = gmsh::model::geo::addLine(2,7);

  // curve loops
  int c1 = gmsh::model::geo::addCurveLoop({l1, l2, l2222, -l51, l4, l4444});
  int p1 = gmsh::model::geo::addPlaneSurface({c1});
  int c2 = gmsh::model::geo::addCurveLoop({l51, l222, l3, l124});
  int p2 = gmsh::model::geo::addPlaneSurface({c2});

  // curve loops
  int c22 = gmsh::model::geo::addCurveLoop({l18, l8, -l52, -l2222, -l2});
  int p7 = gmsh::model::geo::addPlaneSurface({c22});
  int c3 = gmsh::model::geo::addCurveLoop({l6, l7, -l18, l17});
  int p3 = gmsh::model::geo::addPlaneSurface({c3});

  int c4 = gmsh::model::geo::addCurveLoop({l19, l9, l10, -l20});
  int p4 = gmsh::model::geo::addPlaneSurface({c4});
  int c5 = gmsh::model::geo::addCurveLoop({l20, l11, -l21, -l3});
  int p5 = gmsh::model::geo::addPlaneSurface({c5});
  
  int c6 = gmsh::model::geo::addCurveLoop({l21, l12, l13, -l22});
  int p6 = gmsh::model::geo::addPlaneSurface({c6});
  int c8 = gmsh::model::geo::addCurveLoop({-l24, l16, l15, l23});
  int p8 = gmsh::model::geo::addPlaneSurface({c8});
  int c9 = gmsh::model::geo::addCurveLoop({l5, -l17, -l1, l24});
  int p9 = gmsh::model::geo::addPlaneSurface({c9});
  
  int c23 = gmsh::model::geo::addCurveLoop({l52, l888, -l19, -l222});
  int p10 = gmsh::model::geo::addPlaneSurface({c23});
  int c24 = gmsh::model::geo::addCurveLoop({-l23, -l14, -l50, -l4, -l4444});
  int p11 = gmsh::model::geo::addPlaneSurface({c24});
  int c25 = gmsh::model::geo::addCurveLoop({l50, -l124, l22, -l114});
  int p12 = gmsh::model::geo::addPlaneSurface({c25});

  // Transfinite curves
  gmsh::model::geo::synchronize();
  std::vector<int> lines1{l1, l5, l11, l3};
  for (unsigned int j=0; j<=lines1.size(); j++) {
    gmsh::model::geo::mesh::setTransfiniteCurve(lines1[j], Nx, "Progression", 1);
  }
  
  std::vector<int> lines2{l114, l124, l222, l888};
  for (unsigned int j=0; j<=lines2.size(); j++) {
    gmsh::model::geo::mesh::setTransfiniteCurve(lines2[j], Ny_short+1, "Progression", 1);
  }
  
  std::vector<int> lines3{l17, l18, l7, l6, l19, l9, l10, l20, l21, l12, l13, l22, l23, l15, l16, l24, l52, l50};
  for (unsigned int j=0; j<=lines3.size(); j++) {
    gmsh::model::geo::mesh::setTransfiniteCurve(lines3[j], NrOfLayers+1, "Progression", 1);
  }
  
  

  // embed points
  gmsh::model::geo::synchronize();
  gmsh::model::mesh::embed(0, {1,523,524}, 2, p1);
  for (int j=501; j<=521; j++) {
    gmsh::model::mesh::embed(0, {j}, 2, p1);
  }

  // translation
  std::vector<std::pair<int, int> > out{{2,p1}, {2,p2}, {2,p3}, {2,p4}, {2,p5}, {2,p6}, {2,p7}, {2,p8}, {2,p9}, {2,p10}, {2,p11}, {2,p12}};
  gmsh::model::geo::translate(out, 2.*sigma, 3.*sigma, 0.);
  
  // physical groups
  gmsh::model::addPhysicalGroup(1, {l1, l2, l2222, l222, l3, l124, l4, l4444}, 1);
  gmsh::model::setPhysicalName(1, 1, "gamma");
  gmsh::model::addPhysicalGroup(1, {l51}, 2);
  gmsh::model::setPhysicalName(1, 2, "control_line");

  gmsh::model::addPhysicalGroup(2, {p1,p2}, 1);
  gmsh::model::setPhysicalName(2, 1, "omega");
  gmsh::model::addPhysicalGroup(2, {p3, p4, p5, p6, p7, p8, p9, p10, p11, p12}, 2);
  gmsh::model::setPhysicalName(2, 2, "omega_pml");

  // generate mesh
  gmsh::model::mesh::generate(2);
  gmsh::model::mesh::setOrder(MeshElemOrder);
  gmsh::write("hotjet.msh");

  // pass physical domain boundary coordinates for gmshfem pml
  ylim.first = -Cy;
  ylim.second = Ly+Cy;
  xlim.first = -5*sigma;
  xlim.second = 18*sigma;
  // gmsh::finalize();
}