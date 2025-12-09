#include <gmsh.h>
#include <gmshfem/Formulation.h>
#include <gmshfem/GmshFem.h>
#include <gmshfem/Mesh.h>

using namespace gmshfem;
using namespace gmshfem::common;
using namespace gmshfem::equation;
using namespace gmshfem::domain;
using namespace gmshfem::field;
using namespace gmshfem::function;
using namespace gmshfem::post;
using namespace gmshfem::problem;

int main(int argc, char **argv)
{
  GmshFem gmshFem(argc, argv);

  double lc = 0.02;
  gmshFem.userDefinedParameter(lc, "lc");

  //****
  //* Gmsh part
  //****

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

  l[9] = gmsh::model::geo::addLine(p[3], p[6]);
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

  //****
  //* GmshFEM part
  //****

  Domain omega[4] = {Domain("omega_1"), Domain("omega_2"), Domain("omega_3"), Domain("omega_4")};
  Domain bottom[4] = {Domain("gammaBottom_1"), Domain("gammaBottom_2"), Domain("gammaBottom_3"), Domain("gammaBottom_4")};
  Domain right[4] = {Domain("gammaRight_1"), Domain("gammaRight_2"), Domain("gammaRight_3"), Domain("gammaRight_4")};
  Domain top[4] = {Domain("gammaTop_1"), Domain("gammaTop_2"), Domain("gammaTop_3"), Domain("gammaTop_4")};
  Domain left[4] = {Domain("gammaLeft_1"), Domain("gammaLeft_2"), Domain("gammaLeft_3"), Domain("gammaLeft_4")};

  Domain omegaTot = omega[0] | omega[1] | omega[2] | omega[3];

  // Non periodic reference formulation

  Formulation< double > formulationNP("non-periodic");

  Field< double, Form::Form0 > uNP("uNP", omegaTot, FunctionSpaceTypeForm0::Lagrange);
  uNP.addConstraint(bottom[2] | right[3] | top[1] | left[0], 1.);

  formulationNP.integral(grad(dof(uNP)), grad(tf(uNP)), omegaTot, "Gauss6");
  formulationNP.integral(1., tf(uNP), omegaTot, "Gauss6");

  formulationNP.pre();

  formulationNP.assemble();
  formulationNP.solve();

  save(uNP);

  // Periodic reference formulation (the lower right square)

  PeriodicLink link("gammaTop_3");

  Formulation< double > formulationP("periodic");

  Field< double, Form::Form0 > uP("uP", omega[2], FunctionSpaceTypeForm0::Lagrange);
  uP.addConstraint(bottom[2], 1.);
  uP.addPeriodicConstraint(link, 1.);

  formulationP.integral(grad(dof(uP)), grad(tf(uP)), omega[2], "Gauss6");
  formulationP.integral(1., tf(uP), omega[2], "Gauss6");

  formulationP.pre();

  formulationP.assemble();
  formulationP.solve();

  save(uP);

  // Check solution

  double diff = integrate(pow(norm(uNP - uP), 2), omega[2], "Gauss6");
  msg::info << "Difference between the periodic solution and the reference one: " << diff << msg::endl;
  save(uNP - uP, omega[2], "diff");

  return 0;
}
