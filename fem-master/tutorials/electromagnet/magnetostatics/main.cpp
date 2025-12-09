#include <gmsh.h>
#include <gmshfem/Formulation.h>
#include <gmshfem/GmshFem.h>
#include <gmshfem/Mesh.h>
#include <gmshfem/Message.h>

using namespace gmshfem;
using namespace gmshfem::common;
using namespace gmshfem::problem;
using namespace gmshfem::domain;
using namespace gmshfem::field;
using namespace gmshfem::function;
using namespace gmshfem::post;
using namespace gmshfem::equation;

int main(int argc, char **argv)
{
  GmshFem gmshFem(argc, argv);

  double ms = 1.;
  gmshFem.userDefinedParameter(ms, "meshSize");
  double murCore = 100.;
  gmshFem.userDefinedParameter(murCore, "murCore");
  double current = 0.01;
  gmshFem.userDefinedParameter(current, "i");
  std::string gauss = "Gauss6";
  gmshFem.userDefinedParameter(gauss, "gauss");
  double rInt = 200.e-3;
  double rExt = 250.e-3;
  double dxCore = 50.e-3;
  double dyCore = 100.e-3;
  double xInd = 75.e-3;
  double dxInd = 25.e-3;
  double dyInd = 50.e-3;
  bool axi = false;
  gmshFem.userDefinedParameter(axi, "axi");

  // Geometry
  const double p0 = 12.e-3 * ms;
  const double pCorex = 4.e-3 * ms;
  const double pCorey0 = 8.e-3 * ms;
  const double pCorey = 4.e-3 * ms;
  const double pIndx = 5.e-3 * ms;
  const double pIndy = 5.e-3 * ms;
  const double pInt = 12.5e-3 * ms;
  const double pExt = 12.5e-3 * ms;

  int p[12];
  p[0] = gmsh::model::geo::addPoint(0., 0., 0., p0);
  p[1] = gmsh::model::geo::addPoint(dxCore, 0., 0., pCorex);
  p[2] = gmsh::model::geo::addPoint(dxCore, dyCore, 0., pCorey);
  p[3] = gmsh::model::geo::addPoint(0., dyCore, 0., pCorey0);
  p[4] = gmsh::model::geo::addPoint(xInd, 0., 0., pIndx);
  p[5] = gmsh::model::geo::addPoint(xInd + dxInd, 0., 0., pIndx);
  p[6] = gmsh::model::geo::addPoint(xInd + dxInd, dyInd, 0., pIndy);
  p[7] = gmsh::model::geo::addPoint(xInd, dyInd, 0., pIndy);
  p[8] = gmsh::model::geo::addPoint(rInt, 0., 0., pInt);
  p[9] = gmsh::model::geo::addPoint(rExt, 0., 0., pExt);
  p[10] = gmsh::model::geo::addPoint(0., rInt, 0., pInt);
  p[11] = gmsh::model::geo::addPoint(0., rExt, 0., pExt);

  int l[15];
  l[0] = gmsh::model::geo::addLine(p[0], p[1]);
  l[1] = gmsh::model::geo::addLine(p[1], p[4]);
  l[2] = gmsh::model::geo::addLine(p[4], p[5]);
  l[3] = gmsh::model::geo::addLine(p[5], p[8]);
  l[4] = gmsh::model::geo::addLine(p[8], p[9]);
  l[5] = gmsh::model::geo::addLine(p[0], p[3]);
  l[6] = gmsh::model::geo::addLine(p[3], p[10]);
  l[7] = gmsh::model::geo::addLine(p[10], p[11]);
  l[8] = gmsh::model::geo::addLine(p[1], p[2]);
  l[9] = gmsh::model::geo::addLine(p[2], p[3]);
  l[10] = gmsh::model::geo::addLine(p[5], p[6]);
  l[11] = gmsh::model::geo::addLine(p[6], p[7]);
  l[12] = gmsh::model::geo::addLine(p[7], p[4]);
  l[13] = gmsh::model::geo::addCircleArc(p[8], p[0], p[10]);
  l[14] = gmsh::model::geo::addCircleArc(p[9], p[0], p[11]);

  int cl[4];
  cl[0] = gmsh::model::geo::addCurveLoop({-l[5], l[0], l[8], l[9]});
  cl[1] = gmsh::model::geo::addCurveLoop({l[10], l[11], l[12], l[2]});
  cl[2] = gmsh::model::geo::addCurveLoop({l[6], -l[13], -l[3], l[10], l[11], l[12], -l[1], l[8], l[9]});
  cl[3] = gmsh::model::geo::addCurveLoop({l[7], -l[14], -l[4], l[13]});

  int s[4];
  s[0] = gmsh::model::geo::addPlaneSurface({cl[0]});
  s[1] = gmsh::model::geo::addPlaneSurface({cl[1]});
  s[2] = gmsh::model::geo::addPlaneSurface({-cl[2]});
  s[3] = gmsh::model::geo::addPlaneSurface({-cl[3]});

  gmsh::model::geo::synchronize();

  gmsh::model::addPhysicalGroup(2, {s[2]}, 1);
  gmsh::model::setPhysicalName(2, 1, "Air");
  gmsh::model::addPhysicalGroup(2, {s[0]}, 2);
  gmsh::model::setPhysicalName(2, 2, "Core");
  gmsh::model::addPhysicalGroup(2, {s[1]}, 3);
  gmsh::model::setPhysicalName(2, 3, "Ind");
  gmsh::model::addPhysicalGroup(2, {s[3]}, 4);
  gmsh::model::setPhysicalName(2, 4, "AirInf");

  gmsh::model::addPhysicalGroup(1, {l[0], l[1], l[2], l[3], l[4]}, 1);
  gmsh::model::setPhysicalName(1, 1, "Left");
  gmsh::model::addPhysicalGroup(1, {l[5], l[6], l[7]}, 2);
  gmsh::model::setPhysicalName(1, 2, "Bottom");
  gmsh::model::addPhysicalGroup(1, {l[14]}, 3);
  gmsh::model::setPhysicalName(1, 3, "Inf");

  gmsh::model::geo::synchronize();
  gmsh::model::mesh::generate();

  // Formulation
  Domain air;
  Domain airInf;
  Domain core;
  Domain ind;

  if(axi) {
    air = Domain("Air", Axisymmetry< double >());
    airInf = Domain("AirInf", AxisymmetryShell< double >(rInt, rExt));
    core = Domain("Core", Axisymmetry< double >());
    ind = Domain("Ind", Axisymmetry< double >());
  }
  else {
    air = Domain("Air");
    airInf = Domain("AirInf", PolarShell< double >(rInt, rExt));
    core = Domain("Core");
    ind = Domain("Ind");
  }

  Domain inf("Inf");
  Domain left("Left");
  Domain bottom("Bottom");

  Formulation< double > formulation("Electromagnet");

  Field< double, Form::Form1 > a("a", air | airInf | core | ind, functionSpaceHCurl::P_Lagrange);
  a.addConstraint(bottom | inf, vector< double >(0., 0., 0.));

  const double mu0 = 4.e-7 * 3.14159265359;
  ScalarPiecewiseFunction< double > nu;
  nu.addFunction(1. / mu0, air | ind | airInf);
  nu.addFunction(1. / (murCore * mu0), core);

  const double nbTurns = 1000.;
  const double surface = integrate(ScalarFunction< double >(1.), ind, "Gauss1");
  const double js_z = -nbTurns * current / surface;
  const VectorFunction< double > js = vector< double >(0., 0., js_z);

  formulation.integral(nu * curl(dof(a)), curl(tf(a)), air | airInf | core | ind, gauss);
  formulation.integral(-js, tf(a), ind, gauss);

  formulation.pre();
  formulation.assemble();
  formulation.solve();

  save(a, air | airInf | core | ind, "a");
  save(zComp(a), air | airInf | core | ind, "az");
  save(curl(a), air | airInf | core | ind, "b");
  save(nu * curl(a), air | airInf | core | ind, "h");
  save(js, ind, "js");

  Line line("line", 0., 0.02, 0., rExt - 0.001, 0.02, 0., 30);
  save(curl(a), line, "by");

  return 0;
}
