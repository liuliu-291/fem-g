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

  double lc = 1e-1;
  gmshFem.userDefinedParameter(lc, "lc");

  //****
  //* Gmsh part
  //****

  gmsh::model::add("poissonNeumann");

  int pC = gmsh::model::geo::addPoint(0., 0., 0., lc); // Center

  // small circle
  int pS[4];
  pS[0] = gmsh::model::geo::addPoint(1., 0., 0., lc);
  pS[1] = gmsh::model::geo::addPoint(0, 1., 0., lc);
  pS[2] = gmsh::model::geo::addPoint(-1., 0., 0., lc);
  pS[3] = gmsh::model::geo::addPoint(0, -1., 0., lc);

  int lS[4];
  lS[0] = gmsh::model::geo::addCircleArc(pS[0], pC, pS[1]);
  lS[1] = gmsh::model::geo::addCircleArc(pS[1], pC, pS[2]);
  lS[2] = gmsh::model::geo::addCircleArc(pS[2], pC, pS[3]);
  lS[3] = gmsh::model::geo::addCircleArc(pS[3], pC, pS[0]);

  int llS = gmsh::model::geo::addCurveLoop({lS[0], lS[1], lS[2], lS[3]});

  // big circle
  int pB[4];
  pB[0] = gmsh::model::geo::addPoint(5., 0., 0., lc);
  pB[1] = gmsh::model::geo::addPoint(0, 5., 0., lc);
  pB[2] = gmsh::model::geo::addPoint(-5., 0., 0., lc);
  pB[3] = gmsh::model::geo::addPoint(0, -5., 0., lc);

  int lB[4];
  lB[0] = gmsh::model::geo::addCircleArc(pB[0], pC, pB[1]);
  lB[1] = gmsh::model::geo::addCircleArc(pB[1], pC, pB[2]);
  lB[2] = gmsh::model::geo::addCircleArc(pB[2], pC, pB[3]);
  lB[3] = gmsh::model::geo::addCircleArc(pB[3], pC, pB[0]);

  int llB = gmsh::model::geo::addCurveLoop({lB[0], lB[1], lB[2], lB[3]});

  // surface
  int s = gmsh::model::geo::addPlaneSurface({llB, llS});

  gmsh::model::geo::synchronize();

  // physicals
  gmsh::model::addPhysicalGroup(2, {s}, 1);
  gmsh::model::setPhysicalName(2, 1, "omega");
  gmsh::model::addPhysicalGroup(1, {lS[0], lS[1], lS[2], lS[3]}, 1);
  gmsh::model::setPhysicalName(1, 1, "gammaSmall");
  gmsh::model::addPhysicalGroup(1, {lB[0], lB[1], lB[2], lB[3]}, 2);
  gmsh::model::setPhysicalName(1, 2, "gammaBig");

  gmsh::model::geo::synchronize();
  gmsh::model::mesh::generate();

  //****
  //* GmshFEM part
  //****

  Formulation< double > formulation("poissonNeumann");

  Domain omega("omega");
  Domain gammaSmall("gammaSmall");
  Domain gammaBig("gammaBig");

  Field< double, Form::Form0 > v("v", omega, FunctionSpaceTypeForm0::Lagrange);
  v.addConstraint(gammaSmall, 1.);

  formulation.integral(grad(dof(v)), grad(tf(v)), omega, "Gauss6");
  formulation.integral(1., tf(v), omega, "Gauss6");
  formulation.integral(-1., tf(v), gammaBig, "Gauss6");

  formulation.pre();

  formulation.assemble();
  formulation.solve();

  save(v);

  return 0;
}
