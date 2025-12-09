#include <gmsh.h>
#include <gmshfem/Formulation.h>
#include <gmshfem/GmshFem.h>
#include <gmshfem/Message.h>
#include <gmshfem/PiecewiseFunction.h>

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

  bool withoutShell = false;
  gmshFem.userDefinedParameter(withoutShell, "withoutShell");
  double lc = 2e-1;
  gmshFem.userDefinedParameter(lc, "lc");
  double U = 2.;
  gmshFem.userDefinedParameter(U, "u");

  //
  // Gmsh part
  //

  double R = 2.5;
  double Rinf = 3.5; // used if withShell = 1

  gmsh::model::add("demoShell");

  int CC = gmsh::model::geo::addPoint(0., 0., 0., lc / 2.); // Center

  // left rectangle
  int pL[4];
  pL[0] = gmsh::model::geo::addPoint(-1., -1., 0., lc / 2.);
  pL[1] = gmsh::model::geo::addPoint(-0.3, -1., 0., lc / 4.);
  pL[2] = gmsh::model::geo::addPoint(-0.3, 1., 0., lc / 4.);
  pL[3] = gmsh::model::geo::addPoint(-1., 1., 0., lc / 2.);

  int lL[4];
  lL[0] = gmsh::model::geo::addLine(pL[0], pL[1]);
  lL[1] = gmsh::model::geo::addLine(pL[1], pL[2]);
  lL[2] = gmsh::model::geo::addLine(pL[2], pL[3]);
  lL[3] = gmsh::model::geo::addLine(pL[3], pL[0]);

  int clL = gmsh::model::geo::addCurveLoop({lL[0], lL[1], lL[2], lL[3]});
  int sL = gmsh::model::geo::addPlaneSurface({clL});

  // right rectangle
  int pR[4];
  pR[0] = gmsh::model::geo::addPoint(0.3, -1., 0., lc / 4.);
  pR[1] = gmsh::model::geo::addPoint(1., -1., 0., lc / 2.);
  pR[2] = gmsh::model::geo::addPoint(1., 1., 0., lc / 2.);
  pR[3] = gmsh::model::geo::addPoint(0.3, 1., 0., lc / 4.);

  int lR[4];
  lR[0] = gmsh::model::geo::addLine(pR[0], pR[1]);
  lR[1] = gmsh::model::geo::addLine(pR[1], pR[2]);
  lR[2] = gmsh::model::geo::addLine(pR[2], pR[3]);
  lR[3] = gmsh::model::geo::addLine(pR[3], pR[0]);

  int clR = gmsh::model::geo::addCurveLoop({lR[0], lR[1], lR[2], lR[3]});
  int sR = gmsh::model::geo::addPlaneSurface({clR});

  // dielectric
  int lD[4];
  lD[0] = gmsh::model::geo::addLine(pL[1], pR[0]);
  lD[1] = -lR[3];
  lD[2] = gmsh::model::geo::addLine(pR[3], pL[2]);
  lD[3] = -lL[1];

  int clD = gmsh::model::geo::addCurveLoop({lD[0], lD[1], lD[2], lD[3]});
  int sD = gmsh::model::geo::addPlaneSurface({clD});

  //circle
  int pC[4];
  pC[0] = gmsh::model::geo::addPoint(R, 0., 0., 2. * lc);
  pC[1] = gmsh::model::geo::addPoint(0, R, 0., 2. * lc);
  pC[2] = gmsh::model::geo::addPoint(-R, 0., 0., 2. * lc);
  pC[3] = gmsh::model::geo::addPoint(0, -R, 0., 2. * lc);

  int lC[4];
  lC[0] = gmsh::model::geo::addCircleArc(pC[0], CC, pC[1]);
  lC[1] = gmsh::model::geo::addCircleArc(pC[1], CC, pC[2]);
  lC[2] = gmsh::model::geo::addCircleArc(pC[2], CC, pC[3]);
  lC[3] = gmsh::model::geo::addCircleArc(pC[3], CC, pC[0]);

  int clC = gmsh::model::geo::addCurveLoop({lC[0], lC[1], lC[2], lC[3]});
  int sC = gmsh::model::geo::addPlaneSurface({clC, clR, clL, clD});

  int pS[4];
  int lS[4];
  int clS;
  int sS;
  if(!withoutShell) {
    pS[0] = gmsh::model::geo::addPoint(Rinf, 0., 0., 4. * lc);
    pS[1] = gmsh::model::geo::addPoint(0, Rinf, 0., 4. * lc);
    pS[2] = gmsh::model::geo::addPoint(-Rinf, 0., 0., 4. * lc);
    pS[3] = gmsh::model::geo::addPoint(0, -Rinf, 0., 4. * lc);

    lS[0] = gmsh::model::geo::addCircleArc(pS[0], CC, pS[1]);
    lS[1] = gmsh::model::geo::addCircleArc(pS[1], CC, pS[2]);
    lS[2] = gmsh::model::geo::addCircleArc(pS[2], CC, pS[3]);
    lS[3] = gmsh::model::geo::addCircleArc(pS[3], CC, pS[0]);

    clS = gmsh::model::geo::addCurveLoop({lS[0], lS[1], lS[2], lS[3]});
    sS = gmsh::model::geo::addPlaneSurface({clS, clC});
  }

  gmsh::model::geo::synchronize();

  // physicals
  gmsh::model::addPhysicalGroup(2, {sC}, 1);
  gmsh::model::setPhysicalName(2, 1, "omega");
  gmsh::model::addPhysicalGroup(2, {sL}, 2);
  gmsh::model::setPhysicalName(2, 2, "left");
  gmsh::model::addPhysicalGroup(2, {sR}, 3);
  gmsh::model::setPhysicalName(2, 3, "right");
  gmsh::model::addPhysicalGroup(2, {sD}, 4);
  gmsh::model::setPhysicalName(2, 4, "dielectric");

  if(!withoutShell) {
    gmsh::model::addPhysicalGroup(1, {lS[0], lS[1], lS[2], lS[3]}, 1);
    gmsh::model::setPhysicalName(1, 1, "gamma");

    gmsh::model::addPhysicalGroup(2, {sS}, 5);
    gmsh::model::setPhysicalName(2, 5, "omegaInf");
  }
  else {
    gmsh::model::addPhysicalGroup(1, {lC[0], lC[1], lC[2], lC[3]}, 1);
    gmsh::model::setPhysicalName(1, 1, "gamma");
  }

  gmsh::model::geo::synchronize();
  gmsh::model::mesh::generate();

  gmsh::model::mesh::setOrder(2);

  //
  // GmshFEM part
  //

  //
  // Problem declaration
  //

  // Allocate the formulation object.
  Formulation< double > formulation("demoShell");

  // Define physical regions (dim, tag).
  Domain omega("omega");
  Domain left("left");
  Domain right("right");
  Domain dielectric("dielectric");
  Domain gamma("gamma");
  if(!withoutShell) {
    Domain omegaInf("omegaInf", PolarShell< double >(R, Rinf));
    omega |= omegaInf;
  }

  const double eps0 = 8.85418782e-12;
  ScalarPiecewiseFunction< double > eps;
  eps.addFunction(eps0, omega);
  eps.addFunction(/*4.8*/ eps0, dielectric);

  // Allocate field object.
  Field< double, Form::Form0 > v("v", omega | dielectric, FunctionSpaceTypeForm0::Lagrange);
  // Apply Dirichlet BC.
  v.addConstraint(gamma, 0.);
  // Define the global quantity associated to the potential of the left electrode
  GlobalQuantity< double > vLeft("vLeftElectrode", left);
  v.assignGlobalQuantity(vLeft);
  // Define the global quantity associated to the potential of the right electrode
  GlobalQuantity< double > vRight("vRightElectrode", right);
  v.assignGlobalQuantity(vRight);

  // Write the corresponding weak formulation terms by terms.
  formulation.integral(eps * grad(dof(v)), grad(tf(v)), omega | dielectric, "Gauss4");
  formulation.globalTerm(vLeft, FixedComponent::Primal, U / 2.);
  formulation.globalTerm(vRight, FixedComponent::Primal, -U / 2.);

  // Prepro
  formulation.pre();

  // Pro
  formulation.assemble();
  formulation.solve();

  // Define and run post-processing operations.
  save(v);

  save(-grad(v), omega | dielectric | left | right, "e");

  Circle circle("circle", 0., 0., 0., 2., 0., 0., 1., 1000);
  double result = integrate(-grad(v) * tangent< double >(), circle, "Gauss20");
  save(-grad(v) * tangent< double >(), circle, "circulation");
  msg::info << "The circulation of the electric field 'e' on around the capacitor is about " << result << msg::endl;

  msg::info << "The charge on the left electrode is " << vLeft.getDualValue() << "[C]" << msg::endl;
  msg::info << "The charge on the right electrode is " << vRight.getDualValue() << "[C]" << msg::endl;

  return 0;
}
