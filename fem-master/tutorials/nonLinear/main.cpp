#include <gmsh.h>
#include <gmshfem/Formulation.h>
#include <gmshfem/Function.h>
#include <gmshfem/GmshFem.h>
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

  double lc = 0.01;
  gmshFem.userDefinedParameter(lc, "lc");


  //*************************
  // Mesh
  //*************************

  // circle
  gmsh::model::geo::addPoint(0., 0., 0., lc, 10);

  // big circle
  gmsh::model::geo::addPoint(0.5, 0., 0., lc, 1);
  gmsh::model::geo::addPoint(0., 0.5, 0., lc, 2);
  gmsh::model::geo::addPoint(-0.5, 0., 0., lc, 3);
  gmsh::model::geo::addPoint(0., -0.5, 0., lc, 4);

  gmsh::model::geo::addCircleArc(1, 10, 2, 1);
  gmsh::model::geo::addCircleArc(2, 10, 3, 2);
  gmsh::model::geo::addCircleArc(3, 10, 4, 3);
  gmsh::model::geo::addCircleArc(4, 10, 1, 4);

  gmsh::model::geo::addCurveLoop({1, 2, 3, 4}, 1);

  // small circle
  const double R = 0.2;
  gmsh::model::geo::addPoint(R, 0., 0., lc, 11);
  gmsh::model::geo::addPoint(0., R, 0., lc, 12);
  gmsh::model::geo::addPoint(-R, 0., 0., lc, 13);
  gmsh::model::geo::addPoint(0., -R, 0., lc, 14);

  gmsh::model::geo::addCircleArc(11, 10, 12, 10);
  gmsh::model::geo::addCircleArc(12, 10, 13, 11);
  gmsh::model::geo::addCircleArc(13, 10, 14, 12);
  gmsh::model::geo::addCircleArc(14, 10, 11, 13);

  gmsh::model::geo::addCurveLoop({10, 11, 12, 13}, 2);

  gmsh::model::geo::addPlaneSurface({1, 2}, 10);
  gmsh::model::geo::addPlaneSurface({2}, 11);

  gmsh::model::geo::synchronize();

  gmsh::model::addPhysicalGroup(2, {10}, 1);
  gmsh::model::setPhysicalName(2, 1, "omegaOut");
  gmsh::model::addPhysicalGroup(2, {11}, 2);
  gmsh::model::setPhysicalName(2, 2, "omegaIn");
  gmsh::model::addPhysicalGroup(1, {1, 2, 3, 4}, 3);
  gmsh::model::setPhysicalName(1, 3, "gamma");

  gmsh::model::geo::synchronize();
  gmsh::model::mesh::generate();
  gmsh::write("mesh.msh");


  //*************************
  // Problem
  //*************************

  std::string gauss = "Gauss6";
  gmshFem.userDefinedParameter(gauss, "gauss");
  double threshold = 1e-10;
  gmshFem.userDefinedParameter(threshold, "th");
  bool useNewtonRaphson = false;
  gmshFem.userDefinedParameter(useNewtonRaphson, "nr");
  double relaxationFactor = 1.0;
  gmshFem.userDefinedParameter(relaxationFactor, "relax");
  double absTol = 1e-6, relTol = 1e-6;
  gmshFem.userDefinedParameter(absTol, "absTol");
  gmshFem.userDefinedParameter(relTol, "relTol");


  Formulation< double > formulation("Poisson");

  Domain omegaOut(2, 1);
  Domain omegaIn(2, 2);
  Domain gamma(1, 3);
  Domain omegaTot = omegaIn | omegaOut;

  Field< double, Form::Form0 > u("u", omegaTot, FunctionSpaceTypeForm0::Lagrange);
  u.addConstraint(gamma, 0.);

  ScalarFunction< double > Dx = exp(100. * u) + 0.1;
  ScalarFunction< double > d_Dx = 100. * exp(100. * u);
  ScalarPiecewiseFunction< double > D;
  D.addFunction(Dx, omegaIn);
  D.addFunction(1., omegaOut);
  ScalarPiecewiseFunction< double > d_D;
  d_D.addFunction(d_Dx, omegaIn);
  d_D.addFunction(0., omegaOut);

  algebra::Vector< double > dx;
  algebra::Vector< double > x;

  // u
  formulation.integral(-D * grad(dof(u)), grad(tf(u)), omegaTot, gauss);
  formulation.integral(+1., tf(u), omegaIn, gauss);

  if(useNewtonRaphson) {
    formulation.integral(-d_D * grad(u) * dof(u), grad(tf(u)), omegaTot, gauss);
    formulation.integral(d_D * grad(u) * u, grad(tf(u)), omegaTot, gauss);
  }

  formulation.pre();
  formulation.assemble();
  double residual = 0.;
  double residual0 = formulation.getResidual();
  unsigned int i = 0;
  do {
    formulation.solve();
    formulation.setSystemToZero();
    formulation.assemble();
    residual = formulation.getResidual();
    msg::info << "********************** Iteration " << ++i << ": abs = " << residual << ", rel = " << residual / residual0 << msg::endl;
  } while(residual > absTol && residual / residual0 > relTol);

  save(u);
  save(D, omegaTot, "D");

  return 0;
}
