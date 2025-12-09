#include <gmsh.h>
#include <gmshfem/Formulation.h>
#include <gmshfem/GmshFem.h>

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

  double lc = 1e-1;
  gmshFem.userDefinedParameter(lc, "lc");

  /****
   * Gmsh part
   ****/

  gmsh::model::add("poissonFunctions");

  gmsh::model::geo::addPoint(0., 0., 0., lc, 1); // Center

  // small circle
  gmsh::model::geo::addPoint(1., 0., 0., lc, 2);
  gmsh::model::geo::addPoint(0, 1., 0., lc, 3);
  gmsh::model::geo::addPoint(-1., 0., 0., lc, 4);
  gmsh::model::geo::addPoint(0, -1., 0., lc, 5);

  gmsh::model::geo::addCircleArc(2, 1, 3, 1);
  gmsh::model::geo::addCircleArc(3, 1, 4, 2);
  gmsh::model::geo::addCircleArc(4, 1, 5, 3);
  gmsh::model::geo::addCircleArc(5, 1, 2, 4);

  // big circle
  gmsh::model::geo::addPoint(5., 0., 0., lc, 6);
  gmsh::model::geo::addPoint(0, 5., 0., lc, 7);
  gmsh::model::geo::addPoint(-5., 0., 0., lc, 8);
  gmsh::model::geo::addPoint(0, -5., 0., lc, 9);

  gmsh::model::geo::addCircleArc(6, 1, 7, 5);
  gmsh::model::geo::addCircleArc(7, 1, 8, 6);
  gmsh::model::geo::addCircleArc(8, 1, 9, 7);
  gmsh::model::geo::addCircleArc(9, 1, 6, 8);

  // radius
  gmsh::model::geo::addLine(2, 6, 9);
  gmsh::model::geo::addLine(3, 7, 10);
  gmsh::model::geo::addLine(4, 8, 11);
  gmsh::model::geo::addLine(5, 9, 12);

  // CurveLoop
  gmsh::model::geo::addCurveLoop({9, 5, -10, -1}, 1);
  gmsh::model::geo::addCurveLoop({10, 6, -11, -2}, 2);
  gmsh::model::geo::addCurveLoop({11, 7, -12, -3}, 3);
  gmsh::model::geo::addCurveLoop({12, 8, -9, -4}, 4);

  // surface
  gmsh::model::geo::addPlaneSurface({1}, 10);
  gmsh::model::geo::addPlaneSurface({2}, 11);
  gmsh::model::geo::addPlaneSurface({3}, 12);
  gmsh::model::geo::addPlaneSurface({4}, 13);

  gmsh::model::geo::synchronize();

  // physicals
  gmsh::model::addPhysicalGroup(2, {10}, 1);
  gmsh::model::setPhysicalName(2, 1, "omega1");
  gmsh::model::addPhysicalGroup(2, {11}, 2);
  gmsh::model::setPhysicalName(2, 2, "omega2");
  gmsh::model::addPhysicalGroup(2, {12}, 3);
  gmsh::model::setPhysicalName(2, 3, "omega3");
  gmsh::model::addPhysicalGroup(2, {13}, 4);
  gmsh::model::setPhysicalName(2, 4, "omega4");
  gmsh::model::addPhysicalGroup(1, {1, 2, 3, 4}, 1);
  gmsh::model::setPhysicalName(1, 1, "gammaSmall");
  gmsh::model::addPhysicalGroup(1, {5, 6, 7, 8}, 2);
  gmsh::model::setPhysicalName(1, 2, "gammaBig");

  gmsh::model::occ::synchronize();
  gmsh::model::mesh::generate();

  /****
   * GmshFem part
   ******/

  //*****
  // Problem declaration
  //*****

  // I. Allocate the formulation object.
  Formulation< double > formulation("poissonFunctions");

  // II. Define physical regions (dim, tag).
  std::vector< Domain > omega(4);
  Domain omegaTot;
  for(unsigned int i = 0; i < 4; ++i) {
    omega[i] = Domain(2, i + 1);
    omegaTot |= omega[i];
  }
  Domain gammaSmall(1, 1);
  Domain gammaBig(1, 2);

  // III. Allocate field object.
  Field< double, Form::Form0 > v("v", omegaTot, FunctionSpaceTypeForm0::Lagrange);
  // Apply Dirichlet BC.
  ScalarFunction< double > fSmall = x< double >() + y< double >(); // fSmall(x,y) = x + y;
  ScalarFunction< double > fBig = sqrt((x< double >() + 10.) * (y< double >() + 10.)); // fSmall(x,y) = sqrt((x + 10) * (y + 10));
  v.addConstraint(gammaSmall, fSmall); // Apply fSmall on gammaSmall
  v.addConstraint(gammaBig, fBig); // Apply fBig on gammaBig

  ScalarPiecewiseFunction< double > pf;
  pf.addFunction(-1., omega[0]);
  pf.addFunction(cos(x< double >()), omega[1]);
  pf.addFunction(pow(y< double >(), 2) / 25., omega[2]);
  pf.addFunction(-1. / r2d< double >(), omega[3]);

  // IV. Write the corresponding weak formulation terms by terms.
  formulation.integral(grad(dof(v)), grad(tf(v)), omegaTot, "Gauss6");
  formulation.integral(10. * pf, tf(v), omegaTot, "Gauss6");

  // Prepro
  formulation.pre();

  // Pro
  formulation.assemble();
  formulation.solve();

  // V. Define and run post-processing operations.
  save(v);
  save(pf, omegaTot, "pf");

  return 0;
}
