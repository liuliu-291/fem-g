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

  gmsh::model::add("coupling");

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

  gmsh::model::geo::addCurveLoop({1, 2, 3, 4}, 1);

  // big circle
  gmsh::model::geo::addPoint(5., 0., 0., lc, 6);
  gmsh::model::geo::addPoint(0, 5., 0., lc, 7);
  gmsh::model::geo::addPoint(-5., 0., 0., lc, 8);
  gmsh::model::geo::addPoint(0, -5., 0., lc, 9);

  gmsh::model::geo::addCircleArc(6, 1, 7, 5);
  gmsh::model::geo::addCircleArc(7, 1, 8, 6);
  gmsh::model::geo::addCircleArc(8, 1, 9, 7);
  gmsh::model::geo::addCircleArc(9, 1, 6, 8);

  gmsh::model::geo::addCurveLoop({5, 6, 7, 8}, 2);

  // surface
  gmsh::model::geo::addPlaneSurface({2, 1}, 3);

  gmsh::model::geo::synchronize();

  // physicals
  gmsh::model::addPhysicalGroup(2, {3}, 1);
  gmsh::model::setPhysicalName(2, 1, "omega");
  gmsh::model::addPhysicalGroup(1, {1, 2, 3, 4}, 1);
  gmsh::model::setPhysicalName(1, 1, "gammaSmall");
  gmsh::model::addPhysicalGroup(1, {5, 6, 7, 8}, 2);
  gmsh::model::setPhysicalName(1, 2, "gammaBig");

  gmsh::model::geo::synchronize();
  gmsh::model::mesh::generate();

  /****
   * GmshFem part
   *****/

  //*****
  // Problem declaration
  //*****

  // I. Allocate the formulation object.
  Formulation< double > formulation("multiRHS");

  // II. Define physical regions (dim, tag).
  Domain omega(2, 1);
  Domain gammaSmall(1, 1);
  Domain gammaBig(1, 2);

  // III. Allocate field object.
  Field< double, Form::Form0 > v("v", omega, FunctionSpaceTypeForm0::Lagrange);
  // Apply Dirichlet BC.
  v.addConstraint(gammaSmall, 0.); // Apply 0 on gammaSmall
  v.addConstraint(gammaBig, 0.); // Apply 0 on gammaBig

  // IV. Write the corresponding weak formulation terms by terms.
  formulation.integral(grad(dof(v)), grad(tf(v)), omega, "Gauss6");
  formulation.integral(1., tf(v), omega, "Gauss6");

  // Prepro
  formulation.pre();

  // Pro
  formulation.assemble();
  formulation.solve();

  // V. Define and run post-processing operations.
  save(v);

  formulation.removeTerms();
  formulation.setRHSToZero();
  formulation.integral(-1., tf(v), omega, "Gauss6");

  formulation.assemble();
  formulation.solve(true);

  // Postpro
  save(v, "msh", "v2");

  return 0;
}
