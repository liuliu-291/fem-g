#include <gmsh.h>
#include <gmshfem/Formulation.h>
#include <gmshfem/GmshFem.h>

using namespace gmshfem;
using namespace gmshfem::algebra;
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
  std::string gauss = "Gauss6";
  gmshFem.userDefinedParameter(gauss, "gauss");
  double L = 1; // [m]
  gmshFem.userDefinedParameter(L, "L");
  unsigned int nMode = 5;
  gmshFem.userDefinedParameter(nMode, "nMode");
  bool withDamping = false;
  gmshFem.userDefinedParameter(withDamping, "withDamping");
  unsigned int basisOrder = 2;
  gmshFem.userDefinedParameter(basisOrder, "order");
  bool alsoSolve = false;
  gmshFem.userDefinedParameter(alsoSolve, "solve");

  {
    gmsh::model::add("square");

    std::vector< int > p(4);
    p[0] = gmsh::model::geo::addPoint(0., 0., 0., lc);
    p[1] = gmsh::model::geo::addPoint(L, 0., 0., lc);
    p[2] = gmsh::model::geo::addPoint(L, L, 0., lc);
    p[3] = gmsh::model::geo::addPoint(0., L, 0., lc);

    std::vector< int > l(4);
    for(unsigned int i = 0; i < p.size(); ++i) {
      l[i] = gmsh::model::geo::addLine(p[i], p[(i + 1) % p.size()]);
    }

    std::vector< int > ll(1);
    ll[0] = gmsh::model::geo::addCurveLoop(l);
    std::vector< int > s(1);
    s[0] = gmsh::model::geo::addPlaneSurface(ll);

    gmsh::model::geo::mesh::setTransfiniteSurface(s[0]);
    gmsh::model::geo::mesh::setRecombine(2, s[0]);

    gmsh::model::geo::synchronize();

    // physicals
    gmsh::model::addPhysicalGroup(2, s, 1);
    gmsh::model::setPhysicalName(2, 1, "square");
    gmsh::model::addPhysicalGroup(1, l, 1);
    gmsh::model::setPhysicalName(1, 1, "boundary");

    gmsh::model::geo::synchronize();
    gmsh::model::mesh::generate();
  }

  //*****
  // Problem declaration
  //*****

  Formulation< std::complex< double > > formulation("Laplacian modes");

  Domain omega("square");
  Domain boundary("boundary");

  Field< std::complex< double >, Form::Form0 > u("u", omega, FunctionSpaceTypeForm0::HierarchicalH1, basisOrder);
  u.addConstraint(boundary, 0.);

  // Mass
  formulation.integral(dt2_dof(u), tf(u), omega, gauss);
  // Damping
  if(withDamping) {
    formulation.integral(dt_dof(u), tf(u), omega, gauss);
  }
  // Stiffness
  formulation.integral(grad(dof(u)), grad(tf(u)), omega, gauss);

  formulation.pre();
  formulation.assemble(true);

  Vector< std::complex< double > > eigenvalues;
  formulation.eigensolve(eigenvalues, true, nMode);

  eigenvalues.save("eigenvalues");
  for(unsigned int i = 0; i < eigenvalues.size(); ++i) {
    msg::info << "Mode " << i << " of eigenvalue " << eigenvalues[i] << msg::endl;
    save(eigenfunction(u, i), omega, "mode_" + std::to_string(i));
  }

  if(alsoSolve) {
    double frequency = 10.;
    gmshFem.userDefinedParameter(frequency, "frequency");

    formulation.integral(-1., tf(u), omega, gauss);
    formulation.setAngularFrequency(frequency);

    formulation.removeSystem();
    formulation.initSystem();
    formulation.pre();
    formulation.assemble();
    formulation.solve();
    save(u, omega, "u");
  }

  return 0;
}
