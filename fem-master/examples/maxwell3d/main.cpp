#include <gmsh.h>
#include <gmshfem/Formulation.h>
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

  gmsh::model::add("waveguide");

  double L = 2;
  gmshFem.userDefinedParameter(L, "L");
  double a = 1;
  gmshFem.userDefinedParameter(a, "a");
  double b = 1;
  gmshFem.userDefinedParameter(b, "b");
  double lc = 0.1;
  gmshFem.userDefinedParameter(lc, "lc");

  // simple extruded geometry
  {
    std::vector< int > border(1);
    std::vector< int > omegas(1);
    int p[4] = {gmsh::model::geo::addPoint(0., 0, 0., lc),
                gmsh::model::geo::addPoint(0, a, 0., lc),
                gmsh::model::geo::addPoint(0, a, b, lc),
                gmsh::model::geo::addPoint(0, 0., b, lc)};
    int l[4] = {gmsh::model::geo::addLine(p[0], p[1]),
                gmsh::model::geo::addLine(p[2], p[1]),
                gmsh::model::geo::addLine(p[2], p[3]),
                gmsh::model::geo::addLine(p[3], p[0])};
    int cl = gmsh::model::geo::addCurveLoop({l[3], l[0], -l[1], l[2]});
    int s = gmsh::model::geo::addPlaneSurface({cl});
    gmsh::model::geo::mesh::setTransfiniteSurface(s);
    gmsh::model::geo::mesh::setRecombine(2, s);
    std::vector< std::pair< int, int > > e;
    gmsh::model::geo::extrude({std::make_pair(2, s)}, L, 0., 0., e,
                              {static_cast< int >((L) / lc)},
                              std::vector< double >(), true);
    gmsh::model::geo::synchronize();
    gmsh::model::addPhysicalGroup(2, {s}, 1);
    gmsh::model::setPhysicalName(2, 1, "gammaDir");
    gmsh::model::addPhysicalGroup(2, {e[0].second}, 2);
    gmsh::model::setPhysicalName(2, 2, "gammaInf");
    gmsh::model::addPhysicalGroup(2, {e[2].second, e[3].second, e[4].second, e[5].second}, 100);
    gmsh::model::setPhysicalName(2, 100, "border");
    gmsh::model::addPhysicalGroup(3, {e[1].second}, 1);
    gmsh::model::setPhysicalName(3, 1, "omega");
    gmsh::model::mesh::generate();
  }

  double k = 10;
  gmshFem.userDefinedParameter(k, "k");
  int order = 0;
  gmshFem.userDefinedParameter(order, "order");
  bool useLagMult = false;
  gmshFem.userDefinedParameter(useLagMult, "useLagMult");

  std::string gauss = "Gauss" + std::to_string(2 * (order + 1));
  double pi = 3.14159265359;
  double mode_m = 2.;
  double mode_n = 1.;
  double ky = mode_m * pi / a;
  double kz = mode_n * pi / b;
  double kc = sqrt(ky * ky + kz * kz);
  std::complex< double > im(0., 1.);
  std::complex< double > beta;
  if(-kc * kc + k * k >= 0) {
    beta = sqrt(-kc * kc + k * k);
  }
  else {
    beta = -im * sqrt(kc * kc - k * k);
  }
  VectorFunction< std::complex< double > > Einc = vector< std::complex< double > >(sin< std::complex< double > >(ky * y< std::complex< double > >()) *
                                                                                     sin< std::complex< double > >(kz * z< std::complex< double > >()),
                                                                                   im * beta * ky / (kc * kc) *
                                                                                     cos< std::complex< double > >(ky * y< std::complex< double > >()) *
                                                                                     sin< std::complex< double > >(kz * z< std::complex< double > >()),
                                                                                   im * beta * kz / (kc * kc) *
                                                                                     cos< std::complex< double > >(kz * z< std::complex< double > >()) *
                                                                                     sin< std::complex< double > >(ky * y< std::complex< double > >()));

  Domain omega(3, 1);
  Domain gammaInf(2, 2);
  Domain gammaDir(2, 1);
  Domain border(2, 100);

  gmshfem::problem::Formulation< std::complex< double > > formulation("maxwell");

  // E field
  Field< std::complex< double >, Form::Form1 >
    E("E", omega | gammaDir | gammaInf,
      FunctionSpaceTypeForm1::HierarchicalHCurl, order);
  E.addConstraint(border, vector< std::complex< double > >(0., 0., 0.));
  if(!useLagMult) // set Einc as Dirichlet BC using a pre-resolution on gammaDir
    E.addConstraint(gammaDir, Einc);

  // Volume terms
  formulation.integral(curl(dof(E)), curl(tf(E)), omega, gauss);
  formulation.integral(-k * k * dof(E), tf(E), omega, gauss);

  // Silver-Muller ABC on gammaInf
  formulation.integral(-im * k * dof(E), tf(E), gammaInf, gauss);

  // Impose Einc using Lagrange mutiplier?
  Field< std::complex< double >, Form::Form1 > lambda;
  if(useLagMult) {
    lambda = Field< std::complex< double >, Form::Form1 >("lambda", gammaDir,
                                                          FunctionSpaceTypeForm1::HierarchicalHCurl, order);
    // lagrange mult should be zero where E is fixed!
    lambda.addConstraint(border, vector< std::complex< double > >(0., 0., 0.));

    formulation.integral(dof(lambda), tf(E), gammaDir, gauss);
    formulation.integral(dof(E), tf(lambda), gammaDir, gauss);
    formulation.integral(Einc, tf(lambda), gammaDir, gauss);
  }

  formulation.pre();
  formulation.assemble();
  formulation.solve();
  save(E, omega, "E");

  return 0;
}
