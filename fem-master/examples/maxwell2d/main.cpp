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

  double lc = 1e-1;
  gmshFem.userDefinedParameter(lc, "lc");
  std::string gauss = "Gauss4";
  gmshFem.userDefinedParameter(gauss, "gauss");
  int dim = 2;
  gmshFem.userDefinedParameter(dim, "dim");

  gmsh::model::add("maxwellMesh");

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
  if(dim == 2 || dim == 0) {
    gmsh::model::addPhysicalGroup(2, {3}, 1);
    gmsh::model::setPhysicalName(2, 1, "omega");
    gmsh::model::addPhysicalGroup(1, {1, 2, 3, 4}, 1);
    gmsh::model::setPhysicalName(1, 1, "gammaSmall");
    gmsh::model::addPhysicalGroup(1, {5, 6, 7, 8}, 2);
    gmsh::model::setPhysicalName(1, 2, "gammaBig");
  }
  else if(dim == 3) {
    std::vector< std::pair< int, int > > outDimTags;
    gmsh::model::geo::extrude({std::make_pair(2, 3)}, 0., 0., 1., outDimTags);
    gmsh::model::addPhysicalGroup(3, {outDimTags[1].second}, 1);
    gmsh::model::setPhysicalName(3, 1, "omega");
    gmsh::model::addPhysicalGroup(2, {outDimTags[6].second, outDimTags[7].second, outDimTags[8].second, outDimTags[9].second}, 1);
    gmsh::model::setPhysicalName(2, 1, "gammaSmall");
    gmsh::model::addPhysicalGroup(2, {outDimTags[2].second, outDimTags[3].second, outDimTags[4].second, outDimTags[5].second}, 2);
    gmsh::model::setPhysicalName(2, 2, "gammaBig");
  }

  gmsh::model::geo::synchronize();
  gmsh::model::mesh::generate();

  //*****
  // Problem declaration
  //*****

  const double pi = 3.14159265359;
  const std::complex< double > im(0., 1.);

  double k = 2 * pi;
  gmshFem.userDefinedParameter(k, "k");

  bool TM = false;
  gmshFem.userDefinedParameter(TM, "tm");

  Domain omega(2, 1);
  Domain gammaSmall(1, 1);
  Domain gammaBig(1, 2);

  Formulation< std::complex< double > > formulation("maxwell");

  Field< std::complex< double >, Form::Form1 > e("e", omega | gammaBig | gammaSmall, (TM ? FunctionSpaceTypeForm1::HierarchicalHCurl : FunctionSpaceTypeForm1::P_Lagrange), 0);
  Field< std::complex< double >, Form::Form1 > lambda("lambda", gammaSmall, (TM ? FunctionSpaceTypeForm1::HierarchicalHCurl : FunctionSpaceTypeForm1::P_Lagrange), 0);

  formulation.integral(dof(lambda), tf(e), gammaSmall, gauss);
  formulation.integral(dof(e), tf(lambda), gammaSmall, gauss);
  VectorFunction< std::complex< double > > eInc = vector< std::complex< double > >(0., (TM ? 1. : 0.), (TM ? 0. : 1.)) * (cos< std::complex< double > >(k * x< std::complex< double > >()) + im * sin< std::complex< double > >(k * x< std::complex< double > >()));
  formulation.integral(eInc, tf(lambda), gammaSmall, gauss);

  formulation.integral(curl(dof(e)), curl(tf(e)), omega, gauss);
  formulation.integral(-k * k * dof(e), tf(e), omega, gauss);

  // Silver Muller
  VectorFunction< std::complex< double > > n = normal< std::complex< double > >();
  formulation.integral(-im * k * n % (dof(e) % n), tf(e), gammaBig, gauss);

  formulation.pre();

  formulation.assemble();
  formulation.solve();

  save(e, omega, "e");
  save(curl(e), omega, "curl_e");
  save(norm(e), omega, "norm_e");
  save(n % e, gammaBig, "trace_e");

  return 0;
}
