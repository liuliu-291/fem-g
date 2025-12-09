#include <gmsh.h>
#include <gmshfem/Formulation.h>
#include <gmshfem/GmshFem.h>
#include <gmshfem/Mesh.h>

using namespace std::literals::complex_literals;

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

  double lc = 5e-1;
  gmshFem.userDefinedParameter(lc, "lc");
  double pi = 3.14159265358979323846264338327950288419716939937510;
  double k = 2. * pi;
  gmshFem.userDefinedParameter(k, "k");
  std::string gauss = "Gauss8";
  gmshFem.userDefinedParameter(gauss, "gauss");

  //****
  //* Gmsh part
  //****

  gmsh::model::add("sphere");

  int s0 = gmsh::model::occ::addSphere(0., 0., 0., 3.);
  int s1 = gmsh::model::occ::addSphere(0., 0., 0., 5.);

  std::vector< std::pair< int, int > > outDimTags;
  std::vector< std::vector< std::pair< int, int > > > outDimTagsMap;
  gmsh::model::occ::cut({std::make_pair(3, s1)}, {std::make_pair(3, s0)}, outDimTags, outDimTagsMap);

  gmsh::option::setNumber("Mesh.MeshSizeFactor", lc);

  gmsh::model::occ::synchronize();

  gmsh::model::addPhysicalGroup(2, {s0}, 1);
  gmsh::model::setPhysicalName(2, 1, "gammaScat");
  gmsh::model::addPhysicalGroup(2, {s1}, 2);
  gmsh::model::setPhysicalName(2, 2, "gammaExt");
  for(unsigned int i = 0; i < outDimTags.size(); ++i) {
    if(outDimTags[i].first == 3) {
      gmsh::model::addPhysicalGroup(3, {outDimTags[i].second}, 3);
    }
  }
  gmsh::model::setPhysicalName(3, 3, "omega");

  gmsh::model::occ::synchronize();
  gmsh::model::mesh::generate();

  //****
  //* GmshFEM part
  //****

  Formulation< std::complex< double > > formulation("helmholtz");

  Domain omega("omega");
  Domain gammaScat("gammaScat");
  Domain gammaExt("gammaExt");

  ScalarFunction< std::complex< double > > planeWave = cos(k * x< std::complex< double > >()) - 1.i * sin(k * x< std::complex< double > >());

  Field< std::complex< double >, Form::Form0 > u("u", omega, FunctionSpaceTypeForm0::HierarchicalH1, 4);
  u.addConstraint(gammaScat, -planeWave);

  formulation.integral(-grad(dof(u)), grad(tf(u)), omega, gauss);
  formulation.integral(k * k * dof(u), tf(u), omega, gauss);
  formulation.integral(-1.i * k * dof(u), tf(u), gammaExt, gauss);

  formulation.pre();
  formulation.assemble();
  formulation.solve();

  const double eps = 0.05;
  save(u, omega, "u", "msh");
  Disk planeXY("plane XY", 0., 0., 0., 5. - eps, 0., 0., 1., 500);
  save(u, planeXY, "uXY", "pos");
  Disk planeXZ("plane XZ", 0., 0., 0., 5. - eps, 0., -1., 0., 500);
  save(u, planeXZ, "uXZ", "pos");
  Disk planeYZ("plane YZ", 0., 0., 0., 5. - eps, 1., 0., 0., 500);
  save(u, planeYZ, "uYZ", "pos");

  Sphere scatter("scatter", 0., 0., 0., 3. + eps, 300);
  save(u, scatter, "uScat", "pos");

  return 0;
}
