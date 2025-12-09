#include <complex>
#include <gmsh.h>
#include <gmshfem/AnalyticalFunction.h>
#include <gmshfem/Domain.h>
#include <gmshfem/Exception.h>
#include <gmshfem/FieldInterface.h>
#include <gmshfem/Formulation.h>
#include <gmshfem/Function.h>
#include <gmshfem/GmshFem.h>
#include <gmshfem/Message.h>
#include <gmshfem/Navier2D.h>
#include <gmshfem/io.h>

using namespace gmshfem;
using namespace gmshfem::analytics;
using namespace gmshfem::common;
using namespace gmshfem::term;
using namespace gmshfem::problem;
using namespace gmshfem::domain;
using namespace gmshfem::field;
using namespace gmshfem::function;
using namespace gmshfem::post;
using namespace gmshfem::equation;

void mesh(const double lc)
{
  gmsh::model::add("ground");

  const double width = 17. * 1e3; // 17 [km]
  const double height = 3.5 * 1e3; // 3.5 [km]
  int p[4];
  p[0] = gmsh::model::geo::addPoint(0., -800., 0., lc);
  p[1] = gmsh::model::geo::addPoint(0., -height, 0., lc);
  p[2] = gmsh::model::geo::addPoint(width, -height, 0., lc);
  p[3] = gmsh::model::geo::addPoint(width, -800., 0., lc);

  int pS = gmsh::model::geo::addPoint(width / 2., -height / 2., 0., lc);

  int l[4];
  l[0] = gmsh::model::geo::addLine(p[0], p[1]);
  l[1] = gmsh::model::geo::addLine(p[1], p[2]);
  l[2] = gmsh::model::geo::addLine(p[2], p[3]);
  l[3] = gmsh::model::geo::addLine(p[3], p[0]);

  int ll = gmsh::model::geo::addCurveLoop({l[0], l[1], l[2], l[3]});

  int s = gmsh::model::geo::addPlaneSurface({ll});

  gmsh::model::geo::synchronize();
  gmsh::model::mesh::embed(0, {pS}, 2, s);

  // physical groups
  gmsh::model::addPhysicalGroup(2, {s}, 1);
  gmsh::model::setPhysicalName(2, 1, "omega");
  gmsh::model::addPhysicalGroup(1, {l[3]}, 2);
  gmsh::model::setPhysicalName(1, 2, "gammaTop");
  gmsh::model::addPhysicalGroup(1, {l[0]}, 3);
  gmsh::model::setPhysicalName(1, 3, "gammaLeft");
  gmsh::model::addPhysicalGroup(1, {l[1]}, 4);
  gmsh::model::setPhysicalName(1, 4, "gammaBottom");
  gmsh::model::addPhysicalGroup(1, {l[2]}, 5);
  gmsh::model::setPhysicalName(1, 5, "gammaRight");
  gmsh::model::addPhysicalGroup(0, {pS}, 6);
  gmsh::model::setPhysicalName(0, 6, "source");

  // generate mesh
  gmsh::model::geo::synchronize();
  gmsh::model::mesh::generate();
}

void readData(const std::string &path, std::vector< double > &x, std::vector< double > &y, std::vector< std::vector< std::complex< double > > > &data)
{
  std::ifstream file(path, std::ios::in | std::ios::binary);

  unsigned int Nx = 0, Ny = 0;
  file.read((char *)&Nx, sizeof(unsigned int));
  file.read((char *)&Ny, sizeof(unsigned int));
  msg::info << "Reading: Nx = " << Nx << ", Ny = " << Ny << "." << msg::endl;

  x.resize(Nx);
  y.resize(Ny);
  file.read((char *)&x[0], Nx * sizeof(double));
  file.read((char *)&y[0], Ny * sizeof(double));

  data.resize(Nx, std::vector< std::complex< double > >(Ny));
  std::vector< double > rawData(Ny, 0.);
  for(unsigned int i = 0; i < Nx; ++i) {
    file.read((char *)&rawData[0], data[i].size() * sizeof(double));
    for(unsigned int j = 0; j < Ny; ++j) {
      data[i][j] = rawData[j];
    }
  }
  file.close();
}

int main(int argc, char **argv)
{
  GmshFem gmshFem(argc, argv);

  double lc = 30;
  gmshFem.userDefinedParameter(lc, "lc");
  std::string gauss = "Gauss10";
  gmshFem.userDefinedParameter(gauss, "gauss");
  double frequency = 0.002;
  gmshFem.userDefinedParameter(frequency, "frequency");
  bool compound = false;
  gmshFem.userDefinedParameter(compound, "compound");

  mesh(lc);

  Domain omega("omega");
  Domain gamma = Domain("gammaLeft") | Domain("gammaBottom") | Domain("gammaRight") | Domain("gammaTop");
  Domain source("source");

  std::vector< double > xD, yD, xP, yP, xS, yS;
  std::vector< std::vector< std::complex< double > > > dataD, dataP, dataS;

  readData("../data/density.dat", xD, yD, dataD);
  ScalarFunction< std::complex< double > > density = bilinearInterpolation< std::complex< double > >(&xD, &yD, &dataD) * 1000;
  save(density, omega, "density");
  readData("../data/vp.dat", xP, yP, dataP);
  ScalarFunction< std::complex< double > > vp = bilinearInterpolation< std::complex< double > >(&xP, &yP, &dataP);
  save(vp, omega, "vp");
  readData("../data/vs.dat", xS, yS, dataS);
  ScalarFunction< std::complex< double > > vs = bilinearInterpolation< std::complex< double > >(&xS, &yS, &dataS);
  save(vs, omega, "vs");

  const std::complex< double > im(0., 1.);
  const double pi = 3.14159265358979323846264338327950288;
  const double w = 2. * pi * frequency;
  ScalarFunction< std::complex< double > > nX = xComp(normal< std::complex< double > >());
  ScalarFunction< std::complex< double > > nY = yComp(normal< std::complex< double > >());

  ScalarFunction< std::complex< double > > mu = vs * vs * density;
  ScalarFunction< std::complex< double > > lambda = vp * vp * density - 2. * mu;
  ScalarFunction< std::complex< double > > lambda2mu = vp * vp * density;

  if(compound) {
    Formulation< std::complex< double > > formulation("2D-Navier (compound)");
    
    TensorFunction< std::complex< double > > C_xx = tensor< std::complex< double > >(lambda2mu, 0., 0.,
                                                                                      0., lambda, 0.,
                                                                                      0., 0., 0.);
    TensorFunction< std::complex< double > > C_xy = tensor< std::complex< double > >(0., mu, 0.,
                                                                                     mu, 0., 0.,
                                                                                     0., 0., 0.);
    TensorFunction< std::complex< double > > C_yx = tensor< std::complex< double > >(0., mu, 0.,
                                                                                     mu, 0., 0.,
                                                                                     0., 0., 0.);
    TensorFunction< std::complex< double > > C_yy = tensor< std::complex< double > >(lambda, 0., 0.,
                                                                                     0., lambda2mu, 0.,
                                                                                     0., 0., 0.);
                                                                                     
    TensorFunction< std::complex< double > > null = tensor< std::complex< double > >(0., 0., 0.,
                                                                                     0., 0., 0.,
                                                                                     0., 0., 0.);

    TensorFunction< std::complex< double >, 4 > C = tensor< std::complex< double >, 4 >(C_xx, C_xy, null,
                                                                                        C_yx, C_yy, null,
                                                                                        null, null, null);

    TensorFunction< std::complex< double > > tenVP = im * w * vp * density * tensor< std::complex< double > >(nX * nX, nX * nY, 0., nX * nY, nY * nY, 0., 0., 0., 0.);
    TensorFunction< std::complex< double > > tenVS = im * w * vs * density * tensor< std::complex< double > >(nY * nY, - (nX * nY), 0., - (nY * nX), nX * nX, 0., 0., 0., 0.);

    CompoundField< std::complex< double >, Form::Form0, 2 > u("u", omega, FunctionSpaceTypeForm0::HierarchicalH1, 4);
    u.addConstraint(source, vector< std::complex< double > >(1., 0., 0.));

    formulation.integral(C * grad(dof(u)), grad(tf(u)), omega, gauss);
    formulation.integral(-w * w * density * dof(u), tf(u), omega, gauss);

    // Lysmer Kuhlemeyer ABC
    formulation.integral(-tenVP * dof(u), tf(u), gamma, gauss);
    formulation.integral(-tenVS * dof(u), tf(u), gamma, gauss);

    formulation.pre();
    formulation.assemble();
    formulation.solve();

    save(xComp(u), omega, "ux_c", "pos");
    save(yComp(u), omega, "uy_c", "pos");
  }
  else {
    TensorFunction< std::complex< double > > C_xx = tensor< std::complex< double > >(lambda2mu, 0., 0.,
                                                                                     0., mu, 0.,
                                                                                     0., 0., 0.);
    TensorFunction< std::complex< double > > C_xy = tensor< std::complex< double > >(0., lambda, 0.,
                                                                                     mu, 0., 0.,
                                                                                     0., 0., 0.);
    TensorFunction< std::complex< double > > C_yx = tensor< std::complex< double > >(0., mu, 0.,
                                                                                     lambda, 0., 0.,
                                                                                     0., 0., 0.);
    TensorFunction< std::complex< double > > C_yy = tensor< std::complex< double > >(mu, 0., 0.,
                                                                                     0., lambda2mu, 0.,
                                                                                     0., 0., 0.);
                                                                                   
    Formulation< std::complex< double > > formulation("2D-Navier");

    Field< std::complex< double >, Form::Form0 > ux("ux", omega, FunctionSpaceTypeForm0::HierarchicalH1, 4);
    Field< std::complex< double >, Form::Form0 > uy("uy", omega, FunctionSpaceTypeForm0::HierarchicalH1, 4);
    ux.addConstraint(source, 1.);
    uy.addConstraint(source, 0.);

    formulation.integral(C_xx * grad(dof(ux)), grad(tf(ux)), omega, gauss);
    formulation.integral(C_xy * grad(dof(uy)), grad(tf(ux)), omega, gauss);
    formulation.integral(C_yx * grad(dof(ux)), grad(tf(uy)), omega, gauss);
    formulation.integral(C_yy * grad(dof(uy)), grad(tf(uy)), omega, gauss);

    formulation.integral(-w * w * density * dof(ux), tf(ux), omega, gauss);
    formulation.integral(-w * w * density * dof(uy), tf(uy), omega, gauss);

    // Lysmer Kuhlemeyer ABC
    formulation.integral(-im * w * vp * density * nX * nX * dof(ux), tf(ux), gamma, gauss);
    formulation.integral(-im * w * vp * density * nX * nY * dof(uy), tf(ux), gamma, gauss);
    formulation.integral(-im * w * vp * density * nX * nY * dof(ux), tf(uy), gamma, gauss);
    formulation.integral(-im * w * vp * density * nY * nY * dof(uy), tf(uy), gamma, gauss);

    formulation.integral(-im * w * vs * density * nY * nY * dof(ux), tf(ux), gamma, gauss);
    formulation.integral(im * w * vs * density * nX * nY * dof(uy), tf(ux), gamma, gauss);
    formulation.integral(-im * w * vs * density * nX * nX * dof(uy), tf(uy), gamma, gauss);
    formulation.integral(im * w * vs * density * nX * nY * dof(ux), tf(uy), gamma, gauss);

    formulation.pre();
    formulation.assemble();
    formulation.solve();

    save(ux, omega, "ux", "pos");
    save(uy, omega, "uy", "pos");
  }

  return 0;
}
