#include <fstream>
#include <gmsh.h>
#include <gmshfem/Formulation.h>
#include <gmshfem/GmshFem.h>
#include <gmshfem/Message.h>
#include <iostream>

using namespace gmshfem;
using namespace gmshfem::common;
using namespace gmshfem::problem;
using namespace gmshfem::domain;
using namespace gmshfem::field;
using namespace gmshfem::function;
using namespace gmshfem::post;
using namespace gmshfem::equation;

// void getInputMode(int, int, double, double, std::complex<double>, double, ScalarFunction< std::complex<double> >, ScalarFunction< std::complex<double> >, ScalarFunction< std::complex<double> >);
void getInputMode(int m, int n, double Rout, double w, double M, std::complex< double > &kzmn, double &krmn, ScalarFunction< std::complex< double > > &psi, ScalarFunction< std::complex< double > > &dx_psi, ScalarFunction< std::complex< double > > &dy_psi)
{
  ScalarFunction< std::complex< double > > Umn, dUmndr;
  std::complex< double > im(0., 1.);
  double pi = 3.14159265358979323846264338327950288;
  // get Laplace-Beltrami eigenvalue associated to the mode (m,n)
  double data[50][30];
  std::ifstream table("../Bessel_circular_zeros_m50_n30.txt");
  if(Rout != 1.2) { // tabulated values for a given Rin, Rout !!
    msg::error << "The Laplace-Beltrami eigenvalues must be recomputed !" << msg::endl;
  }
  if(!table) {
    std::cout << "Error, file couldn't be opened" << std::endl;
  }
  else {
    for(int row = 0; row < 50; row++) { // stop loops if nothing to read
      for(int column = 0; column < 30; column++) {
        table >> data[row][column];
        // std::cout << "test: row: " << row << "column: " << column << "data: " << data[row][column] << std::endl;
        if(!table) {
          std::cout << "Error reading file for element " << row << "," << column << std::endl;
        }
      }
    }
  }
  krmn = data[m][n];

  double beta = sqrt(1 - M * M);
  double Argsqrt = w * w - beta * beta * krmn * krmn;
  if(Argsqrt >= 0) {
    kzmn = (1 / (beta * beta)) * (-M * w + sqrt(Argsqrt));
    msg::info << " - Input propagative mode " << msg::endl;
  }
  else {
    kzmn = (1 / (beta * beta)) * (-M * w - im * sqrt(abs(Argsqrt)));
    msg::info << " - Input evanescent mode " << msg::endl;
  }
  // print modeal info
  if(w > (beta * krmn) && w < krmn) {
    msg::info << " - Input inverse upstream mode ! " << msg::endl;
  }
  if(Argsqrt < 5.) {
    msg::warning << " - The mode is close to cut-off !" << msg::endl;
  }

  ScalarFunction< std::complex< double > > kr = krmn * y< std::complex< double > >();
  if(m == 0) {
    dUmndr = -krmn * cylBesselJ(1, kr);
  }
  else {
    dUmndr = krmn * 0.5 * (cylBesselJ(m - 1, kr) - cylBesselJ(m + 1, kr));
  }
  Umn = cylBesselJ(m, kr);

  // mode amplitude for validation
  double a = 0.020991195871490; // (8,2)
  // double a = 0.030601838437059; // (8,2) no flow M=0
  // double a = 0.012923264535242; // (8,0)
  // double a = 0.013887559096146; // (8,1)
  std::complex< double > A(a * cos(pi / 2), a * sin(pi / 2)); // optional modal phase shift (normalization)
  psi = A * exp(-im * kzmn * x< std::complex< double > >()) * Umn;
  dx_psi = -A * im * kzmn * exp(-im * kzmn * x< std::complex< double > >()) * Umn;
  dy_psi = A * exp(-im * kzmn * x< std::complex< double > >()) * dUmndr;
}

void meshDuct(const double L, const double Lliner, const double Rout, const double h, const int Rfactor, const int Npml, const int NpmlActive, const int MeshOrder)
{
  gmsh::model::add("geometry");
  // duct definition
  gmsh::model::geo::addPoint(0, 0., 0., h, 1);
  gmsh::model::geo::addPoint(L, 0., 0., h, 2);
  gmsh::model::geo::addPoint(L, Rout, 0., h, 3);
  gmsh::model::geo::addPoint(0, Rout, 0., h, 4);
  gmsh::model::geo::addPoint(0.5 * (L - Lliner), Rout, 0., h, 5);
  gmsh::model::geo::addPoint(0.5 * (L + Lliner), Rout, 0., h, 6);

  // axis
  gmsh::model::geo::addLine(1, 2, 1);
  // outlet
  gmsh::model::geo::addLine(3, 2, 2);
  // wall
  gmsh::model::geo::addLine(3, 6, 3);
  // liner
  gmsh::model::geo::addLine(6, 5, 4);
  // wall
  gmsh::model::geo::addLine(5, 4, 5);
  // inlet
  gmsh::model::geo::addLine(4, 1, 6);

  gmsh::model::geo::addCurveLoop({6, 1, -2, 3, 4, 5}, 1);
  gmsh::model::geo::addPlaneSurface({1}, 1);

  // refine mesh along the upper wall
  int N = std::round(Lliner / h);
  int Nliner = std::round(0.5 * (L - Lliner) / h);
  gmsh::model::geo::mesh::setTransfiniteCurve(4, Rfactor * N);
  gmsh::model::geo::mesh::setTransfiniteCurve(3, Rfactor * Nliner);
  gmsh::model::geo::mesh::setTransfiniteCurve(5, Rfactor * Nliner);

  // mesh the pml accordingly
  double growth = 1.1;
  double hLiner = h / Rfactor;
  double temp = (1 / std::log10(growth)) * std::log10(1 - ((1 - growth) * Rout / hLiner));
  double N2 = std::round(temp) + 1;

  if(Npml != 0) {
    // pml extension
    gmsh::model::geo::addPoint(L + Npml * h, Rout, 0., h, 7);
    gmsh::model::geo::addPoint(L + Npml * h, 0., 0., h, 8);
    gmsh::model::geo::addLine(2, 8, 8);
    gmsh::model::geo::addLine(7, 8, 9);
    gmsh::model::geo::addLine(7, 3, 10);

    gmsh::model::geo::addCurveLoop({8, -9, 10, 2}, 2);
    gmsh::model::geo::addPlaneSurface({2}, 2);
    gmsh::model::geo::mesh::setTransfiniteCurve(2, N2, "Progression", growth);
    gmsh::model::geo::mesh::setTransfiniteCurve(9, N2, "Progression", growth);

    gmsh::model::geo::mesh::setTransfiniteSurface(2);
    gmsh::model::geo::mesh::setRecombine(2, 2);

    gmsh::model::addPhysicalGroup(1, {8}, 8);
    gmsh::model::setPhysicalName(1, 8, "PmlBottom");
    gmsh::model::addPhysicalGroup(2, {2}, 2);
    gmsh::model::setPhysicalName(2, 2, "omegaPml");
  }
  if(NpmlActive != 0) {
    // active pml extension
    gmsh::model::geo::addPoint(-NpmlActive * h, Rout, 0., h, 9);
    gmsh::model::geo::addPoint(-NpmlActive * h, 0., 0., h, 10);
    gmsh::model::geo::addLine(1, 10, 11);
    gmsh::model::geo::addLine(9, 10, 12);
    gmsh::model::geo::addLine(9, 4, 13);

    gmsh::model::geo::addCurveLoop({11, -12, 13, 6}, 3);
    gmsh::model::geo::addPlaneSurface({3}, 3);
    gmsh::model::geo::mesh::setTransfiniteCurve(6, N2, "Progression", growth);
    gmsh::model::geo::mesh::setTransfiniteCurve(12, N2, "Progression", growth);

    gmsh::model::geo::mesh::setTransfiniteSurface(3);
    gmsh::model::geo::mesh::setRecombine(3, 3);

    gmsh::model::addPhysicalGroup(1, {11}, 11);
    gmsh::model::setPhysicalName(1, 11, "PmlBottomA");
    gmsh::model::addPhysicalGroup(2, {3}, 3);
    gmsh::model::setPhysicalName(2, 3, "omegaPmlActive");
  }

  // physical groups
  gmsh::model::addPhysicalGroup(2, {1}, 1);
  gmsh::model::setPhysicalName(2, 1, "omega");
  gmsh::model::addPhysicalGroup(1, {6}, 1);
  gmsh::model::setPhysicalName(1, 1, "Inlet");
  gmsh::model::addPhysicalGroup(1, {2}, 2);
  gmsh::model::setPhysicalName(1, 2, "Outlet");
  gmsh::model::addPhysicalGroup(1, {4}, 3);
  gmsh::model::setPhysicalName(1, 3, "Liner");
  gmsh::model::addPhysicalGroup(1, {1}, 4);
  gmsh::model::setPhysicalName(1, 4, "Bottom");

  // generate mesh
  gmsh::model::geo::synchronize();
  gmsh::model::mesh::generate(2);
  gmsh::model::mesh::setOrder(MeshOrder);
  gmsh::option::setNumber("Mesh.ElementOrder", MeshOrder);
}

void getPmlFunctions(double sigma0, double Lmax, double Lmin, double Xmin, double Xmax, double k0, ScalarFunction< std::complex< double > > &gamma)
{
  std::complex< double > im(0., 1.);
  ScalarFunction< std::complex< double > > Sigma = heaviside(x< std::complex< double > >() - Xmax) * sigma0 / (Lmax - abs(x< std::complex< double > >())) + heaviside(-x< std::complex< double > >() + Xmin) * sigma0 / (abs(Lmin) - abs(x< std::complex< double > >()));
  gamma = 1 - (im / k0) * Sigma;
}

int main(int argc, char **argv)
{
  GmshFem gmshFem(argc, argv);
  double pi = 3.14159265358979323846264338327950288;
  const std::complex< double > im(0., 1.);

  // meshing
  double L = 2.08;
  double Rout = 1.2;
  double Lliner = 0.9 * Rout;
  double h = 0.1 * Rout; // meshsize
  int Rfactor = 5.0; // refinement factor along the liner
  int MeshOrder = 1;
  int Npml = 4; // number of PML layers
  int NpmlActive = Npml; // number of active PML layers
  meshDuct(L, Lliner, Rout, h, Rfactor, Npml, NpmlActive, MeshOrder);

  // Problem declaration
  int FEMorder = 5;
  std::string gauss = "Gauss" + std::to_string(2 * FEMorder + 1);
  double M = -0.3; // Mean flow
  double beta = sqrt(1 - M * M);
  double rho0 = 1.225; // density
  double k0 = 18.0 / Rout; // normalized wavenumber
  int m = 8;
  int n = 2; // azimuthal and radial mode order

  // Domains declaration
  Domain omega;
  Domain bottom;
  if(NpmlActive != 0) {
    omega = Domain(2, 1) | Domain(2, 2) | Domain(2, 3);
    bottom = Domain(1, 4) | Domain(1, 8) | Domain(1, 11);
  }
  else {
    omega = Domain(2, 1) | Domain(2, 2);
    bottom = Domain(1, 4) | Domain(1, 8);
  }
  Domain inlet = Domain(1, 1);
  Domain outlet = Domain(1, 2);
  Domain liner = Domain(1, 3);

  // Get Input Modal function and wavenumber
  std::complex< double > kx; // propagative wavenumber
  double ky; // Laplace-Beltrami eigenvalue
  ScalarFunction< std::complex< double > > psi;
  ScalarFunction< std::complex< double > > dx_psi;
  ScalarFunction< std::complex< double > > dy_psi;
  getInputMode(m, n, Rout, k0, M, kx, ky, psi, dx_psi, dy_psi);
  msg::info << " - propagative wavenumber kx = " << kx << ", transverse wavenumber ky = " << ky << msg::endl;

  float pointsByWl = (2 * pi * FEMorder * (1 - abs(M))) / (h * k0);
  msg::info << "********************************" << msg::endl;
  msg::info << " - convected Helmholtz duct problem " << msg::endl;
  msg::info << " - FEM basis order = " << FEMorder << "" << msg::endl;
  msg::info << " - max. dofs by wavelength = " << pointsByWl << msg::endl;
  msg::info << " - normalized wave number = " << k0 << msg::endl;
  msg::info << "********************************" << msg::endl;

  Formulation< std::complex< double > > formulation("helmholtzflow");
  Field< std::complex< double >, Form::Form0 > v("v", omega | inlet | outlet | liner | bottom, FunctionSpaceTypeForm0::HierarchicalH1, FEMorder);

  // modified density for cylindrical formulation (y <-> r)
  ScalarFunction< std::complex< double > > rho0S = rho0 * y< std::complex< double > >();
  // Dirichlet constraint on central axis for non-zero azimuthal mode order
  if(m != 0) {
    v.addConstraint(bottom, 0.);
  }

  // Active PML - volume terms
  const double LpmlA = -NpmlActive * h;
  const double Lpml = L + Npml * h;
  double sigma0 = 1; //4*beta*beta;
  msg::info << " - use hyperbolic stable PML with sigma0 = " << sigma0 << msg::endl;
  ScalarFunction< std::complex< double > > gamma;
  getPmlFunctions(sigma0, Lpml, LpmlA, 0., L, k0, gamma);

  // Input BC
  if(NpmlActive != 0) {
    // Active PML - incident field
    formulation.integral(-rho0S * beta * beta / gamma * dx_psi, vector< std::complex< double > >(1., 0., 0.) * grad(tf(v)), Domain(2, 3), gauss);
    formulation.integral(-rho0S * gamma * dy_psi, vector< std::complex< double > >(0., 1., 0.) * grad(tf(v)), Domain(2, 3), gauss);
    formulation.integral(-rho0S * im * k0 * M * dx_psi / gamma, tf(v), Domain(2, 3), gauss);
    formulation.integral(+rho0S * im * k0 * M / gamma * psi, vector< std::complex< double > >(1., 0., 0.) * grad(tf(v)), Domain(2, 3), gauss);
    formulation.integral(+rho0S * k0 * k0 * (gamma - M * M / gamma) / (beta * beta) * psi, tf(v), Domain(2, 3), gauss);
    formulation.integral(-rho0S * pow(m / y< std::complex< double > >(), 2) * gamma * psi, tf(v), Domain(2, 3), gauss);
    // Flux term
    formulation.integral(rho0S * beta * beta * dx_psi, tf(v), inlet, gauss);
    formulation.integral(-rho0S * im * k0 * M * psi, tf(v), inlet, gauss);
  }
  else {
    // Usual Input BC
    formulation.integral(rho0S * beta * beta * dx_psi, tf(v), inlet, gauss);
    formulation.integral(-rho0S * im * k0 * M * dof(v), tf(v), inlet, gauss);
  }

  // weak form for all domains - stable PML for the 2.5D convected Helmholtz equation
  formulation.integral(rho0S * vector< std::complex< double > >(beta * beta / gamma, 0., 0.) * grad(dof(v)), vector< std::complex< double > >(1., 0., 0.) * grad(tf(v)), omega, gauss);
  formulation.integral(rho0S * vector< std::complex< double > >(0., gamma, 0.) * grad(dof(v)), vector< std::complex< double > >(0., 1., 0.) * grad(tf(v)), omega, gauss);
  formulation.integral(rho0S * vector< std::complex< double > >(im * k0 * M / gamma, 0., 0.) * grad(dof(v)), tf(v), omega, gauss);
  formulation.integral(-rho0S * im * k0 * M / gamma * dof(v), vector< std::complex< double > >(1., 0., 0.) * grad(tf(v)), omega, gauss);
  formulation.integral(-rho0S * k0 * k0 * (gamma - (M * M / gamma)) / (beta * beta) * dof(v), tf(v), omega, gauss);
  formulation.integral(rho0S * pow(m / y< std::complex< double > >(), 2) * gamma * dof(v), tf(v), omega, gauss);

  // Myer's BC on Liner
  double R = 2; // resistance
  double Mf = 0.020; // inertance
  double hh = 0.02286; // liner depth
  std::complex< double > Z(R, -(1 / std::tan(k0 * hh)) + k0 * Mf); // normalized impedance (by rho0*c0)
  formulation.integral(im * rho0S * k0 / Z * dof(v), tf(v), liner, gauss);
  formulation.integral(-rho0S / Z * dof(v), vector< std::complex< double > >(M, 0., 0.) * grad(tf(v)), liner, gauss);
  formulation.integral(+rho0S / Z * vector< std::complex< double > >(M, 0., 0.) * grad(dof(v)), tf(v), liner, gauss);
  formulation.integral(im * rho0S / (Z * k0) * vector< std::complex< double > >(M, 0., 0.) * grad(dof(v)), vector< std::complex< double > >(M, 0., 0.) * grad(tf(v)), liner, gauss);
  msg::info << " - Applying Myer's boundary condition on Liner with impedance Re(Z)=" << std::real(Z) << ", Im(Z)=" << std::imag(Z) << msg::endl;

  formulation.pre();
  formulation.assemble();
  formulation.solve();

  save(+v, omega, "v");

  Line myline("lineDuct", 0., Rout, 0., L, Rout, 0., 500);
  save(v, myline, "vLine", "csv");
  return 0;
}
