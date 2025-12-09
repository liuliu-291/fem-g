#include <gmsh.h>
#include <gmshfem/Formulation.h>
#include <gmshfem/GmshFem.h>
#include <gmshfem/Message.h>
#include "myFunctions.h"
#include <fstream>
#include <iostream>

using namespace gmshfem;
using namespace gmshfem::common;
using namespace gmshfem::problem;
using namespace gmshfem::domain;
using namespace gmshfem::field;
using namespace gmshfem::function;
using namespace gmshfem::analytics;
using namespace gmshfem::post;
using namespace gmshfem::equation;


int main(int argc, char **argv)
{
  GmshFem gmshFem(argc, argv);

  int ElemOrder = 4; // mesh order
  double h_near = 0.015;
  double h_far = 0.03;
  double Lx = 3.;
  double Ly = 2.5;
  double L_duct = -0.5; // fan face position
  int Npml = 4;
  // Build geometry and mesh using gmsh API
  // h_near, h_far, Lx, Ly, elemOrder, pml
  Nacelle(h_near, h_far, Lx, Ly, L_duct, ElemOrder, Npml);

  int FEMorder = 4; // FEM basis order
  std::string gauss = "Gauss" + std::to_string(2 * FEMorder + 1);
  double pi = 3.14159265358979323846264338327950288;
  std::complex< double > im(0., 1.);

  // ***************************************************************************************
  // 1 - Compute mean flow properties, Poisson problem
  // Mean flow parameters
  double gamma0 = 1.4;
  double rho_inf = 1.23;
  double c_inf = 340;
  double v_inf = -85;
  double v_face = 187; // 74.8,153,187,181.02: approach, cutback, sideline, cruise
  double c_face = sqrt(c_inf * c_inf * (1 + (gamma0 - 1) / 2. * (v_inf * v_inf - v_face * v_face) / (c_inf * c_inf)));
  double Mf = -v_face / c_face;
  msg::info << "Fan face Mach number = " << Mf << msg::endl;

  // Define FEM domains from CAD
  Domain fanFace = Domain(1, 5000);
  Domain omegaPhy = Domain(2, 2000);
  Domain omegaPml = Domain(2, 3000);
  Domain GammaInfR = Domain(1, 24000);
  Domain GammaInfL = Domain(1, 22000);
  Domain GammaExt = Domain(1, 18000) | Domain(1, 19000) | Domain(1, 20000);
  Domain SymAxis = Domain(1, 17000);

  Field< double, Form::Form0 > phi("phi", omegaPhy | omegaPml | fanFace | GammaInfR | GammaInfL, FunctionSpaceTypeForm0::HierarchicalH1, FEMorder);
  Formulation< double > MeanFlow("poisson");

  // Apply Dirichlet BC
  phi.addConstraint(GammaInfL, 0.);
  // Weak formulation
  MeanFlow.integral(-rho_inf * grad(dof(phi)), grad(tf(phi)), omegaPhy | omegaPml, gauss);
  MeanFlow.integral(v_face, tf(phi), fanFace, gauss);
  // Apply Neumann BC
  MeanFlow.integral(v_inf, tf(phi), GammaInfR, gauss); // Apply Neumann BCl

  MeanFlow.pre();
  MeanFlow.assemble();
  MeanFlow.solve();
  // Compute mean flow properties
  ScalarFunction< double > v, c2, c0d, rho0d, Mx, My;
  ScalarFunction< std::complex< double > > c0;
  v = norm(grad(phi), 2);

  // density
  rho0d = rho_inf * pow((1 + (gamma0 - 1) / 2. * (v_inf * v_inf - v * v) / (c_inf * c_inf)), 1 / (gamma0 - 1));
  save(rho0d, omegaPhy | omegaPml, "rho0");
  // speed of sound
  c0d = sqrt(c_inf * c_inf * (1 + (gamma0 - 1) / 2. * (v_inf * v_inf - v * v) / (c_inf * c_inf)));
  save(c0d, omegaPhy | omegaPml, "c0");
  c0 = complex(c0d);
  // Mach-velocity field
  Mx = xComp(grad(phi)) / c0d;
  My = yComp(grad(phi)) / c0d;
  VectorFunction< double > MMd = vector< double >(Mx, My, 0.);
  VectorFunction< std::complex< double > > MM = complex(MMd);

  save(Mx, omegaPhy | omegaPml, "Mx");
  save(My, omegaPhy | omegaPml, "My");
  save(norm(MM, 2), omegaPhy | omegaPml, "Mach");
  save(phi, omegaPhy | omegaPml, "phi");
  // ***************************************************************************************

  // -------------------------------------------------------------------------
  // acoustic problem
  double w = 2 * pi * 1300 * 2; // frequency
  gmshFem.userDefinedParameter(w, "omega");

  // Get Input Modal function and wavenumber
  std::complex< double > kx; // propagative wavenumber
  double ky; // Laplace-Beltrami eigenvalue
  int m = 24;
  int n = 0;
  ScalarFunction< std::complex< double > > psi, dy_psi;
  getInputMode(m, n, w / c_face, Mf, kx, ky, psi, dy_psi);
  msg::info << "propagative wavenumber kx = " << kx << ", transverse wavenumber ky = " << ky << msg::endl;

  float pointsByWl = (2 * pi * c_inf * FEMorder * (1 - abs(Mf))) / (h_far * w);
  msg::info << " - FEM basis order = " << FEMorder << "" << msg::endl;
  msg::info << " - Approx. dofs by wavelength at fan face = " << pointsByWl << msg::endl;

  // -------------------------------------------------------------------------
  // Parameters for Inverse Lorentz transformation
  ScalarFunction< double > betad = sqrt(1 - pow(norm(MMd, 2), 2)); // Jacobian
  ScalarFunction< std::complex< double > > beta = complex(betad);
  ScalarFunction< double > Alphax = 1 + pow(Mx, 2) / (betad * (1 + betad));
  ScalarFunction< double > Alphay = 1 + pow(My, 2) / (betad * (1 + betad));
  ScalarFunction< double > K = Mx * My / (betad * (1 + betad));
  // Inverse Lorentz transformation tensor
  TensorFunction< double > Linv = tensor< double >(betad * Alphay, -betad * K, 0., -betad * K, betad * Alphax, 0., 0., 0., 0.);

  // Mach-velocity normal and tangent vectors
  ScalarFunction< std::complex< double > > Mn = MM * normal< std::complex< double > >();
  ScalarFunction< std::complex< double > > Mt = MM * tangent< std::complex< double > >();

  // Cartesian PML definition
  ScalarFunction< double > Sigma0 = betad; // use heterogeneous parameter
  double L_ductpml = L_duct - Npml * h_far;
  double Lxpml = Lx + Npml * h_far;
  double Lypml = Ly + Npml * h_far;
  ScalarFunction< double > SigmaX = heaviside(x< double >() - Lx) * Sigma0 / (Lxpml - abs(x< double >())) + heaviside(-x< double >() + L_duct) * Sigma0 / (abs(L_ductpml) - abs(x< double >()));
  ScalarFunction< double > SigmaY = heaviside(abs(y< double >()) - Ly) * Sigma0 / (Lypml - abs(y< double >()));
  ScalarFunction< std::complex< double > > gammaX = 1 - (c0 * im / w) * complex(SigmaX);
  ScalarFunction< std::complex< double > > gammaY = 1 - (c0 * im / w) * complex(SigmaY);

  ScalarPiecewiseFunction< std::complex< double > > det_J;
  det_J.addFunction(1., omegaPhy);
  det_J.addFunction(gammaX * gammaY, omegaPml);

  // define Lorentz-PML coupling parameters for the weak form
  TensorFunction< std::complex< double > > J_PML_inv_T = tensor< std::complex< double > >(1. / gammaX, 0., 0., 0., 1. / gammaY, 0., 0., 0., 0.);
  VectorPiecewiseFunction< std::complex< double > > J_PML_inv_T_M;
  J_PML_inv_T_M.addFunction(MM, omegaPhy);
  J_PML_inv_T_M.addFunction(J_PML_inv_T * MM, omegaPml);

  TensorPiecewiseFunction< std::complex< double > > J_PML_Linv;
  J_PML_Linv.addFunction(complex(Linv), omegaPhy);
  J_PML_Linv.addFunction(J_PML_inv_T * complex(Linv), omegaPml);


  Domain Omega = omegaPhy | omegaPml;
  Field< std::complex< double >, Form::Form0 > u("u", Omega | GammaExt | SymAxis | fanFace, FunctionSpaceTypeForm0::HierarchicalH1, FEMorder);
  Formulation< std::complex< double > > formulation("HelmholtzFlow");

  if(m != 0) {
    u.addConstraint(SymAxis, 0.);
  }

  ScalarFunction< double > rho0Sd = rho0d * y< double >();
  ScalarFunction< std::complex< double > > rho0S = complex(rho0Sd);
  // use Piecewise functions to keep a single weak form in omega_phy | omega_pml
  formulation.integral(rho0S * det_J * J_PML_Linv*grad(dof(u)) , J_PML_Linv*grad(tf(u)) , Omega, gauss, term::ProductType::Scalar);
  formulation.integral(+ rho0S * (w*w)/(c0*c0*beta*beta) * det_J * J_PML_inv_T_M * dof(u) , J_PML_inv_T_M * tf(u), Omega, gauss, term::ProductType::Scalar);
  formulation.integral(+ rho0S * im*w/(c0*beta) * det_J * J_PML_Linv * grad(dof(u)), J_PML_inv_T_M * tf(u), Omega, gauss, term::ProductType::Scalar);
  formulation.integral(- rho0S * im*w/(c0*beta) * det_J * J_PML_inv_T_M * dof(u), J_PML_Linv * grad(tf(u)), Omega, gauss, term::ProductType::Scalar);
  formulation.integral(- rho0S * (w*w)/(c0*c0*beta*beta) * det_J * dof(u), tf(u), Omega, gauss);
  formulation.integral(+ rho0S *(m*m) / ( y< std::complex< double > >()*y< std::complex< double > >() ) * det_J * dof(u), tf(u), Omega, gauss);
  // Input BC fan face
  formulation.integral(-rho0S * im * kx * (betaf * betaf) * psi, tf(u), fanFace, gauss);
  formulation.integral(-rho0S * im * (w / c0) * Mf * dof(u), tf(u), fanFace, gauss);

  formulation.pre();
  formulation.assemble();
  formulation.solve();

  save(+u, Omega, "u");

  return 0;
}
