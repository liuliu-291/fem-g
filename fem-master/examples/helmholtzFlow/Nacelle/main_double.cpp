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
  double Lx = 3.; //3
  double Ly = 2.5; //2.5
  double L_duct = -0.5; // fan face position

  int Npml = 0, Npml_active = 0;
  bool withPml = false;
  bool withActivePml = false;
  gmshFem.userDefinedParameter(withPml, "Pml");
  gmshFem.userDefinedParameter(withActivePml, "aPml");
  if(withPml) {
    Npml = 4;
  }
  if(withActivePml) {
    Npml_active = 4;
  }
  // Build geometry and mesh using gmsh API
  // h_near, h_far, Lx, Ly, elemOrder, pml
  Nacelle(h_near, h_far, Lx, Ly, L_duct, ElemOrder, Npml, Npml_active);
  int FEMorder = 4; // FEM basis order
  std::string gauss = "Gauss" + std::to_string(2 * FEMorder + 1);
  double pi = 3.14159265358979323846264338327950288;
  std::complex< double > im(0., 1.);

  // ***************************************************************************************
  // Mean flow parameters
  double gamma0 = 1.4;
  double rho_inf = 1.23;
  double c_inf = 340;
  double v_inf = -85;
  double v_face = 187; // 74.8,153,187,181.02: approach, cutback, sideline, cruise
  double c_face = sqrt(c_inf * c_inf * (1 + (gamma0 - 1) / 2. * (v_inf * v_inf - v_face * v_face) / (c_inf * c_inf)));
  double M_face = -v_face / c_face;
  msg::info << "Fan face Mach number = " << M_face << msg::endl;

  ScalarFunction< double > c0, rho0, Mx, My, v;
  // Define FEM domains from CAD
  Domain fanFace = Domain(1, 5000);
  Domain SymAxis = Domain(1, 17000);
  Domain Liner = Domain(1, 2000);

  Domain omegaPhy, omegaPml, omegaPmlA, GammaExt, GammaInfR, GammaInfL, fanFace_flow;
  omegaPhy = Domain(2, 2000);
  Domain Omega;
  if(withPml) {
    GammaInfR = Domain(1, 24000);
    GammaInfL = Domain(1, 22000);
    fanFace_flow = Domain(1, 5000);
    omegaPml = Domain(2, 3000);
    Omega = omegaPhy | omegaPml;
  }
  else {
    fanFace_flow = Domain(1, 5000);
    GammaInfL = Domain(1, 20000);
    GammaInfR = Domain(1, 18000);
    GammaExt = GammaInfR | Domain(1, 19000) | GammaInfL;
    Omega = omegaPhy;
  }
  if(withActivePml) {
    omegaPmlA = Domain(2, 4000);
    fanFace_flow = Domain(1, 25000);
    if(withPml) {
      omegaPml = Domain(2, 3000) | omegaPmlA;
      Omega = omegaPhy | omegaPml | omegaPmlA;
    }
    else {
      Omega = omegaPhy | omegaPmlA;
    }
  }

  // compute mean flow...
  Field< double, Form::Form0 > phi("phi", Omega | fanFace_flow | GammaInfR | GammaInfL, FunctionSpaceTypeForm0::HierarchicalH1, FEMorder);
  Formulation< double > MeanFlow("poisson");

  phi.addConstraint(GammaInfL, 0.); // Dirichlet BC
  // Weak formulation
  MeanFlow.integral(-y<double>() * rho_inf * grad(dof(phi)), grad(tf(phi)), Omega , gauss);
  MeanFlow.integral(y<double>() * rho_inf * v_face, tf(phi), fanFace_flow, gauss); // Neumann BC fan face
  MeanFlow.integral(y<double>() * rho_inf * v_inf, tf(phi), GammaInfR, gauss); // Neumann BC
  MeanFlow.pre();
  MeanFlow.assemble();
  MeanFlow.solve();
  v = norm(grad(phi), 2);
  c0 = sqrt(c_inf * c_inf * (1 + (gamma0 - 1) / 2. * (v_inf * v_inf - v * v) / (c_inf * c_inf))); // speed of sound
  Mx = xComp(grad(phi)) / c0;
  My = yComp(grad(phi)) / c0; // Mach number components
  rho0 = rho_inf * pow((1 + (gamma0 - 1) / 2. * (v_inf * v_inf - v * v) / (c_inf * c_inf)), 1 / (gamma0 - 1)); // density

  bool save_flow = false;
  if(save_flow) {
    save(c0, Omega, "c0_ref");
    save(rho0, Omega, "rho0_ref");
    save(Mx, Omega, "Mx_ref");
    save(My, Omega, "My_ref");
  }

  // -------------------------------------------------------------------------
  // acoustic problem
  double w = 2 * pi * 1300 * 2; // frequency - 1BPF
  gmshFem.userDefinedParameter(w, "omega");

  // Get Input Modal function and wavenumber
  std::complex< double > kx; // propagative wavenumber
  double ky; // Laplace-Beltrami eigenvalue
  int m = 24; // azimuthal mode number
  int n = 0; // radial mode number
  ScalarFunction< std::complex< double > > psi, dy_psi;

  getInputMode(m, n, w / c_face, M_face, kx, ky, psi, dy_psi);
  VectorFunction< std::complex< double > > grad_psi = vector< std::complex< double > >(-im * kx * psi, dy_psi, 0.); // for active PML

  msg::info << "propagative wavenumber kx = " << kx << ", transverse wavenumber ky = " << ky << msg::endl;
  float pointsByWl_face = (2 * pi * c_face * FEMorder * (1 - abs(M_face))) / (h_near * w);
  float pointsByWl_far = (2 * pi * c_inf * FEMorder * (1 - abs(v_inf / c_inf))) / (h_far * w);
  msg::info << " - FEM basis order = " << FEMorder << "" << msg::endl;
  msg::info << " - Approx. dofs by wavelength at fan face = " << pointsByWl_face << msg::endl;
  msg::info << " - Approx. dofs by wavelength at outer boundary = " << pointsByWl_far << msg::endl;

  // Parameters for Inverse Lorentz transformation
  VectorFunction< double > MM = vector< double >(Mx, My, 0.);
  ScalarFunction< double > beta = sqrt(1 - pow(norm(MM, 2), 2)); // Jacobian
  ScalarFunction< double > Alphax = 1 + pow(Mx, 2) / (beta * (1 + beta));
  ScalarFunction< double > Alphay = 1 + pow(My, 2) / (beta * (1 + beta));
  ScalarFunction< double > K = Mx * My / (beta * (1 + beta));
  // Inverse Lorentz transformation tensor
  TensorFunction< double > Linv = tensor< double >(beta * Alphay, -beta * K, 0., -beta * K, beta * Alphax, 0., 0., 0., 0.);

  // Mach-velocity normal and tangent vectors
  ScalarFunction< double > Mn = MM * normal< double >();
  ScalarFunction< double > Mt = MM * tangent< double >();

  // Cartesian PML definition
  ScalarFunction< double > Sigma0 = beta; // use heterogeneous parameter
  double Lxpml = Lx + Npml * h_far;
  double L_ductpml = L_duct - Npml * h_far;
  double Lypml = Ly + Npml * h_far;

  ScalarFunction< double > SigmaX = heaviside(x< double >() - Lx) * Sigma0 / (Lxpml - abs(x< double >())) + heaviside(-x< double >() + L_duct) * Sigma0 / (abs(L_ductpml) - abs(x< double >()));
  ScalarFunction< double > SigmaY = heaviside(abs(y< double >()) - Ly) * Sigma0 / (Lypml - abs(y< double >()));
  ScalarFunction< std::complex< double > > gammaX = 1 - (complex(c0 * SigmaX) * im / w);
  ScalarFunction< std::complex< double > > gammaY = 1 - (complex(c0 * SigmaY) * im / w);

  ScalarFunction< std::complex< double > > det_J = gammaX * gammaY;
  TensorFunction< std::complex< double > > J_PML_inv_T = tensor< std::complex< double > >(1. / gammaX, 0., 0., 0., 1. / gammaY, 0., 0., 0., 0.);
  // define Lorentz-PML coupling parameters for the weak form
  VectorFunction< std::complex< double > > J_PML_inv_T_M = J_PML_inv_T * complex(MM);
  TensorFunction< std::complex< double > > J_PML_Linv = J_PML_inv_T * complex(Linv);

  Field< std::complex< double >, Form::Form0 > u("u", Omega | GammaExt | SymAxis | fanFace | Liner, FunctionSpaceTypeForm0::HierarchicalH1, FEMorder);
  Formulation< std::complex< double > > formulation("HelmholtzFlow");

  if(m != 0) {
    u.addConstraint(SymAxis, 0.);
  }

  ScalarFunction< double > rho0S = rho0 * y< double >(); // rho0*jacobian in 2D cylindrical coordinates
  ScalarFunction< std::complex< double > > rho0Sc = complex(rho0S);
  ScalarFunction< std::complex< double > > c0c = complex(c0);
  ScalarFunction< std::complex< double > > betac = complex(beta);

  // convected Helmholz weak form
  formulation.integral(rho0Sc * complex(vector< double >(1-Mx*Mx, -Mx*My, 0.)) * grad(dof(u)), vector< std::complex< double > >(1., 0., 0.) * grad(tf(u)), omegaPhy, gauss);
  formulation.integral(rho0Sc * complex(vector< double >(-Mx*My, 1-My*My, 0.)) * grad(dof(u)), vector< std::complex< double > >(0., 1., 0.) * grad(tf(u)), omegaPhy, gauss);
  formulation.integral(+ im * rho0Sc * (w/c0c) * complex(vector< double >(Mx, My, 0.)) * grad(dof(u)), tf(u), omegaPhy, gauss);
  formulation.integral(- im * rho0Sc * (w/c0c) * dof(u), complex(vector< double >(Mx, My, 0.)) * grad(tf(u)), omegaPhy, gauss);
  formulation.integral(- rho0Sc *(w*w) / (c0c*c0c) * dof(u), tf(u), omegaPhy, gauss);
  formulation.integral(+ rho0Sc *(m*m) / complex( y< double >()*y< double >() ) * dof(u), tf(u), omegaPhy, gauss);

  if (withPml || withActivePml) {
    // stable pml formulation
    msg::info << " - Use stable pml formulation "  << msg::endl;
    formulation.integral(+ rho0Sc * det_J * J_PML_Linv*grad(dof(u)) , J_PML_Linv*grad(tf(u)) , omegaPml, gauss, term::ProductType::Scalar);
    formulation.integral(+ rho0Sc * (w*w)/(c0c*c0c*betac*betac) * det_J * J_PML_inv_T_M * dof(u) , J_PML_inv_T_M * tf(u), omegaPml, gauss, term::ProductType::Scalar);
    formulation.integral(+ im * rho0Sc * w /(c0c*betac) * det_J * J_PML_Linv * grad(dof(u)), J_PML_inv_T_M * tf(u), omegaPml, gauss, term::ProductType::Scalar);
    formulation.integral(- im * rho0Sc * w /(c0c*betac) * det_J * J_PML_inv_T_M * dof(u), J_PML_Linv * grad(tf(u)), omegaPml, gauss, term::ProductType::Scalar);
    formulation.integral(- rho0Sc * (w*w)/(c0c*c0c*betac*betac) * det_J * dof(u), tf(u), omegaPml, gauss);
    formulation.integral(+ rho0Sc *(m*m) / complex(y< double >()*y< double >()) * det_J * dof(u), tf(u), omegaPml, gauss);
  }
  else {
    // Sommerfeld abc
    msg::info << " - Use Sommerfeld ABC "  << msg::endl;
    formulation.integral( im * rho0Sc * (w/c0c) * dof(u), tf(u), GammaExt, gauss);
  }

  if(withActivePml) {
    // Active PML - incident field
    msg::info << " - Use active pml "  << msg::endl;
    formulation.integral(- rho0Sc * det_J * J_PML_Linv*grad_psi , J_PML_Linv*grad(tf(u)) , omegaPmlA, gauss, term::ProductType::Scalar);
    formulation.integral(- rho0Sc * (w*w)/(c0c*c0c*betac*betac) * det_J * J_PML_inv_T_M * psi , J_PML_inv_T_M * tf(u), omegaPmlA, gauss, term::ProductType::Scalar);
    formulation.integral(- im * rho0Sc * w/(c0c*betac) * det_J * J_PML_Linv * grad_psi, J_PML_inv_T_M * tf(u), omegaPmlA, gauss, term::ProductType::Scalar);
    formulation.integral(+ im * rho0Sc * w/(c0c*betac) * det_J * J_PML_inv_T_M * psi, J_PML_Linv * grad(tf(u)), omegaPmlA, gauss, term::ProductType::Scalar);
    formulation.integral(+ rho0Sc * (w*w)/(c0c*c0c*betac*betac) * det_J * psi, tf(u), omegaPmlA, gauss);
    formulation.integral(- rho0Sc *(m*m) / complex( y< double >()*y< double >() ) * det_J * psi, tf(u), omegaPmlA, gauss);
    // Flux term
    formulation.integral(- im * kx * rho0Sc * (betac*betac) * psi, tf(u), fanFace, gauss);
    formulation.integral(- im * rho0Sc * (w/c0c) * complex(Mx) * psi, tf(u), fanFace, gauss);
  }
  else {
    formulation.integral(- im * kx * rho0Sc * (betac*betac) * psi, tf(u), fanFace, gauss);
    formulation.integral(- im* rho0Sc * (w/c0c) * complex(Mx) * dof(u), tf(u), fanFace, gauss);
  }
  // Myer's BC on Liner
  std::complex< double > Z(0.4, -0.2); // normalized impedance (by rho0*c0)
  formulation.integral( im*rho0Sc*w/c0c / Z * dof(u), tf(u), Liner, gauss);
  formulation.integral( -rho0Sc / Z * dof(u), complex(Mt * tangent< double >()) * grad(tf(u)), Liner, gauss);
  formulation.integral( +rho0Sc / Z * complex(Mt * tangent< double >()) * grad(dof(u)), tf(u), Liner, gauss);
  formulation.integral( im*rho0Sc / Z * (c0c/w) * complex(Mt * tangent< double >()) * grad(dof(u)),  complex(Mt * tangent< double >()) * grad(tf(u)), Liner, gauss);
  msg::info << " - Applying Myer's boundary condition on Liner with impedance Re(Z)=" << std::real(Z) <<", Im(Z)=" << std::imag(Z) << msg::endl;

  formulation.pre();
  formulation.assemble();
  formulation.solve();

  save(+u, Omega, "u");
  return 0;
}
