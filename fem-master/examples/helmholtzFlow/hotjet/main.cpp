//#include <gmshfem/Formulation.h>
//#include <gmshfem/GmshFem.h>
#include "mesh.h"
#include "flowPmlFunctions.h"

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
  
  // 
  double pi = 3.14159265358979323846264338327950288;
  std::complex< double > im(0., 1.);
  std::string abcType = "pml"; // available choices: sommerfeld, pml
  gmshFem.userDefinedParameter(abcType, "abcType");

  // 
  double Sigma = 1.0018;
  gmshFem.userDefinedParameter(Sigma, "Sigma");
  int ElemOrder = 2;
  int Npml = 4;
  gmshFem.userDefinedParameter(Npml, "Npml");
  int RefSource = 18;

  std::pair<double, double> xlim;
  std::pair<double, double> ylim;
  double lc;
  meshJet(Sigma, RefSource, Npml, ElemOrder, xlim, ylim, lc);
  
  msg::info << "mesh size: " << lc << msg::endl;
  Domain omegaPhy = Domain("omega");
  Domain gammaABC = Domain("gamma");
  Domain omegaPml = Domain("omega_pml");
  Domain omega;
  if (abcType == "pml") {
    omega = omegaPhy | omegaPml;
    msg::info << "x-pml-limit-left: " << xlim.first << msg::endl;
    msg::info << "x-pml-limit-right: " << xlim.second << msg::endl;
    msg::info << "y-pml-limit-down: " << ylim.first << msg::endl;
    msg::info << "y-pml-limit-up: " << ylim.second << msg::endl;
  }
  else {
    omega = omegaPhy;
  }

  // flow acoustic problem
  double w = 200*pi;
  int FEMorder = 5;
  std::string gauss = "Gauss" + std::to_string(2 * FEMorder + 1);
  double SigmaS = 0.136*Sigma;
  
  // define source
  ScalarFunction< std::complex< double > > source = SigmaS*exp(- r2d<std::complex< double > >() / (2.*SigmaS*SigmaS) );
  // save(source, omega, "source");
  // define mean flow
  double Tj = 600; // Jet temperature [Kelvin]
  double gamma = 1.4;
  double R = 287.;
  double Mj = 0.9;
  double cj = sqrt(gamma*R*Tj);
  double uj = Mj*cj;
  ScalarFunction< std::complex< double > > u0y = uj*exp(- y<std::complex< double > >()*y<std::complex< double > >() / (2.*Sigma*Sigma) );
  
  double Tinf = 300; // Reference Temperature [Kelvin]
  double p0 = 103330;
  double rhoj = gamma*p0/(cj*cj);
  ScalarFunction< std::complex< double > > rho0_inv = Tinf/Tj - (Tinf/Tj - 1.) * (u0y/uj) + ((gamma-1)/2) * (Mj*Mj) * (u0y/uj) * (1-u0y/uj);
  ScalarFunction< std::complex< double > > rho0 = rhoj/rho0_inv;
  ScalarFunction< std::complex< double > > c0 = sqrt(gamma*p0/rho0);
   
  ScalarFunction< std::complex< double > > Mx = u0y/c0;
  ScalarFunction< std::complex< double > > My = 0.;
  ScalarFunction< std::complex< double > > k = w/c0;

  save(rho0, omega, "rho0");
  save(u0y, omega, "u0");
  save(u0y/c0, omega, "Mach");

  // set variational formulation
  Formulation< std::complex< double > > formulation("helmholtzflow");
  Field< std::complex< double >, Form::Form0 > u("u", omega | gammaABC, FunctionSpaceTypeForm0::HierarchicalH1, FEMorder);

  // convected Helmholz weak form
  formulation.integral((1./rho0) * vector< std::complex< double > >(1 - Mx * Mx, -Mx * My, 0.) * grad(dof(u)), vector< std::complex< double > >(1., 0., 0.) * grad(tf(u)), omegaPhy, gauss);
  formulation.integral((1./rho0) * vector< std::complex< double > >(-Mx * My, 1 - My * My, 0.) * grad(dof(u)), vector< std::complex< double > >(0., 1., 0.) * grad(tf(u)), omegaPhy, gauss);
  formulation.integral(- (1./rho0) * (k*k) * dof(u), tf(u), omegaPhy, gauss);
  formulation.integral((1./rho0) * dof(u), vector< std::complex< double > >(-im * k * Mx, -im * k * My, 0.) * grad(tf(u)), omegaPhy, gauss, term::ProductType::Scalar);
  formulation.integral((1./rho0) * vector< std::complex< double > >(im * k * Mx, im * k * My, 0.) * grad(dof(u)), tf(u), omegaPhy, gauss);

  formulation.integral(source, tf(u), omegaPhy, gauss);

  
  if (abcType == "sommerfeld") { // Zeroth-order ABC
    msg::info << "use sommerfeld abc " << msg::endl;
    formulation.integral((1./rho0) * im * k * dof(u), tf(u), gammaABC, gauss);
  }
  else if(abcType == "pml") { // PML
    msg::info << "use a pml " << msg::endl;
    ScalarFunction< std::complex< double > > beta = sqrt(1-Mx*Mx-My*My);
    TensorFunction< std::complex< double > > Linv;
    Linv = getInverseLorentzTensor(Mx, My, beta);

    ScalarFunction< std::complex< double > > detJ;
    TensorFunction< std::complex< double > > JPml_invT;
    JPml_invT = getInversePMLJacobian(k, xlim, ylim, Npml*lc, beta, detJ);

    VectorFunction< std::complex< double > > MM = vector< std::complex< double > >(Mx, My, 0.);
    VectorFunction< std::complex< double > > JPml_invT_M = JPml_invT * MM;
    TensorFunction< std::complex< double > > JPml_Linv = JPml_invT * Linv;

    // stabilized PML variational formulation
    formulation.integral((1./rho0) * detJ * JPml_Linv*grad(dof(u)) , JPml_Linv*grad(tf(u)) , omegaPml, gauss, term::ProductType::Scalar);
    formulation.integral((1./rho0) * k * k/(beta * beta) * detJ * JPml_invT_M * dof(u) , JPml_invT_M * tf(u), omegaPml, gauss, term::ProductType::Scalar);
    formulation.integral((1./rho0) * im * k / beta * detJ * JPml_Linv * grad(dof(u)), JPml_invT_M * tf(u), omegaPml, gauss, term::ProductType::Scalar);
    formulation.integral(- (1./rho0) * im * k /beta * detJ * JPml_invT_M * dof(u), JPml_Linv * grad(tf(u)), omegaPml, gauss, term::ProductType::Scalar);
    formulation.integral(- (1./rho0) * k * k/(beta * beta) * detJ * dof(u), tf(u), omegaPml, gauss);
  }
  else {
    msg::error << "ABC Type not available !" << msg::endl;
    exit(0);
  }


  // solve
  formulation.pre();
  formulation.assemble();
  formulation.solve();

  // save
  save(+u, omega, "u");


  return 0;
}