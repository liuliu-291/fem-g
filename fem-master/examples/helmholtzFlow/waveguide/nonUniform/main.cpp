#include <gmsh.h>
#include <gmshfem/AnalyticalFunction.h>
#include <gmshfem/Formulation.h>
#include <gmshfem/GmshFem.h>

using namespace gmshfem;
using namespace gmshfem::common;
using namespace gmshfem::problem;
using namespace gmshfem::domain;
using namespace gmshfem::field;
using namespace gmshfem::function;
using namespace gmshfem::analytics;
using namespace gmshfem::post;
using namespace gmshfem::equation;

std::complex< double > GetInputKx(double, double, double, double, double, double);

void meshDuct(const double L, const double Lext, const double H, const double h, const int Npml, const int MeshOrder)
{
  msg::info << " - Duct of length L = " << L << msg::endl;
  msg::info << " - Duct of heigth H = " << H << msg::endl;
  gmsh::model::add("geometry");

  // duct definition
  gmsh::model::geo::addPoint(L, H, 0., h, 1);
  gmsh::model::geo::addPoint(0., H, 0., h, 2);
  gmsh::model::geo::addPoint(0., 0., 0., h, 3);
  gmsh::model::geo::addPoint(L, 0., 0., h, 4);

  gmsh::model::geo::addLine(1, 2, 1);
  gmsh::model::geo::addLine(2, 3, 2);
  gmsh::model::geo::addLine(3, 4, 3);
  gmsh::model::geo::addLine(4, 1, 4);
  gmsh::model::geo::addCurveLoop({1, 2, 3, 4}, 1);
  gmsh::model::geo::addPlaneSurface({1}, 1);
  gmsh::model::geo::mesh::setTransfiniteSurface(1);
  gmsh::model::geo::mesh::setRecombine(1, 1);

  // extended duct
  gmsh::model::geo::addPoint(Lext, H, 0., h, 5);
  gmsh::model::geo::addPoint(Lext, 0, 0., h, 6);
  gmsh::model::geo::addLine(4, 6, 5);
  gmsh::model::geo::addLine(6, 5, 6);
  gmsh::model::geo::addLine(5, 1, 7);
  gmsh::model::geo::addCurveLoop({5, 6, 7, -4}, 2);
  gmsh::model::geo::addPlaneSurface({2}, 2);
  gmsh::model::geo::mesh::setTransfiniteSurface(2);
  gmsh::model::geo::mesh::setRecombine(2, 2);

  // pml
  gmsh::model::geo::addPoint(Lext + Npml * h, H, 0., h, 7);
  gmsh::model::geo::addPoint(Lext + Npml * h, 0, 0., h, 8);
  gmsh::model::geo::addLine(6, 8, 8);
  gmsh::model::geo::addLine(8, 7, 9);
  gmsh::model::geo::addLine(7, 5, 10);

  gmsh::model::geo::addCurveLoop({8, 9, 10, -6}, 3);
  gmsh::model::geo::addPlaneSurface({3}, 3);
  gmsh::model::geo::mesh::setTransfiniteSurface(3);
  gmsh::model::geo::mesh::setRecombine(3, 3);

  // physical groups - surfaces
  gmsh::model::addPhysicalGroup(2, {1}, 1);
  gmsh::model::setPhysicalName(2, 1, "omega");
  gmsh::model::addPhysicalGroup(2, {2}, 2);
  gmsh::model::setPhysicalName(2, 2, "omegaExt");
  gmsh::model::addPhysicalGroup(2, {3}, 3);
  gmsh::model::setPhysicalName(2, 3, "omegaPml");

  // boundaries - physical domain
  gmsh::model::addPhysicalGroup(1, {1}, 1);
  gmsh::model::setPhysicalName(1, 1, "gammaTop");
  gmsh::model::addPhysicalGroup(1, {2}, 2);
  gmsh::model::setPhysicalName(1, 2, "gammaLeft");
  gmsh::model::addPhysicalGroup(1, {3}, 3);
  gmsh::model::setPhysicalName(1, 3, "gammaBottom");
  gmsh::model::addPhysicalGroup(1, {4}, 4);
  gmsh::model::setPhysicalName(1, 4, "gammaRight");
  // boundaries - extended domain
  gmsh::model::addPhysicalGroup(1, {5}, 5);
  gmsh::model::setPhysicalName(1, 5, "gammaBottomExt");
  gmsh::model::addPhysicalGroup(1, {6}, 6);
  gmsh::model::setPhysicalName(1, 6, "gammaRightExt");
  gmsh::model::addPhysicalGroup(1, {7}, 7);
  gmsh::model::setPhysicalName(1, 7, "gammaTopExt");
  // boundaries - pml domain
  gmsh::model::addPhysicalGroup(1, {8}, 8);
  gmsh::model::setPhysicalName(1, 8, "gammaBottomPml");
  gmsh::model::addPhysicalGroup(1, {9}, 9);
  gmsh::model::setPhysicalName(1, 9, "gammaRightPml");
  gmsh::model::addPhysicalGroup(1, {10}, 10);
  gmsh::model::setPhysicalName(1, 10, "gammaTopPml");

  // generate mesh
  gmsh::model::geo::synchronize();
  gmsh::model::mesh::generate(2);
  gmsh::model::mesh::setOrder(MeshOrder);
  gmsh::option::setNumber("Mesh.ElementOrder", MeshOrder);
  gmsh::model::mesh::setOrder(MeshOrder);
  gmsh::write("m.msh");
}

int main(int argc, char **argv)
{
  //*****
  // Problem declaration
  //*****
  GmshFem gmshFem(argc, argv);
  const double pi = 3.14159265358979323846264338327950288;
  const std::complex< double > im(0., 1.);

  // numerical parameters
  int FEMorder = 8;
  gmshFem.userDefinedParameter(FEMorder, "FEMorder");
  std::string gauss = "Gauss" + std::to_string(2 * FEMorder + 1);

  // meshing
  double L = 2.; // duct length
  gmshFem.userDefinedParameter(L, "L");
  double Lext = 4.; // reference duct length
  double H = 0.5; // duct height
  gmshFem.userDefinedParameter(H, "H");
  double h = 1. / 40; // meshsize
  gmshFem.userDefinedParameter(h, "h");
  int Npml = 10;
  meshDuct(L, Lext, H, h, Npml, 1);

  // duct problem specifications
  double kinf = 40;
  gmshFem.userDefinedParameter(kinf, "kinf");
  int mode = 4; // transverse mode
  gmshFem.userDefinedParameter(mode, "mode");
  double ky = mode * pi / H;

  // Choose ABC type
  // Available conditions: ABC-0, ABC-2, Pade, Pade-S2
  std::string abcName = "Pade";
  gmshFem.userDefinedParameter(abcName, "abcName");

  // Mean flow profile - define vx, Mx, c0, rho0
  const double cinf = 340;
  double k = cinf * kinf;
  const double rhoinf = 1.2;
  const double deltav = 0.3 * cinf;
  const double vinf = -0.1 * cinf;
  const double v0 = -0.4 * cinf;
  const double s = 5;
  const double xc = 1;
  const double gamm = 1.4;
  ScalarFunction< std::complex< double > > vx = v0 + deltav * tanh(s * (x< std::complex< double > >() - xc));
  ScalarFunction< std::complex< double > > c0 = cinf * sqrt(1 + (vinf * vinf - vx * vx) * (gamm - 1) / (cinf * cinf * 2.));
  ScalarFunction< std::complex< double > > rho0 = rhoinf * abs(vinf / vx);
  ScalarFunction< std::complex< double > > Mx = vx / c0;
  ScalarFunction< std::complex< double > > betax = sqrt(1 - Mx * Mx);
  // derivatives
  ScalarFunction< std::complex< double > > dvx = deltav * s * (1 - pow(tanh(s * (x< std::complex< double > >() - xc)), 2));
  ScalarFunction< std::complex< double > > dc0 = vx * dvx * (1 - gamm) / (2. * c0);
  // wavenumber
  ScalarFunction< std::complex< double > > k0 = k / c0;

  std::complex< double > kx = GetInputKx(k, vinf, v0, gamm, cinf, ky);

  float MinPtsByWl = 2 * pi * FEMorder * (1 - abs(v0 / cinf)) / kinf / h; // mean flow is max at x=0
  msg::info << "********************************" << msg::endl;
  msg::info << " - convected Helmholtz duct problem " << msg::endl;
  msg::info << " - reference wave number = " << kinf << msg::endl;
  msg::info << " - FEM basis order = " << FEMorder << "" << msg::endl;
  msg::info << " - Min dofs by wavelength = " << MinPtsByWl << "" << msg::endl;
  msg::info << "********************************" << msg::endl;

  // declare formulation and FEM domains
  std::vector< FieldInterface< std::complex< double > > * > fieldBucket;
  Formulation< std::complex< double > > formulation("helmholtz-duct-flow");
  Formulation< std::complex< double > > formulationRef("helmholtz-duct-flow-refPML");

  // physical domains
  Domain omega(2, 1);
  Domain omegaExt(2, 2);
  Domain omegaPml(2, 3);
  // boundaries
  Domain gammaTop(1, 1);
  Domain gammaLeft(1, 2);
  Domain gammaBottom(1, 3);
  Domain gammaRight(1, 4);
  Domain gammaBottomExt(1, 5);
  Domain gammaRightExt(1, 6);
  Domain gammaTopExt(1, 7);
  Domain gammaBottomPml(1, 8);
  Domain gammaRightPml(1, 9);
  Domain gammaTopPml(1, 10);

  Field< std::complex< double >, Form::Form0 > v("v", omega | gammaLeft | gammaRight | gammaTop | gammaBottom, FunctionSpaceTypeForm0::HierarchicalH1, FEMorder);

  // weak formulation
  formulation.integral(rho0 * vector< std::complex< double > >(betax * betax, 0., 0.) * grad(dof(v)), vector< std::complex< double > >(1., 0., 0.) * grad(tf(v)), omega, gauss);
  formulation.integral(rho0 * vector< std::complex< double > >(0., 1., 0.) * grad(dof(v)), vector< std::complex< double > >(0., 1., 0.) * grad(tf(v)), omega, gauss);
  formulation.integral(rho0 * vector< std::complex< double > >(im * k0 * Mx, 0., 0.) * grad(dof(v)), tf(v), omega, gauss);
  formulation.integral(-rho0 * im * k0 * Mx * dof(v), vector< std::complex< double > >(1., 0., 0.) * grad(tf(v)), omega, gauss);
  formulation.integral(-rho0 * k0 * k0 * dof(v), tf(v), omega, gauss);

  // input BC
  formulation.integral(-betax * betax * im * kx * cos(ky * y< std::complex< double > >()), tf(v), gammaLeft, gauss);
  formulation.integral(-im * k0 * Mx * dof(v), tf(v), gammaLeft, gauss);

  // output BC
  if(abcName == "ABC-0") {
    msg::info << "Use Sommerfeld ABC." << msg::endl;
    formulation.integral(im * k0 * (1 - Mx) * dof(v), tf(v), gammaRight, gauss);
    formulation.integral(im * k0 * Mx * dof(v), tf(v), gammaRight, gauss);
  }
  if(abcName == "ABC-2") {
    msg::info << "Use second order ABC." << msg::endl;
    formulation.integral(im * k0 * (1 - Mx) * dof(v), tf(v), gammaRight, gauss);
    formulation.integral(im * k0 * Mx * dof(v), tf(v), gammaRight, gauss);
    formulation.integral(im * betax * betax / (2. * k) * grad(dof(v)), grad(tf(v)), gammaRight, gauss);
  }
  else if(abcName.std::string::find("Pade") != std::string::npos) {
    int padeOrder = 4;
    gmshFem.userDefinedParameter(padeOrder, "padeOrder");
    double angle = -pi / 4.;
    gmshFem.userDefinedParameter(angle, "angle");
    msg::info << "Use Pade localization for the principal symbol of order " << padeOrder << " with angle " << angle << " rad" << msg::endl;

    const double Np = 2. * padeOrder + 1.;
    const std::complex< double > exp1 = std::complex< double >(std::cos(angle), std::sin(angle));
    const std::complex< double > exp2 = std::complex< double >(std::cos(angle / 2.), std::sin(angle / 2.));
    const std::complex< double > coef = 2. / Np;
    std::vector< std::complex< double > > c(padeOrder, 0.);
    for(int i = 0; i < padeOrder; ++i) {
      c[i] = std::tan((i + 1) * pi / Np);
      c[i] *= c[i];
    }

    // define the auxiliary fields
    std::vector< Field< std::complex< double >, Form::Form0 > * > phi;
    for(int i = 0; i < padeOrder; ++i) {
      phi.push_back(new Field< std::complex< double >, Form::Form0 >("phi_" + std::to_string(i), gammaRight, FunctionSpaceTypeForm0::HierarchicalH1, FEMorder));
      fieldBucket.push_back(phi.back());
    }

    // write the augmented weak form
    formulation.integral(rho0 * im * k0 * exp2 * dof(v), tf(v), gammaRight, gauss);
    for(int i = 0; i < padeOrder; ++i) {
      // boundary integral terms relating the auxiliary fields
      formulation.integral(rho0 * im * k0 * exp2 * coef * c[i] * dof(*phi[i]), tf(v), gammaRight, gauss);
      formulation.integral(rho0 * im * k0 * exp2 * coef * c[i] * dof(v), tf(v), gammaRight, gauss);

      // coupling of the auxiliary equations
      formulation.integral(grad(dof(*phi[i])), grad(tf(*phi[i])), gammaRight, gauss);
      formulation.integral(-(k0 * k0) / (betax * betax) * (exp1 * c[i] + 1.) * dof(*phi[i]), tf(*phi[i]), gammaRight, gauss);
      formulation.integral(-(k0 * k0) / (betax * betax) * exp1 * (c[i] + 1.) * dof(v), tf(*phi[i]), gammaRight, gauss);
    }
    if(abcName == "Pade-S2") {
      msg::info << "Add second symbol to Pade ABC " << msg::endl;
      // define a new auxiliary field
      std::vector< Field< std::complex< double >, Form::Form0 > * > psi;
      psi.push_back(new Field< std::complex< double >, Form::Form0 >("psi", gammaRight, FunctionSpaceTypeForm0::HierarchicalH1, FEMorder));
      fieldBucket.push_back(psi.back());
      // boundary integral terms relating the auxiliary fields
      formulation.integral(rho0 * 0.5 * (Mx * Mx - 1) * (dvx / vx + dc0 / c0) * k0 * k0 * dof(*psi[0]), tf(v), gammaRight, gauss);
      formulation.integral(-rho0 * 0.5 * (Mx * Mx - 1) * (dvx / vx - Mx * Mx * (dc0 / c0)) * grad(dof(*psi[0])), grad(tf(v)), gammaRight, gauss);
      // coupling of the auxiliary equations
      formulation.integral((1 - Mx * Mx) * grad(dof(*psi[0])), grad(tf(*psi[0])), gammaRight, gauss);
      formulation.integral(-k0 * k0 * dof(*psi[0]), tf(*psi[0]), gammaRight, gauss);
      formulation.integral(dof(v), tf(*psi[0]), gammaRight, gauss);
    }
  }

  formulation.pre();
  formulation.assemble();
  formulation.solve();
  save(+v, omega, "v");

  for(unsigned int i = 0; i < fieldBucket.size(); ++i) {
    delete fieldBucket[i];
  }

  // Compute reference solution with large PML
  msg::info << "Use a PML with N = " << Npml << " layers" << msg::endl;
  const double Lpml = Lext + Npml * h;
  ScalarFunction< std::complex< double > > Sigma0 = 4;
  // Bermudez function
  ScalarFunction< std::complex< double > > SigmaX = Sigma0 / (Lpml - x< std::complex< double > >());
  ScalarFunction< std::complex< double > > gammaX = 1 - (im / k0) * SigmaX;
  ScalarFunction< std::complex< double > > gammaXinv = 1. / gammaX;

  Field< std::complex< double >, Form::Form0 > vRef("vRef", omega | omegaExt | omegaPml | gammaLeft | gammaTop | gammaBottom | gammaTopExt | gammaBottomExt | gammaTopPml | gammaBottomPml, FunctionSpaceTypeForm0::HierarchicalH1, FEMorder);
  formulationRef.integral(rho0 * vector< std::complex< double > >(betax * betax, 0., 0.) * grad(dof(vRef)), vector< std::complex< double > >(1., 0., 0.) * grad(tf(vRef)), omega | omegaExt, gauss);
  formulationRef.integral(rho0 * vector< std::complex< double > >(0., 1., 0.) * grad(dof(vRef)), vector< std::complex< double > >(0., 1., 0.) * grad(tf(vRef)), omega | omegaExt, gauss);
  formulationRef.integral(rho0 * vector< std::complex< double > >(im * k0 * Mx, 0., 0.) * grad(dof(vRef)), tf(vRef), omega | omegaExt, gauss);
  formulationRef.integral(-rho0 * im * k0 * Mx * dof(vRef), vector< std::complex< double > >(1., 0., 0.) * grad(tf(vRef)), omega | omegaExt, gauss);
  formulationRef.integral(-rho0 * k0 * k0 * dof(vRef), tf(vRef), omega | omegaExt, gauss);

  // input BC
  formulationRef.integral(-betax * betax * im * kx * cos(ky * y< std::complex< double > >()), tf(vRef), gammaLeft, gauss);
  formulationRef.integral(-im * k0 * Mx * dof(vRef), tf(vRef), gammaLeft, gauss);

  // pml
  formulationRef.integral(rho0 * vector< std::complex< double > >(gammaXinv * betax * betax, 0., 0.) * grad(dof(vRef)), vector< std::complex< double > >(1., 0., 0.) * grad(tf(vRef)), omegaPml, gauss);
  formulationRef.integral(rho0 * vector< std::complex< double > >(0., gammaX, 0.) * grad(dof(vRef)), vector< std::complex< double > >(0., 1., 0.) * grad(tf(vRef)), omegaPml, gauss);
  formulationRef.integral(rho0 * vector< std::complex< double > >(gammaXinv * im * k0 * Mx, 0., 0.) * grad(dof(vRef)), tf(vRef), omegaPml, gauss);
  formulationRef.integral(-rho0 * gammaXinv * im * k0 * Mx * dof(vRef), vector< std::complex< double > >(1., 0., 0.) * grad(tf(vRef)), omegaPml, gauss);
  formulationRef.integral(-rho0 * k0 * k0 * (gammaX - gammaXinv * Mx * Mx) / (betax * betax) * dof(vRef), tf(vRef), omegaPml, gauss);

  formulationRef.pre();
  formulationRef.assemble();
  formulationRef.solve();
  save(+vRef, omega | omegaExt | omegaPml, "vpml");
  save(vRef - v, omega, "error");

  std::complex< double > num = integrate(pow(abs(vRef - v), 2), omega, gauss);
  std::complex< double > den = integrate(pow(abs(vRef), 2), omega, gauss);
  msg::info << "L_2 error w.r.t pml solution = " << 100. * sqrt(num / den) << " %" << msg::endl;

  return 0;
}

std::complex< double > GetInputKx(double k, double vinf, double v0, double gamm, double cinf, double ky)
{
  std::complex< double > kx;
  const std::complex< double > im(0., 1.);
  double cL0 = cinf * sqrt(1 + (vinf * vinf - v0 * v0) * (gamm - 1) / (cinf * cinf * 2.));
  double beta0 = sqrt(1 - v0 * v0 / (cL0 * cL0));
  double sigma = (k * k) / (cL0 * cL0) - beta0 * beta0 * ky * ky;

  if(sigma >= 0) {
    kx = (-v0 * k / (cL0 * cL0) + sqrt(sigma)) / (beta0 * beta0);
    msg::info << "propagative input mode " << msg::endl;
    if(k > (cL0 * beta0 * beta0 * ky) && k < ky) {
      msg::info << "inverse upstream input mode " << msg::endl;
    }
  }
  else {
    kx = (-v0 * k / (cL0 * cL0) - im * sqrt(abs(sigma))) / (beta0 * beta0);
    msg::info << "evanescent input mode " << msg::endl;
  }

  return kx;
}
