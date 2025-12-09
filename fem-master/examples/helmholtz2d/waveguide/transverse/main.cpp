#include <gmsh.h>
#include <gmshfem/Formulation.h>
#include <gmshfem/GmshFem.h>

using namespace gmshfem;
using namespace gmshfem::common;
using namespace gmshfem::problem;
using namespace gmshfem::domain;
using namespace gmshfem::field;
using namespace gmshfem::function;
using namespace gmshfem::post;
using namespace gmshfem::equation;

void meshDuct(const double L, const double Lext, const double H, const double h, const int Npml, const int MeshOrder)
{
  msg::info << " - Duct of length L = " << L << msg::endl;
  msg::info << " - Duct of heigth H = " << H << msg::endl;
  gmsh::model::add("geometry");

  // duct definition
  gmsh::model::geo::addPoint(L, H, 0., h, 6);
  gmsh::model::geo::addPoint(0, H, 0., h, 7);
  gmsh::model::geo::addPoint(0, 0, 0., h, 8);
  gmsh::model::geo::addPoint(L, 0, 0., h, 9);

  gmsh::model::geo::addLine(6, 7, 5);
  gmsh::model::geo::addLine(7, 8, 6);
  gmsh::model::geo::addLine(8, 9, 7);
  gmsh::model::geo::addLine(9, 6, 8);
  gmsh::model::geo::addCurveLoop({5, 6, 7, 8}, 1);
  gmsh::model::geo::addPlaneSurface({1}, 1);
  gmsh::model::geo::mesh::setTransfiniteSurface(1);
  gmsh::model::geo::mesh::setRecombine(1, 1);

  // pml extension
  gmsh::model::geo::addPoint(L + Lext + Npml * h, H, 0., h, 10);
  gmsh::model::geo::addPoint(L + Lext + Npml * h, 0, 0., h, 11);
  gmsh::model::geo::addLine(6, 10, 12);
  gmsh::model::geo::addLine(10, 11, 13);
  gmsh::model::geo::addLine(11, 9, 14);
  gmsh::model::geo::addCurveLoop({12, 13, 14, 8}, 2);
  gmsh::model::geo::addPlaneSurface({2}, 2);
  gmsh::model::geo::mesh::setTransfiniteSurface(2);
  gmsh::model::geo::mesh::setRecombine(2, 2);

  // physical groups
  gmsh::model::addPhysicalGroup(2, {1}, 1);
  gmsh::model::setPhysicalName(2, 1, "omega");
  gmsh::model::addPhysicalGroup(1, {5}, 2);
  gmsh::model::setPhysicalName(1, 2, "gammaTop");
  gmsh::model::addPhysicalGroup(1, {6}, 3);
  gmsh::model::setPhysicalName(1, 3, "gammaLeft");
  gmsh::model::addPhysicalGroup(1, {7}, 4);
  gmsh::model::setPhysicalName(1, 4, "gammaBottom");
  gmsh::model::addPhysicalGroup(1, {8}, 5);
  gmsh::model::setPhysicalName(1, 5, "gammaRight");

  // physical corners
  gmsh::model::addPhysicalGroup(0, {6}, 6);
  gmsh::model::setPhysicalName(0, 6, "gammaTopCorner");
  gmsh::model::addPhysicalGroup(0, {9}, 7);
  gmsh::model::setPhysicalName(0, 7, "gammaBottomCorner");

  // pml
  gmsh::model::addPhysicalGroup(2, {2}, 8);
  gmsh::model::setPhysicalName(2, 8, "omegaPml");
  gmsh::model::addPhysicalGroup(1, {12}, 9);
  gmsh::model::setPhysicalName(1, 9, "gammaTopPml");
  gmsh::model::addPhysicalGroup(1, {13}, 10);
  gmsh::model::setPhysicalName(1, 10, "gammaRightPml");
  gmsh::model::addPhysicalGroup(1, {14}, 11);
  gmsh::model::setPhysicalName(1, 11, "gammaBottomPml");

  // generate mesh
  gmsh::model::geo::synchronize();
  gmsh::model::mesh::generate(2);
  gmsh::model::mesh::setOrder(MeshOrder);
  gmsh::option::setNumber("Mesh.ElementOrder", MeshOrder);
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
  int FEMorder = 6;
  gmshFem.userDefinedParameter(FEMorder, "FEMorder");
  std::string gauss = "Gauss" + std::to_string(2 * FEMorder + 1);

  // meshing
  double L = 2; // duct length
  gmshFem.userDefinedParameter(L, "L");
  double H = 1; // duct height
  gmshFem.userDefinedParameter(H, "H");
  double h = 1. / 40;
  gmshFem.userDefinedParameter(h, "h");
  int Npml = 40;
  gmshFem.userDefinedParameter(Npml, "Npml");
  double Lext = 0;
  gmshFem.userDefinedParameter(Lext, "Lext");
  meshDuct(L, Lext, H, h, Npml, 4);

  // duct problem specifications
  std::string problem = "hard";
  gmshFem.userDefinedParameter(problem, "problem");
  double w = 30;
  gmshFem.userDefinedParameter(w, "omega");
  int mode = 4; // transverse mode for the input BC
  gmshFem.userDefinedParameter(mode, "mode");
  double ky = mode * pi / H;

  // choose ABC type
  // Available conditions:
  // rho0y: ABC-0, ABC-2, ABC-22, Pade, Pade-S, Pade-2
  // c0: ABC-0, ABC-2, ABC-22, Pade, Pade-C
  std::string abcName = "Pade";
  gmshFem.userDefinedParameter(abcName, "abcName");

  // transverse speed of sound c_0(y) and density rho_0(y) profiles
  double Ampli = 1.25;
  double A = -0.4;
  double delta = 32;
  double deltaH = H / 4;
  ScalarFunction< std::complex< double > > c0, dyc0, rho0, dyrho0;
  double minc0;
  int runcase = 3;
  gmshFem.userDefinedParameter(runcase, "runcase");
  switch(runcase) {
  case 1:
    msg::info << " - variable density, constant speed of sound " << msg::endl;
    rho0 = pow(Ampli * (1 + A * exp(-delta * pow(y< std::complex< double > >() - H / 2., 2))), 2);
    dyrho0 = 2 * sqrt(rho0) * (-2 * Ampli * A * delta * (y< std::complex< double > >() - H / 2) * exp(-delta * pow(y< std::complex< double > >() - H / 2., 2)));
    c0 = 1;
    dyc0 = 0;
    minc0 = 1;
    break;
  case 2:
    msg::info << " - variable speed of sound, constant density " << msg::endl;
    rho0 = 1;
    dyrho0 = 0;
    c0 = Ampli * (1 + A * exp(-delta * pow(y< std::complex< double > >() - H / 2., 2)));
    dyc0 = -2 * Ampli * A * delta * (y< std::complex< double > >() - H / 2) * exp(-delta * pow(y< std::complex< double > >() - H / 2., 2));
    minc0 = 0.75;
    break;
  case 3:
    msg::info << " - piecewise constant speed of sound, constant density " << msg::endl;
    rho0 = 1;
    dyrho0 = 0;
    c0 = 1.0 - 0.75 * heaviside(y< std::complex< double > >() - H / 2. - deltaH) + 0.75 * heaviside(y< std::complex< double > >() - H / 2. + deltaH);
    //c0 = 1.0 + 1.0*heaviside(y< std::complex<double> >() - H/8.) - 1.75*heaviside(y< std::complex<double> >() - 0.4)
    //+ 1.75*heaviside(y< std::complex<double> >() - 0.6) - 1.0*heaviside(y< std::complex<double> >() - 7*H/8.);
    dyc0 = 0;
    minc0 = 0.25;
  }
  // wavenumber
  ScalarFunction< std::complex< double > > k0 = w / c0;

  double MinPtsByWl = 2 * pi * FEMorder * minc0 / w / h; // min(c0) = 1, 0.75 or 0.25
  msg::info << "********************************" << msg::endl;
  msg::info << " - heterogeneous Helmholtz duct problem in the transverse direction " << msg::endl;
  msg::info << " - Input wavenumber = " << w << msg::endl;
  msg::info << " - FEM basis order = " << FEMorder << "" << msg::endl;
  msg::info << " - Min dofs by wavelength = " << MinPtsByWl << "" << msg::endl;
  msg::info << "********************************" << msg::endl;

  // declare formulation and FEM domains
  std::vector< FieldInterface< std::complex< double > > * > fieldBucket;
  Formulation< std::complex< double > > formulation("helmholtz-duct-y");
  Formulation< std::complex< double > > formulationRef("helmholtz-duct-y-refPML");

  // domains
  Domain omega(2, 1);
  Domain omegaPml(2, 8);
  // boundaries
  Domain gammaTop(1, 2);
  Domain gammaLeft(1, 3);
  Domain gammaBottom(1, 4);
  Domain gammaRight(1, 5);
  // pml boundaries
  Domain gammaTopPml(1, 9);
  Domain gammaRightPml(1, 10);
  Domain gammaBottomPml(1, 11);
  // corners
  Domain gammaTopCorner(0, 6);
  Domain gammaBottomCorner(0, 7);

  Field< std::complex< double >, Form::Form0 > v("v", omega | gammaLeft | gammaRight | gammaTop | gammaBottom, FunctionSpaceTypeForm0::HierarchicalH1, FEMorder);
  if(problem == "soft") {
    v.addConstraint(gammaTop, 0.);
    v.addConstraint(gammaBottom, 0.);
  }

  // weak formulation
  formulation.integral(grad(dof(v)), grad(tf(v)), omega, gauss);
  formulation.integral(-1 / rho0 * vector< std::complex< double > >(0., dyrho0, 0.) * grad(dof(v)), tf(v), omega, gauss);
  formulation.integral(-(w * w) / (c0 * c0) * dof(v), tf(v), omega, gauss);

  // input boundary condition
  if(problem == "soft") {
    formulation.integral(-sin(ky * y< std::complex< double > >()), tf(v), gammaLeft, gauss);
  }
  else {
    formulation.integral(-cos(ky * y< std::complex< double > >()), tf(v), gammaLeft, gauss);
  }

  // output BC
  if(abcName == "ABC-0") {
    msg::info << "Use Sommerfeld ABC." << msg::endl;
    formulation.integral(im * k0 * dof(v), tf(v), gammaRight, gauss);
  }
  else if(abcName == "ABC-2") {
    msg::info << "Use second order ABC." << msg::endl;
    formulation.integral(im * k0 * dof(v), tf(v), gammaRight, gauss);
    formulation.integral(-im / (2. * k0) * grad(dof(v)), grad(tf(v)), gammaRight, gauss);
  }
  else if(abcName == "ABC-22") {
    msg::info << "Use second order ABC with 2 symbols " << msg::endl;
    formulation.integral(im * k0 * dof(v), tf(v), gammaRight, gauss);
    formulation.integral(-im / (2. * k0) * grad(dof(v)), grad(tf(v)), gammaRight, gauss);
    formulation.integral(-im * dyc0 / (2. * w) * tangent< std::complex< double > >() * grad(dof(v)), tf(v), gammaRight, gauss);
    formulation.integral(im * dyrho0 / (2. * k0 * rho0) * tangent< std::complex< double > >() * grad(dof(v)), tf(v), gammaRight, gauss);
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
    Field< std::complex< double >, Form::Form0 > *phis;
    for(int i = 0; i < padeOrder; ++i) {
      phis = new Field< std::complex< double >, Form::Form0 >("phi_" + std::to_string(i), gammaRight | gammaTopCorner | gammaBottomCorner, FunctionSpaceTypeForm0::HierarchicalH1, FEMorder);
      if(problem == "soft") { // impose Dirichlet constraints on the auxiliary fields
        phis->addConstraint(gammaTopCorner, 0.);
        phis->addConstraint(gammaBottomCorner, 0.);
      }
      phi.push_back(phis);
      fieldBucket.push_back(phi.back());
    }

    if(abcName == "Pade-C") {
      msg::info << "Use modified square-root operator for variable speed of sound" << msg::endl;
      // write the augmented weak form
      formulation.integral(im * w * exp2 * dof(v), tf(v), gammaRight, gauss);
      for(int i = 0; i < padeOrder; ++i) {
        // boundary integral terms relating the auxiliary fields
        formulation.integral(im * w * exp2 * coef * c[i] * dof(*phi[i]), tf(v), gammaRight, gauss);
        formulation.integral(im * w * exp2 * coef * c[i] * dof(v), tf(v), gammaRight, gauss);
        // coupling of the auxiliary equations
        formulation.integral(grad(dof(*phi[i])), grad(tf(*phi[i])), gammaRight, gauss);
        formulation.integral((w * w) * (1. - 1. / (c0 * c0)) * dof(*phi[i]), tf(*phi[i]), gammaRight, gauss);
        formulation.integral(-(w * w) * (exp1 * c[i] + 1.) * dof(*phi[i]), tf(*phi[i]), gammaRight, gauss);
        formulation.integral(-(w * w) * exp1 * (c[i] + 1.) * dof(v), tf(*phi[i]), gammaRight, gauss);
      }
    }
    else {
      // write the augmented weak form
      formulation.integral(im * k0 * exp2 * dof(v), tf(v), gammaRight, gauss);
      for(int i = 0; i < padeOrder; ++i) {
        // boundary integral terms relating the auxiliary fields
        formulation.integral(im * k0 * exp2 * coef * c[i] * dof(*phi[i]), tf(v), gammaRight, gauss);
        formulation.integral(im * k0 * exp2 * coef * c[i] * dof(v), tf(v), gammaRight, gauss);
        // coupling of the auxiliary equations
        formulation.integral(grad(dof(*phi[i])), grad(tf(*phi[i])), gammaRight, gauss);
        if(abcName == "Pade-S") {
          msg::info << "Enrich square-root operator with density variation" << msg::endl;
          //          formulation.integral( -dyrho0/rho0 * tangent< std::complex< double > >() * grad(dof(*phi[i])), tf(*phi[i]), gammaRight, gauss);
        }
        formulation.integral(-(k0 * k0) * (exp1 * c[i] + 1.) * dof(*phi[i]), tf(*phi[i]), gammaRight, gauss);
        formulation.integral(-(k0 * k0) * exp1 * (c[i] + 1.) * dof(v), tf(*phi[i]), gammaRight, gauss);
      }
    }
    if(abcName == "Pade-2") {
      msg::info << "Use a second Pade localization for the zeroth order symbol of order " << padeOrder << " with angle " << angle << " rad" << msg::endl;
      // define coefficients for inverse square root approximation
      std::vector< std::complex< double > > d(padeOrder, 0.);
      std::vector< std::complex< double > > R(padeOrder, 0.);
      std::vector< std::complex< double > > S(padeOrder, 0.);
      for(int l = 0; l < padeOrder; ++l) {
        d[l] = 1. + pow(std::tan((l + 0.5) * 0.5 * pi / (double)padeOrder), 2.);
        R[l] = exp2 * d[l] / (double)padeOrder;
        S[l] = 1. + exp1 * (d[l] - 1.);
      }
      // define a new auxiliary field
      std::vector< Field< std::complex< double >, Form::Form0 > * > psi;
      Field< std::complex< double >, Form::Form0 > *psis;
      for(int i = 0; i < padeOrder; ++i) {
        psis = new Field< std::complex< double >, Form::Form0 >("psi_" + std::to_string(i), gammaRight | gammaTopCorner | gammaBottomCorner, FunctionSpaceTypeForm0::HierarchicalH1, FEMorder);
        if(problem == "soft") { // impose Dirichlet constraints on the auxiliary fields
          psis->addConstraint(gammaTopCorner, 0.);
          psis->addConstraint(gammaBottomCorner, 0.);
        }
        psi.push_back(psis);
        fieldBucket.push_back(psi.back());
        // boundary integral terms relating the auxiliary fields
        formulation.integral(im * dyrho0 / (2 * k0 * rho0) * R[i] * tangent< std::complex< double > >() * grad(dof(*psi[i])), tf(v), gammaRight, gauss);
        // coupling of the auxiliary equations
        formulation.integral(grad(dof(*psi[i])), grad(tf(*psi[i])), gammaRight, gauss);
        formulation.integral(-(k0 * k0) * S[i] * dof(*psi[i]), tf(*psi[i]), gammaRight, gauss);
        formulation.integral((k0 * k0) * dof(v), tf(*psi[i]), gammaRight, gauss);
      }
    }
  }
  else {
    msg::error << "ABC not defined !" << msg::endl;
    exit(0);
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
  // Bermudez PML function
  double Lpml = L + Lext + Npml * h;
  double Sigma0 = 40;
  ScalarFunction< std::complex< double > > SigmaX = heaviside(x< std::complex< double > >() - (L + Lext)) * Sigma0 / (Lpml - x< std::complex< double > >());
  // Cubic PML function
  // double R = 1000;
  // ScalarFunction< std::complex< double > > SigmaX = 3*R*pow( (x< std::complex< double > >() - L)/(Npml*h) ,2);
  ScalarFunction< std::complex< double > > gammaX = 1 - (im / w) * SigmaX;
  ScalarFunction< std::complex< double > > gammaXinv = 1. / gammaX;

  Field< std::complex< double >, Form::Form0 > vRef("vRef", omega | omegaPml | gammaLeft | gammaRight | gammaTop | gammaBottom | gammaTopPml | gammaBottomPml, FunctionSpaceTypeForm0::HierarchicalH1, FEMorder);
  if(problem == "soft") {
    vRef.addConstraint(gammaTop, 0.);
    vRef.addConstraint(gammaBottom, 0.);
    vRef.addConstraint(gammaTopPml, 0.);
    vRef.addConstraint(gammaBottomPml, 0.);
  }
  // weak formulation
  formulationRef.integral(grad(dof(vRef)), grad(tf(vRef)), omega, gauss);
  formulationRef.integral(-1 / rho0 * vector< std::complex< double > >(0., dyrho0, 0.) * grad(dof(vRef)), tf(vRef), omega, gauss);
  formulationRef.integral(-(w * w) / (c0 * c0) * dof(vRef), tf(vRef), omega, gauss);
  // weak formulation in the pml domain
  formulationRef.integral(vector< std::complex< double > >(gammaXinv, 0., 0.) * grad(dof(vRef)), vector< std::complex< double > >(1., 0., 0.) * grad(tf(vRef)), omegaPml, gauss);
  formulationRef.integral(vector< std::complex< double > >(0., gammaX, 0.) * grad(dof(vRef)), vector< std::complex< double > >(0., 1., 0.) * grad(tf(vRef)), omegaPml, gauss);
  formulationRef.integral(-1 / rho0 * vector< std::complex< double > >(0., gammaX * dyrho0, 0.) * grad(dof(vRef)), tf(vRef), omegaPml, gauss);
  formulationRef.integral(-gammaX * (w * w) / (c0 * c0) * dof(vRef), tf(vRef), omegaPml, gauss);

  // input BC
  if(problem == "soft") {
    formulationRef.integral(-sin(ky * y< std::complex< double > >()), tf(vRef), gammaLeft, gauss);
  }
  else {
    formulationRef.integral(-cos(ky * y< std::complex< double > >()), tf(vRef), gammaLeft, gauss);
  }
  formulationRef.pre();
  formulationRef.assemble();
  formulationRef.solve();

  save(+vRef, omega | omegaPml, "vRef");
  save(vRef - v, omega, "error");
  std::complex< double > num = integrate(pow(abs(vRef - v), 2), omega, gauss);
  std::complex< double > den = integrate(pow(abs(vRef), 2), omega, gauss);
  msg::info << "L_2 error = " << 100. * sqrt(num / den) << " %" << msg::endl;

  return 0;
}
