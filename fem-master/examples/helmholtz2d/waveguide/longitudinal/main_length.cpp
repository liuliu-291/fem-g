#include <fstream>
#include <gmsh.h>
#include <gmshfem/AnalyticalFunction.h>
#include <gmshfem/Bessel.h>
#include <gmshfem/Formulation.h>
#include <gmshfem/Function.h>
#include <gmshfem/GmshFem.h>
#include <gmshfem/Message.h>
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

std::complex< double > EvalAiryArg(double, double, double, double, double, std::complex< double >);
std::complex< double > airy(std::complex< double >, bool);
double GetepsilonOpt(double, double, double, double, double);

void meshDuct(const double L, const double H, const double h, const int MeshOrder)
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

  // physical groups
  gmsh::model::addPhysicalGroup(2, {1}, 1);
  gmsh::model::setPhysicalName(2, 3, "omega");
  gmsh::model::addPhysicalGroup(1, {5}, 2);
  gmsh::model::setPhysicalName(1, 2, "gammaTop");
  gmsh::model::addPhysicalGroup(1, {6}, 3);
  gmsh::model::setPhysicalName(1, 3, "gammaLeft");
  gmsh::model::addPhysicalGroup(1, {7}, 4);
  gmsh::model::setPhysicalName(1, 4, "gammaBottom");
  gmsh::model::addPhysicalGroup(1, {8}, 5);
  gmsh::model::setPhysicalName(1, 5, "gammaRight");

  gmsh::model::addPhysicalGroup(0, {6}, 6);
  gmsh::model::setPhysicalName(0, 6, "gammaTopCorner");
  gmsh::model::addPhysicalGroup(0, {9}, 7);
  gmsh::model::setPhysicalName(0, 7, "gammaBottomCorner");
  // generate mesh
  gmsh::model::geo::synchronize();
  gmsh::model::mesh::generate(2);
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
  const std::complex< double > exp_T = std::complex< double >(std::cos(2 * pi / 3.), -std::sin(2 * pi / 3.));

  // numerical parameters
  int FEMorder = 6;
  gmshFem.userDefinedParameter(FEMorder, "FEMorder");
  std::string gauss = "Gauss" + std::to_string(2 * FEMorder + 1);

  // meshing
  double H = 0.5; // duct height
  gmshFem.userDefinedParameter(H, "H");
  double h = 1. / 40.; // mesh size
  gmshFem.userDefinedParameter(h, "h");

  // duct problem specifications
  std::string problem = "hard";
  gmshFem.userDefinedParameter(problem, "problem");
  double w = 19.5; // frequency
  gmshFem.userDefinedParameter(w, "omega");
  int mode = 5; // transverse mode
  gmshFem.userDefinedParameter(mode, "mode");
  double ky = mode * pi / H;

  // speed of sound c_0(x)^(-2)=nbar*x+ni, linear profile prescribed analytically
  double nbar = 5.0;
  double ni = 0.1;
  ScalarFunction< std::complex< double > > c0_m2 = nbar * x< std::complex< double > >() + ni;
  ScalarFunction< std::complex< double > > dx_c0_m2 = nbar; // x-derivative

  std::vector< std::string > abcName;
  std::vector< int > NPade;
  std::vector< double > angles, eps_fact;
  int runcase = 1;
  gmshFem.userDefinedParameter(runcase, "runcase");
  switch(runcase) {
  case 1:
    abcName.push_back("Pade-S2");
    abcName.push_back("Pade-S2");
    abcName.push_back("Pade-S2");
    abcName.push_back("Pade-S2");
    abcName.push_back("DtN");
    NPade = {3, 3, 3, 3, 3};
    angles = {0, -pi / 4., -pi / 2., -3 * pi / 4., 0};
    eps_fact = {0, 0, 0, 0, 0};
    break;
  case 2:
    abcName.push_back("Pade-S2");
    abcName.push_back("Pade-S2");
    abcName.push_back("Pade-S2");
    abcName.push_back("Pade-S2");
    abcName.push_back("DtN");
    NPade = {5, 5, 5, 5, 0};
    angles = {-pi / 2., -pi / 2., -pi / 2., -pi / 2., 0};
    eps_fact = {0, 0.5, 1, 1.5, 0};
    break;
  case 3:
    abcName.push_back("Pade-S2");
    abcName.push_back("Pade-S2");
    abcName.push_back("Pade-S2");
    abcName.push_back("Pade-S2");
    abcName.push_back("Pade-S2");
    abcName.push_back("DtN");
    NPade = {1, 2, 3, 5, 7, 0};
    angles = {-pi / 2., -pi / 2., -pi / 2., -pi / 2., -pi / 2., 0};
    eps_fact = {0, 0, 0, 0, 0, 0};
    break;
  }

  double sigma = w * w * ni - ky * ky; // radicand of principal symbol at x=0
  if(sigma >= 0) {
    msg::info << " Input propagative mode " << msg::endl;
  }
  else {
    msg::info << " Input evanescent mode " << msg::endl;
  }
  ScalarFunction< std::complex< double > > k0 = w * sqrt(c0_m2);
  double Xt = ky * ky / (nbar * w * w) - ni / nbar; // evalaute location of the turning point, i.e where the radicand of the principal symbol vanishes
  // use optimal epsilon for cut-on/cut-off transistion
  double epsopt = GetepsilonOpt(ky, ni, nbar, Xt, pi);

  msg::info << "********************************" << msg::endl;
  msg::info << " - heterogeneous Helmholtz duct problem, variable speed of sound along the wavefront " << msg::endl;
  msg::info << " - Input wavenumber = " << w << msg::endl;
  msg::info << " - Input mode = " << mode << msg::endl;
  msg::info << " - FEM basis order = " << FEMorder << "" << msg::endl;
  msg::info << "********************************" << msg::endl;

  // store L2 errors
  std::vector< std::vector< double > > L2err;
  // frequency loop
  for(int j = 0; j <= 39; j++) {
    double L = h * (1 + j);
    msg::info << " - ABC position = " << L << msg::endl;
    std::vector< double > L2errABC;
    meshDuct(L, H, h, 1);
    L2errABC.push_back(L);

    if((Xt >= 0) && (Xt <= L)) {
      msg::info << " Turning point within the duct at Xt = " << Xt << msg::endl;
    }
    else {
      msg::info << " Turning point outside of the physical domain at Xt = " << Xt << msg::endl;
    }

    // domains
    Domain omega(2, 1);
    Domain gammaTop(1, 2);
    Domain gammaLeft(1, 3);
    Domain gammaBottom(1, 4);
    Domain gammaRight(1, 5);
    Domain gammaTopCorner(0, 6);
    Domain gammaBottomCorner(0, 7);

    // declare formulation
    Formulation< std::complex< double > > formulation("helmholtz-duct-c0x");
    Field< std::complex< double >, Form::Form0 > v("v", omega | gammaLeft | gammaRight | gammaTop | gammaBottom, FunctionSpaceTypeForm0::HierarchicalH1, FEMorder);
    // define analytical solution
    Function< std::complex< double >, Degree::Degree0 > *solution = nullptr;

    if(problem == "soft") {
      v.addConstraint(gammaTop, 0.);
      v.addConstraint(gammaBottom, 0.);
      solution = new AnalyticalFunction< helmholtz2D::DuctModeSolutionAiry< std::complex< double > > >(w, H, nbar, ni, 0., 0., mode, 1);
    }
    else {
      solution = new AnalyticalFunction< helmholtz2D::DuctModeSolutionAiry< std::complex< double > > >(w, H, nbar, ni, 0., 0., mode, 0);
    }

    // ABC loop
    std::vector< Field< std::complex< double >, Form::Form0 > * > vNum;
    for(unsigned int j = 0; j <= (abcName.size() - 1); j++) {
      vNum.push_back(new Field< std::complex< double >, Form::Form0 >("vNum_" + std::to_string(j), omega | gammaLeft | gammaRight | gammaTop | gammaBottom, FunctionSpaceTypeForm0::HierarchicalH1, FEMorder));
    }
    std::vector< FieldInterface< std::complex< double > > * > fieldBucket;

    for(unsigned int j = 0; j <= (abcName.size() - 1); j++) {
      formulation.integral(grad(dof(*vNum[j])), grad(tf(*vNum[j])), omega, gauss);
      formulation.integral(-w * w * c0_m2 * dof(*vNum[j]), tf(*vNum[j]), omega, gauss);

      // input BC
      std::complex< double > z0 = EvalAiryArg(0., w, nbar, ni, ky, exp_T);
      std::complex< double > g = exp_T * pow(nbar * w * w, 1. / 3.) * airy(z0, 1);
      if(problem == "soft") {
        formulation.integral(-g * sin(ky * y< std::complex< double > >()), tf(*vNum[j]), gammaLeft, gauss);
      }
      else {
        formulation.integral(-g * cos(ky * y< std::complex< double > >()), tf(*vNum[j]), gammaLeft, gauss);
      }

      // output BC
      if(abcName[j] == "ABC-0") {
        msg::info << "Use Sommerfeld ABC." << msg::endl;
        formulation.integral(im * k0 * dof(*vNum[j]), tf(*vNum[j]), gammaRight, gauss);
      }
      else if(abcName[j] == "ABC-2") {
        msg::info << "Use second order ABC." << msg::endl;
        formulation.integral(im * k0 * dof(*vNum[j]), tf(*vNum[j]), gammaRight, gauss);
        formulation.integral(-im / (2. * k0) * grad(dof(*vNum[j])), grad(tf(*vNum[j])), gammaRight, gauss);
      }
      else if(abcName[j] == "ABC-22") {
        msg::info << "Use second order ABC with 2 symbols." << msg::endl;
        formulation.integral(im * k0 * dof(*vNum[j]), tf(*vNum[j]), gammaRight, gauss);
        formulation.integral(-im / (2. * k0) * grad(dof(*vNum[j])), grad(tf(*vNum[j])), gammaRight, gauss);
        formulation.integral(dx_c0_m2 / (4. * c0_m2) * dof(*vNum[j]), tf(*vNum[j]), gammaRight, gauss);
        formulation.integral(dx_c0_m2 / (4. * k0 * k0) * grad(dof(*vNum[j])), grad(tf(*vNum[j])), gammaRight, gauss);
      }
      else if(abcName[j].std::string::find("Pade") != std::string::npos) { // if the string contains "Pade"
        int padeOrder = NPade[j];
        gmshFem.userDefinedParameter(padeOrder, "padeOrder");
        double angle = angles[j];
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
        // write the augmented weak form
        formulation.integral(im * k0 * exp2 * dof(*vNum[j]), tf(*vNum[j]), gammaRight, gauss);
        for(int i = 0; i < padeOrder; ++i) {
          // boundary integral terms relating the auxiliary fields
          formulation.integral(im * k0 * exp2 * coef * c[i] * dof(*phi[i]), tf(*vNum[j]), gammaRight, gauss);
          formulation.integral(im * k0 * exp2 * coef * c[i] * dof(*vNum[j]), tf(*vNum[j]), gammaRight, gauss);

          // coupling of the auxiliary equations
          formulation.integral(grad(dof(*phi[i])), grad(tf(*phi[i])), gammaRight, gauss);
          formulation.integral(-(k0 * k0) * (exp1 * c[i] + 1.) * dof(*phi[i]), tf(*phi[i]), gammaRight, gauss);
          formulation.integral(-(k0 * k0) * exp1 * (c[i] + 1.) * dof(*vNum[j]), tf(*phi[i]), gammaRight, gauss);
        }

        if(abcName[j] == "Pade-S2") {
          msg::info << "Use Pade ABC with 2 symbols." << msg::endl;
          // define a new auxiliary field
          std::vector< Field< std::complex< double >, Form::Form0 > * > psi;
          Field< std::complex< double >, Form::Form0 > *psis;
          psis = new Field< std::complex< double >, Form::Form0 >("psi", gammaRight | gammaTopCorner | gammaBottomCorner, FunctionSpaceTypeForm0::HierarchicalH1, FEMorder);
          if(problem == "soft") { // impose Dirichlet constraints on the auxiliary fields
            psis->addConstraint(gammaTopCorner, 0.);
            psis->addConstraint(gammaBottomCorner, 0.);
          }
          psi.push_back(psis);
          fieldBucket.push_back(psi.back());

          // boundary integral terms relating the auxiliary fields
          formulation.integral(dx_c0_m2 / (4.0 * c0_m2) * dof(*psi[0]), tf(*vNum[j]), gammaRight, gauss);

          // define complexified frequency
          std::complex< double > weps = w - im * epsopt * eps_fact[j];
          if(eps_fact[j] != 0) {
            msg::info << "Use regularization parameter for grazing waves, epsilon = " << epsopt * eps_fact[j] << msg::endl;
          }
          // coupling of the auxiliary equations
          formulation.integral(grad(dof(*psi[0])), grad(tf(*psi[0])), gammaRight, gauss);
          formulation.integral(-(weps * weps) * (c0_m2)*dof(*psi[0]), tf(*psi[0]), gammaRight, gauss);
          formulation.integral((weps * weps) * (c0_m2)*dof(*vNum[j]), tf(*psi[0]), gammaRight, gauss);
        }
      }
      else if(abcName[j] == "DtN") {
        msg::info << "Use analytical DtN." << msg::endl;
        std::complex< double > zL = EvalAiryArg(L, w, nbar, ni, ky, exp_T);
        std::complex< double > kx = -im * exp_T * pow(nbar * w * w, 1. / 3.) * airy(zL, 1) / airy(zL, 0);
        formulation.integral(im * kx * dof(*vNum[j]), tf(*vNum[j]), gammaRight, gauss);
      }
      else {
        msg::error << "ABC not defined !" << msg::endl;
        exit(0);
      }
    } // end ABC loop
    formulation.pre();
    formulation.assemble();
    formulation.solve();
    // save errors
    for(unsigned int j = 0; j <= (abcName.size() - 1); j++) {
      std::complex< double > num = integrate(pow(abs(*vNum[j] - *solution), 2), omega, gauss);
      std::complex< double > den = integrate(pow(abs(*solution), 2), omega, gauss);
      double err = std::real(100. * sqrt(num / den));
      msg::info << "L_2 error w.r.t pml solution = " << err << " with ABC " << abcName[j] << msg::endl;
      L2errABC.push_back(err);
    }

    for(unsigned int i = 0; i < fieldBucket.size(); ++i) {
      delete fieldBucket[i];
    }

    L2err.push_back(L2errABC);
  } // end frequency loop

  std::ofstream Errfile;
  Errfile.open("L2-err-case" + std::to_string(runcase) + ".txt");
  // write results
  for(unsigned int i = 0; i < L2err.size(); i++) {
    for(unsigned int j = 0; j < L2err[i].size(); j++) {
      Errfile << L2err[i][j] << " ";
    }
    Errfile << std::endl;
  }
  Errfile.close();

  return 0;
}


std::complex< double > EvalAiryArg(double x, double w, double a, double b, double ky, std::complex< double > angle)
{
  return angle * (ky * ky - w * w * (a * x + b)) / pow(w * w * a, 2. / 3.);
}

std::complex< double > airy(std::complex< double > z, bool derivative)
{

  double valr, vali;
  int err = AiryComplex(std::real(z), std::imag(z), derivative, &valr, &vali);

  if(err != 0) {
    msg::warning << "issue with complex airy function, error output " << err << msg::endl;
  }

  return std::complex< double >(valr, vali);
}

double GetepsilonOpt(double ky, double ni, double nbar, double Xt, double pi)
{
  double wt = ky / std::sqrt(nbar * Xt + ni); // turning frequency
  double s = std::sin(2. * pi / 3.);

  // optimal epsilon for cut-on/cut-off transition
  double eps_opt;
  eps_opt = 2. * wt * s * (1 - std::sqrt(1 + std::real(airy(std::complex< double >(0, 0), 0)) * pow(nbar * wt * wt, 2. / 3.) / (8. * pow(s * ky, 2) * std::real(airy(std::complex< double >(0, 0), 1)))));
  msg::info << " Optimal regularization parameter at turning point is " << eps_opt << msg::endl;
  return eps_opt;
}
