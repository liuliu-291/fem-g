#include <fstream>
#include <gmsh.h>
#include <gmshfem/Formulation.h>
#include <gmshfem/GmshFem.h>
#include <iostream>

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
  GmshFem gmshFem(argc, argv);
  double pi = 3.14159265358979323846264338327950288;
  const std::complex< double > im(0., 1.);

  int FEMorder = 6;
  std::string gauss = "Gauss" + std::to_string(2 * FEMorder + 1);
  double Lref = 2.; // duct length
  msg::info << " - Reference length Lref = " << Lref << msg::endl;
  double H = 1.; // duct height
  msg::info << " - Duct of heihgt H = " << H << msg::endl;
  double h = 0.025;
  int Npml = 40;
  // transverse mode
  int mode = 4;
  double ky = mode * pi / H;
  // frequency
  double w = 30;
  gmshFem.userDefinedParameter(w, "omega");

  // For abc loop
  std::vector< std::string > abcName;
  std::vector< int > NPade;
  std::vector< double > PadeAngle;

  // transverse speed of sound c_0(y) and density rho_0(y) profiles
  double Ampli = 1.25;
  double A = -0.4;
  double delta = 32;
  ScalarFunction< std::complex< double > > c0, dyc0, rho0, dyrho0;
  int runcase = 1;
  gmshFem.userDefinedParameter(runcase, "runcase");
  switch(runcase) {
  case 1:
    msg::info << " - variable density, constant speed of sound " << msg::endl;
    rho0 = pow(Ampli * (1 + A * exp(-delta * pow(y< std::complex< double > >() - H / 2., 2))), 2);
    dyrho0 = 2 * sqrt(rho0) * (-2 * Ampli * A * delta * (y< std::complex< double > >() - H / 2) * exp(-delta * pow(y< std::complex< double > >() - H / 2., 2)));
    c0 = 1;
    dyc0 = 0;
    abcName.push_back("Pade-S");
    abcName.push_back("Pade-S");
    abcName.push_back("Pade-S");
    abcName.push_back("Pade-S");
    abcName.push_back("Pade-S");
    abcName.push_back("Pade-S");
    NPade = {4, 4, 4, 8, 8, 8};
    PadeAngle = {0., -pi / 4., -pi / 2., 0, -pi / 4., -pi / 2.};
    break;
  case 2:
    msg::info << " - variable speed of sound, constant density " << msg::endl;
    rho0 = 1;
    dyrho0 = 0;
    c0 = Ampli * (1 + A * exp(-delta * pow(y< std::complex< double > >() - H / 2., 2)));
    dyc0 = -2 * Ampli * A * delta * (y< std::complex< double > >() - H / 2) * exp(-delta * pow(y< std::complex< double > >() - H / 2., 2));
    // ABCs tests
    abcName.push_back("Pade");
    abcName.push_back("Pade");
    abcName.push_back("Pade-C");
    abcName.push_back("Pade-C");
    abcName.push_back("Pade-C");
    NPade = {4, 8, 2, 4, 8};
    PadeAngle = {-pi / 4., -pi / 2., 0., -pi / 4., -pi / 2.};
    break;
  }

  ScalarFunction< std::complex< double > > k0 = w / c0;
  double MinPtsByWl = 2 * pi * FEMorder * 0.75 / w / h; // min(c0) = 0.75
  msg::info << "********************************" << msg::endl;
  msg::info << " - Transverse heterogeneous Helmholtz duct problem " << msg::endl;
  msg::info << " - Input wavenumber = " << w << msg::endl;
  msg::info << " - FEM basis order = " << FEMorder << "" << msg::endl;
  msg::info << " - Approx dofs by wavelength = " << MinPtsByWl << "" << msg::endl;
  msg::info << "********************************" << msg::endl;

  // store L2 errors
  std::vector< std::vector< double > > L2err;
  // frequency loop
  for(int j = 0; j <= 70; j++) // length loop L=[0.025,1.75] with step 0.025
  {
    double L = h * (1 + j);
    msg::info << " - ABC position = " << L << msg::endl;
    std::vector< double > L2errABC;
    meshDuct(L, Lref, H, h, Npml, 1);
    L2errABC.push_back(L);

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

    Formulation< std::complex< double > > formulation("helmholtz");
    Formulation< std::complex< double > > formulationRef("helmholtz-ref");

    Field< std::complex< double >, Form::Form0 > vRef("vRef", omega | omegaExt | omegaPml | gammaLeft | gammaTop | gammaBottom | gammaTopExt | gammaBottomExt | gammaTopPml | gammaBottomPml, FunctionSpaceTypeForm0::HierarchicalH1, FEMorder);
    // Compute reference solution with large PML
    msg::info << "Use a PML with N = " << Npml << " layers at Lref = " << Lref << msg::endl;
    double Lpml = Lref + Npml * h;
    double Sigma0 = 40;
    // Bermudez function
    ScalarFunction< std::complex< double > > SigmaX = Sigma0 / (Lpml - x< std::complex< double > >());
    ScalarFunction< std::complex< double > > gammaX = 1 - (im / w) * SigmaX;
    ScalarFunction< std::complex< double > > gammaXinv = 1. / gammaX;

    formulationRef.integral(grad(dof(vRef)), grad(tf(vRef)), omega | omegaExt, gauss);
    formulationRef.integral(-1 / rho0 * vector< std::complex< double > >(0., dyrho0, 0.) * grad(dof(vRef)), tf(vRef), omega | omegaExt, gauss);
    formulationRef.integral(-(w * w) / (c0 * c0) * dof(vRef), tf(vRef), omega | omegaExt, gauss);
    // input BC
    formulationRef.integral(-1. * cos(ky * y< std::complex< double > >()), tf(vRef), gammaLeft, gauss);
    // pml
    formulationRef.integral(vector< std::complex< double > >(gammaXinv, 0., 0.) * grad(dof(vRef)), vector< std::complex< double > >(1., 0., 0.) * grad(tf(vRef)), omegaPml, gauss);
    formulationRef.integral(vector< std::complex< double > >(0., gammaX, 0.) * grad(dof(vRef)), vector< std::complex< double > >(0., 1., 0.) * grad(tf(vRef)), omegaPml, gauss);
    formulationRef.integral(-1 / rho0 * vector< std::complex< double > >(0., gammaX * dyrho0, 0.) * grad(dof(vRef)), tf(vRef), omegaPml, gauss);
    formulationRef.integral(-gammaX * (w * w) / (c0 * c0) * dof(vRef), tf(vRef), omegaPml, gauss);
    // solve
    formulationRef.pre();
    formulationRef.assemble();
    formulationRef.solve();

    // ABC loop
    std::vector< Field< std::complex< double >, Form::Form0 > * > vNum;
    for(unsigned int j = 0; j <= (abcName.size() - 1); j++) {
      vNum.push_back(new Field< std::complex< double >, Form::Form0 >("vNum_" + std::to_string(j), omega | gammaLeft | gammaRight | gammaTop | gammaBottom, FunctionSpaceTypeForm0::HierarchicalH1, FEMorder));
    }
    std::vector< FieldInterface< std::complex< double > > * > fieldBucket;

    for(unsigned int j = 0; j <= (abcName.size() - 1); j++) {
      formulation.integral(grad(dof(*vNum[j])), grad(tf(*vNum[j])), omega, gauss);
      formulation.integral(-1 / rho0 * vector< std::complex< double > >(0., dyrho0, 0.) * grad(dof(*vNum[j])), tf(*vNum[j]), omega, gauss);
      formulation.integral(-(w * w) / (c0 * c0) * dof(*vNum[j]), tf(*vNum[j]), omega, gauss);
      // input BC
      formulation.integral(-1. * cos(ky * y< std::complex< double > >()), tf(*vNum[j]), gammaLeft, gauss);
      // ABC
      if(abcName[j] == "ABC-0") {
        formulation.integral(im * k0 * dof(*vNum[j]), tf(*vNum[j]), gammaRight, gauss);
      }
      else if(abcName[j] == "ABC-2") {
        formulation.integral(im * k0 * dof(*vNum[j]), tf(*vNum[j]), gammaRight, gauss);
        formulation.integral(-im / (2. * k0) * grad(dof(*vNum[j])), grad(tf(*vNum[j])), gammaRight, gauss);
      }
      else if(abcName[j] == "ABC-22") {
        formulation.integral(im * k0 * dof(*vNum[j]), tf(*vNum[j]), gammaRight, gauss);
        formulation.integral(-im / (2. * k0) * grad(dof(*vNum[j])), grad(tf(*vNum[j])), gammaRight, gauss);
        formulation.integral(-im * dyc0 / (2. * w) * tangent< std::complex< double > >() * grad(dof(*vNum[j])), tf(*vNum[j]), gammaRight, gauss);
        formulation.integral(im * dyrho0 / (2. * k0 * rho0) * tangent< std::complex< double > >() * grad(dof(*vNum[j])), tf(*vNum[j]), gammaRight, gauss);
      }
      else if(abcName[j].std::string::find("Pade") != std::string::npos) {
        // pade parameters
        int padeOrder = NPade[j];
        double angle = PadeAngle[j];
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

        if(abcName[j] == "Pade-C") {
          // write the augmented weak form
          formulation.integral(im * w * exp2 * dof(*vNum[j]), tf(*vNum[j]), gammaRight, gauss);
          for(int i = 0; i < padeOrder; ++i) {
            // boundary integral terms relating the auxiliary fields
            formulation.integral(im * w * exp2 * coef * c[i] * dof(*phi[i]), tf(*vNum[j]), gammaRight, gauss);
            formulation.integral(im * w * exp2 * coef * c[i] * dof(*vNum[j]), tf(*vNum[j]), gammaRight, gauss);
            // coupling of the auxiliary equations
            formulation.integral(grad(dof(*phi[i])), grad(tf(*phi[i])), gammaRight, gauss);
            formulation.integral((w * w) * (1. - 1. / (c0 * c0)) * dof(*phi[i]), tf(*phi[i]), gammaRight, gauss);
            formulation.integral(-(w * w) * (exp1 * c[i] + 1.) * dof(*phi[i]), tf(*phi[i]), gammaRight, gauss);
            formulation.integral(-(w * w) * exp1 * (c[i] + 1.) * dof(*vNum[j]), tf(*phi[i]), gammaRight, gauss);
          }
        }
        else {
          // write the augmented weak form
          formulation.integral(im * k0 * exp2 * dof(*vNum[j]), tf(*vNum[j]), gammaRight, gauss);
          for(int i = 0; i < padeOrder; ++i) {
            // boundary integral terms relating the auxiliary fields
            formulation.integral(im * k0 * exp2 * coef * c[i] * dof(*phi[i]), tf(*vNum[j]), gammaRight, gauss);
            formulation.integral(im * k0 * exp2 * coef * c[i] * dof(*vNum[j]), tf(*vNum[j]), gammaRight, gauss);
            // coupling of the auxiliary equations
            formulation.integral(grad(dof(*phi[i])), grad(tf(*phi[i])), gammaRight, gauss);
            if(abcName[j] == "Pade-S") {
              msg::info << "Enrich square-root operator with density variation" << msg::endl;
              formulation.integral(-dyrho0 / rho0 * tangent< std::complex< double > >() * grad(dof(*phi[i])), tf(*phi[i]), gammaRight, gauss);
            }
            formulation.integral(-(k0 * k0) * (exp1 * c[i] + 1.) * dof(*phi[i]), tf(*phi[i]), gammaRight, gauss);
            formulation.integral(-(k0 * k0) * exp1 * (c[i] + 1.) * dof(*vNum[j]), tf(*phi[i]), gammaRight, gauss);
          }
          if(abcName[j] == "Pade-2") {
            // define coefficients for inverse square root approximation
            std::vector< std::complex< double > > d(padeOrder, 0.);
            std::vector< std::complex< double > > R(padeOrder, 0.);
            std::vector< std::complex< double > > S(padeOrder, 0.);
            for(int l = 0; l < padeOrder; ++l) {
              d[l] = 1. + pow(std::tan((l + 0.5) * 0.5 * pi / (double)padeOrder), 2.);
              R[l] = exp2 * d[l] / (double)padeOrder;
              S[l] = 1. + exp1 * (d[l] - 1.);
            }
            // define new auxiliary fields
            std::vector< Field< std::complex< double >, Form::Form0 > * > psi;
            for(int i = 0; i < padeOrder; ++i) {
              psi.push_back(new Field< std::complex< double >, Form::Form0 >("psi_" + std::to_string(i), gammaRight, FunctionSpaceTypeForm0::HierarchicalH1, FEMorder));
              fieldBucket.push_back(psi.back());
              // boundary integral terms relating the auxiliary fields
              formulation.integral(im * dyrho0 / (2 * k0 * rho0) * R[i] * tangent< std::complex< double > >() * grad(dof(*psi[i])), tf(*vNum[j]), gammaRight, gauss);
              // coupling of the auxiliary equations
              formulation.integral(grad(dof(*psi[i])), grad(tf(*psi[i])), gammaRight, gauss);
              formulation.integral(-(k0 * k0) * S[i] * dof(*psi[i]), tf(*psi[i]), gammaRight, gauss);
              formulation.integral((k0 * k0) * dof(*vNum[j]), tf(*psi[i]), gammaRight, gauss);
            }
          }
        }
      }
      else {
        msg::error << "ABC not defined !" << msg::endl;
        exit(0);
      }

    } // end ABC loop

    // solve all ABCs at once
    formulation.pre();
    formulation.assemble();
    formulation.solve();
    // save errors
    for(unsigned int j = 0; j <= (abcName.size() - 1); j++) {
      std::complex< double > num = integrate(pow(abs(*vNum[j] - vRef), 2), omega, gauss);
      std::complex< double > den = integrate(pow(abs(vRef), 2), omega, gauss);
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
