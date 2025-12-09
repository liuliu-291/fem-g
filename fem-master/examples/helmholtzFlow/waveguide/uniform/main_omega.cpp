#include <fstream>
#include <gmsh.h>
#include <gmshfem/AnalyticalFunction.h>
#include <gmshfem/Formulation.h>
#include <gmshfem/GmshFem.h>
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
  gmsh::option::setNumber("Mesh.ElementOrder", MeshOrder);
  gmsh::model::mesh::setOrder(MeshOrder);
}

int main(int argc, char **argv)
{
  //*****
  // Problem declaration
  //*****
  GmshFem gmshFem(argc, argv);
  double pi = 3.14159265358979323846264338327950288;
  const std::complex< double > im(0., 1.);

  double L; // duct length
  double H; // duct height
  double h = 1 / 40.; // meshsize
  gmshFem.userDefinedParameter(h, "h");

  int FEMorder = 8;
  gmshFem.userDefinedParameter(FEMorder, "FEMorder");
  std::string gauss = "Gauss" + std::to_string(2 * FEMorder + 1);

  std::vector< unsigned int > n; // single or multiple modes
  const double M = -0.8;
  double beta = sqrt(1 - M * M);

  std::vector< std::string > abcName;
  int runcase = 1;
  gmshFem.userDefinedParameter(runcase, "runcase");
  std::vector< int > NPade;
  std::vector< double > alphaPade;
  double kmax, step;
  //
  switch(runcase) {
  case 1:
    L = 1;
    kmax = 100;
    step = 1;
    H = 0.5;
    n = {6};
    abcName.push_back("Pade");
    abcName.push_back("Pade");
    abcName.push_back("Pade");
    abcName.push_back("Pade");
    abcName.push_back("DtN");
    NPade = {0, 1, 4, 8, 0};
    alphaPade = {0, 0, 0, 0, 0};
    break;
  case 2:
    L = 0.2;
    kmax = 100;
    step = 0.25;
    H = 1;
    n = {8};
    abcName.push_back("Pade");
    abcName.push_back("Pade");
    abcName.push_back("Pade");
    abcName.push_back("Pade");
    abcName.push_back("Pade");
    abcName.push_back("DtN");
    NPade = {4, 4, 4, 4, 4, 0};
    alphaPade = {0, -pi / 4, -pi / 2, -3 * pi / 4, -pi, 0};
    break;
  }

  meshDuct(L, H, h, 1);
  std::vector< double > A(n.size(), 1.); // select amplitudes
  int N_modes = n.size();
  msg::info << " There are " << N_modes << " input modes " << msg::endl;

  Domain omega(2, 1);
  Domain gammaTop(1, 2);
  Domain gammaLeft(1, 3);
  Domain gammaBottom(1, 4);
  Domain gammaRight(1, 5);
  Domain gammaTopCorner(0, 6);
  Domain gammaBottomCorner(0, 7);

  // store L2 errors
  std::vector< std::vector< double > > L2err;
  // frequency loop
  for(int j = 0; j <= kmax; j++) // frequency loop
  {
    double k = 2 + j * step;
    std::vector< double > L2errABC;
    L2errABC.push_back(k);

    std::vector< std::complex< double > > KX(N_modes);
    std::vector< double > KY(N_modes), SIGMA(N_modes);
    for(unsigned int ll = 0; ll < n.size(); ++ll) {
      KY[ll] = n[ll] * pi / H;
      SIGMA[ll] = k * k - beta * beta * KY[ll] * KY[ll];
      if(SIGMA[ll] >= 0) {
        KX[ll] = (1 / (beta * beta)) * (-M * k + sqrt(SIGMA[ll]));
        msg::info << " Mode " << n[ll] << " is propagative " << msg::endl;
        if(k > (beta * KY[ll]) && k < KY[ll]) {
          msg::info << " Mode " << n[ll] << " is inverse upstream ! " << msg::endl;
        }
      }
      else {
        KX[ll] = (1 / (beta * beta)) * (-M * k - im * sqrt(abs(SIGMA[ll])));
        msg::info << " Mode " << n[ll] << " is evanescent " << msg::endl;
      }
      msg::info << " - propagating wavenumber = " << KX[ll] << " of mode " << n[ll] << msg::endl;
    }

    float pointsByWl = (2 * pi * FEMorder * (1 + M) / k) / h;
    msg::info << "********************************" << msg::endl;
    msg::info << " - convected Helmholtz duct problem " << msg::endl;
    msg::info << " - wave number = " << k << msg::endl;
    msg::info << " - FEM basis order = " << FEMorder << "" << msg::endl;
    msg::info << " - Approximate dofs by wavelength = " << pointsByWl << "" << msg::endl;
    msg::info << "********************************" << msg::endl;

    Formulation< std::complex< double > > formulation("helmholtzflow");

    // define analytical solution
    Function< std::complex< double >, Degree::Degree0 > *solution = nullptr;
    solution = new AnalyticalFunction< helmholtz2D::DuctModeSolutionMulti< std::complex< double > > >(k, M, H, A, n, 0., 0., 0);

    std::vector< Field< std::complex< double >, Form::Form0 > * > vNum;
    for(unsigned int j = 0; j <= (abcName.size() - 1); j++) {
      vNum.push_back(new Field< std::complex< double >, Form::Form0 >("vNum_" + std::to_string(j), omega | gammaLeft | gammaRight | gammaTop | gammaBottom, FunctionSpaceTypeForm0::HierarchicalH1, FEMorder));
    }
    std::vector< FieldInterface< std::complex< double > > * > fieldBucket;

    for(unsigned int j = 0; j <= (abcName.size() - 1); j++) {
      formulation.integral(vector< std::complex< double > >(beta * beta, 0., 0.) * grad(dof(*vNum[j])), vector< std::complex< double > >(1., 0., 0.) * grad(tf(*vNum[j])), omega, gauss);
      formulation.integral(vector< std::complex< double > >(0., 1., 0.) * grad(dof(*vNum[j])), vector< std::complex< double > >(0., 1., 0.) * grad(tf(*vNum[j])), omega, gauss);
      formulation.integral(vector< std::complex< double > >(im * k * M, 0., 0.) * grad(dof(*vNum[j])), tf(*vNum[j]), omega, gauss);
      formulation.integral(-im * k * M * dof(*vNum[j]), vector< std::complex< double > >(1., 0., 0.) * grad(tf(*vNum[j])), omega, gauss);
      formulation.integral(-k * k * dof(*vNum[j]), tf(*vNum[j]), omega, gauss);

      // input BC
      ScalarFunction< std::complex< double > > Input_modal_sum = 0.;
      for(unsigned int j = 0; j < n.size(); ++j) {
        Input_modal_sum = Input_modal_sum + A[j] * KX[j] * cos(KY[j] * y< std::complex< double > >());
      }
      formulation.integral(-beta * beta * im * Input_modal_sum, tf(*vNum[j]), gammaLeft, gauss);
      formulation.integral(-im * k * M * dof(*vNum[j]), tf(*vNum[j]), gammaLeft, gauss);

      // output BC
      if(abcName[j] == "ABC-0") {
        msg::info << "Use Sommerfeld ABC." << msg::endl;
        const double angle = 0.;
        std::complex< double > kM = k * (M - exp(im * angle / 2.)) / (beta * beta); //angle=0 -> -k/(1+M);
        formulation.integral(-beta * beta * im * kM * dof(*vNum[j]), tf(*vNum[j]), gammaRight, gauss);
        formulation.integral(im * k * M * dof(*vNum[j]), tf(*vNum[j]), gammaRight, gauss);
      }
      else if(abcName[j] == "ABC-2") {
        msg::info << "Use second order ABC." << msg::endl;
        const double angle = 0.;
        // contribution to the mass and rigidity matrices
        std::complex< double > kM = k * (M - cos(angle / 2.)) / (beta * beta);
        std::complex< double > kK = -exp(-im * angle / 2.) / (2 * k);
        formulation.integral(-beta * beta * im * kM * dof(*vNum[j]), tf(*vNum[j]), gammaRight, gauss);
        formulation.integral(beta * beta * im * kK * grad(dof(*vNum[j])), grad(tf(*vNum[j])), gammaRight, gauss);
        formulation.integral(im * k * M * dof(*vNum[j]), tf(*vNum[j]), gammaRight, gauss);
      }
      else if(abcName[j] == "Pade") {
        int padeOrder = NPade[j];
        const double angle = alphaPade[j];
        msg::info << "Use Pade ABC of order " << padeOrder << " with angle " << angle << " rad" << msg::endl;

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
          phi.push_back(phis);
          fieldBucket.push_back(phi.back());
        }

        // write the augmented weak form
        formulation.integral(im * k * exp2 * dof(*vNum[j]), tf(*vNum[j]), gammaRight, gauss);
        for(int i = 0; i < padeOrder; ++i) {
          // boundary integral terms relating the auxiliary fields
          formulation.integral(im * k * exp2 * coef * c[i] * dof(*phi[i]), tf(*vNum[j]), gammaRight, gauss);
          formulation.integral(im * k * exp2 * coef * c[i] * dof(*vNum[j]), tf(*vNum[j]), gammaRight, gauss);

          // coupling of the auxiliary equations
          formulation.integral(grad(dof(*phi[i])), grad(tf(*phi[i])), gammaRight, gauss);
          formulation.integral(-(k * k) / (beta * beta) * (exp1 * c[i] + 1.) * dof(*phi[i]), tf(*phi[i]), gammaRight, gauss);
          formulation.integral(-(k * k) / (beta * beta) * exp1 * (c[i] + 1.) * dof(*vNum[j]), tf(*phi[i]), gammaRight, gauss);
        }
      }
      else if(abcName[j] == "DtN") {
        msg::info << "Use analytical DtN." << msg::endl;
        formulation.integral(beta * beta * im * KX[0] * dof(*vNum[j]), tf(*vNum[j]), gammaRight, gauss);
        formulation.integral(im * k * M * dof(*vNum[j]), tf(*vNum[j]), gammaRight, gauss);
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
