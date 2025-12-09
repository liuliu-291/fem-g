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
  const double pi = 3.14159265358979323846264338327950288;
  const std::complex< double > im(0., 1.);

  // numerical parameters
  int FEMorder = 8;
  gmshFem.userDefinedParameter(FEMorder, "FEMorder");
  std::string gauss = "Gauss" + std::to_string(2 * FEMorder + 1);

  // meshing
  double L = 0.5; // duct length
  gmshFem.userDefinedParameter(L, "L");
  double H = 0.25; // duct height
  gmshFem.userDefinedParameter(H, "H");
  double h = 1. / 40; // meshsize
  gmshFem.userDefinedParameter(h, "h");
  meshDuct(L, H, h, 1);

  // duct problem specifications
  std::string problem = "hard";
  gmshFem.userDefinedParameter(problem, "problem");
  double w = 30; // frequency
  gmshFem.userDefinedParameter(w, "omega");
  std::vector< unsigned int > n = {3}; // single mode
  // std::vector< unsigned int > n = {0,2,4,6,8}; // multiple modes modes
  std::vector< double > A(n.size(), 1.); // select amplitudes
  int N_modes = n.size();
  msg::info << " There are " << N_modes << " input modes " << msg::endl;

  // Choose ABC type
  std::string abcName = "Pade";
  gmshFem.userDefinedParameter(abcName, "abcName");

  // mean flow
  const double M = 0.8;
  double beta = sqrt(1 - M * M);

  std::vector< std::complex< double > > KX(N_modes);
  std::vector< double > KY(N_modes), SIGMA(N_modes);
  for(unsigned int j = 0; j < n.size(); ++j) {
    KY[j] = n[j] * pi / H;
    SIGMA[j] = w * w - beta * beta * KY[j] * KY[j];
    if(SIGMA[j] >= 0) {
      KX[j] = (1 / (beta * beta)) * (-M * w + sqrt(SIGMA[j]));
      msg::info << " Mode " << n[j] << " is propagative " << msg::endl;
      if(w > (beta * KY[j]) && w < KY[j]) {
        msg::info << " Mode " << n[j] << " is inverse upstream ! " << msg::endl;
      }
    }
    else {
      KX[j] = (1 / (beta * beta)) * (-M * w - im * sqrt(abs(SIGMA[j])));
      msg::info << " Mode " << n[j] << " is evanescent " << msg::endl;
    }
    msg::info << " - propagating wavenumber, kx = " << KX[j] << " of mode " << n[j] << msg::endl;
  }

  float PtsByWl = 2 * pi * FEMorder * (1 + M) / w / h;
  msg::info << "********************************" << msg::endl;
  msg::info << " - convected Helmholtz duct problem " << msg::endl;
  msg::info << " - wave number = " << w << msg::endl;
  msg::info << " - FEM basis order = " << FEMorder << "" << msg::endl;
  msg::info << " - Approximate dofs by wavelength = " << PtsByWl << "" << msg::endl;
  msg::info << "********************************" << msg::endl;

  // declare formulation and FEM domains
  std::vector< FieldInterface< std::complex< double > > * > fieldBucket;
  Formulation< std::complex< double > > formulation("helmholtzflow");

  Domain omega(2, 1);
  Domain gammaTop(1, 2);
  Domain gammaLeft(1, 3);
  Domain gammaBottom(1, 4);
  Domain gammaRight(1, 5);
  Domain gammaTopCorner(0, 6);
  Domain gammaBottomCorner(0, 7);

  Field< std::complex< double >, Form::Form0 > v("v", omega | gammaLeft | gammaRight | gammaTop | gammaBottom, FunctionSpaceTypeForm0::HierarchicalH1, FEMorder);
  // define analytical solution
  Function< std::complex< double >, Degree::Degree0 > *solution = nullptr;

  if(problem == "soft") {
    v.addConstraint(gammaTop, 0.);
    v.addConstraint(gammaBottom, 0.);
    solution = new AnalyticalFunction< helmholtz2D::DuctModeSolutionMulti< std::complex< double > > >(w, M, H, A, n, 0., 0., 1);
  }
  else {
    solution = new AnalyticalFunction< helmholtz2D::DuctModeSolutionMulti< std::complex< double > > >(w, M, H, A, n, 0., 0., 0);
  }

  formulation.integral(vector< std::complex< double > >(beta * beta, 0., 0.) * grad(dof(v)), vector< std::complex< double > >(1., 0., 0.) * grad(tf(v)), omega, gauss);
  formulation.integral(vector< std::complex< double > >(0., 1., 0.) * grad(dof(v)), vector< std::complex< double > >(0., 1., 0.) * grad(tf(v)), omega, gauss);
  formulation.integral(vector< std::complex< double > >(im * w * M, 0., 0.) * grad(dof(v)), tf(v), omega, gauss);
  formulation.integral(-im * w * M * dof(v), vector< std::complex< double > >(1., 0., 0.) * grad(tf(v)), omega, gauss);
  formulation.integral(-w * w * dof(v), tf(v), omega, gauss);

  // input BC
  ScalarFunction< std::complex< double > > Input_modal_sum = 0.;
  for(unsigned int j = 0; j < n.size(); ++j) {
    if(problem == "soft") {
      Input_modal_sum = Input_modal_sum + A[j] * KX[j] * sin(KY[j] * y< std::complex< double > >());
    }
    else if(problem == "hard") {
      Input_modal_sum = Input_modal_sum + A[j] * KX[j] * cos(KY[j] * y< std::complex< double > >());
    }
    else {
      msg::error << " this problem is not defined " << msg::endl;
    }
  }
  formulation.integral(-beta * beta * im * Input_modal_sum, tf(v), gammaLeft, gauss);
  formulation.integral(-im * w * M * dof(v), tf(v), gammaLeft, gauss);

  // output BC
  if(abcName == "ABC-0") {
    msg::info << "Use Sommerfeld ABC." << msg::endl;
    const double angle = 0.;
    std::complex< double > wM = w * (M - exp(im * angle / 2.)) / (beta * beta); //angle=0 -> -k/(1+M);
    formulation.integral(-beta * beta * im * wM * dof(v), tf(v), gammaRight, gauss);
    formulation.integral(im * w * M * dof(v), tf(v), gammaRight, gauss);
  }
  else if(abcName == "ABC-2") {
    msg::info << "Use second order ABC." << msg::endl;
    const double angle = 0.;
    // contribution to the mass and rigidity matrices
    std::complex< double > wM = w * (M - cos(angle / 2.)) / (beta * beta);
    std::complex< double > wK = -exp(-im * angle / 2.) / (2 * w);
    formulation.integral(-beta * beta * im * wM * dof(v), tf(v), gammaRight, gauss);
    formulation.integral(beta * beta * im * wK * grad(dof(v)), grad(tf(v)), gammaRight, gauss);
    formulation.integral(im * w * M * dof(v), tf(v), gammaRight, gauss);

    /* use DtN eigenvalues instead of Laplace Beltrami operator
    std::complex< double > kT = -k + KY[0]*KY[0]/(2*k);
    formulation.integral(-beta*beta * im * kT * dof(v), tf(v), gammaRight, gauss);
    */
  }
  else if(abcName == "Pade") {
    int padeOrder = 4;
    gmshFem.userDefinedParameter(padeOrder, "padeOrder");
    double angle = -pi / 4.;
    gmshFem.userDefinedParameter(angle, "angle");
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
      if(problem == "soft") { // impose Dirichlet constraints on the auxiliary fields
        phis->addConstraint(gammaTopCorner, 0.);
        phis->addConstraint(gammaBottomCorner, 0.);
      }
      phi.push_back(phis);
      fieldBucket.push_back(phi.back());
    }

    // write the augmented weak form
    formulation.integral(im * w * exp2 * dof(v), tf(v), gammaRight, gauss);
    for(int i = 0; i < padeOrder; ++i) {
      // boundary integral terms relating the auxiliary fields
      formulation.integral(im * w * exp2 * coef * c[i] * dof(*phi[i]), tf(v), gammaRight, gauss);
      formulation.integral(im * w * exp2 * coef * c[i] * dof(v), tf(v), gammaRight, gauss);

      // coupling of the auxiliary equations
      formulation.integral(grad(dof(*phi[i])), grad(tf(*phi[i])), gammaRight, gauss);
      formulation.integral(-(w * w) / (beta * beta) * (exp1 * c[i] + 1.) * dof(*phi[i]), tf(*phi[i]), gammaRight, gauss);
      formulation.integral(-(w * w) / (beta * beta) * exp1 * (c[i] + 1.) * dof(v), tf(*phi[i]), gammaRight, gauss);
    }
  }
  else if(abcName == "DtN") {
    msg::info << "Use analytical DtN." << msg::endl;
    formulation.integral(beta * beta * im * KX[0] * dof(v), tf(v), gammaRight, gauss);
    formulation.integral(im * w * M * dof(v), tf(v), gammaRight, gauss);
  }

  formulation.pre();
  formulation.assemble();
  formulation.solve();
  save(+v, omega, "v");
  save(*solution, omega, "v_exact");
  save(*solution - v, omega, "error");
  std::complex< double > num = integrate(pow(abs(*solution - v), 2), omega, gauss);
  std::complex< double > den = integrate(pow(abs(*solution), 2), omega, gauss);
  msg::info << "L_2 error = " << 100. * sqrt(num / den) << " %" << msg::endl;

  for(unsigned int i = 0; i < fieldBucket.size(); ++i) {
    delete fieldBucket[i];
  }

  return 0;
}
