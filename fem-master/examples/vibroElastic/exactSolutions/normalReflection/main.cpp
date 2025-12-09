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

using Scalar = std::complex< double >;

// Defined below, after main()
static void setupMesh(double lc, double L = 1, double waterDepth = 1, double groundDepth = 1);

int main(int argc, char **argv)
{
  GmshFem gmshFem(argc, argv);

  double lc = 1e-2;
  double waterDepth = 1;
  double groundDepth = 1;
  double L = 1;
  gmshFem.userDefinedParameter(lc, "lc");
  gmshFem.userDefinedParameter(waterDepth, "waterDepth");
  gmshFem.userDefinedParameter(groundDepth, "groundDepth");
  gmshFem.userDefinedParameter(L, "L");

  double rho1{1}, rho2{1};
  double lambda1{1}, lambda2{1};
  double omega = 10;
  gmshFem.userDefinedParameter(rho1, "rho1");
  gmshFem.userDefinedParameter(rho2, "rho2");
  gmshFem.userDefinedParameter(lambda1, "lambda1");
  gmshFem.userDefinedParameter(lambda2, "lambda2");
  gmshFem.userDefinedParameter(omega, "omega");

  setupMesh(lc, L, waterDepth, groundDepth);

  /****
   * GmshFem part
   *****/

  //*****
  // Problem declaration
  //*****

  // I. Allocate the formulation object.
  Formulation< Scalar > formulation("poissonDirichlet");

  // II. Define physical regions (dim, tag).
  Domain omegaAcoustic("omegaAcoustic");
  Domain omegaElastic("omegaElastic");
  Domain top("top");
  Domain interface("interface");
  Domain bottom("bottom");

  // III. Allocate field object.

  // Acoustic displacement potential (u = grad(phi))
  Field< Scalar, Form::Form0 > phiA("phi", omegaAcoustic | interface, FunctionSpaceTypeForm0::HierarchicalH1, 4);
  // Elastic displacement potential
  Field< Scalar, Form::Form0 > uy("uy", omegaElastic | interface | bottom, FunctionSpaceTypeForm0::HierarchicalH1, 4);

  std::complex< double > incident(5, 0);
  std::complex< double > im(0, 1);

  double k1 = sqrt(rho1 / lambda1 * omega * omega);
  double k2 = sqrt(rho2 / lambda2 * omega * omega);

  Scalar Z1 = sqrt(rho1 * lambda1);
  Scalar Z2 = sqrt(rho2 * lambda2);
  Scalar expectedReflection = (Z2 - Z1) / (Z2 + Z1);
  Scalar expectedTransmission = im * k1 * (1.0 - expectedReflection);
  gmshfem::msg::print << "Expected reflection coefficient : " << expectedReflection << gmshfem::msg::endl;



  // IV. Write the corresponding weak formulation terms by terms.

  // Convention i*omega*t for the time derivative

  // Acoustic region
  {
    formulation.integral(-grad(dof(phiA)), grad(tf(phiA)), omegaAcoustic, "Gauss6");
    formulation.integral(k1 * k1 * dof(phiA), tf(phiA), omegaAcoustic, "Gauss8");

    // Force incident (negative y) wave (inhomogenuous Robin BC) -> remove exp(-iky)
    formulation.integral((-im * k1) * dof(phiA), tf(phiA), top, "Gauss8");
    formulation.integral(2. * im * k1 * incident * exp(im * k1 * y< Scalar >()), tf(phiA), top, "Gauss8"); // Line is y = 0

    // Continuity (grad(phi).n = u.n with elastic displacement)
    // In  weak form, dphi/dn * tf(phi) occurs
    // Negatif sign because d/dn = -d/dy, but uy = d(phi)/dn
    formulation.integral(-dof(uy), tf(phiA), interface, "Gauss8");
  }

  // Elastic region
  if(true) {
    formulation.integral(-lambda2 * grad(dof(uy)), grad(tf(uy)), omegaElastic, "Gauss6");
    formulation.integral(rho2 * omega * omega * dof(uy), tf(uy), omegaElastic, "Gauss8");

    // Continuity of stress (w * sigma * n term with sigma * n the traction and w the TF)
    formulation.integral(-rho1 * omega * omega * dof(phiA), tf(uy), interface, "Gauss8");
    // Non-reflecting BC
    formulation.integral((-im * k2 * lambda2) * dof(uy), tf(uy), bottom, "Gauss8");
  }

  // Prepro
  formulation.pre();

  // Pro
  formulation.assemble();
  formulation.solve();

  // V. Define and run post-processing operations.
  save(phiA, omegaAcoustic, "phiA");
  save(uy, omegaElastic, "uy");

  auto expected = incident * (exp(im * k1 * y< Scalar >()) + expectedReflection * exp(-im * k1 * y< Scalar >()));

  Scalar integralErrorPhia = integrate(norm(expected - phiA), omegaAcoustic, "Gauss15");
  Scalar integralErrorUy = integrate(norm(incident * expectedTransmission * exp(im * k2 * y< Scalar >()) - uy), omegaElastic, "Gauss15");

  if (real(integralErrorPhia) < 1e-4 && real(integralErrorUy) < 1e-4)
    return 0;
  else
    return 1;
}

static void setupMesh(double lc, double L, double waterDepth, double groundDepth)
{
  /****
   * Gmsh part
   ****/

  gmsh::model::add("HelmholtzNavier");

  /*
  Geometry:
  1 -- 2
  |    |
  3 -- 4
  |    |
  5 -- 6
  
  */
  gmsh::model::geo::addPoint(0., waterDepth, 0., lc, 1);
  gmsh::model::geo::addPoint(L, waterDepth, 0., lc, 2);
  gmsh::model::geo::addPoint(0., 0., 0., lc, 3);
  gmsh::model::geo::addPoint(L, 0., 0., lc, 4);
  gmsh::model::geo::addPoint(0., -groundDepth, 0., lc, 5);
  gmsh::model::geo::addPoint(L, -groundDepth, 0., lc, 6);

  int l12 = gmsh::model::geo::addLine(1, 2);
  int l34 = gmsh::model::geo::addLine(3, 4);
  int l56 = gmsh::model::geo::addLine(5, 6);
  int l13 = gmsh::model::geo::addLine(1, 3);
  int l35 = gmsh::model::geo::addLine(3, 5);
  int l24 = gmsh::model::geo::addLine(2, 4);
  int l46 = gmsh::model::geo::addLine(4, 6);

  int waterLoop = gmsh::model::geo::addCurveLoop({l13, l34, -l24, -l12});
  int groundLoop = gmsh::model::geo::addCurveLoop({l35, l56, -l46, -l34});
  int waterSurface = gmsh::model::geo::addPlaneSurface({waterLoop});
  int groundSurface = gmsh::model::geo::addPlaneSurface({groundLoop});

  gmsh::model::geo::synchronize();

  gmsh::model::addPhysicalGroup(2, {waterSurface}, 1);
  gmsh::model::setPhysicalName(2, 1, "omegaAcoustic");

  gmsh::model::addPhysicalGroup(2, {groundSurface}, 2);
  gmsh::model::setPhysicalName(2, 2, "omegaElastic");

  gmsh::model::addPhysicalGroup(1, {l12}, 1);
  gmsh::model::setPhysicalName(1, 1, "top");

  gmsh::model::addPhysicalGroup(1, {l34}, 2);
  gmsh::model::setPhysicalName(1, 2, "interface");

  gmsh::model::addPhysicalGroup(1, {l56}, 3);
  gmsh::model::setPhysicalName(1, 3, "bottom");

  gmsh::model::mesh::generate();
}