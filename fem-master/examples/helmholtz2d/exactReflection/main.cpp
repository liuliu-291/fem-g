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
static void setupMesh(double lc, double L = 1, double waterDepth = 1);

int main(int argc, char **argv)
{
  GmshFem gmshFem(argc, argv);

  // Define the ADD_PARAMETER macro
#define ADD_PARAMETER(name) gmshFem.userDefinedParameter(name, #name)

  // Declare the parameters
  double lc{1e-2};
  double waterDepth{1};
  double groundDepth{1};
  double L{1};
  double k{35};
  double thetaInc{0.3};

  // Add the parameters to gmshFem using the ADD_PARAMETER macro
  ADD_PARAMETER(lc);
  ADD_PARAMETER(waterDepth);
  ADD_PARAMETER(groundDepth);
  ADD_PARAMETER(L);
  ADD_PARAMETER(k);
  ADD_PARAMETER(thetaInc);


  setupMesh(lc, L, waterDepth);

  Formulation< Scalar > formulation("Helmholtz");

  Domain omegaAcoustic("omegaAcoustic");
  Domain top("top");
  Domain interface("interface");
  PeriodicLink link("gammaLeft");

  std::complex< double > incident(5, 0);
  std::complex< double > im(0, 1);

  Field< Scalar, Form::Form0 > phiA("phi", omegaAcoustic | interface, FunctionSpaceTypeForm0::Lagrange);

  phiA.addPeriodicConstraint(link, 1.0 * exp(im * k * L * sin(thetaInc)));

  // Convention i*omega*t for the time derivative

  // The analytical form of the incident wave as well as its analytical gradient are needed.
  auto phiInc = incident * exp(-im * k * (x< Scalar >() * sin(thetaInc) - y< Scalar >() * cos(thetaInc)));
  auto dPhiInc_dy = incident * im * k * cos(thetaInc) * exp(-im * k * (x< Scalar >() * sin(thetaInc) - y< Scalar >() * cos(thetaInc)));


  // Acoustic region
  {
    formulation.integral(-grad(dof(phiA)), grad(tf(phiA)), omegaAcoustic, "Gauss6");
    formulation.integral(k * k * dof(phiA), tf(phiA), omegaAcoustic, "Gauss8");

    // Absorbing BC on top
    formulation.integral(-(im * k) * dof(phiA), tf(phiA), top, "Gauss8");
    formulation.integral(im / (2. * k) * grad(dof(phiA)), grad(tf(phiA)), top, "Gauss8");

    // Neumann : du_scat/dn = - du_inc / dn
    // The total derivative is zero but we compute the scattered wave only.
    formulation.integral(yComp(normal< Scalar >()) * dPhiInc_dy, tf(phiA), interface, "Gauss8");
  }

  // Prepro
  formulation.pre();

  // Pro
  formulation.assemble();
  formulation.solve();

  // V. Define and run post-processing operations.
  save(phiA, omegaAcoustic, "phiA");
  auto expected = incident * exp(-im * k * (x< Scalar >() * sin(thetaInc) + y< Scalar >() * cos(thetaInc)));
  auto err = phiA - expected;
  save(expected, omegaAcoustic, "phiA_ref");
  save(err, omegaAcoustic, "phiA_err");
  Scalar integralErrorPhia = integrate(norm(err), omegaAcoustic, "Gauss15");
  Scalar integralPhia = integrate(norm(phiA), omegaAcoustic, "Gauss15");
  gmshfem::msg::print << "Relative error : " << integralErrorPhia / integralPhia << gmshfem::msg::endl;

  auto normalDir = integrate(yComp(normal< Scalar >()), interface, "Gauss1");

  const double tol = 1e-2;
  return real(integralErrorPhia / integralPhia) < tol ? 0 : 1;
}

static void setupMesh(double lc, double L, double waterDepth)
{
  gmsh::model::add("Helmholtz");


  // Create a new 2D geometry
  std::vector< int > points = {
    gmsh::model::geo::addPoint(0, 0, 0, lc, 1),
    gmsh::model::geo::addPoint(L, 0, 0, lc, 2),
    gmsh::model::geo::addPoint(L, waterDepth, 0, lc, 3),
    gmsh::model::geo::addPoint(0, waterDepth, 0, lc, 4),
  };
  std::vector< int > lines = {
    gmsh::model::geo::addLine(points[1], points[0]),
    gmsh::model::geo::addLine(points[2], points[1]),
    gmsh::model::geo::addLine(points[3], points[2]),
    gmsh::model::geo::addLine(points[0], points[3]),
  };

  int line_loop = gmsh::model::geo::addCurveLoop({lines[0], lines[1], lines[2], lines[3]}, 7);
  int waterSurface = gmsh::model::geo::addPlaneSurface({line_loop});


  gmsh::model::geo::synchronize();

  std::vector< double > translation({1, 0, 0, -L,
                                     0, 1, 0, 0,
                                     0, 0, 1, 0,
                                     0, 0, 0, 1});
  gmsh::model::mesh::setPeriodic(1, {lines[3]}, {lines[1]}, translation);

  // Create a physical surface from the line loop and assign it the name "omegaAcoustic"
  gmsh::model::addPhysicalGroup(2, {waterSurface}, 1);
  gmsh::model::setPhysicalName(2, 1, "omegaAcoustic");

  gmsh::model::addPhysicalGroup(1, {lines[0]}, 2);
  gmsh::model::setPhysicalName(1, 2, "interface");

  gmsh::model::addPhysicalGroup(1, {lines[2]}, 3);
  gmsh::model::setPhysicalName(1, 3, "top");

  gmsh::model::addPhysicalGroup(1, {lines[3]}, 5);
  gmsh::model::setPhysicalName(1, 5, "gammaLeft");

  gmsh::model::mesh::generate();
  gmsh::model::mesh::setOrder(2);
}