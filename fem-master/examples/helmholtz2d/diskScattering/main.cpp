#include <gmsh.h>
#include <gmshfem/AnalyticalFunction.h>
#include <gmshfem/CSVio.h>
#include <gmshfem/Formulation.h>
#include <gmshfem/Function.h>
#include <gmshfem/GmshFem.h>
#include <gmshfem/Message.h>

using namespace gmshfem;
using namespace gmshfem::common;
using namespace gmshfem::problem;
using namespace gmshfem::domain;
using namespace gmshfem::field;
using namespace gmshfem::function;
using namespace gmshfem::analytics;
using namespace gmshfem::post;
using namespace gmshfem::equation;

struct Edge {
  Domain gamma;
  Domain corner[2];
};

struct HelmholtzDomain {
  Domain omega;
  Domain gammaInf;
  Domain gammaDir;
  Domain omegaPml;
  Domain corners;
  std::vector< Domain > corner;
  std::vector< Edge > edge;
};

//***************************************
// GEOMETRY
//***************************************

void meshCircle(GmshFem &gmshFem, HelmholtzDomain &domains, const double R, const double RInf, const double RPml, bool withPml)
{
  double lc = 0.05;
  gmshFem.userDefinedParameter(lc, "lc");
  int order = 1;
  gmshFem.userDefinedParameter(order, "order");

  gmsh::model::add("geometry");

  gmsh::model::geo::addPoint(0., 0., 0., lc, 1); // Center

  // small circle
  gmsh::model::geo::addPoint(R, 0., 0., lc, 2);
  gmsh::model::geo::addPoint(0, R, 0., lc, 3);
  gmsh::model::geo::addPoint(-R, 0., 0., lc, 4);
  gmsh::model::geo::addPoint(0, -R, 0., lc, 5);

  gmsh::model::geo::addCircleArc(2, 1, 3, 1);
  gmsh::model::geo::addCircleArc(3, 1, 4, 2);
  gmsh::model::geo::addCircleArc(4, 1, 5, 3);
  gmsh::model::geo::addCircleArc(5, 1, 2, 4);

  gmsh::model::geo::addCurveLoop({1, 2, 3, 4}, 1);

  // big circle
  gmsh::model::geo::addPoint(RInf, 0., 0., lc, 6);
  gmsh::model::geo::addPoint(0, RInf, 0., lc, 7);
  gmsh::model::geo::addPoint(-RInf, 0., 0., lc, 8);
  gmsh::model::geo::addPoint(0, -RInf, 0., lc, 9);

  gmsh::model::geo::addCircleArc(6, 1, 7, 5);
  gmsh::model::geo::addCircleArc(7, 1, 8, 6);
  gmsh::model::geo::addCircleArc(8, 1, 9, 7);
  gmsh::model::geo::addCircleArc(9, 1, 6, 8);

  gmsh::model::geo::addCurveLoop({5, 6, 7, 8}, 2);

  // surface
  gmsh::model::geo::addPlaneSurface({2, 1}, 3);

  if(withPml) {
    gmsh::model::geo::addPoint(RPml, 0., 0., lc, 10);
    gmsh::model::geo::addPoint(0, RPml, 0., lc, 11);
    gmsh::model::geo::addPoint(-RPml, 0., 0., lc, 12);
    gmsh::model::geo::addPoint(0, -RPml, 0., lc, 13);

    gmsh::model::geo::addCircleArc(10, 1, 11, 9);
    gmsh::model::geo::addCircleArc(11, 1, 12, 10);
    gmsh::model::geo::addCircleArc(12, 1, 13, 11);
    gmsh::model::geo::addCircleArc(13, 1, 10, 12);

    gmsh::model::geo::addCurveLoop({9, 10, 11, 12}, 3);
    gmsh::model::geo::addPlaneSurface({3, 2}, 4);
  }

  gmsh::model::geo::synchronize();

  // physicals
  gmsh::model::addPhysicalGroup(2, {3}, 3);
  gmsh::model::setPhysicalName(2, 3, "omega");
  gmsh::model::addPhysicalGroup(1, {1, 2, 3, 4}, 1);
  gmsh::model::setPhysicalName(1, 1, "gammaDir");
  gmsh::model::addPhysicalGroup(1, {5, 6, 7, 8}, 2);
  gmsh::model::setPhysicalName(1, 2, "gammaInf");

  if(withPml) {
    gmsh::model::addPhysicalGroup(2, {4}, 4);
    gmsh::model::setPhysicalName(2, 4, "pml");
  }

  gmsh::model::geo::synchronize();
  gmsh::model::mesh::generate();

  gmsh::model::mesh::setOrder(order);

  // Domains

  domains.omega = Domain(2, 3);
  domains.gammaInf = Domain(1, 2);
  domains.gammaDir = Domain(1, 1);
  if(withPml) {
    domains.omegaPml = Domain(2, 4);
  }
}

void meshSquare(GmshFem &gmshFem, HelmholtzDomain &domains, const double R, const double RInf, const double RPml, bool withPml)
{
  double lc = 0.05;
  gmshFem.userDefinedParameter(lc, "lc");
  int order = 1;
  gmshFem.userDefinedParameter(order, "order");

  gmsh::model::add("geometry");

  gmsh::model::geo::addPoint(0., 0., 0., lc, 1); // Center

  // small circle
  gmsh::model::geo::addPoint(R, 0., 0., lc, 2);
  gmsh::model::geo::addPoint(0, R, 0., lc, 3);
  gmsh::model::geo::addPoint(-R, 0., 0., lc, 4);
  gmsh::model::geo::addPoint(0, -R, 0., lc, 5);

  gmsh::model::geo::addCircleArc(2, 1, 3, 1);
  gmsh::model::geo::addCircleArc(3, 1, 4, 2);
  gmsh::model::geo::addCircleArc(4, 1, 5, 3);
  gmsh::model::geo::addCircleArc(5, 1, 2, 4);

  gmsh::model::geo::addCurveLoop({1, 2, 3, 4}, 1);

  // polygone
  gmsh::model::geo::addPoint(RInf, RInf, 0., lc, 6);
  gmsh::model::geo::addPoint(-RInf, RInf, 0., lc, 7);
  gmsh::model::geo::addPoint(-RInf, -RInf, 0., lc, 8);
  gmsh::model::geo::addPoint(RInf, -RInf, 0., lc, 9);

  gmsh::model::geo::addLine(6, 7, 5);
  gmsh::model::geo::addLine(7, 8, 6);
  gmsh::model::geo::addLine(8, 9, 7);
  gmsh::model::geo::addLine(9, 6, 8);

  gmsh::model::geo::addCurveLoop({5, 6, 7, 8}, 2);

  gmsh::model::geo::addPlaneSurface({2, 1}, 3);

  if(withPml) {
    msg::warning << "PML not yet implemented for square domain." << msg::endl;
  }

  gmsh::model::geo::synchronize();

  // physicals
  gmsh::model::addPhysicalGroup(2, {3}, 3);
  gmsh::model::setPhysicalName(2, 3, "omega");
  gmsh::model::addPhysicalGroup(1, {1, 2, 3, 4}, 1);
  gmsh::model::setPhysicalName(1, 1, "gammaDir");

  gmsh::model::addPhysicalGroup(1, {5}, 2);
  gmsh::model::setPhysicalName(1, 2, "gammaInf_top");
  gmsh::model::addPhysicalGroup(1, {6}, 3);
  gmsh::model::setPhysicalName(1, 3, "gammaInf_left");
  gmsh::model::addPhysicalGroup(1, {7}, 4);
  gmsh::model::setPhysicalName(1, 4, "gammaInf_bottom");
  gmsh::model::addPhysicalGroup(1, {8}, 5);
  gmsh::model::setPhysicalName(1, 5, "gammaInf_right");

  gmsh::model::addPhysicalGroup(0, {6}, 1);
  gmsh::model::setPhysicalName(0, 1, "corner_north_east");
  gmsh::model::addPhysicalGroup(0, {7}, 2);
  gmsh::model::setPhysicalName(0, 2, "corner_north_west");
  gmsh::model::addPhysicalGroup(0, {8}, 3);
  gmsh::model::setPhysicalName(0, 3, "corner_south_west");
  gmsh::model::addPhysicalGroup(0, {9}, 4);
  gmsh::model::setPhysicalName(0, 4, "corner_south_east");

  if(withPml) {
  }

  gmsh::model::geo::synchronize();
  gmsh::model::mesh::generate();

  gmsh::model::mesh::setOrder(order);

  // Domains

  domains.omega = Domain(2, 3);
  domains.gammaDir = Domain(1, 1);
  if(withPml) {
    domains.omegaPml = Domain(2, 4);
  }
  domains.edge.resize(4);
  domains.corner.resize(4);
  for(int i = 0; i < 4; ++i) {
    domains.edge[i].gamma = Domain(1, 2 + i);
    domains.edge[i].corner[0] = Domain(0, i + 1);
    domains.edge[i].corner[1] = Domain(0, (i + 1) % 4 + 1);

    domains.corner[i] = Domain(0, i + 1);

    domains.corners |= domains.edge[i].corner[0];
    domains.gammaInf |= domains.edge[i].gamma;
  }
}

//***************************************
// FORMULATION
//***************************************

int main(int argc, char **argv)
{
  GmshFem gmshFem(argc, argv);

  std::string problem = "hard"; // hard
  gmshFem.userDefinedParameter(problem, "problem");
  std::string gauss = "Gauss4";
  gmshFem.userDefinedParameter(gauss, "gauss");

  double pi = 3.14159265359;
  double k = 2. * pi;
  gmshFem.userDefinedParameter(k, "k");
  double theta = 0.;
  gmshFem.userDefinedParameter(theta, "theta");

  std::string abcName = "sommerfeld";
  gmshFem.userDefinedParameter(abcName, "abcName");

  bool saveUExact = false;
  gmshFem.userDefinedParameter(saveUExact, "saveUExact");
  bool saveError = false;
  gmshFem.userDefinedParameter(saveError, "saveError");
  bool computeError = false;
  gmshFem.userDefinedParameter(computeError, "computeError");

  int geometry = 0;
  gmshFem.userDefinedParameter(geometry, "geometry");

  double lc = 0.05;
  gmshFem.userDefinedParameter(lc, "lc");
  int order = 1;
  gmshFem.userDefinedParameter(order, "order");

  double R = 1.;
  double RInf = 5.;
  double RPml = 6.;

  HelmholtzDomain domains;

  switch(geometry) {
  case 0:
    meshCircle(gmshFem, domains, R, RInf, RPml, abcName == "pml");
    break;
  case 1:
    meshSquare(gmshFem, domains, R, RInf, RPml, abcName == "pml");
    break;
  case 2:
    RInf = 1.1;
    meshSquare(gmshFem, domains, R, RInf, RPml, abcName == "pml");
    break;
  default:
    break;
  }

  msg::info << "********************************" << msg::endl;
  msg::info << " Wave properties:" << msg::endl;
  msg::info << " - wave number = " << k << " [1/m]." << msg::endl;
  msg::info << " - wavelength = " << 2 * pi / k << " [m]." << msg::endl;
  msg::info << " Problem properties:" << msg::endl;
  msg::info << " - cylinder scattering by a plan wave." << msg::endl;
  msg::info << " - using " << abcName << " as boundary condition." << msg::endl;
  msg::info << " Mesh properties:" << msg::endl;
  unsigned int pointsByWl = (2 * pi / k) / lc;
  msg::info << " - Approximate number of points by wavelength = " << pointsByWl << "" << msg::endl;
  msg::info << " - Order = " << order << "" << msg::endl;
  msg::info << "********************************" << msg::endl;

  std::vector< FieldInterface< std::complex< double > > * > fieldBucket;

  //*****
  // Problem declaration
  //*****

  Formulation< std::complex< double > > formulation("2D-Helmholtz");

  std::complex< double > im = std::complex< double >(0., 1.);

  Field< std::complex< double >, Form::Form0 > u("u", domains.omega | domains.omegaPml | domains.gammaInf | domains.gammaDir | domains.corners, FunctionSpaceTypeForm0::Lagrange);

  Function< std::complex< double >, Degree::Degree0 > *solution = nullptr;
  if(problem == "soft") {
    solution = new AnalyticalFunction< helmholtz2D::ScatteringByASoftCylinder< std::complex< double > > >(k, R, 0., 0., 2 * k, theta);
    u.addConstraint(domains.gammaDir, *solution);
  }
  else if(problem == "hard") {
    solution = new AnalyticalFunction< helmholtz2D::ScatteringByAHardCylinder< std::complex< double > > >(k, R, 0., 0., 2 * k, theta);
    formulation.integral(AnalyticalFunction< helmholtz2D::dr_ScatteringByAHardCylinder< std::complex< double > > >(k, R, 0., 0., 2 * k, theta), tf(u), domains.gammaDir, gauss);
  }
  else {
    msg::error << "Unknown problem: " << problem << msg::endl;
  }

  formulation.integral(grad(dof(u)), grad(tf(u)), domains.omega, gauss);
  formulation.integral(-k * k * dof(u), tf(u), domains.omega, gauss);

  if(abcName == "sommerfeld") {
    msg::info << "Use Sommerfeld ABC." << msg::endl;
    formulation.integral(-im * k * dof(u), tf(u), domains.gammaInf, gauss);
  }
  else if(abcName == "bayliss-turkel") {
    if(geometry == 0) {
      msg::info << "Use Bayliss-Turkel ABC." << msg::endl;
      const std::complex< double > alpha = 1. / (2. * RInf) - im / (8. * k * RInf * RInf * (1. + im / (k * RInf)));
      const std::complex< double > beta = -1. / (2. * im * k * (1. + im / (k * RInf)));

      formulation.integral(-im * k * dof(u), tf(u), domains.gammaInf, gauss);
      formulation.integral(alpha * dof(u), tf(u), domains.gammaInf, gauss);
      formulation.integral(beta * grad(dof(u)), grad(tf(u)), domains.gammaInf, gauss);
    }
    else {
      msg::warning << "Bayliss-Turkel can only be appied on cylindrical border." << msg::endl;
    }
  }
  else if(abcName == "pade") {
    int padeOrder = 2;
    gmshFem.userDefinedParameter(padeOrder, "padeOrder");

    msg::info << "Use order " << padeOrder << " Pade ABC." << msg::endl;

    const double angle = pi / 3.;
    const std::complex< double > exp1 = std::complex< double >(std::cos(angle), std::sin(angle));
    const std::complex< double > exp2 = std::complex< double >(std::cos(angle / 2.), std::sin(angle / 2.));
    const std::complex< double > expM = std::complex< double >(std::cos(-angle), std::sin(-angle));
    const double M = 2. * padeOrder + 1.;
    const std::complex< double > coef = 2. / M;
    std::vector< std::complex< double > > c(padeOrder, 0.);
    for(int i = 0; i < padeOrder; ++i) {
      c[i] = std::tan((i + 1) * pi / M);
      c[i] *= c[i];
    }

    if(geometry == 0) {
      // Auxilary fields
      std::vector< Field< std::complex< double >, Form::Form0 > * > phi;
      for(int i = 0; i < padeOrder; ++i) {
        phi.push_back(new Field< std::complex< double >, Form::Form0 >("phi_" + std::to_string(i), domains.gammaInf, FunctionSpaceTypeForm0::Lagrange));
        fieldBucket.push_back(phi.back());
      }

      formulation.integral(-im * k * exp2 * dof(u), tf(u), domains.gammaInf, gauss);
      for(int i = 0; i < padeOrder; ++i) {
        // Boundary conditions
        formulation.integral(-im * k * exp2 * coef * c[i] * dof(u), tf(u), domains.gammaInf, gauss);
        formulation.integral(-im * k * exp2 * coef * c[i] * dof(*phi[i]), tf(u), domains.gammaInf, gauss);

        // Auxiliary equations
        formulation.integral(grad(dof(*phi[i])), grad(tf(*phi[i])), domains.gammaInf, gauss);
        formulation.integral(-k * k * (exp1 * c[i] + 1.) * dof(*phi[i]), tf(*phi[i]), domains.gammaInf, gauss);
        formulation.integral(-k * k * exp1 * (c[i] + 1.) * dof(u), tf(*phi[i]), domains.gammaInf, gauss);
      }
    }
    else if(geometry == 1 || geometry == 2 || geometry == 3) {
      //
      //  (1)****** 0 ******(0)
      //   *                 *
      //   *                 *
      //   *                 *
      //   1                 3
      //   *                 *
      //   *                 *
      //   *                 *
      //  (2)****** 2 ******(3)
      //

      // Auxilary fields
      std::vector< Field< std::complex< double >, Form::Form0 > * > phi[4];
      for(int x = 0; x < 4; ++x) {
        for(int i = 0; i < padeOrder; ++i) {
          phi[x].push_back(new Field< std::complex< double >, Form::Form0 >("phi_" + std::to_string(i) + "^" + std::to_string(x), domains.edge[x].gamma | domains.edge[x].corner[0] | domains.edge[x].corner[1], FunctionSpaceTypeForm0::Lagrange));
          fieldBucket.push_back(phi[x].back());
        }
      }

      formulation.integral(-im * k * exp2 * dof(u), tf(u), domains.gammaInf, gauss);
      for(int x = 0; x < 4; ++x) {
        for(int i = 0; i < padeOrder; ++i) {
          // Boundary conditions
          formulation.integral(-im * k * exp2 * coef * c[i] * dof(u), tf(u), domains.edge[x].gamma, gauss);
          formulation.integral(-im * k * exp2 * coef * c[i] * dof(*phi[x][i]), tf(u), domains.edge[x].gamma, gauss);

          // Auxiliary equations
          formulation.integral(grad(dof(*phi[x][i])), grad(tf(*phi[x][i])), domains.edge[x].gamma, gauss);
          formulation.integral(-k * k * (exp1 * c[i] + 1.) * dof(*phi[x][i]), tf(*phi[x][i]), domains.edge[x].gamma, gauss);
          formulation.integral(-k * k * exp1 * (c[i] + 1.) * dof(u), tf(*phi[x][i]), domains.edge[x].gamma, gauss);
        }
      }

      bool withCorner = true;
      gmshFem.userDefinedParameter(withCorner, "withCorner");

      if(withCorner) {
        // Corner fields
        std::vector< std::vector< std::vector< Field< std::complex< double >, Form::Form0 > * > > > corner_phi(4);
        for(int k = 0; k < 4; ++k) {
          corner_phi[k].resize(padeOrder);
          for(int i = 0; i < padeOrder; ++i) {
            for(int j = 0; j < padeOrder; ++j) {
              corner_phi[k][i].push_back(new Field< std::complex< double >, Form::Form0 >("corner_phi_(" + std::to_string(i) + ", " + std::to_string(j) + ")^(" + std::to_string(k) + ")", domains.corner[k], FunctionSpaceTypeForm0::Lagrange));
              fieldBucket.push_back(corner_phi[k][i].back());
            }
          }
        }

        for(int ci = 0; ci < 4; ++ci) {
          int x = ci;
          int y = (ci - 1 < 0 ? 4 - 1 : ci - 1);
          for(int i = 0; i < padeOrder; ++i) {
            formulation.integral(-im * k * exp2 * dof(*phi[x][i]), tf(*phi[x][i]), domains.corner[ci], gauss);
            formulation.integral(-im * k * exp2 * dof(*phi[y][i]), tf(*phi[y][i]), domains.corner[ci], gauss);
            for(int j = 0; j < padeOrder; ++j) {
              // Auxiliary equations for corner
              formulation.integral(dof(*corner_phi[ci][i][j]), tf(*corner_phi[ci][i][j]), domains.corner[ci], gauss);
              formulation.integral((c[j] + 1.) / (c[i] + c[j] + expM) * dof(*phi[x][i]), tf(*corner_phi[ci][i][j]), domains.corner[ci], gauss);
              formulation.integral((c[i] + 1.) / (c[i] + c[j] + expM) * dof(*phi[y][j]), tf(*corner_phi[ci][i][j]), domains.corner[ci], gauss);

              // Corner conditions
              formulation.integral(-im * k * exp2 * coef * c[j] * dof(*phi[x][i]), tf(*phi[x][i]), domains.corner[ci], gauss);
              formulation.integral(-im * k * exp2 * coef * c[j] * dof(*phi[y][i]), tf(*phi[y][i]), domains.corner[ci], gauss);
              formulation.integral(-im * k * exp2 * coef * c[j] * dof(*corner_phi[ci][i][j]), tf(*phi[x][i]), domains.corner[ci], gauss);
              formulation.integral(-im * k * exp2 * coef * c[j] * dof(*corner_phi[ci][j][i]), tf(*phi[y][i]), domains.corner[ci], gauss);
            }
          }
        }
      }
    }
  }
  else if(abcName == "pml") {
    msg::info << "Use a PML." << msg::endl;
    const double WPml = RPml - RInf;
    ScalarFunction< std::complex< double > > cosT = x< std::complex< double > >() / r2d< std::complex< double > >();
    ScalarFunction< std::complex< double > > sinT = y< std::complex< double > >() / r2d< std::complex< double > >();
    ScalarFunction< std::complex< double > > dampingProfileR = 1. / (WPml - (r2d< std::complex< double > >() - RInf));
    ScalarFunction< std::complex< double > > dampingProfileInt = -ln((WPml - (r2d< std::complex< double > >() - RInf)) / WPml);
    ScalarFunction< std::complex< double > > cR = 1. + im * dampingProfileR / k;
    ScalarFunction< std::complex< double > > cStretch = 1. + im * (1. / r2d< std::complex< double > >()) * dampingProfileInt / k;
    ScalarFunction< std::complex< double > > S_PML = cR * cStretch;
    TensorFunction< std::complex< double > > D = tensor< std::complex< double > >(cStretch / cR * cosT * cosT + cR / cStretch * sinT * sinT,
                                                                                  cStretch / cR * cosT * sinT - cR / cStretch * cosT * sinT,
                                                                                  0.,
                                                                                  cStretch / cR * cosT * sinT - cR / cStretch * cosT * sinT,
                                                                                  cStretch / cR * sinT * sinT + cR / cStretch * cosT * cosT,
                                                                                  0.,
                                                                                  0., 0., 0.);

    formulation.integral(D * grad(dof(u)), grad(tf(u)), domains.omegaPml, gauss);
    formulation.integral(-k * k * S_PML * dof(u), tf(u), domains.omegaPml, gauss);
  }

  // Prepro
  formulation.pre();

  // Pro
  formulation.assemble();
  formulation.solve();

  // Postpro
  save(u);
  if(saveUExact) {
    save(*solution, domains.omega, "u_exact");
  }
  if(saveError) {
    save(*solution - u, domains.omega, "e");
  }
  if(computeError) {
    std::complex< double > num = integrate(pow(abs(*solution - u), 2), domains.omega, gauss);
    std::complex< double > den = integrate(pow(abs(*solution), 2), domains.omega, gauss);
    msg::info << "L_2 error = " << sqrt(num / den) << msg::endl;
    CSVio file("convergence", ';', OpeningMode::Append);
    file << order << formulation.getTotalNumberOfDof() << std::sqrt((double)formulation.getTotalNumberOfDof()) << pointsByWl << std::real(sqrt(num / den)) << csv::endl;
  }

  for(unsigned int i = 0; i < fieldBucket.size(); ++i) {
    delete fieldBucket[i];
  }

  delete solution;

  return 0;
}
