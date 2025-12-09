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

struct Edge {
  Domain gamma;
  Domain corner[2];
};

struct HelmholtzDomain {
  Domain omega;
  Domain gammaInf;
  Domain gammaPml;
  Domain omegaPml;
  Domain gammaDir;
  Domain source;
  Domain corners;
  std::vector< Domain > corner;
  std::vector< Edge > edge;
};

//***************************************
// GEOMETRY
//***************************************

void meshCircle(GmshFem &gmshFem, HelmholtzDomain &domains, const double R, const double Rpml, const double Ro, const double lc, bool withPml)
{
  gmsh::model::add("geometry");

  gmsh::model::geo::addPoint(0., 0., 0., lc, 1); // Center
  gmsh::model::geo::addPoint(-Ro - 0.4, 0., 0., lc / 2., 100); // Source point
  // physical circle
  gmsh::model::geo::addPoint(R, 0., 0., lc, 2);
  gmsh::model::geo::addPoint(0, R, 0., lc, 3);
  gmsh::model::geo::addPoint(-R, 0., 0., lc, 4);
  gmsh::model::geo::addPoint(0, -R, 0., lc, 5);

  gmsh::model::geo::addCircleArc(2, 1, 3, 1);
  gmsh::model::geo::addCircleArc(3, 1, 4, 2);
  gmsh::model::geo::addCircleArc(4, 1, 5, 3);
  gmsh::model::geo::addCircleArc(5, 1, 2, 4);

  // scattered disk
  gmsh::model::geo::addPoint(Ro, 0., 0., lc / 2., 22);
  gmsh::model::geo::addPoint(0, Ro, 0., lc / 2., 33);
  gmsh::model::geo::addPoint(-Ro, 0., 0., lc / 2., 44);
  gmsh::model::geo::addPoint(0, -Ro, 0., lc / 2., 55);

  gmsh::model::geo::addCircleArc(22, 1, 33, 11);
  gmsh::model::geo::addCircleArc(33, 1, 44, 22);
  gmsh::model::geo::addCircleArc(44, 1, 55, 33);
  gmsh::model::geo::addCircleArc(55, 1, 22, 44);

  // surface
  gmsh::model::geo::addCurveLoop({1, 2, 3, 4}, 1);
  gmsh::model::geo::addCurveLoop({11, 22, 33, 44}, 11);
  gmsh::model::geo::addPlaneSurface({1, 11}, 1);

  if(withPml) {
    gmsh::model::geo::addPoint(Rpml, 0., 0., lc, 6);
    gmsh::model::geo::addPoint(0, Rpml, 0., lc, 7);
    gmsh::model::geo::addPoint(-Rpml, 0., 0., lc, 8);
    gmsh::model::geo::addPoint(0, -Rpml, 0., lc, 9);

    gmsh::model::geo::addCircleArc(6, 1, 7, 5);
    gmsh::model::geo::addCircleArc(7, 1, 8, 6);
    gmsh::model::geo::addCircleArc(8, 1, 9, 7);
    gmsh::model::geo::addCircleArc(9, 1, 6, 8);
    gmsh::model::geo::addCurveLoop({5, 6, 7, 8}, 2);
    gmsh::model::geo::addPlaneSurface({1, 2}, 2);
  }
  gmsh::model::geo::synchronize();
  gmsh::model::mesh::embed(0, {100}, 2, 1);

  // physicals
  gmsh::model::addPhysicalGroup(2, {1}, 1);
  gmsh::model::setPhysicalName(2, 1, "omega");
  gmsh::model::addPhysicalGroup(1, {1, 2, 3, 4}, 3);
  gmsh::model::setPhysicalName(1, 3, "gammaInf");
  gmsh::model::addPhysicalGroup(1, {11, 22, 33, 44}, 55);
  gmsh::model::setPhysicalName(1, 55, "gammaDir");
  gmsh::model::addPhysicalGroup(0, {100}, 100);
  gmsh::model::setPhysicalName(0, 100, "source");
  if(withPml) {
    gmsh::model::addPhysicalGroup(2, {2}, 2);
    gmsh::model::setPhysicalName(2, 2, "omegaPml");
    gmsh::model::addPhysicalGroup(1, {5, 6, 7, 8}, 4);
    gmsh::model::setPhysicalName(1, 4, "gammaPml");
  }

  // generate mesh
  gmsh::model::geo::synchronize();
  gmsh::model::mesh::generate(2);
  // gmsh::model::mesh::setOrder(2); // T6 elem
  // gmsh::write("m.msh");

  // Domains
  domains.omega = Domain(2, 1);
  domains.gammaInf = Domain(1, 3);
  domains.gammaDir = Domain(1, 55);
  domains.source = Domain(0, 100);

  if(withPml) {
    domains.omegaPml = Domain(2, 2);
    domains.gammaPml = Domain(1, 4);
  }
}

void meshSquare(GmshFem &gmshFem, HelmholtzDomain &domains, const double L, const double Lpml, const double Ro, const double lc, bool withPml)
{
  gmsh::model::add("geometry");

  gmsh::model::geo::addPoint(0., 0., 0., lc, 1); // Center
  gmsh::model::geo::addPoint(-Ro - 0.4, 0., 0., lc / 2., 100); // Source point
  // polygone
  gmsh::model::geo::addPoint(L, L, 0., lc, 2);
  gmsh::model::geo::addPoint(-L, L, 0., lc, 3);
  gmsh::model::geo::addPoint(-L, -L, 0., lc, 4);
  gmsh::model::geo::addPoint(L, -L, 0., lc, 5);

  gmsh::model::geo::addLine(2, 3, 1);
  gmsh::model::geo::addLine(3, 4, 2);
  gmsh::model::geo::addLine(4, 5, 3);
  gmsh::model::geo::addLine(5, 2, 4);

  // scattered disk
  gmsh::model::geo::addPoint(Ro, 0., 0., lc / 2., 22);
  gmsh::model::geo::addPoint(0, Ro, 0., lc / 2., 33);
  gmsh::model::geo::addPoint(-Ro, 0., 0., lc / 2., 44);
  gmsh::model::geo::addPoint(0, -Ro, 0., lc / 2., 55);

  gmsh::model::geo::addCircleArc(22, 1, 33, 11);
  gmsh::model::geo::addCircleArc(33, 1, 44, 22);
  gmsh::model::geo::addCircleArc(44, 1, 55, 33);
  gmsh::model::geo::addCircleArc(55, 1, 22, 44);

  // surface
  gmsh::model::geo::addCurveLoop({1, 2, 3, 4}, 1);
  gmsh::model::geo::addCurveLoop({11, 22, 33, 44}, 11);
  gmsh::model::geo::addPlaneSurface({11, 1}, 1);

  if(withPml) {
    gmsh::model::geo::addPoint(Lpml, Lpml, 0., lc, 6);
    gmsh::model::geo::addPoint(-Lpml, Lpml, 0., lc, 7);
    gmsh::model::geo::addPoint(-Lpml, -Lpml, 0., lc, 8);
    gmsh::model::geo::addPoint(Lpml, -Lpml, 0., lc, 9);

    gmsh::model::geo::addLine(6, 7, 5);
    gmsh::model::geo::addLine(7, 8, 6);
    gmsh::model::geo::addLine(8, 9, 7);
    gmsh::model::geo::addLine(9, 6, 8);

    gmsh::model::geo::addCurveLoop({5, 6, 7, 8}, 2);
    gmsh::model::geo::addPlaneSurface({1, 2}, 2);
  }
  gmsh::model::geo::synchronize();
  gmsh::model::mesh::embed(0, {100}, 2, 1);

  // physicals
  gmsh::model::addPhysicalGroup(2, {1}, 1);
  gmsh::model::setPhysicalName(2, 1, "omega");
  gmsh::model::addPhysicalGroup(1, {1, 2, 3, 4}, 3);
  gmsh::model::setPhysicalName(1, 3, "gammaInf");
  gmsh::model::addPhysicalGroup(1, {11, 22, 33, 44}, 55);
  gmsh::model::setPhysicalName(1, 55, "gammaDir");
  gmsh::model::addPhysicalGroup(0, {100}, 100);
  gmsh::model::setPhysicalName(0, 100, "source");

  // Physical edges
  gmsh::model::addPhysicalGroup(1, {1}, 20);
  gmsh::model::setPhysicalName(1, 20, "gammaInf_top");
  gmsh::model::addPhysicalGroup(1, {2}, 30);
  gmsh::model::setPhysicalName(1, 30, "gammaInf_left");
  gmsh::model::addPhysicalGroup(1, {3}, 40);
  gmsh::model::setPhysicalName(1, 40, "gammaInf_bottom");
  gmsh::model::addPhysicalGroup(1, {4}, 50);
  gmsh::model::setPhysicalName(1, 50, "gammaInf_right");

  // Physical domain corners
  gmsh::model::addPhysicalGroup(0, {2}, 2);
  gmsh::model::setPhysicalName(0, 2, "UpRightCorner");
  gmsh::model::addPhysicalGroup(0, {3}, 3);
  gmsh::model::setPhysicalName(0, 3, "UpLeftCorner");
  gmsh::model::addPhysicalGroup(0, {4}, 4);
  gmsh::model::setPhysicalName(0, 4, "BottomLeftCorner");
  gmsh::model::addPhysicalGroup(0, {5}, 5);
  gmsh::model::setPhysicalName(0, 5, "BottomRightCorner");
  if(withPml) {
    gmsh::model::addPhysicalGroup(2, {2}, 2);
    gmsh::model::setPhysicalName(2, 2, "omegaPml");
    gmsh::model::addPhysicalGroup(1, {5, 6, 7, 8}, 4);
    gmsh::model::setPhysicalName(1, 4, "gammaPml");
  }

  // generate mesh
  gmsh::model::geo::synchronize();
  gmsh::model::mesh::generate(2);
  // gmsh::model::mesh::setOrder(2); // T6 elem
  // gmsh::write("m.msh");

  // Domains
  domains.omega = Domain(2, 1);
  domains.gammaInf = Domain(1, 3);
  domains.gammaDir = Domain(1, 55);
  domains.source = Domain(0, 100);

  if(withPml) {
    domains.omegaPml = Domain(2, 2);
    domains.gammaPml = Domain(1, 4);
  }
  domains.edge.resize(4);
  domains.corner.resize(4);
  for(int i = 0; i < 4; ++i) {
    domains.edge[i].gamma = Domain(1, 10 * (i + 2));
    domains.edge[i].corner[0] = Domain(0, i + 2);
    domains.edge[i].corner[1] = Domain(0, (i + 1) % 4 + 2);

    domains.corner[i] = Domain(0, i + 2);

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

  // geometry parameters
  int geometry = 0; // 0-circle, 1-square
  gmshFem.userDefinedParameter(geometry, "geometry");
  double R = 2.0; // Domain size: circle radius or square of size [-R,R] x [-R,R]
  double Ro = 1.0; // inner disk radius
  gmshFem.userDefinedParameter(R, "R");
  gmshFem.userDefinedParameter(Ro, "Ro");

  // physical parameters
  const double pi = 3.14159265358979323846264338327950288;
  const double k = 25; // free field wavenumber
  const std::complex< double > im(0., 1.);
  double M = -0.4; // Mach number at infinity, |M| < 0.5 to avoid supersonic flow
  gmshFem.userDefinedParameter(M, "M");
  if(std::abs(M) >= 0.5) {
    msg::info << "Supersonic mean flow !" << msg::endl;
    exit(0);
  }

  // numerical parameters
  int order = 4; // FEM shape function order
  gmshFem.userDefinedParameter(order, "order");
  std::string gauss = "Gauss10";
  double lc = 0.03; // mesh size
  gmshFem.userDefinedParameter(lc, "lc");

  // Choose ABC - exterior boundary condition
  std::string abcName = "Pade"; // available choices: ABC-0, ABC-2, Pade, pml
  gmshFem.userDefinedParameter(abcName, "abcName");
  int Npml;
  double Rpml;
  if(abcName == "pml") { // PML parameters
    Npml = 4; // number of layers
    gmshFem.userDefinedParameter(Npml, "Npml");
    Rpml = R + Npml * lc; // extended domain
  }

  HelmholtzDomain domains;

  switch(geometry) {
  case 0:
    meshCircle(gmshFem, domains, R, Rpml, Ro, lc, abcName == "pml");
    break;
  case 1:
    meshSquare(gmshFem, domains, R, Rpml, Ro, lc, abcName == "pml");
    break;
  default:
    msg::error << "geometry not available ! " << msg::endl;
    exit(0);
    break;
  }

  // Potential flow past a cylinder in cartesian coordinates
  ScalarFunction< std::complex< double > > Mx = M + M * (pow(Ro, 2) / pow(r2d< std::complex< double > >(), 2)) * (1 - 2. * pow(x< std::complex< double > >(), 2) / pow(r2d< std::complex< double > >(), 2));
  ScalarFunction< std::complex< double > > My = -2. * M * pow(Ro, 2) * x< std::complex< double > >() * y< std::complex< double > >() / pow(r2d< std::complex< double > >(), 4);

  float pointsByWl = (2 * pi * order * (1 - abs(2. * M))) / (lc * k);
  msg::info << " - FEM basis order = " << order << "" << msg::endl;
  msg::info << " - Smallest number of dofs by wavelength = " << pointsByWl << "" << msg::endl;
  if(pointsByWl < 6) {
    msg::warning << " - Less than 6 points per wavelength ! " << msg::endl;
  }

  // heterogeneous Mach-velocity vector
  VectorFunction< std::complex< double > > MM = vector< std::complex< double > >(Mx, My, 0.);
  // Normal and tangential projections
  ScalarFunction< std::complex< double > > Mn = vector< std::complex< double > >(Mx, My, 0.) * normal< std::complex< double > >();
  ScalarFunction< std::complex< double > > Mt = vector< std::complex< double > >(Mx, My, 0.) * tangent< std::complex< double > >();
  // Parameters for Inverse Lorentz transformation
  ScalarFunction< std::complex< double > > beta = sqrt(1 - Mx * Mx - My * My); // Jacobian
  ScalarFunction< std::complex< double > > Alphax = 1 + Mx * Mx / (beta * (1 + beta));
  ScalarFunction< std::complex< double > > Alphay = 1 + My * My / (beta * (1 + beta));
  ScalarFunction< std::complex< double > > K = Mx * My / (beta * (1 + beta));
  // Tensor of the Inverse Lorentz transformation
  TensorFunction< std::complex< double > > Linv = tensor< std::complex< double > >(beta * Alphay, -beta * K, 0., -beta * K, beta * Alphax, 0., 0., 0., 0.);

  std::vector< FieldInterface< std::complex< double > > * > fieldBucket;
  Formulation< std::complex< double > > formulation("helmholtzflow");
  Field< std::complex< double >, Form::Form0 > u("u", domains.omega | domains.omegaPml | domains.gammaInf | domains.gammaPml | domains.corners | domains.gammaDir | domains.source, FunctionSpaceTypeForm0::HierarchicalH1, order);

  // convected Helmholz weak form
  formulation.integral(vector< std::complex< double > >(1 - Mx * Mx, -Mx * My, 0.) * grad(dof(u)), vector< std::complex< double > >(1., 0., 0.) * grad(tf(u)), domains.omega, gauss);
  formulation.integral(vector< std::complex< double > >(-Mx * My, 1 - My * My, 0.) * grad(dof(u)), vector< std::complex< double > >(0., 1., 0.) * grad(tf(u)), domains.omega, gauss);
  formulation.integral(-k * k * dof(u), tf(u), domains.omega, gauss);
  formulation.integral(dof(u), vector< std::complex< double > >(-im * k * Mx, -im * k * My, 0.) * grad(tf(u)), domains.omega, gauss, term::ProductType::Scalar);
  formulation.integral(vector< std::complex< double > >(im * k * Mx, im * k * My, 0.) * grad(dof(u)), tf(u), domains.omega, gauss);

  formulation.integral(-1., tf(u), domains.source, gauss);
  // Use plane wave instead of monopole
  // formulation.integral(- (cos(k * x< std::complex< double > >()) + im * sin(k * x< std::complex< double > >())), tf(u), domains.gammaDir, gauss);

  // Exact DtN: beta_n**2 * kx = - Mn * (k0 - Mt*ky) + sqrt( k0**2 -2*k0*Mt*ky - (1-abs(M))^2 ky**2 )
  if(abcName == "ABC-0") { // zeroth order Taylor approximation
    // boundary contributions
    formulation.integral(im * k * Mn * dof(u), tf(u), domains.gammaInf, gauss);
    formulation.integral(Mn * Mt * tangent< std::complex< double > >() * grad(dof(u)), dof(u), domains.gammaInf, gauss);
    // ABC
    msg::info << "Use zeroth order ABC" << msg::endl;
    formulation.integral(im * k * (1 - Mn) * dof(u), tf(u), domains.gammaInf, gauss);
    formulation.integral(-Mn * Mt * tangent< std::complex< double > >() * grad(dof(u)), dof(u), domains.gammaInf, gauss);
  }
  else if(abcName == "ABC-2") {
    // boundary contributions
    formulation.integral(im * k * Mn * dof(u), tf(u), domains.gammaInf, gauss);
    formulation.integral(Mn * Mt * tangent< std::complex< double > >() * grad(dof(u)), dof(u), domains.gammaInf, gauss);
    // ABC
    msg::info << "Use second order ABC" << msg::endl;
    // zeroth order contribution
    formulation.integral(im * k * (1 - Mn) * dof(u), tf(u), domains.gammaInf, gauss);
    // first order contribution
    formulation.integral((1 - Mn) * Mt * tangent< std::complex< double > >() * grad(dof(u)), dof(u), domains.gammaInf, gauss);
    // second order contribution
    formulation.integral(-im * (beta * beta) / (2 * k) * grad(dof(u)), grad(tf(u)), domains.gammaInf, gauss);
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

    if(geometry == 0) {
      // define the auxiliary fields
      std::vector< Field< std::complex< double >, Form::Form0 > * > phi;
      for(int i = 0; i < padeOrder; ++i) {
        phi.push_back(new Field< std::complex< double >, Form::Form0 >("phi_" + std::to_string(i), domains.gammaInf, FunctionSpaceTypeForm0::HierarchicalH1, order));
        fieldBucket.push_back(phi.back());
      }

      // write the augmented weak form - approximation of the square-root
      formulation.integral(im * k * exp2 * dof(u), tf(u), domains.gammaInf, gauss);
      for(int i = 0; i < padeOrder; ++i) {
        // boundary integral terms relating the auxiliary fields
        formulation.integral(im * k * exp2 * coef * c[i] * dof(*phi[i]), tf(u), domains.gammaInf, gauss);
        formulation.integral(im * k * exp2 * coef * c[i] * dof(u), tf(u), domains.gammaInf, gauss);

        // coupling of the auxiliary equations
        formulation.integral(-(beta * beta) * grad(dof(*phi[i])), grad(tf(*phi[i])), domains.gammaInf, gauss);
        formulation.integral(-2. * im * k * Mt * tangent< std::complex< double > >() * grad(dof(*phi[i])), tf(*phi[i]), domains.gammaInf, gauss);
        formulation.integral((k * k) * (exp1 * c[i] + 1.) * dof(*phi[i]), tf(*phi[i]), domains.gammaInf, gauss);
        formulation.integral((k * k) * exp1 * (c[i] + 1.) * dof(u), tf(*phi[i]), domains.gammaInf, gauss);
      }
    }
    else if(geometry == 1) {
      // define the auxilary fields
      std::vector< Field< std::complex< double >, Form::Form0 > * > phi[4];
      for(int x = 0; x < 4; ++x) {
        for(int i = 0; i < padeOrder; ++i) {
          phi[x].push_back(new Field< std::complex< double >, Form::Form0 >("phi_" + std::to_string(i) + "^" + std::to_string(x), domains.edge[x].gamma | domains.edge[x].corner[0] | domains.edge[x].corner[1], FunctionSpaceTypeForm0::HierarchicalH1, order));
          fieldBucket.push_back(phi[x].back());
        }
      }

      formulation.integral(im * k * exp2 * dof(u), tf(u), domains.gammaInf, gauss);
      for(int x = 0; x < 4; ++x) {
        for(int i = 0; i < padeOrder; ++i) {
          // Boundary conditions
          formulation.integral(im * k * exp2 * coef * c[i] * dof(*phi[x][i]), tf(u), domains.edge[x].gamma, gauss);
          formulation.integral(im * k * exp2 * coef * c[i] * dof(u), tf(u), domains.edge[x].gamma, gauss);

          // Auxiliary equations
          formulation.integral(-(beta * beta) * grad(dof(*phi[x][i])), grad(tf(*phi[x][i])), domains.edge[x].gamma, gauss);
          formulation.integral(-2. * im * k * Mt * tangent< std::complex< double > >() * grad(dof(*phi[x][i])), tf(*phi[x][i]), domains.edge[x].gamma, gauss);
          formulation.integral((k * k) * (exp1 * c[i] + 1.) * dof(*phi[x][i]), tf(*phi[x][i]), domains.edge[x].gamma, gauss);
          formulation.integral((k * k) * exp1 * (c[i] + 1.) * dof(u), tf(*phi[x][i]), domains.edge[x].gamma, gauss);
        }
      }

      bool withCorner = true;
      gmshFem.userDefinedParameter(withCorner, "withCorner");
      std::string CornerStrategy = "Sommerfeld"; // Sommerfeld or HABC
      gmshFem.userDefinedParameter(CornerStrategy, "CornerStrategy");

      if(withCorner) {
        msg::info << "with corner treatment" << msg::endl;
        msg::info << "use " << CornerStrategy << " condition at corners" << msg::endl;
        // define the corner fields
        std::vector< std::vector< std::vector< Field< std::complex< double >, Form::Form0 > * > > > corner_phi(4);
        for(int k = 0; k < 4; ++k) {
          corner_phi[k].resize(padeOrder);
          for(int i = 0; i < padeOrder; ++i) {
            for(int j = 0; j < padeOrder; ++j) {
              corner_phi[k][i].push_back(new Field< std::complex< double >, Form::Form0 >("corner_phi_(" + std::to_string(i) + ", " + std::to_string(j) + ")^(" + std::to_string(k) + ")", domains.corner[k], FunctionSpaceTypeForm0::HierarchicalH1, order));
              fieldBucket.push_back(corner_phi[k][i].back());
            }
          }
        }

        for(int ci = 0; ci < 4; ++ci) {
          int x = ci;
          int y = (ci - 1 < 0 ? 4 - 1 : ci - 1);
          for(int i = 0; i < padeOrder; ++i) {
            if(CornerStrategy == "Sommerfeld") {
              // strategy 1 - Sommerfeld condition at corners
              formulation.integral(-im * k * dof(*phi[x][i]), tf(*phi[x][i]), domains.corner[ci], gauss);
              formulation.integral(-im * k * dof(*phi[y][i]), tf(*phi[y][i]), domains.corner[ci], gauss);
            }

            else if(CornerStrategy == "HABC") {
              // strategy 2 - HABC at corners
              formulation.integral(-im * k * exp2 * dof(*phi[x][i]), tf(*phi[x][i]), domains.corner[ci], gauss);
              formulation.integral(-im * k * exp2 * dof(*phi[y][i]), tf(*phi[y][i]), domains.corner[ci], gauss);
              for(int j = 0; j < padeOrder; ++j) {
                // Auxiliary equations for corner
                formulation.integral((exp1 * c[i] + exp1 * c[j] + 1.) * dof(*corner_phi[ci][i][j]), tf(*corner_phi[ci][i][j]), domains.corner[ci], gauss);
                formulation.integral(exp1 * (c[j] + 1.) * dof(*phi[x][i]), tf(*corner_phi[ci][i][j]), domains.corner[ci], gauss);
                formulation.integral(exp1 * (c[i] + 1.) * dof(*phi[y][j]), tf(*corner_phi[ci][i][j]), domains.corner[ci], gauss);

                // Corner conditions
                formulation.integral(-im * k * exp2 * coef * c[j] * dof(*phi[x][i]), tf(*phi[x][i]), domains.corner[ci], gauss);
                formulation.integral(-im * k * exp2 * coef * c[j] * dof(*phi[y][i]), tf(*phi[y][i]), domains.corner[ci], gauss);
                formulation.integral(-im * k * exp2 * coef * c[j] * dof(*corner_phi[ci][i][j]), tf(*phi[x][i]), domains.corner[ci], gauss);
                formulation.integral(-im * k * exp2 * coef * c[j] * dof(*corner_phi[ci][j][i]), tf(*phi[y][i]), domains.corner[ci], gauss);
              }
            }
            else {
              msg::error << "Corner strategy not available !" << msg::endl;
              exit(0);
            }
          }
        }
      }
    }
  }
  else if(abcName == "pml") {
    msg::info << "Use a PML with " << Npml << " layers" << msg::endl;
    ScalarFunction< std::complex< double > > Sigma0 = beta; // use heterogeneous parameter
    ScalarFunction< std::complex< double > > det_J;
    TensorFunction< std::complex< double > > J_PML_inv_T;
    const double Wpml = Rpml - R;
    if(geometry == 0) {
      msg::info << "stabilized PML - circular domain" << msg::endl;
      ScalarFunction< std::complex< double > > cosT = x< std::complex< double > >() / r2d< std::complex< double > >();
      ScalarFunction< std::complex< double > > sinT = y< std::complex< double > >() / r2d< std::complex< double > >();
      ScalarFunction< std::complex< double > > dampingProfileR = Sigma0 / (Wpml - (r2d< std::complex< double > >() - R));
      ScalarFunction< std::complex< double > > dampingProfileInt = -Sigma0 * ln((Wpml - (r2d< std::complex< double > >() - R)) / Wpml);
      ScalarFunction< std::complex< double > > gamma = 1. - im * dampingProfileR / k;
      ScalarFunction< std::complex< double > > gamma_hat = 1. - im * (1. / r2d< std::complex< double > >()) * dampingProfileInt / k;
      det_J = gamma * gamma_hat;
      J_PML_inv_T = tensor< std::complex< double > >(cosT / gamma, sinT / gamma, 0., -sinT / gamma_hat, cosT / gamma_hat, 0., 0., 0., 0.);
    }
    else if(geometry == 1) {
      msg::info << "stabilized PML - square domain" << msg::endl;
      ScalarFunction< std::complex< double > > SigmaX = heaviside(abs(x< std::complex< double > >()) - R) * Sigma0 / (Rpml - abs(x< std::complex< double > >()));
      ScalarFunction< std::complex< double > > gammaX = 1 - (im / k) * SigmaX;
      ScalarFunction< std::complex< double > > SigmaY = heaviside(abs(y< std::complex< double > >()) - R) * Sigma0 / (Rpml - abs(y< std::complex< double > >()));
      ScalarFunction< std::complex< double > > gammaY = 1 - (im / k) * SigmaY;
      det_J = gammaX * gammaY;
      J_PML_inv_T = tensor< std::complex< double > >(1. / gammaX, 0., 0., 0., 1. / gammaY, 0., 0., 0., 0.);
    }
    // build vector and matrix for the weak form
    VectorFunction< std::complex< double > > J_PML_inv_T_M = J_PML_inv_T * MM;
    TensorFunction< std::complex< double > > J_PML_Linv = J_PML_inv_T * Linv;
    // stabilized PML weak form - valid for a general domain
    formulation.integral( det_J * J_PML_Linv*grad(dof(u)) , J_PML_Linv*grad(tf(u)) , domains.omegaPml, gauss, term::ProductType::Scalar);

    formulation.integral(+ k * k/(beta * beta) * det_J * J_PML_inv_T_M * dof(u) , J_PML_inv_T_M * tf(u), domains.omegaPml, gauss, term::ProductType::Scalar);
    formulation.integral(+ im * k / beta * det_J * J_PML_Linv * grad(dof(u)), J_PML_inv_T_M * tf(u), domains.omegaPml, gauss, term::ProductType::Scalar);

    formulation.integral(- im * k /beta * det_J * J_PML_inv_T_M * dof(u), J_PML_Linv * grad(tf(u)), domains.omegaPml, gauss, term::ProductType::Scalar);
    formulation.integral(- k * k/(beta * beta) * det_J * dof(u), tf(u), domains.omegaPml, gauss);
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
  save(+u, domains.omega, "u");
  save(+(sqrt(Mx * Mx + My * My)), domains.omega, "Mach");
  if(abcName == "pml") {
    save(+u, domains.omegaPml, "u_pml");
  }

  for(unsigned int i = 0; i < fieldBucket.size(); ++i) {
    delete fieldBucket[i];
  }

  return 0;
}
