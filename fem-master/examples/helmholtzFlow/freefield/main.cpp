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

struct Edge {
  Domain gamma;
  Domain corner[2];
};

struct HelmholtzDomain {
  Domain omega;
  Domain omegaErr;
  Domain gammaPhy;
  Domain gammaPml;
  Domain omegaPml;
  Domain source;
  Domain corners;
  std::vector< Domain > corner;
  std::vector< Edge > edge;
};

//***************************************
// GEOMETRY
//***************************************

void meshCircle(GmshFem &gmshFem, HelmholtzDomain &domains, const double R, const double Rpml, const double xs, const double ys, const double lc, const int ElemOrder, const bool withPml)
{
  gmsh::model::add("geometry");

  gmsh::model::geo::addPoint(0., 0., 0., lc, 1); // Center
  gmsh::model::geo::addPoint(xs, ys, 0., lc / 5., 100); // Source point
  // physical circle
  gmsh::model::geo::addPoint(R, 0., 0., lc, 2);
  gmsh::model::geo::addPoint(0, R, 0., lc, 3);
  gmsh::model::geo::addPoint(-R, 0., 0., lc, 4);
  gmsh::model::geo::addPoint(0, -R, 0., lc, 5);

  gmsh::model::geo::addCircleArc(2, 1, 3, 1);
  gmsh::model::geo::addCircleArc(3, 1, 4, 2);
  gmsh::model::geo::addCircleArc(4, 1, 5, 3);
  gmsh::model::geo::addCircleArc(5, 1, 2, 4);
  gmsh::model::geo::addCurveLoop({1, 2, 3, 4}, 1);

  double Rs = 2. * lc;
  // mini circle to surround the point source - to compute the L2error
  gmsh::model::geo::addPoint(xs + Rs, ys, 0., lc / 2., 6);
  gmsh::model::geo::addPoint(xs, Rs + ys, 0., lc / 2., 7);
  gmsh::model::geo::addPoint(xs - Rs, ys, 0., lc / 2., 8);
  gmsh::model::geo::addPoint(xs, ys - Rs, 0., lc / 2., 9);

  gmsh::model::geo::addCircleArc(6, 100, 7, 5);
  gmsh::model::geo::addCircleArc(7, 100, 8, 6);
  gmsh::model::geo::addCircleArc(8, 100, 9, 7);
  gmsh::model::geo::addCircleArc(9, 100, 6, 8);

  gmsh::model::geo::addCurveLoop({5, 6, 7, 8}, 2);

  // surfaces
  gmsh::model::geo::addPlaneSurface({2}, 2);
  gmsh::model::geo::addPlaneSurface({1, 2}, 1);

  if(withPml) {
    gmsh::model::geo::addPoint(Rpml, 0., 0., lc, 10);
    gmsh::model::geo::addPoint(0, Rpml, 0., lc, 11);
    gmsh::model::geo::addPoint(-Rpml, 0., 0., lc, 12);
    gmsh::model::geo::addPoint(0, -Rpml, 0., lc, 13);

    gmsh::model::geo::addCircleArc(10, 1, 11, 9);
    gmsh::model::geo::addCircleArc(11, 1, 12, 10);
    gmsh::model::geo::addCircleArc(12, 1, 13, 11);
    gmsh::model::geo::addCircleArc(13, 1, 10, 12);
    gmsh::model::geo::addCurveLoop({9, 10, 11, 12}, 3);
    gmsh::model::geo::addPlaneSurface({3, 1}, 3);
  }
  gmsh::model::geo::synchronize();
  gmsh::model::mesh::embed(0, {100}, 2, 2);

  // physicals
  gmsh::model::addPhysicalGroup(2, {1}, 1);
  gmsh::model::setPhysicalName(2, 1, "omega");
  gmsh::model::addPhysicalGroup(2, {2}, 2);
  gmsh::model::setPhysicalName(2, 2, "omegaSource");
  gmsh::model::addPhysicalGroup(1, {1, 2, 3, 4}, 3);
  gmsh::model::setPhysicalName(1, 3, "gammaPhy");
  gmsh::model::addPhysicalGroup(0, {100}, 100);
  gmsh::model::setPhysicalName(0, 100, "source");
  if(withPml) {
    gmsh::model::addPhysicalGroup(2, {3}, 3);
    gmsh::model::setPhysicalName(2, 3, "omegaPml");
    gmsh::model::addPhysicalGroup(1, {9, 10, 11, 12}, 4);
    gmsh::model::setPhysicalName(1, 4, "gammaPml");
  }

  // generate mesh
  gmsh::model::geo::synchronize();
  gmsh::model::mesh::generate(2);
  gmsh::model::mesh::setOrder(ElemOrder);
  gmsh::write("m_circle.msh");

  // Domains
  domains.omega = Domain(2, 1) | Domain(2, 2);
  domains.omegaErr = Domain(2, 1);
  domains.gammaPhy = Domain(1, 3);
  domains.source = Domain(0, 100);

  if(withPml) {
    domains.omegaPml = Domain(2, 3);
    domains.gammaPml = Domain(1, 4);
  }
}

void meshSquare(GmshFem &gmshFem, HelmholtzDomain &domains, const double L, const double Lpml, const double xs, const double ys, const double lc, const int ElemOrder, const bool withPml)
{
  gmsh::model::add("geometry");

  gmsh::model::geo::addPoint(0., 0., 0., lc, 1); // Center
  gmsh::model::geo::addPoint(xs, ys, 0., lc / 5., 100); // Source point

  // polygone
  gmsh::model::geo::addPoint(L, L, 0., lc, 2);
  gmsh::model::geo::addPoint(-L, L, 0., lc, 3);
  gmsh::model::geo::addPoint(-L, -L, 0., lc, 4);
  gmsh::model::geo::addPoint(L, -L, 0., lc, 5);

  gmsh::model::geo::addLine(2, 3, 1);
  gmsh::model::geo::addLine(3, 4, 2);
  gmsh::model::geo::addLine(4, 5, 3);
  gmsh::model::geo::addLine(5, 2, 4);
  gmsh::model::geo::addCurveLoop({1, 2, 3, 4}, 1);

  double Rs = 2. * lc;
  // mini circle to surround the point source - to compute the L2error
  gmsh::model::geo::addPoint(xs + Rs, ys, 0., lc / 2., 6);
  gmsh::model::geo::addPoint(xs, Rs + ys, 0., lc / 2., 7);
  gmsh::model::geo::addPoint(xs - Rs, ys, 0., lc / 2., 8);
  gmsh::model::geo::addPoint(xs, ys - Rs, 0., lc / 2., 9);

  gmsh::model::geo::addCircleArc(6, 100, 7, 5);
  gmsh::model::geo::addCircleArc(7, 100, 8, 6);
  gmsh::model::geo::addCircleArc(8, 100, 9, 7);
  gmsh::model::geo::addCircleArc(9, 100, 6, 8);
  gmsh::model::geo::addCurveLoop({5, 6, 7, 8}, 2);

  // surfaces
  gmsh::model::geo::addPlaneSurface({2}, 2);
  gmsh::model::geo::addPlaneSurface({1, 2}, 1);

  if(withPml) {
    gmsh::model::geo::addPoint(Lpml, Lpml, 0., lc, 10);
    gmsh::model::geo::addPoint(-Lpml, Lpml, 0., lc, 11);
    gmsh::model::geo::addPoint(-Lpml, -Lpml, 0., lc, 12);
    gmsh::model::geo::addPoint(Lpml, -Lpml, 0., lc, 13);

    gmsh::model::geo::addLine(10, 11, 9);
    gmsh::model::geo::addLine(11, 12, 10);
    gmsh::model::geo::addLine(12, 13, 11);
    gmsh::model::geo::addLine(13, 10, 12);
    gmsh::model::geo::addCurveLoop({9, 10, 11, 12}, 3);
    gmsh::model::geo::addPlaneSurface({3, 1}, 3);
  }
  gmsh::model::geo::synchronize();
  gmsh::model::mesh::embed(0, {100}, 2, 2);
  // physicals
  gmsh::model::addPhysicalGroup(2, {1}, 1);
  gmsh::model::setPhysicalName(2, 1, "omega");
  gmsh::model::addPhysicalGroup(2, {2}, 2);
  gmsh::model::setPhysicalName(2, 2, "omegaSource");
  gmsh::model::addPhysicalGroup(1, {1, 2, 3, 4}, 3);
  gmsh::model::setPhysicalName(1, 3, "gammaPhy");
  gmsh::model::addPhysicalGroup(0, {100}, 100);
  gmsh::model::setPhysicalName(0, 100, "source");
  // Physical edges
  gmsh::model::addPhysicalGroup(1, {1}, 20);
  gmsh::model::setPhysicalName(1, 20, "gammaPhy_top");
  gmsh::model::addPhysicalGroup(1, {2}, 30);
  gmsh::model::setPhysicalName(1, 30, "gammaPhy_left");
  gmsh::model::addPhysicalGroup(1, {3}, 40);
  gmsh::model::setPhysicalName(1, 40, "gammaPhy_bottom");
  gmsh::model::addPhysicalGroup(1, {4}, 50);
  gmsh::model::setPhysicalName(1, 50, "gammaPhy_right");
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
    gmsh::model::addPhysicalGroup(2, {3}, 3);
    gmsh::model::setPhysicalName(2, 3, "omegaPml");
    gmsh::model::addPhysicalGroup(1, {9, 10, 11, 12}, 4);
    gmsh::model::setPhysicalName(1, 4, "gammaPml");
  }

  // generate mesh
  gmsh::model::geo::synchronize();
  gmsh::model::mesh::generate(2);
  gmsh::model::mesh::setOrder(ElemOrder);
  gmsh::write("m_square.msh");

  // Domains
  domains.omega = Domain(2, 1) | Domain(2, 2);
  domains.omegaErr = Domain(2, 1);
  domains.gammaPhy = Domain(1, 3);
  domains.source = Domain(0, 100);

  if(withPml) {
    domains.omegaPml = Domain(2, 3);
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
    domains.gammaPhy |= domains.edge[i].gamma;
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
  double R = 1.; // Circle radius or square of size [-R,R] x [-R,R]
  gmshFem.userDefinedParameter(R, "R");
  // Mesh order
  int ElemOrder = 4;
  // Point source location
  double xs = 0.;
  double ys = 0.;
  gmshFem.userDefinedParameter(xs, "xs");
  gmshFem.userDefinedParameter(ys, "ys");
  msg::info << "Point source located at xs=" << xs << ", ys=" << ys << msg::endl;

  // physical parameters
  const double pi = 3.14159265358979323846264338327950288;
  double k = 6 * pi; // free field wavenumber
  gmshFem.userDefinedParameter(k, "k");
  const std::complex< double > im(0., 1.);
  double M = 0.8; // Mach number
  gmshFem.userDefinedParameter(M, "M");
  double theta = pi / 4.; // mean flow orientation
  gmshFem.userDefinedParameter(theta, "theta");
  msg::info << "Mach number " << M << " oriented at theta=" << theta << " rad" << msg::endl;
  
  // numerical parameters
  int order; // FEM shape function order
  if (M >= 0.9) {
    order = 9;
  }
  else if (M < 0.9 && M >= 0.7) {
    order = 6;
  }
  else {
    order = 4;
  }
  
  std::string gauss = "Gauss" + std::to_string(2 * order + 1);
  double pointsByWl = 8; // pointsByWl = (2 * pi * order * (1 - abs(M))) / (lc * k);
  double lc = (2 * pi * order * (1 - abs(M))) / (pointsByWl * k); // choose lc based on pointsByWl
  msg::info << " - pointsByWl =  " << (2 * pi * order * (1 - abs(M))) / (lc * k) << msg::endl;
  msg::info << " - lc =  " << lc << msg::endl;
  msg::info << " - FEM basis order = " << order << "" << msg::endl;
 
  // Choose ABC - exterior boundary condition
  std::string abcName = "Pade"; // available choices: ABC-0, ABC-2, Pade, pml
  gmshFem.userDefinedParameter(abcName, "abcName");
  bool L0 = false; // activate zeroth-order symbol
  gmshFem.userDefinedParameter(L0, "L0");
  bool withCorner = true; // try corner treatment for square domain 
  gmshFem.userDefinedParameter(withCorner, "withCorner");
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
    meshCircle(gmshFem, domains, R, Rpml, xs, ys, lc, ElemOrder, abcName == "pml");
    break;
  case 1:
    meshSquare(gmshFem, domains, R, Rpml, xs, ys, lc, ElemOrder, abcName == "pml");
    break;
  default:
    msg::error << "geometry not available ! " << msg::endl;
    exit(0);
    break;
  }

  // Mach number projections
  double Mx = M * std::cos(theta);
  double My = M * std::sin(theta);
  // Mach-velocity vector
  VectorFunction< std::complex< double > > MM = vector< std::complex< double > >(Mx, My, 0.);
  // Normal and tangential projections
  ScalarFunction< std::complex< double > > Mn = vector< std::complex< double > >(Mx, My, 0.) * normal< std::complex< double > >();
  ScalarFunction< std::complex< double > > Mt = vector< std::complex< double > >(Mx, My, 0.) * tangent< std::complex< double > >();
  // Parameters for Inverse Lorentz transformation
  double beta = sqrt(1 - M * M); // Jacobian
  double Alphax = 1 + Mx * Mx / (beta * (1 + beta));
  double Alphay = 1 + My * My / (beta * (1 + beta));
  double K = Mx * My / (beta * (1 + beta));
  // Tensor of the Inverse Lorentz transformation
  TensorFunction< std::complex< double > > Linv = tensor< std::complex< double > >(beta * Alphay, -beta * K, 0., -beta * K, beta * Alphax, 0., 0., 0., 0.);

  std::vector< FieldInterface< std::complex< double > > * > fieldBucket;
  Formulation< std::complex< double > > formulation("helmholtzflow");
  Field< std::complex< double >, Form::Form0 > u("u", domains.omega | domains.omegaPml | domains.gammaPhy | domains.gammaPml | domains.corners | domains.source, FunctionSpaceTypeForm0::HierarchicalH1, order);

  // Define analytic solution
  Function< std::complex< double >, Degree::Degree0 > *solution = nullptr;
  solution = new AnalyticalFunction< helmholtz2D::MonopoleFreeField< std::complex< double > > >(k, M, theta, 1., xs, ys, 0., 0.);

  // convected Helmholz weak form
  formulation.integral(vector< std::complex< double > >(1 - Mx * Mx, -Mx * My, 0.) * grad(dof(u)), vector< std::complex< double > >(1., 0., 0.) * grad(tf(u)), domains.omega, gauss);
  formulation.integral(vector< std::complex< double > >(-Mx * My, 1 - My * My, 0.) * grad(dof(u)), vector< std::complex< double > >(0., 1., 0.) * grad(tf(u)), domains.omega, gauss);
  formulation.integral(-k * k * dof(u), tf(u), domains.omega, gauss);
  formulation.integral(dof(u), vector< std::complex< double > >(-im * k * Mx, -im * k * My, 0.) * grad(tf(u)), domains.omega, gauss, term::ProductType::Scalar);
  formulation.integral(vector< std::complex< double > >(im * k * Mx, im * k * My, 0.) * grad(dof(u)), tf(u), domains.omega, gauss);

  // Source point
  formulation.integral(-1., tf(u), domains.source, gauss);

  // Exact DtN: beta_n**2 * kx = - Mn * (k0 - Mt*ky) + sqrt( k0**2 -2*k0*Mt*ky - (1-abs(M))^2 ky**2 )
  if(abcName == "ABC-0") { // zeroth order Taylor approximation
    // boundary contributions
    formulation.integral(im * k * Mn * dof(u), tf(u), domains.gammaPhy, gauss);
    formulation.integral(Mn * Mt * tangent< std::complex< double > >() * grad(dof(u)), dof(u), domains.gammaPhy, gauss);
    // ABC
    msg::info << "Use zeroth order ABC" << msg::endl;
    formulation.integral(im * k * (1 - Mn) * dof(u), tf(u), domains.gammaPhy, gauss);
    formulation.integral(-Mn * Mt * tangent< std::complex< double > >() * grad(dof(u)), dof(u), domains.gammaPhy, gauss);
    if ( L0 && (geometry==0) ) {
      msg::info << "Use curvature correction" << msg::endl;
      formulation.integral( + (beta * 0.5 / R) * dof(u), tf(u), domains.gammaPhy, gauss);
    }
  }
  else if(abcName == "ABC-2") {
    // boundary contributions
    formulation.integral(im * k * Mn * dof(u), tf(u), domains.gammaPhy, gauss);
    formulation.integral(Mn * Mt * tangent< std::complex< double > >() * grad(dof(u)), dof(u), domains.gammaPhy, gauss);
    // ABC
    msg::info << "Use second order ABC" << msg::endl;
    // zeroth order contribution
    formulation.integral(im * k * (1 - Mn) * dof(u), tf(u), domains.gammaPhy, gauss);
    // first order contribution
    formulation.integral((1 - Mn) * Mt * tangent< std::complex< double > >() * grad(dof(u)), dof(u), domains.gammaPhy, gauss);
    // second order contribution
    formulation.integral(- im * (beta * beta) /(2 * k) * grad(dof(u)), grad(tf(u)), domains.gammaPhy, gauss);
    if ( L0 && (geometry==0) ) {
      formulation.integral( + (beta * 0.5 / R) * dof(u), tf(u), domains.gammaPhy, gauss);
    }
    if ( withCorner && (geometry==1) ) {
      msg::info << "with corner treatment" << msg::endl;
      // try a corner condition for 2nd order abc, M = 0
      for(int ci = 0; ci < 4; ++ci) {
        formulation.integral( 0.75 * dof(u), tf(u), domains.corner[ci], gauss);
      }
    }
  }
  else if(abcName == "Pade") {
    int padeOrder = 4;
    gmshFem.userDefinedParameter(padeOrder, "padeOrder");
    double angle = 0.;
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
        phi.push_back(new Field< std::complex< double >, Form::Form0 >("phi_" + std::to_string(i), domains.gammaPhy, FunctionSpaceTypeForm0::HierarchicalH1, order));
        fieldBucket.push_back(phi.back());
      }

      // write the augmented weak form - approximation of the square-root
      formulation.integral(im * k * exp2 * dof(u), tf(u), domains.gammaPhy, gauss);
      if ( L0 && (geometry==0) ) {
        formulation.integral( + (beta * 0.5 / R) * dof(u), tf(u), domains.gammaPhy, gauss);
      }
      for(int i = 0; i < padeOrder; ++i) {
        // boundary integral terms relating the auxiliary fields
        formulation.integral(im * k * exp2 * coef * c[i] * dof(*phi[i]), tf(u), domains.gammaPhy, gauss);
        formulation.integral(im * k * exp2 * coef * c[i] * dof(u), tf(u), domains.gammaPhy, gauss);

        // coupling of the auxiliary equations
        formulation.integral(- (beta * beta) * grad(dof(*phi[i])), grad(tf(*phi[i])), domains.gammaPhy, gauss);
        formulation.integral(-2. * im * k * Mt * tangent< std::complex< double > >() * grad(dof(*phi[i])), tf(*phi[i]), domains.gammaPhy, gauss);
        // regularization
        std::complex<double> keps = k;//-0.*im*0.4*pow(k,(1/3.))*pow((beta/R),2/3.);
        formulation.integral((keps * keps) * (exp1 * c[i] + 1.) * dof(*phi[i]), tf(*phi[i]), domains.gammaPhy, gauss);
        formulation.integral((keps * keps) * exp1 * (c[i] + 1.) * dof(u), tf(*phi[i]), domains.gammaPhy, gauss);
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

      formulation.integral(im * k * exp2 * dof(u), tf(u), domains.gammaPhy, gauss);
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
          msg::info << "Corners: x=" << x << ", y=" << y << msg::endl;
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
    double Sigma0 = beta;
    gmshFem.userDefinedParameter(Sigma0, "Sigma0");
    ScalarFunction< std::complex< double > > det_J;
    TensorFunction< std::complex< double > > J_PML_inv_T;
    const double Wpml = Rpml - R;
    if(geometry == 0) {
      msg::info << "stabilized PML - circular domain" << msg::endl;
      msg::info << "use unbounded absorbing function with sigma_0 = " << Sigma0 << msg::endl;
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
      msg::info << "use unbounded absorbing function with sigma_0 = " << Sigma0 << msg::endl;
      ScalarFunction< std::complex< double > > SigmaX = heaviside(abs(x< std::complex< double > >()) - R) * Sigma0 / (Rpml - abs(x< std::complex< double > >()));
      ScalarFunction< std::complex< double > > gammaX = 1 - (im / k) * SigmaX;
      ScalarFunction< std::complex< double > > SigmaY = heaviside(abs(y< std::complex< double > >()) - R) * Sigma0 / (Rpml - abs(y< std::complex< double > >()));
      ScalarFunction< std::complex< double > > gammaY = 1 - (im / k) * SigmaY;
      det_J = gammaX * gammaY;
      J_PML_inv_T = tensor< std::complex< double > >(1. / gammaX, 0., 0., 0., 1. / gammaY, 0., 0., 0., 0.);

      /*
      msg::info << "classical PML - square domain - unstable ! " << msg::endl;
      formulation.integral(vector< std::complex< double > >((gammaY/gammaX)*(1-Mx*Mx), -Mx*My, 0.) * grad(dof(u)), vector< std::complex< double > >(1., 0., 0.) * grad(tf(u)), domains.omegaPml, gauss);
      formulation.integral(vector< std::complex< double > >(-Mx*My, (gammaX/gammaY)*(1-My*My), 0.) * grad(dof(u)), vector< std::complex< double > >(0., 1., 0.) * grad(tf(u)), domains.omegaPml, gauss);
      formulation.integral(dof(u), vector< std::complex< double > >(-gammaY*im*k*Mx,-gammaX*im*k*My, 0.) * grad(tf(u)), domains.omegaPml, gauss);
      formulation.integral(vector< std::complex< double > >(gammaY*im*k*Mx, gammaX*im*k*My, 0.) * grad(dof(u)), tf(u), domains.omegaPml, gauss);
      formulation.integral(- k*k* gammaX*gammaY * dof(u), tf(u), domains.omegaPml, gauss);
      */

      /*
      msg::info << "Alternative stable PML - square domain only " << msg::endl;
      ScalarFunction< std::complex< double > > Ax = (gammaY/gammaX)*(1-Mx*Mx);
      std::complex< double > Cx =-im*k*Mx*( 1 - (My*My)/(1-Mx*Mx) ) / (beta*beta);
      double Ka = -Mx*My/(1-Mx*Mx);
      ScalarFunction< std::complex< double > > Ay = (gammaX/gammaY)*( 1 - (My*My)/(1-Mx*Mx) );
      std::complex< double > Cy = -im*k*My/(beta*beta);

      // modified bilinear term in dx
      formulation.integral( Ax * vector< std::complex< double > >(1. , Ka, 0.) * grad(dof(u)), vector< std::complex< double > >(1., Ka, 0.) * grad(tf(u)), domains.omegaPml, gauss);
      formulation.integral( Ax * Cx * dof(u), vector< std::complex< double > >(1., Ka, 0.) * grad(tf(u)), domains.omegaPml, gauss);
      formulation.integral( Ax * (-Cx) * vector< std::complex< double > >(1. , Ka, 0.) * grad(dof(u)), tf(u), domains.omegaPml, gauss);
      formulation.integral( Ax * Cx * (-Cx) *dof(u), tf(u), domains.omegaPml, gauss);
      // modified bilinear term in dy
      formulation.integral( Ay * vector< std::complex< double > >(0. , 1., 0.) * grad(dof(u)), vector< std::complex< double > >(0., 1., 0.) * grad(tf(u)), domains.omegaPml, gauss);
      formulation.integral( Ay * Cy * dof(u), vector< std::complex< double > >(0., 1., 0.) * grad(tf(u)), domains.omegaPml, gauss);
      formulation.integral( Ay * vector< std::complex< double > >(0. , -Cy, 0.) * grad(dof(u)), tf(u), domains.omegaPml, gauss);
      formulation.integral( Ay * Cy * (-Cy) *dof(u), tf(u), domains.omegaPml, gauss);
      //
      formulation.integral(- k*k/(beta*beta) * gammaX*gammaY * dof(u), tf(u), domains.omegaPml, gauss);
      */
    }
    // build vector and matrix for the weak form
    VectorFunction< std::complex< double > > J_PML_inv_T_M = J_PML_inv_T * MM;
    TensorFunction< std::complex< double > > J_PML_Linv = J_PML_inv_T * Linv;
    // stabilized PML weak form - valid for a general domain
    formulation.integral(det_J * J_PML_Linv*grad(dof(u)) , J_PML_Linv*grad(tf(u)) , domains.omegaPml, gauss, term::ProductType::Scalar);

    formulation.integral(+k * k / (beta * beta) * det_J * J_PML_inv_T_M * dof(u) , J_PML_inv_T_M * tf(u), domains.omegaPml, gauss, term::ProductType::Scalar);
    formulation.integral(+im * k / beta * det_J * J_PML_Linv * grad(dof(u)), J_PML_inv_T_M * tf(u), domains.omegaPml, gauss, term::ProductType::Scalar);

    formulation.integral(-im * k / beta * det_J * J_PML_inv_T_M * dof(u), J_PML_Linv * grad(tf(u)), domains.omegaPml, gauss, term::ProductType::Scalar);
    formulation.integral(-k * k /(beta * beta) * det_J * dof(u), tf(u), domains.omegaPml, gauss);
  }
  else {
    msg::error << "ABC Type not available !" << msg::endl;
    exit(0);
  }

  // solve
  formulation.pre();
  formulation.assemble();
  formulation.solve();

  // compute the projection of the analytic solution on the FEM basis - best FEM approximation
  Formulation< std::complex< double > > projection("helmholtzflow");
  Field< std::complex< double >, Form::Form0 > uP("uP", domains.omega, FunctionSpaceTypeForm0::HierarchicalH1, order);
  projection.integral(-dof(uP), tf(uP), domains.omega, gauss);
  projection.integral(*solution, tf(uP), domains.omega, gauss);

  msg::info << "Compute best approximation... " << msg::endl;
  projection.pre();
  projection.assemble();
  projection.solve();

  // save results and compute errors
  save(+u, domains.omega, "u");
  save(+uP, domains.omega, "uP");
  if(abcName == "pml") {
    save(+u, domains.omegaPml, "u_pml");
  }
  save(*solution, domains.omega, "u_exact");
  save(*solution - u, domains.omega, "error");

  std::complex< double > num = integrate(pow(abs(*solution - u), 2), domains.omegaErr, gauss);
  std::complex< double > numP = integrate(pow(abs(*solution - uP), 2), domains.omegaErr, gauss);
  std::complex< double > den = integrate(pow(abs(*solution), 2), domains.omegaErr, gauss);
  msg::info << "Relative L_2 error = " << 100. * sqrt(num / den) << " %" << msg::endl;
  msg::info << "Best L_2 error = " << 100. * sqrt(numP / den) << " %" << msg::endl;

  std::complex< double > numGamma = integrate(pow(abs(*solution - u), 2), domains.gammaPhy, gauss);
  std::complex< double > denGamma = integrate(pow(abs(*solution), 2), domains.gammaPhy, gauss);
  msg::info << "Interface relative L_2 error = " << 100. * sqrt(numGamma / denGamma) << " %" << msg::endl;

  for(unsigned int i = 0; i < fieldBucket.size(); ++i) {
    delete fieldBucket[i];
  }

  return 0;
}
