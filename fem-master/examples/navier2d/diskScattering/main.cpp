#include <complex>
#include <gmsh.h>
#include <gmshfem/AnalyticalFunction.h>
#include <gmshfem/Domain.h>
#include <gmshfem/Exception.h>
#include <gmshfem/FieldInterface.h>
#include <gmshfem/Formulation.h>
#include <gmshfem/Function.h>
#include <gmshfem/GmshFem.h>
#include <gmshfem/Message.h>
#include <gmshfem/Navier2D.h>
#include <gmshfem/io.h>

using namespace gmshfem;
using namespace gmshfem::analytics;
using namespace gmshfem::common;
using namespace gmshfem::term;
using namespace gmshfem::problem;
using namespace gmshfem::domain;
using namespace gmshfem::field;
using namespace gmshfem::function;
using namespace gmshfem::post;
using namespace gmshfem::equation;

int main(int argc, char **argv)
{
  GmshFem gmshFem(argc, argv);

  gmsh::model::add("geometry");

  const double pi = 3.14159265358979323846264338327950288;
  const double rho = 1.; //.............[kg/m3]    masse volumique
  double f = 1.; //.............[Hz]       fréquence
  gmshFem.userDefinedParameter(f, "f");
  double w = 2 * pi * f; //.............[rad/s]    pulsation
  const double lambda = 1.; //.............[Pa]       premier  coefficient de Lamé
  const double mu = 1.; //.............[Pa]       deuxième coefficient de Lamé
  const double cP = sqrt((lambda + 2.0 * mu) / rho); //.............[m/s]      celerite de l'onde P
  const double cS = sqrt(mu / rho); //.............[m/s]      celerite de l'onde S
  double kP = w / cP; //.............[1/m]      nombre d'onde p
  double kS = w / cS; //.............[1/m]      nombre d'onde s

  int pointsPerSWavelength = 8; //.............[-]        number of points per S-wavelength
  gmshFem.userDefinedParameter(pointsPerSWavelength, "pointsPerSWavelength");
  double lc = 2 * pi / (pointsPerSWavelength * kS); //.............[m]        characteristic length

  const double rExt = 2.; //.............[m]        external radius of the annulus geometry
  const double rInt = 1.; //.............[m]        internal radius of the annulus geometry

  int ABCOrder = 2;
  gmshFem.userDefinedParameter(ABCOrder, "abcOrder");
  std::string meshPath = "";
  gmshFem.userDefinedParameter(meshPath, "mesh");

  if(meshPath == "") {
    gmsh::model::geo::addPoint(0., 0., 0., lc, 1); // Center

    // small circle
    gmsh::model::geo::addPoint(rInt, 0., 0., lc, 2);
    gmsh::model::geo::addPoint(0, rInt, 0., lc, 3);
    gmsh::model::geo::addPoint(-rInt, 0., 0., lc, 4);
    gmsh::model::geo::addPoint(0, -rInt, 0., lc, 5);

    gmsh::model::geo::addCircleArc(2, 1, 3, 1);
    gmsh::model::geo::addCircleArc(3, 1, 4, 2);
    gmsh::model::geo::addCircleArc(4, 1, 5, 3);
    gmsh::model::geo::addCircleArc(5, 1, 2, 4);

    gmsh::model::geo::addCurveLoop({1, 2, 3, 4}, 1);

    // big circle
    gmsh::model::geo::addPoint(rExt, 0., 0., lc, 6);
    gmsh::model::geo::addPoint(0, rExt, 0., lc, 7);
    gmsh::model::geo::addPoint(-rExt, 0., 0., lc, 8);
    gmsh::model::geo::addPoint(0, -rExt, 0., lc, 9);

    gmsh::model::geo::addCircleArc(6, 1, 7, 5);
    gmsh::model::geo::addCircleArc(7, 1, 8, 6);
    gmsh::model::geo::addCircleArc(8, 1, 9, 7);
    gmsh::model::geo::addCircleArc(9, 1, 6, 8);

    gmsh::model::geo::addCurveLoop({5, 6, 7, 8}, 2);

    // surface
    gmsh::model::geo::addPlaneSurface({2, 1}, 3);

    gmsh::model::geo::synchronize();

    // physicals
    gmsh::model::addPhysicalGroup(2, {3}, 3);
    gmsh::model::setPhysicalName(2, 3, "omega");
    gmsh::model::addPhysicalGroup(1, {1, 2, 3, 4}, 1);
    gmsh::model::setPhysicalName(1, 1, "gammaDir");
    gmsh::model::addPhysicalGroup(1, {5, 6, 7, 8}, 2);
    gmsh::model::setPhysicalName(1, 2, "gammaInf");

    gmsh::model::geo::synchronize();
    gmsh::model::mesh::generate();
  }
  else {
    gmsh::open(meshPath);
  }

  // modif Vanessa mai 2019
  Formulation< std::complex< double > > formulation("2D-Navier");

  Domain surface(2, 3);
  Domain gammaInf(1, 2);
  Domain gammaDir(1, 1);

  const std::complex< double > im = std::complex< double >(0., 1.);
  const double thetaInc = 0.;

  msg::info << " kp : " << kP << msg::endl;
  msg::info << " ks : " << kS << msg::endl;

  ScalarFunction< std::complex< double > > nX = xComp(normal< std::complex< double > >()); //x< std::complex< double > >() / r2d< std::complex< double > >();//
  ScalarFunction< std::complex< double > > nY = yComp(normal< std::complex< double > >()); //y< std::complex< double > >() / r2d< std::complex< double > >();//

  const unsigned int padeOrder = 4;
  const double padeAngle = pi / 4.;

  // Fonctions analytiques disponibles pour Navier
  AnalyticalFunction< navier2D::SoftPWavesScatteringByACylinder< std::complex< double > > > solution(kP, kS, rInt, 2. * floor(w) + 1, thetaInc);
  //AnalyticalFunction< navier2D::SoftSWavesScatteringByACylinder< std::complex< double > > > solution(kP, kS, rInt, 2. * floor(w) + 1, thetaInc);
  //AnalyticalFunction< navier2D::SoftPWavesScatteringByACylinderWithLOABC< std::complex< double > > > solution(kP, kS, rInt, rExt, lambda, mu, 2. * floor(w) + 1, thetaInc);
  //AnalyticalFunction< navier2D::SoftSWavesScatteringByACylinderWithLOABC< std::complex< double > > > solution(kP, kS, rInt, rExt, lambda, mu, 2. * floor(w) + 1, thetaInc);
  //AnalyticalFunction< navier2D::SoftPWavesScatteringByACylinderWithHOABC< std::complex< double > > > solution(f, rInt, rExt, 2. * floor(w) + 1);
  //AnalyticalFunction< navier2D::SoftSWavesScatteringByACylinderWithHOABC< std::complex< double > > > solution(f, rInt, rExt, 2. * floor(w) + 1);

  Field< std::complex< double >, Form::Form0 > ux("ux", surface | gammaInf, FunctionSpaceTypeForm0::Lagrange);
  Field< std::complex< double >, Form::Form0 > uy("uy", surface | gammaInf, FunctionSpaceTypeForm0::Lagrange);

  Field< std::complex< double >, Form::Form0 > vx("vx", gammaInf, FunctionSpaceTypeForm0::Lagrange);
  Field< std::complex< double >, Form::Form0 > vy("vy", gammaInf, FunctionSpaceTypeForm0::Lagrange);

  Field< std::complex< double >, Form::Form0 > v0("v0", gammaInf, FunctionSpaceTypeForm0::Lagrange);
  Field< std::complex< double >, Form::Form0 > v1("v1", gammaInf, FunctionSpaceTypeForm0::Lagrange);

  Field< std::complex< double >, Form::Form0 > qx("qx", gammaInf, FunctionSpaceTypeForm0::Lagrange);
  Field< std::complex< double >, Form::Form0 > qy("qy", gammaInf, FunctionSpaceTypeForm0::Lagrange);

  Field< std::complex< double >, Form::Form0 > q0("q0", gammaInf, FunctionSpaceTypeForm0::Lagrange);
  Field< std::complex< double >, Form::Form0 > q1("q1", gammaInf, FunctionSpaceTypeForm0::Lagrange);

  Field< std::complex< double >, Form::Form0 > tx("tx", gammaInf, FunctionSpaceTypeForm0::Lagrange);
  Field< std::complex< double >, Form::Form0 > ty("ty", gammaInf, FunctionSpaceTypeForm0::Lagrange);

  Field< std::complex< double >, Form::Form0 > un("un", gammaInf, FunctionSpaceTypeForm0::Lagrange);
  Field< std::complex< double >, Form::Form0 > ut("ut", gammaInf, FunctionSpaceTypeForm0::Lagrange);

  std::vector< Field< std::complex< double >, Form::Form0 > > h;
  std::vector< Field< std::complex< double >, Form::Form0 > > ix;
  std::vector< Field< std::complex< double >, Form::Form0 > > j;
  std::vector< Field< std::complex< double >, Form::Form0 > > k;

  for(unsigned int i = 0; i < padeOrder; ++i) {
    h.push_back(Field< std::complex< double >, Form::Form0 >("h_" + std::to_string(i), gammaInf, FunctionSpaceTypeForm0::Lagrange));
    ix.push_back(Field< std::complex< double >, Form::Form0 >("ix_" + std::to_string(i), gammaInf, FunctionSpaceTypeForm0::Lagrange));
    j.push_back(Field< std::complex< double >, Form::Form0 >("j_" + std::to_string(i), gammaInf, FunctionSpaceTypeForm0::Lagrange));
    k.push_back(Field< std::complex< double >, Form::Form0 >("k_" + std::to_string(i), gammaInf, FunctionSpaceTypeForm0::Lagrange));
  }

  const std::complex< double > kPEps = kP + im * 0.4 * pow(kP, 1. / 3) * pow(rExt, -2. / 3);
  const std::complex< double > kSEps = kS + im * 0.4 * pow(kS, 1. / 3) * pow(rExt, -2. / 3);

  //msg::info <<  " kpEps : " <<  kPEps << msg::endl << "\n";
  //msg::info <<  " ksEps : " <<  kSEps << msg::endl << "\n"; // ok vérifié 25/06

  std::complex< double > padeD[padeOrder];
  std::complex< double > padeC[padeOrder];
  std::complex< double > padeR[padeOrder];
  std::complex< double > padeS[padeOrder];

  for(unsigned int i = 0; i < padeOrder; ++i) {
    padeD[i] = 1.0 + pow(std::tan((3.14159265359 / (2.0 * (double)padeOrder)) * (0.5 + (double)i)), 2);
    padeC[i] = padeD[i] / (double)padeOrder;
    padeR[i] = padeC[i] * std::complex< double >(std::cos(padeAngle / 2.0), std::sin(padeAngle / 2.0));
    padeS[i] = 1.0 + std::complex< double >(std::cos(padeAngle), std::sin(padeAngle)) * (-1.0 + padeD[i]);

    //msg::info <<  " padeD : " <<  padeD[i] << msg::endl;
    //msg::info <<  " padeC : " <<  padeC[i] << msg::endl;
    msg::info << i << " padeR : " << padeR[i] << msg::endl;
    msg::info << i << " padeS : " << padeS[i] << msg::endl; //OK vérifié 25/06
  }

  // Choisir l'onde incidente (P ou S)

  ux.addConstraint(gammaDir, im * kP * cos(-kP * x< std::complex< double > >()) - kP * sin(-kP * x< std::complex< double > >()));
  uy.addConstraint(gammaDir, 0.);

  //ux.addConstraint(gammaDir, 0.);
  //uy.addConstraint(gammaDir, -im*kS*cos(-kS * x< std::complex< double > >()) + kS * sin(-kS * x< std::complex< double > >()));

  TensorFunction< std::complex< double > > C_xx = tensor< std::complex< double > >((lambda + 2. * mu), 0., 0., 0., mu, 0., 0., 0., 0.);
  TensorFunction< std::complex< double > > C_xy = tensor< std::complex< double > >(0., lambda, 0., mu, 0., 0., 0., 0., 0.);
  TensorFunction< std::complex< double > > C_yx = tensor< std::complex< double > >(0., mu, 0., lambda, 0., 0., 0., 0., 0.);
  TensorFunction< std::complex< double > > C_yy = tensor< std::complex< double > >(mu, 0., 0., 0., (lambda + 2. * mu), 0., 0., 0., 0.);

  //-------------------
  //- Navier equation -
  //-------------------
  // In surface int sigma(u) : epsilon(u')
  formulation.integral(C_xx * grad(dof(ux)), grad(tf(ux)), surface, "Gauss2");
  formulation.integral(C_xy * grad(dof(uy)), grad(tf(ux)), surface, "Gauss2");
  formulation.integral(C_yx * grad(dof(ux)), grad(tf(uy)), surface, "Gauss2");
  formulation.integral(C_yy * grad(dof(uy)), grad(tf(uy)), surface, "Gauss2");

  // In surface int -w^2 rho u . u'
  formulation.integral(-w * w * rho * dof(ux), tf(ux), surface, "Gauss2");
  formulation.integral(-w * w * rho * dof(uy), tf(uy), surface, "Gauss2");

  //---------------------------------
  //- Absorbing Boundary Conditions -
  //---------------------------------
  if(ABCOrder == 0) {
    msg::info << "Use zero order Lysmer Kuhlemeyer ABC." << msg::endl;

    // On ABC boundary int -i w cp rho In u.u' -i w cs rho It u.u'
    formulation.integral(-im * w * cP * rho * nX * nX * dof(ux), tf(ux), gammaInf, "Gauss2");
    formulation.integral(-im * w * cP * rho * nX * nY * dof(uy), tf(ux), gammaInf, "Gauss2");
    formulation.integral(-im * w * cP * rho * nX * nY * dof(ux), tf(uy), gammaInf, "Gauss2");
    formulation.integral(-im * w * cP * rho * nY * nY * dof(uy), tf(uy), gammaInf, "Gauss2");

    formulation.integral(-im * w * cS * rho * nY * nY * dof(ux), tf(ux), gammaInf, "Gauss2");
    formulation.integral(im * w * cS * rho * nX * nY * dof(uy), tf(ux), gammaInf, "Gauss2");
    formulation.integral(-im * w * cS * rho * nX * nX * dof(uy), tf(uy), gammaInf, "Gauss2");
    formulation.integral(im * w * cS * rho * nX * nY * dof(ux), tf(uy), gammaInf, "Gauss2");
  }
  else if(ABCOrder == 2) {
    msg::info << "Use high order ABC." << msg::endl;
    //------------------------------------------------------------------------------------------------------------------------------
    //- See p 1713 of "A high-order ABC for 2D time-harmonic elastodynamic scattering problems", V. Mattesi, M. Darbas C. Geuzaine -
    //------------------------------------------------------------------------------------------------------------------------------


    //-------------------------------------------------------------------
    //- STEP 1 : Find v in H^1/2(Gamma) such that v =  Lambda_{1,eps} u -
    //-------------------------------------------------------------------

    // On gammaInf int v.v' - i rho w^2 [  1/kPEps (n vO,v') + 1/kSEps (t v1,v') ] (=0)
    formulation.integral(dof(vx), tf(vx), gammaInf, "Gauss15");
    formulation.integral(dof(vy), tf(vy), gammaInf, "Gauss15");

    formulation.integral(-im * rho * pow(w, 2) * (1.0 / kP) * nX * dof(v0), tf(vx), gammaInf, "Gauss15");
    formulation.integral(-im * rho * pow(w, 2) * (1.0 / kP) * nY * dof(v0), tf(vy), gammaInf, "Gauss15");

    formulation.integral(-im * rho * pow(w, 2) * (1.0 / kS) * (-nY) * dof(v1), tf(vx), gammaInf, "Gauss15");
    formulation.integral(-im * rho * pow(w, 2) * (1.0 / kS) * (nX)*dof(v1), tf(vy), gammaInf, "Gauss15");

    // (v0,v0') - sum R_i (h_i,v0') =0
    formulation.integral(dof(v0), tf(v0), gammaInf, "Gauss15");

    for(unsigned int i = 0; i < padeOrder; ++i) {
      formulation.integral(-padeR[i] * dof(h[i]), tf(v0), gammaInf, "Gauss15");

      // Sl(hl,hl') - (grad_G/kPEps hl,grad_G/kPEps hl') -(n.u,hl') (=0) forall l=1,...,L
      formulation.integral(padeS[i] * dof(h[i]), tf(h[i]), gammaInf, "Gauss15");

      formulation.integral(-pow(1.0 / kPEps, 2) * grad(dof(h[i])), grad(tf(h[i])), gammaInf, "Gauss15");

      formulation.integral(-nX * dof(ux), tf(h[i]), gammaInf, "Gauss15");
      formulation.integral(-nY * dof(uy), tf(h[i]), gammaInf, "Gauss15");
    }
    // (v1,v1') - sum Rl (il,v1') =0
    formulation.integral(dof(v1), tf(v1), gammaInf, "Gauss15");
    for(unsigned int i = 0; i < padeOrder; ++i) {
      formulation.integral(-padeR[i] * dof(ix[i]), tf(v1), gammaInf, "Gauss15");

      // Sl( il,il') -(grad_G/kSEps il,grad_G/kSEps il') - (t. u,il') (=0) forall l=1,...L
      formulation.integral(padeS[i] * dof(ix[i]), tf(ix[i]), gammaInf, "Gauss15");

      formulation.integral(-pow(1.0 / kSEps, 2) * grad(dof(ix[i])), grad(tf(ix[i])), gammaInf, "Gauss15");

      formulation.integral(-(-nY) * dof(ux), tf(ix[i]), gammaInf, "Gauss15");
      formulation.integral(-(nX)*dof(uy), tf(ix[i]), gammaInf, "Gauss15");
    }

    //-------------------------------------------------------------------
    //- Step 2 : Find q in H^-1/2(Gamma) st : (I+ Lambda_{2,eps}) q = v -
    //-------------------------------------------------------------------

    // (q,q') -i[1/ks(grad_G q0,q') - 1/kp (q1 n,q')] - (v,q') (=0)
    formulation.integral(dof(qx), tf(qx), gammaInf, "Gauss15");
    formulation.integral(dof(qy), tf(qy), gammaInf, "Gauss15");

    formulation.integral(-im * (1.0 / kS) * vector< std::complex< double > >(1., 0., 0.) * grad(dof(q0)), tf(qx), gammaInf, "Gauss15");
    formulation.integral(-im * (1.0 / kS) * vector< std::complex< double > >(0., 1., 0.) * grad(dof(q0)), tf(qy), gammaInf, "Gauss15");

    formulation.integral(im * (1.0 / kP) * nX * vector< std::complex< double > >(-nY, nX, 0.) * grad(dof(q1)), tf(qx), gammaInf, "Gauss15");
    formulation.integral(im * (1.0 / kP) * nY * vector< std::complex< double > >(-nY, nX, 0.) * grad(dof(q1)), tf(qy), gammaInf, "Gauss15");

    formulation.integral(-dof(vx), tf(qx), gammaInf, "Gauss15");
    formulation.integral(-dof(vy), tf(qy), gammaInf, "Gauss15");

    //(q0,q0')=sum Rl (jl,q0')
    formulation.integral(dof(q0), tf(q0), gammaInf, "Gauss15");
    for(unsigned int i = 0; i < padeOrder; ++i) {
      formulation.integral(-padeR[i] * dof(j[i]), tf(q0), gammaInf, "Gauss15");

      //Sl(jl,jl')-(grad/kSEps jl,grad/kSEps jl') - (q.n,jl') =0
      formulation.integral(padeS[i] * dof(j[i]), tf(j[i]), gammaInf, "Gauss15");

      formulation.integral(-pow(1.0 / kSEps, 2) * grad(dof(j[i])), grad(tf(j[i])), gammaInf, "Gauss15");

      formulation.integral(-nX * dof(qx), tf(j[i]), gammaInf, "Gauss15");
      formulation.integral(-nY * dof(qy), tf(j[i]), gammaInf, "Gauss15");
    }
    // (q1,q1')= sum (Rl kl , q1')
    formulation.integral(dof(q1), tf(q1), gammaInf, "Gauss15");

    for(unsigned int i = 0; i < padeOrder; ++i) {
      formulation.integral(-padeR[i] * dof(k[i]), tf(q1), gammaInf, "Gauss15");

      // (Sl kl,kl') -(grad_G/kPEps kl,grad_G/kPEps kl') - (t.q, kl') (=0)
      formulation.integral(padeS[i] * dof(k[i]), tf(k[i]), gammaInf, "Gauss15");

      formulation.integral(-pow(1.0 / kPEps, 2) * grad(dof(k[i])), grad(tf(k[i])), gammaInf, "Gauss15");

      formulation.integral(-(-nY) * dof(qx), tf(k[i]), gammaInf, "Gauss15");
      formulation.integral(-(nX)*dof(qy), tf(k[i]), gammaInf, "Gauss15");
    }

    //----------------------------------
    //- Step 3: -tu.u' with t=q+2mu Mu -
    //----------------------------------
    formulation.integral(dof(tx), tf(tx), gammaInf, "Gauss15");
    formulation.integral(dof(ty), tf(ty), gammaInf, "Gauss15");

    // - (q,u') - 2 mu (Mu,u')
    formulation.integral(-dof(qx), tf(tx), gammaInf, "Gauss15");
    formulation.integral(-dof(qy), tf(ty), gammaInf, "Gauss15");

    formulation.integral(-2.0 * mu * vector< std::complex< double > >(1., 0., 0.) * grad(dof(un)), tf(tx), gammaInf, "Gauss15");
    formulation.integral(-2.0 * mu * vector< std::complex< double > >(0., 1., 0.) * grad(dof(un)), tf(ty), gammaInf, "Gauss15");

    formulation.integral(2.0 * mu * nX * (-nY) * vector< std::complex< double > >(1., 0., 0.) * grad(dof(ut)), tf(tx), gammaInf, "Gauss15");
    formulation.integral(2.0 * mu * nX * nX * vector< std::complex< double > >(0., 1., 0.) * grad(dof(ut)), tf(tx), gammaInf, "Gauss15");
    formulation.integral(2.0 * mu * nY * (-nY) * vector< std::complex< double > >(1., 0., 0.) * grad(dof(ut)), tf(ty), gammaInf, "Gauss15");
    formulation.integral(2.0 * mu * nY * nX * vector< std::complex< double > >(0., 1., 0.) * grad(dof(ut)), tf(ty), gammaInf, "Gauss15");

    // (un,un')=(u.n,un')
    formulation.integral(dof(un), tf(un), gammaInf, "Gauss15");

    formulation.integral(-nX * dof(ux), tf(un), gammaInf, "Gauss15");
    formulation.integral(-nY * dof(uy), tf(un), gammaInf, "Gauss15");

    // (ut,ut')=(u.t,ut')
    formulation.integral(dof(ut), tf(ut), gammaInf, "Gauss15");

    formulation.integral(-(-nY) * dof(ux), tf(ut), gammaInf, "Gauss15");
    formulation.integral(-(nX)*dof(uy), tf(ut), gammaInf, "Gauss15");

    formulation.integral(-dof(tx), tf(ux), gammaInf, "Gauss15");
    formulation.integral(-dof(ty), tf(uy), gammaInf, "Gauss15");
  }

  // Prepro
  formulation.pre();

  // Pro
  formulation.assemble();
  formulation.solve();

  save(ux, surface, "ux", "pos");
  save(uy);

  //save(grad(ux), gammaInf, "grad_ux","pos");

  //save(xComp(solution), surface, "ux_exact");
  //save(yComp(solution), surface, "uy_exact");

  //std::complex< double > num = integrate(pow(abs(xComp(solution) - ux), 2) + pow(abs(yComp(solution) - uy), 2), surface, "Gauss2");
  //std::complex< double > den = integrate(pow(abs(xComp(solution) ), 2) + pow(abs(yComp(solution) ), 2), surface, "Gauss2");

  //std::complex<double> L2error = sqrt(num / den);
  //msg::info << "L_2 error = " << L2error << msg::endl;

  //CSVio file("convergence.txt", ' ', OpeningMode::Append);
  //file << 1. / pointsPerSWavelength << std::real(sqrt(num / den)) << csv::endl;

  return 0;
}
