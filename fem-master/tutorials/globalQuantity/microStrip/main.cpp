#include <gmsh.h>
#include <gmshfem/Formulation.h>
#include <gmshfem/GmshFem.h>
#include <gmshfem/Message.h>

using namespace gmshfem;
using namespace gmshfem::common;
using namespace gmshfem::problem;
using namespace gmshfem::domain;
using namespace gmshfem::field;
using namespace gmshfem::function;
using namespace gmshfem::post;
using namespace gmshfem::equation;

int main(int argc, char **argv)
{
  GmshFem gmshFem(argc, argv);

  double Q = 0.001;
  gmshFem.userDefinedParameter(Q, "Q");
  double V = 0.001;
  gmshFem.userDefinedParameter(V, "V");
  double epsilonRelDiel = 9.8;
  gmshFem.userDefinedParameter(epsilonRelDiel, "eps");
  bool enforceCharge = false;
  gmshFem.userDefinedParameter(enforceCharge, "enforceCharge");

  gmsh::open("../floating.geo");
  gmsh::model::mesh::generate();

  //*****
  // Problem declaration
  //*****

  // Allocate the formulation object
  Formulation< double > formulation("demoFloatingPotential");

  Domain air(2, 101);
  Domain dielectric(2, 111);
  Domain ground(1, 120);
  Domain electrode(1, 121);
  Domain surfaceInfinity(1, 130);

  const double eps0 = 8.854187818e-12;
  ScalarPiecewiseFunction< double > epsilon;
  epsilon.addFunction(eps0, air);
  epsilon.addFunction(eps0 * epsilonRelDiel, dielectric);

  // Allocate field object
  Field< double, Form::Form0 > v("v", air | dielectric, FunctionSpaceTypeForm0::Lagrange);
  // Apply Dirichlet BC on surfaceInfinity
  v.addConstraint(surfaceInfinity, 0.);

  // Define the global quantity associated to the potential of the electrode
  GlobalQuantity< double > vElectrode("vElectrode", electrode);
  v.assignGlobalQuantity(vElectrode);

  // Define the global quantity associated to the potential of the ground
  GlobalQuantity< double > vGround("vGround", ground);
  v.assignGlobalQuantity(vGround);

  // Write the corresponding weak formulation terms by terms
  formulation.integral(epsilon * grad(dof(v)), grad(tf(v)), air | dielectric, "Gauss6");

  // Write the global terms
  if(enforceCharge) {
    formulation.globalTerm(vElectrode, FixedComponent::Dual, Q);
  }
  else {
    formulation.globalTerm(vElectrode, FixedComponent::Primal, V);
  }
  formulation.globalTerm(vGround, FixedComponent::Primal, 0.);

  // Prepro
  formulation.pre();

  // Pro
  formulation.assemble();
  formulation.solve();

  // Compute the capacity : C = Q / V
  const double C = std::abs(vElectrode.getDualValue() / vElectrode.getPrimalValue());

  // Compute the energy : epsilon / 2 * ||e||^2
  const double energy = integrate(epsilon / 2. * pow(norm(-grad(v)), 2), air | dielectric, "Gauss6");

  msg::info << "Global:" << msg::endl;
  msg::info << "* The capacitance is " << C << "[F]" << msg::endl;
  msg::info << "* The energy is " << energy << "[J]" << msg::endl;
  msg::info << "Ground:" << msg::endl;
  msg::info << "* The charge is " << vGround.getDualValue() << "[C]" << msg::endl;
  msg::info << "Microstrip:" << msg::endl;
  msg::info << "* The charge is " << vElectrode.getDualValue() << "[C]" << msg::endl;
  msg::info << "* The potential is " << vElectrode.getPrimalValue() << "[V]" << msg::endl;

  //  Do not check this since Gmsh is compiled without mesh module in CI dockers...
  //
  //  if(C < 2.54063e-10 - 1e-14 || C > 2.54063e-10 + 1e-14) {
  //    throw common::Exception("It seems that there is an error with the global quantities implementation. The value of the capacitance should be 2.54063e-10 [F].");
  //  }

  save(v);

  return 0;
}
