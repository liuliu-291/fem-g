// GmshFEM - Copyright (C) 2019-2021 Université de Liège
//
// See the LICENSE.txt file for license information. Please report all
// issues on https://gitlab.onelab.info/gmsh/fem/issues

#include <Formulation.h>
#include <Function.h>
#include <Post.h>

#include "formulations.h"
#include "geo.h"


using gmshfem::equation::dof;
using gmshfem::equation::dt_dof;
using gmshfem::equation::dt2_dof;
using gmshfem::equation::tf;
using gmshfem::function::operator-;

void formulations(gmshfem::common::GmshFem &gmshfem, unsigned int &numTest)
{
  gmshfem::msg::info << ++numTest << ") Formulation" << gmshfem::msg::endl;

  Geo2D::triangle();
  
  try {
    gmshfem::msg::info << "Global quantity:" << gmshfem::msg::endl;
    globalQuantity< std::complex< double > >();
    globalQuantity< std::complex< float > >();
    globalQuantity< double >();
    globalQuantity< float >();
  } catch(...) { throw; }
  
  try {
    gmshfem::msg::info << "Eigensolve" << gmshfem::msg::endl;
    eigensolve< std::complex< double > >();
    eigensolve< std::complex< float > >();
    eigensolve< double >();
    eigensolve< float >();
  } catch(...) { throw; }
    
  removeGeo();
  
  Geo1D::line();
  
  try {
    gmshfem::msg::info << "Non-symmetric matrix" << gmshfem::msg::endl;
    nonsymmetricMatrix< std::complex< double > >();
    nonsymmetricMatrix< std::complex< float > >();
    nonsymmetricMatrix< double >();
    nonsymmetricMatrix< float >();
  } catch(...) { throw; }
  
  removeGeo();
  
  Geo2D::periodic();
  
  try {
    gmshfem::msg::info << "Periodic" << gmshfem::msg::endl;
    periodic< std::complex< double > >();
    periodic< std::complex< float > >();
    periodic< double >();
    periodic< float >();
  } catch(...) { throw; }
  
  removeGeo();
  
  return;
}

template< class T_Scalar >
void globalQuantity()
{
  // Primal, Form0, Lagrange
  {
    // play with Domain class
    gmshfem::domain::Domain omega("omega");
    gmshfem::domain::Domain dirichlet(1, 1);
    gmshfem::domain::Domain neumann;
    neumann = gmshfem::domain::Domain(std::make_pair(1, 2));
    gmshfem::domain::Domain dir(std::move(dirichlet));
    gmshfem::domain::Domain neu(neumann);
        
    gmshfem::problem::Formulation< T_Scalar > formulation("globalQuantity");

    gmshfem::field::Field< T_Scalar, gmshfem::field::Form::Form0 > u("u", omega | neu, gmshfem::field::FunctionSpaceTypeForm0::Lagrange);
    gmshfem::field::GlobalQuantity< T_Scalar > Ubis("U", dir);
    gmshfem::field::GlobalQuantity< T_Scalar > U(Ubis);
    u.assignGlobalQuantity(U);
    
    U.name("myU");
    
    formulation.integral(grad(dof(u)), d(tf(u)), omega, "Gauss5");
    formulation.integral(-1., tf(u), omega, "Gauss5");
    formulation.integral(1., tf(u), neu, "Gauss5");
    
    formulation.globalTerm(U, gmshfem::field::FixedComponent::Primal, 1.);

    formulation.pre(gmshfem::problem::DofsSort::Algorithm::Hilbert);
    formulation.assemble(false, gmshfem::problem::ElementsSort::Algorithm::Hilbert);
    formulation.solve();
    
    // play with it
    formulation.setRHSToZero();
    formulation.setSystemToZero();
    formulation.removeSystem();
    formulation.removeTerms();
    
    gmshfem::post::save(u);
    
    const T_Scalar check1 = gmshfem::post::integrate(u, neu, "Gauss5");
    const T_Scalar check2 = U.getDualValue();
    const T_Scalar expected1 = 0.5;
    const T_Scalar expected2 = 0.;
    gmshfem::msg::info << "\t- Check value = (" << check1 << ", " << check2 << "), expected = (" << expected1 << ", " << expected2 << ")" << gmshfem::msg::endl;
    if(std::abs(check1-expected1) > 100. * gmshfem::scalar::Epsilon< T_Scalar >::value || std::abs(check2-expected2) > 100. * gmshfem::scalar::Epsilon< T_Scalar >::value || U.name() != "myU") {
      throw gmshfem::common::Exception("Error in " + std::string(__FILE__) + ", line " + std::to_string(__LINE__) + " -> function: globalQuantity< " + gmshfem::scalar::Name< T_Scalar >::name + " >");
    }
  }
  
  // Primal, Form0, HierarchicalH1
  {
    // play with Domain class
    gmshfem::domain::Domain omega("omega");
    gmshfem::domain::Domain dirichlet(1, 1);
    gmshfem::domain::Domain neumann;
    neumann = gmshfem::domain::Domain(std::make_pair(1, 2));
    gmshfem::domain::Domain dir(std::move(dirichlet));
    gmshfem::domain::Domain neu(neumann);
        
    gmshfem::problem::Formulation< T_Scalar > formulation("globalQuantity");

    gmshfem::field::Field< T_Scalar, gmshfem::field::Form::Form0 > u("u", omega | neu, gmshfem::field::FunctionSpaceTypeForm0::HierarchicalH1, 1);
    gmshfem::field::GlobalQuantity< T_Scalar > U;
    U = gmshfem::field::GlobalQuantity< T_Scalar >("U", dir);
    u.assignGlobalQuantity(U);
    
    formulation.integral(grad(dof(u)), d(tf(u)), omega, "Gauss5");
    formulation.integral(-1., tf(u), omega, "Gauss5");
    formulation.integral(1., tf(u), neu, "Gauss5");
    
    formulation.globalTerm(U, gmshfem::field::FixedComponent::Primal, 1.);

    formulation.pre(gmshfem::problem::DofsSort::Algorithm::RCM);
    formulation.assemble();
    formulation.solve();
    
    // play with it
    formulation.setRHSToZero();
    formulation.setSystemToZero();
    formulation.removeSystem();
    formulation.removeTerms();
    
    gmshfem::post::save(u);
    
    const T_Scalar check1 = gmshfem::post::integrate(u, neu, "Gauss5");
    const T_Scalar check2 = U.getDualValue();
    const T_Scalar expected1 = 0.5;
    const T_Scalar expected2 = 0.;
    gmshfem::msg::info << "\t- Check value = (" << check1 << ", " << check2 << "), expected = (" << expected1 << ", " << expected2 << ")" << gmshfem::msg::endl;
    if(std::abs(check1-expected1) > 100. * gmshfem::scalar::Epsilon< T_Scalar >::value || std::abs(check2-expected2) > 100. * gmshfem::scalar::Epsilon< T_Scalar >::value || U.model() != u.model()) {
      throw gmshfem::common::Exception("Error in " + std::string(__FILE__) + ", line " + std::to_string(__LINE__) + " -> function: globalQuantity< " + gmshfem::scalar::Name< T_Scalar >::name + " >");
    }
  }
  
  // Dual, Form0, Lagrange
  {
    // play with Domain class
    gmshfem::domain::Domain omega("omega");
    gmshfem::domain::Domain dirichlet(1, 1);
    gmshfem::domain::Domain neumann;
    neumann = gmshfem::domain::Domain(std::make_pair(1, 2));
    gmshfem::domain::Domain dir(std::move(dirichlet));
    gmshfem::domain::Domain neu(neumann);
        
    gmshfem::problem::Formulation< T_Scalar > formulation("globalQuantity");

    gmshfem::field::Field< T_Scalar, gmshfem::field::Form::Form0 > u("u", omega | dir, gmshfem::field::FunctionSpaceTypeForm0::Lagrange);
    gmshfem::field::GlobalQuantity< T_Scalar > U("U", neu);
    u.assignGlobalQuantity(U);
    u.addConstraint(dir, 1.);
        
    formulation.integral(grad(dof(u)), d(tf(u)), omega, "Gauss5");
    formulation.integral(-1., tf(u), omega, "Gauss5");
    
    formulation.globalTerm(U, gmshfem::field::FixedComponent::Dual, -1.);

    formulation.pre(gmshfem::problem::DofsSort::Algorithm::None);
    formulation.assemble();
    formulation.solve();
    
    // play with it
    formulation.setRHSToZero();
    formulation.setSystemToZero();
    formulation.removeSystem();
    formulation.removeTerms();
    
    gmshfem::post::save(u);
    
    const T_Scalar check1 = gmshfem::post::integrate(u, neu, "Gauss5");
    const T_Scalar check2 = U.getPrimalValue();
    const T_Scalar expected1 = 0.5;
    const T_Scalar expected2 = 0.5;
    gmshfem::msg::info << "\t- Check value = (" << check1 << ", " << check2 << "), expected = (" << expected1 << ", " << expected2 << ")" << gmshfem::msg::endl;
    if(std::abs(check1-expected1) > 100. * gmshfem::scalar::Epsilon< T_Scalar >::value || std::abs(check2-expected2) > 100. * gmshfem::scalar::Epsilon< T_Scalar >::value || U.domain() != neu) {
      throw gmshfem::common::Exception("Error in " + std::string(__FILE__) + ", line " + std::to_string(__LINE__) + " -> function: globalQuantity< " + gmshfem::scalar::Name< T_Scalar >::name + " >");
    }
  }
  
  // Dual, Form0, Hierarchical
  {
    // play with Domain class
    gmshfem::domain::Domain omega("omega");
    gmshfem::domain::Domain dirichlet(1, 1);
    gmshfem::domain::Domain neumann;
    neumann = gmshfem::domain::Domain(std::make_pair(1, 2));
    gmshfem::domain::Domain dir(std::move(dirichlet));
    gmshfem::domain::Domain neu(neumann);
        
    gmshfem::problem::Formulation< T_Scalar > formulation("globalValue");

    gmshfem::field::Field< T_Scalar, gmshfem::field::Form::Form0 > u("u", omega | dir, gmshfem::field::FunctionSpaceTypeForm0::HierarchicalH1, 2);
    gmshfem::field::GlobalQuantity< T_Scalar > U("U", neu);
    u.assignGlobalQuantity(U);
    u.addConstraint(dir, 1.);
        
    formulation.integral(grad(dof(u)), d(tf(u)), omega, "Gauss5");
    formulation.integral(-1., tf(u), omega, "Gauss5");
    
    formulation.globalTerm(U, gmshfem::field::FixedComponent::Dual, -1.);

    formulation.pre();
    formulation.assemble();
    formulation.solve();
    
    // play with it
    formulation.setRHSToZero();
    formulation.setSystemToZero();
    formulation.removeSystem();
    formulation.removeTerms();
    
    gmshfem::post::save(u, omega, "u", "pos");
    
    const T_Scalar check1 = gmshfem::post::integrate(u, neu, "Gauss5");
    const T_Scalar check2 = U.getPrimalValue();
    const T_Scalar expected1 = 0.5;
    const T_Scalar expected2 = 0.5;
    gmshfem::msg::info << "\t- Check value = (" << check1 << ", " << check2 << "), expected = (" << expected1 << ", " << expected2 << ")" << gmshfem::msg::endl;
    if(std::abs(check1-expected1) > 100. * gmshfem::scalar::Epsilon< T_Scalar >::value || std::abs(check2-expected2) > 100. * gmshfem::scalar::Epsilon< T_Scalar >::value || U.domain() != neu) {
      throw gmshfem::common::Exception("Error in " + std::string(__FILE__) + ", line " + std::to_string(__LINE__) + " -> function: globalQuantity< " + gmshfem::scalar::Name< T_Scalar >::name + " >");
    }
  }
}

template void globalQuantity< std::complex< double > >();
template void globalQuantity< std::complex< float > >();
template void globalQuantity< double >();
template void globalQuantity< float >();

template< class T_Scalar >
void eigensolve()
{
  // Form0, MCK
  {
    gmshfem::domain::Domain omega("omega");

    gmshfem::problem::Formulation< T_Scalar > formulation("eigensolve");

    gmshfem::field::Field< T_Scalar, gmshfem::field::Form::Form0 > u("u", omega, gmshfem::field::FunctionSpaceTypeForm0::Lagrange);
    u.addConstraint(omega.getBoundary(), 0.);

    formulation.integral(dt2_dof(u), tf(u), omega, "Gauss5");
    formulation.integral(dt_dof(u), tf(u), omega, "Gauss5");
    formulation.integral(grad(dof(u)), grad(tf(u)), omega, "Gauss5");

    formulation.pre();
    formulation.assemble(true);
    gmshfem::algebra::Vector< gmshfem::scalar::ComplexPrecision< T_Scalar > > eigenvalues;
    formulation.eigensolve(eigenvalues, true, 1);

    gmshfem::post::save(gmshfem::function::eigenfunction(u, 0), omega, "eig");
    eigenvalues.save("eigenvalues");
    const T_Scalar check = std::abs(eigenvalues[0]); // use norm to remove mode degeneracy

    T_Scalar expected; // value measured on commit ead7f0c4620d949ad6756a0a86739bdd78213108
    expected = std::abs(gmshfem::scalar::ComplexPrecision< T_Scalar >(-4.999999999999073e-01,-4.441996170000912e+00));
    gmshfem::msg::info << "\t- Check value = " << gmshfem::msg::precision(gmshfem::scalar::PrecisionDigits< T_Scalar >::value+5) << check << ", expected = " << expected << gmshfem::msg::endl;
    if(std::abs(check-expected) > 10 * std::sqrt(gmshfem::scalar::Epsilon< T_Scalar >::value)) {
      throw gmshfem::common::Exception("Error in " + std::string(__FILE__) + ", line " + std::to_string(__LINE__) + " -> function: eigensolve< " + gmshfem::scalar::Name< T_Scalar >::name + " >");
    }
  }

  // Form0, MK
  {
    gmshfem::domain::Domain omega("omega");

    gmshfem::problem::Formulation< T_Scalar > formulation("eigensolve");

    gmshfem::field::Field< T_Scalar, gmshfem::field::Form::Form0 > u("u", omega, gmshfem::field::FunctionSpaceTypeForm0::Lagrange);
    u.addConstraint(omega.getBoundary(), 0.);

    formulation.integral(dt2_dof(u), tf(u), omega, "Gauss5");
    formulation.integral(grad(dof(u)), grad(tf(u)), omega, "Gauss5");

    formulation.pre();
    formulation.assemble(true);
    gmshfem::algebra::Vector< gmshfem::scalar::ComplexPrecision< T_Scalar > > eigenvalues;
    formulation.eigensolve(eigenvalues, true, 1);
    eigenvalues.save("eigenvalues");
    const T_Scalar check = std::abs(eigenvalues[0]); // use norm to remove mode degeneracy

    T_Scalar expected; // value measured on commit ead7f0c4620d949ad6756a0a86739bdd78213108
    expected = std::abs(gmshfem::scalar::ComplexPrecision< T_Scalar >(1.998132997430292e+01,-1.146783584060433e-15));
    gmshfem::msg::info << "\t- Check value = " << gmshfem::msg::precision(gmshfem::scalar::PrecisionDigits< T_Scalar >::value+5) << check << ", expected = " << expected << gmshfem::msg::endl;
    if(std::abs(check-expected) > 10 * std::sqrt(gmshfem::scalar::Epsilon< T_Scalar >::value)) {
      throw gmshfem::common::Exception("Error in " + std::string(__FILE__) + ", line " + std::to_string(__LINE__) + " -> function: eigensolve< " + gmshfem::scalar::Name< T_Scalar >::name + " >");
    }
  }

  // Form0, CK
  {
    gmshfem::domain::Domain omega("omega");

    gmshfem::problem::Formulation< T_Scalar > formulation("eigensolve");

    gmshfem::field::Field< T_Scalar, gmshfem::field::Form::Form0 > u("u", omega, gmshfem::field::FunctionSpaceTypeForm0::Lagrange);
    u.addConstraint(omega.getBoundary(), 0.);

    formulation.integral(dt_dof(u), tf(u), omega, "Gauss5");
    formulation.integral(grad(dof(u)), grad(tf(u)), omega, "Gauss5");

    formulation.pre();
    formulation.assemble(true);
    gmshfem::algebra::Vector< gmshfem::scalar::ComplexPrecision< T_Scalar > > eigenvalues;
    formulation.eigensolve(eigenvalues, true, 1);
    eigenvalues.save("eigenvalues");
    const T_Scalar check = std::abs(eigenvalues[0]); // use norm to remove mode degeneracy

    T_Scalar expected; // value measured on commit ead7f0c4620d949ad6756a0a86739bdd78213108
    expected = std::abs(gmshfem::scalar::ComplexPrecision< T_Scalar >(1.998132997430292e+01,-1.146783584060433e-15));
    gmshfem::msg::info << "\t- Check value = " << gmshfem::msg::precision(gmshfem::scalar::PrecisionDigits< T_Scalar >::value+5) << check << ", expected = " << expected << gmshfem::msg::endl;
    if(std::abs(check-expected) > 10 * std::sqrt(gmshfem::scalar::Epsilon< T_Scalar >::value)) {
      throw gmshfem::common::Exception("Error in " + std::string(__FILE__) + ", line " + std::to_string(__LINE__) + " -> function: eigensolve< " + gmshfem::scalar::Name< T_Scalar >::name + " >");
    }
  }

  // Form1, MK
  {
    gmshfem::domain::Domain omega("omega");

    gmshfem::problem::Formulation< T_Scalar > formulation("eigensolve");

    gmshfem::field::Field< T_Scalar, gmshfem::field::Form::Form1 > u("u", omega, gmshfem::field::FunctionSpaceTypeForm1::HierarchicalHCurl, 2);
    u.addConstraint(omega.getBoundary(), gmshfem::function::vector< T_Scalar >(1., 1., 0.));

    formulation.integral(dt2_dof(u), tf(u), omega, "Gauss5");
    formulation.integral(curl(dof(u)), d(tf(u)), omega, "Gauss5");

    formulation.pre();
    formulation.assemble(true);
    gmshfem::algebra::Vector< gmshfem::scalar::ComplexPrecision< T_Scalar > > eigenvalues;
    const gmshfem::scalar::Precision< T_Scalar > target = 9.8; // Try to skip the grad kernel
    formulation.eigensolve(eigenvalues, true, 20, target);
    eigenvalues.save("eigenvalues");

    std::vector< gmshfem::scalar::Precision< T_Scalar > > eigenvaluesN(eigenvalues.size());
    for(size_t i = 0; i < eigenvaluesN.size(); ++i)
      eigenvaluesN[i] = std::abs(eigenvalues[i]); // use norm to remove mode degeneracy
    std::sort(eigenvaluesN.begin(), eigenvaluesN.end());
    const gmshfem::scalar::Precision< T_Scalar > check = *std::upper_bound(eigenvaluesN.begin(), eigenvaluesN.end(), target); // Use the eigenvalue whcich is the closest to the target

    T_Scalar expected; // value measured on commit ead7f0c4620d949ad6756a0a86739bdd78213108
    expected = std::abs(gmshfem::scalar::ComplexPrecision< T_Scalar >(9.869644904987688e+00,6.036597792859928e-15));
    gmshfem::msg::info << "\t- Check value = " << gmshfem::msg::precision(gmshfem::scalar::PrecisionDigits< T_Scalar >::value+5) << check << ", expected = " << expected << gmshfem::msg::endl;
    if(std::abs(check-expected) > 10 * std::sqrt(gmshfem::scalar::Epsilon< T_Scalar >::value)) {
      throw gmshfem::common::Exception("Error in " + std::string(__FILE__) + ", line " + std::to_string(__LINE__) + " -> function: eigensolve< " + gmshfem::scalar::Name< T_Scalar >::name + " >");
    }
  }
  
  // Form0, AFrequency
  {
    gmshfem::domain::Domain omega("omega");
        
    gmshfem::problem::Formulation< gmshfem::scalar::ComplexPrecision< T_Scalar > > formulation("eigensolve");

    gmshfem::field::Field< gmshfem::scalar::ComplexPrecision< T_Scalar >, gmshfem::field::Form::Form0 > u("u", omega, gmshfem::field::FunctionSpaceTypeForm0::Lagrange);
    u.addConstraint(omega.getBoundary(), 0.);
        
    formulation.integral(dt2_dof(u), tf(u), omega, "Gauss5");
    formulation.integral(dt_dof(u), tf(u), omega, "Gauss5");
    formulation.integral(grad(dof(u)), grad(tf(u)), omega, "Gauss5");
    
    formulation.integral(-1., tf(u), omega, "Gauss5");
    formulation.setAngularFrequency(1.);

    formulation.pre();
    formulation.assemble();
    formulation.solve();
    gmshfem::post::save(u, omega, "u");
    
    const gmshfem::scalar::ComplexPrecision< T_Scalar > check = gmshfem::post::evaluate(u, 0.8147, 0.1270, 0.);
    gmshfem::scalar::ComplexPrecision< T_Scalar > expected; // value measured on commit ef5338b17e59be542c56f8e815108e910aa9879c
    if(std::is_same< gmshfem::scalar::Precision< T_Scalar >, double >::value) {
      expected = gmshfem::scalar::ComplexPrecision< T_Scalar >(0.024810279029661257144,0.0010343330161330598592);
    }
    else {
      expected = gmshfem::scalar::ComplexPrecision< T_Scalar >(0.0248196646571,0.00103548483457);
    }
    gmshfem::msg::info << "\t- Check value = " << gmshfem::msg::precision(gmshfem::scalar::PrecisionDigits< T_Scalar >::value+5) << check << ", expected = " << expected << gmshfem::msg::endl;
    if(std::abs(check-expected) > 10 * std::sqrt(gmshfem::scalar::Epsilon< T_Scalar >::value)) {
      throw gmshfem::common::Exception("Error in " + std::string(__FILE__) + ", line " + std::to_string(__LINE__) + " -> function: eigensolve< " + gmshfem::scalar::Name< T_Scalar >::name + " >");
    }
  }
}

template void eigensolve< std::complex< double > >();
template void eigensolve< std::complex< float > >();
template void eigensolve< double >();
template void eigensolve< float >();

template< class T_Scalar >
void nonsymmetricMatrix()
{
  // non-symmetric formution
  {
    gmshfem::problem::Formulation< T_Scalar > formulation("testSymmetry");

    gmshfem::domain::Domain line(1, 1);
    gmshfem::field::Field< T_Scalar, gmshfem::field::Form::Form0 > u("u", line, gmshfem::field::FunctionSpaceTypeForm0::Lagrange);
    gmshfem::field::Field< T_Scalar, gmshfem::field::Form::Form0 > v("v", line, gmshfem::field::FunctionSpaceTypeForm0::Lagrange);
    // u = x   in line
    formulation.integral(dof(u), tf(u), line, "Gauss2");
    formulation.integral(- gmshfem::function::x< T_Scalar >(), tf(u), line, "Gauss2");
    // v = u   in line
    formulation.integral(dof(v), tf(v), line, "Gauss2");
    formulation.integral(-dof(u), tf(v), line, "Gauss2");

    formulation.pre();
    formulation.assemble();
    formulation.solve();
    
    const T_Scalar check = gmshfem::post::integrate(v - u, line, "Gauss2");
    const T_Scalar expected = 0.;
    gmshfem::msg::info << "\t- Check value = " << check << ", expected = " << expected << gmshfem::msg::endl;
    if(std::abs(check-expected) > 10. * gmshfem::scalar::Epsilon< T_Scalar >::value) {
      throw gmshfem::common::Exception("Error in " + std::string(__FILE__) + ", line " + std::to_string(__LINE__) + " -> function: nonsymmetricMatrix< " + gmshfem::scalar::Name< T_Scalar >::name + " >");
    }
  }
  
  // non-symmetric elementary matrix
  {
    gmshfem::problem::Formulation< T_Scalar > formulation("testSymmetry");

    gmshfem::domain::Domain line(1, 1);
    gmshfem::field::Field< T_Scalar, gmshfem::field::Form::Form0 > u("u", line, gmshfem::field::FunctionSpaceTypeForm0::Lagrange);
    gmshfem::field::Field< T_Scalar, gmshfem::field::Form::Form0 > v("v", line, gmshfem::field::FunctionSpaceTypeForm0::HierarchicalH1, 2);
    // u = x   in line
    formulation.integral(dof(u), tf(u), line, "Gauss2");
    formulation.integral(- gmshfem::function::x< T_Scalar >(), tf(u), line, "Gauss2");
    // v = u   in line
    formulation.integral(dof(v), tf(v), line, "Gauss2");
    formulation.integral(-dof(u), tf(v), line, "Gauss2");

    formulation.pre();
    formulation.assemble();
    formulation.solve();
    
    const T_Scalar check = gmshfem::post::integrate(v - u, line, "Gauss2");
    const T_Scalar expected = 0.;
    gmshfem::msg::info << "\t- Check value = " << check << ", expected = " << expected << gmshfem::msg::endl;
    if(std::abs(check-expected) > 10. * gmshfem::scalar::Epsilon< T_Scalar >::value) {
      throw gmshfem::common::Exception("Error in " + std::string(__FILE__) + ", line " + std::to_string(__LINE__) + " -> function: nonsymmetricMatrix< " + gmshfem::scalar::Name< T_Scalar >::name + " >");
    }
  }
}

template void nonsymmetricMatrix< std::complex< double > >();
template void nonsymmetricMatrix< std::complex< float > >();
template void nonsymmetricMatrix< double >();
template void nonsymmetricMatrix< float >();

template< class T_Scalar >
void periodic()
{
  gmshfem::domain::Domain omega[4] = { gmshfem::domain::Domain("omega_1"), gmshfem::domain::Domain("omega_2"), gmshfem::domain::Domain("omega_3"), gmshfem::domain::Domain("omega_4") };
  gmshfem::domain::Domain bottom[4] = { gmshfem::domain::Domain("gammaBottom_1"), gmshfem::domain::Domain("gammaBottom_2"), gmshfem::domain::Domain("gammaBottom_3"), gmshfem::domain::Domain("gammaBottom_4") };
  gmshfem::domain::Domain right[4] = { gmshfem::domain::Domain("gammaRight_1"), gmshfem::domain::Domain("gammaRight_2"), gmshfem::domain::Domain("gammaRight_3"), gmshfem::domain::Domain("gammaRight_4") };
  gmshfem::domain::Domain top[4] = { gmshfem::domain::Domain("gammaTop_1"), gmshfem::domain::Domain("gammaTop_2"), gmshfem::domain::Domain("gammaTop_3"), gmshfem::domain::Domain("gammaTop_4") };
  gmshfem::domain::Domain left[4] = { gmshfem::domain::Domain("gammaLeft_1"), gmshfem::domain::Domain("gammaLeft_2"), gmshfem::domain::Domain("gammaLeft_3"), gmshfem::domain::Domain("gammaLeft_4") };
  
  gmshfem::domain::Domain omegaTot = omega[0] | omega[1] | omega[2] | omega[3];
  
  // Non periodic reference formulation
  
  gmshfem::problem::Formulation< T_Scalar > formulationNP("non-periodic");

  gmshfem::field::Field< T_Scalar, gmshfem::field::Form::Form0 > uNP("uNP", omegaTot, gmshfem::field::FunctionSpaceTypeForm0::Lagrange);
  uNP.addConstraint(bottom[2] | right[3] | top[1] | left[0], 1.);
  
  formulationNP.integral(grad(dof(uNP)), grad(tf(uNP)), omegaTot, "Gauss4");
  formulationNP.integral(1., tf(uNP), omegaTot, "Gauss4");

  formulationNP.pre();
  
  formulationNP.assemble();
  formulationNP.solve();
  
  // Periodic reference formulation (the lower right square)
  
  gmshfem::domain::PeriodicLink link("gammaTop_3");
  
  gmshfem::problem::Formulation< T_Scalar > formulationP("periodic");

  gmshfem::field::Field< T_Scalar, gmshfem::field::Form::Form0 > uP("uP", omega[2], gmshfem::field::FunctionSpaceTypeForm0::Lagrange);
  uP.addConstraint(bottom[2], 1.);
  uP.addPeriodicConstraint(link, 1.);
  
  formulationP.integral(grad(dof(uP)), grad(tf(uP)), omega[2], "Gauss4");
  formulationP.integral(1., tf(uP), omega[2], "Gauss4");

  formulationP.pre();
  
  formulationP.assemble();
  formulationP.solve();
  
  // Check solution
  
  const T_Scalar check = gmshfem::post::integrate(gmshfem::function::pow(gmshfem::function::norm(uNP-uP),2), omega[2], "Gauss4");
  const T_Scalar expected = 0.;
  gmshfem::msg::info << "\t- Check value = " << check << ", expected = " << expected << gmshfem::msg::endl;
  if(std::abs(check-expected) > gmshfem::scalar::Epsilon< T_Scalar >::value) {
    throw gmshfem::common::Exception("Error in " + std::string(__FILE__) + ", line " + std::to_string(__LINE__) + " -> function: nonsymmetricMatrix< " + gmshfem::scalar::Name< T_Scalar >::name + " >");
  }
}

template void periodic< std::complex< double > >();
template void periodic< std::complex< float > >();
template void periodic< double >();
template void periodic< float >();
