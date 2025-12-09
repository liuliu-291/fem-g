// GmshFEM - Copyright (C) 2019-2021 Université de Liège
//
// See the LICENSE.txt file for license information. Please report all
// issues on https://gitlab.onelab.info/gmsh/fem/issues

#include <Formulation.h>
#include <Function.h>
#include <Post.h>

#include "functions.h"
#include "geo.h"


using gmshfem::equation::dof;
using gmshfem::equation::tf;

void functions(gmshfem::common::GmshFem &gmshfem, unsigned int &numTest)
{
  gmshfem::msg::info << ++numTest << ") Functions" << gmshfem::msg::endl;

  Geo3D::tetrahedra();
  
  try {
    gmshfem::msg::info << "Trigonometric functions:" << gmshfem::msg::endl;
    trigo< std::complex< double > >();
    trigo< std::complex< float > >();
    trigo< double >();
    trigo< float >();
  } catch(...) { throw; }

  try {
    gmshfem::msg::info << "Hyperbolic functions:" << gmshfem::msg::endl;
    hyperbolic< std::complex< double > >();
    hyperbolic< std::complex< float > >();
    hyperbolic< double >();
    hyperbolic< float >();
  } catch(...) { throw; }

  try {
    gmshfem::msg::info << "Logarithm functions:" << gmshfem::msg::endl;
    logarithm< std::complex< double > >();
    logarithm< std::complex< float > >();
    logarithm< double >();
    logarithm< float >();
  } catch(...) { throw; }

  try {
    gmshfem::msg::info << "Root functions:" << gmshfem::msg::endl;
    root< std::complex< double > >();
    root< std::complex< float > >();
    root< double >();
    root< float >();
  } catch(...) { throw; }
  
  try {
    gmshfem::msg::info << "Bessel functions:" << gmshfem::msg::endl;
    bessel< std::complex< double > >();
    bessel< std::complex< float > >();
    bessel< double >();
    bessel< float >();
  } catch(...) { throw; }

  try {
    gmshfem::msg::info << "Scalar functions:" << gmshfem::msg::endl;
    scalarFunction< std::complex< double > >();
    scalarFunction< std::complex< float > >();
    scalarFunction< double >();
    scalarFunction< float >();
  } catch(...) { throw; }
  
  try {
    gmshfem::msg::info << "Vector functions:" << gmshfem::msg::endl;
    vectorFunction< std::complex< double > >();
    vectorFunction< std::complex< float > >();
    vectorFunction< double >();
    vectorFunction< float >();
  } catch(...) { throw; }

  try {
    gmshfem::msg::info << "Tensor functions:" << gmshfem::msg::endl;
    tensorFunction< std::complex< double > >();
    tensorFunction< std::complex< float > >();
    tensorFunction< double >();
    tensorFunction< float >();
  } catch(...) { throw; }
  
  try {
    gmshfem::msg::info << "Piecewise functions:" << gmshfem::msg::endl;
    piecewiseFunction< std::complex< double > >();
    piecewiseFunction< std::complex< float > >();
    piecewiseFunction< double >();
    piecewiseFunction< float >();
  } catch(...) { throw; }
  
  try {
    gmshfem::msg::info << "Normal and tangent functions:" << gmshfem::msg::endl;
    normalAndTangentFunction< std::complex< double > >();
    normalAndTangentFunction< std::complex< float > >();
    normalAndTangentFunction< double >();
    normalAndTangentFunction< float >();
  } catch(...) { throw; }
  
  try {
    gmshfem::msg::info << "Interpolation functions:" << gmshfem::msg::endl;
    interpolation< std::complex< double > >();
    interpolation< std::complex< float > >();
    interpolation< double >();
    interpolation< float >();
  } catch(...) { throw; }
  
  try {
    gmshfem::msg::info << "Norm functions:" << gmshfem::msg::endl;
    norm< std::complex< double > >();
    norm< std::complex< float > >();
    norm< double >();
    norm< float >();
  } catch(...) { throw; }
  
  try {
    gmshfem::msg::info << "Complex functions:" << gmshfem::msg::endl;
    complexFunction< std::complex< double > >();
    complexFunction< std::complex< float > >();
  } catch(...) { throw; }
  
  try {
    gmshfem::msg::info << "Coordinates functions:" << gmshfem::msg::endl;
    coordinates< std::complex< double > >();
    coordinates< std::complex< float > >();
    coordinates< double >();
    coordinates< float >();
  } catch(...) { throw; }
  
  try {
    gmshfem::msg::info << "Other functions:" << gmshfem::msg::endl;
    others< std::complex< double > >();
    others< std::complex< float > >();
    others< double >();
    others< float >();
  } catch(...) { throw; }
  
  removeGeo();
  
  return;
}

template< class T_Scalar >
void trigo()
{
  // cos^2(x) + sin^2(x)
  {
    gmshfem::function::ScalarFunction< T_Scalar > test = (gmshfem::function::pow(gmshfem::function::cos(gmshfem::function::x< T_Scalar >()), 2) + gmshfem::function::sin(gmshfem::function::x< T_Scalar >()) * gmshfem::function::sin(gmshfem::function::x< T_Scalar >())) / 2.;
    
    const T_Scalar check = gmshfem::post::evaluate(test, 0.8147, 0.9058, 0.2785);
    const T_Scalar expected = 0.5;
    gmshfem::msg::info << "\t- Check value = " << check << ", expected = " << expected << gmshfem::msg::endl;
    if(std::abs(check-expected) > 1. * gmshfem::scalar::Epsilon< T_Scalar >::value) {
      throw gmshfem::common::Exception("Error in " + std::string(__FILE__) + ", line " + std::to_string(__LINE__) + " -> function: trigo< " + gmshfem::scalar::Name< T_Scalar >::name + " >");
    }
  }
  
  // sin(y)/cos(y) - 2 * tan(y)
  {
    gmshfem::function::ScalarFunction< T_Scalar > test = gmshfem::function::sin(gmshfem::function::y< T_Scalar >())/gmshfem::function::cos(gmshfem::function::y< T_Scalar >()) - 2. * gmshfem::function::tan(gmshfem::function::y< T_Scalar >());
    
    const T_Scalar check = gmshfem::post::evaluate(test, 0.1270, 0.9134, 0.9649);
    const T_Scalar expected = -std::tan(0.9134);
    gmshfem::msg::info << "\t- Check value = " << check << ", expected = " << expected << gmshfem::msg::endl;
    if(std::abs(check-expected) > 2. * gmshfem::scalar::Epsilon< T_Scalar >::value) {
      throw gmshfem::common::Exception("Error in " + std::string(__FILE__) + ", line " + std::to_string(__LINE__) + " -> function: trigo< " + gmshfem::scalar::Name< T_Scalar >::name + " >");
    }
  }
}

template void trigo< std::complex< double > >();
template void trigo< std::complex< float > >();
template void trigo< double >();
template void trigo< float >();

template< class T_Scalar >
void hyperbolic()
{
  // cosh^2(x) - sinh^2(x)
  {
    gmshfem::function::ScalarFunction< T_Scalar > test = (gmshfem::function::pow(gmshfem::function::cosh(gmshfem::function::x< T_Scalar >()), 2) - gmshfem::function::sinh(gmshfem::function::x< T_Scalar >()) * gmshfem::function::sinh(gmshfem::function::x< T_Scalar >())) / 2.;
    
    const T_Scalar check = gmshfem::post::evaluate(test, 0.1576, 0.9706, 0.9572);
    const T_Scalar expected = 0.5;
    gmshfem::msg::info << "\t- Check value = " << check << ", expected = " << expected << gmshfem::msg::endl;
    if(std::abs(check-expected) > 1. * gmshfem::scalar::Epsilon< T_Scalar >::value) {
      throw gmshfem::common::Exception("Error in " + std::string(__FILE__) + ", line " + std::to_string(__LINE__) + " -> function: hyperbolic< " + gmshfem::scalar::Name< T_Scalar >::name + " >");
    }
  }
  
  // sinh(y)/cosh(y) - 2 * tanh(y)
  {
    gmshfem::function::ScalarFunction< T_Scalar > test = gmshfem::function::sinh(gmshfem::function::y< T_Scalar >())/gmshfem::function::cosh(gmshfem::function::y< T_Scalar >()) - 2. * gmshfem::function::tanh(gmshfem::function::y< T_Scalar >());
    
    const T_Scalar check = gmshfem::post::evaluate(test, 0.4854, 0.8003, 0.1419);
    const T_Scalar expected = -std::tanh(0.8003);
    gmshfem::msg::info << "\t- Check value = " << check << ", expected = " << expected << gmshfem::msg::endl;
    if(std::abs(check-expected) > 2. * gmshfem::scalar::Epsilon< T_Scalar >::value) {
      throw gmshfem::common::Exception("Error in " + std::string(__FILE__) + ", line " + std::to_string(__LINE__) + " -> function: hyperbolic< " + gmshfem::scalar::Name< T_Scalar >::name + " >");
    }
  }
  
  // asinh(sinh(r2d))
  {
    gmshfem::function::ScalarFunction< T_Scalar > test = gmshfem::function::asinh(gmshfem::function::sinh(gmshfem::function::r2d< T_Scalar >()));
    
    const T_Scalar check = gmshfem::post::evaluate(test, 0.4218, 0.9157, 0.7922);
    const T_Scalar expected = std::sqrt(0.4218 * 0.4218 + 0.9157 * 0.9157);
    gmshfem::msg::info << "\t- Check value = " << check << ", expected = " << expected << gmshfem::msg::endl;
    if(std::abs(check-expected) > 2. * gmshfem::scalar::Epsilon< T_Scalar >::value) {
      throw gmshfem::common::Exception("Error in " + std::string(__FILE__) + ", line " + std::to_string(__LINE__) + " -> function: hyperbolic< " + gmshfem::scalar::Name< T_Scalar >::name + " >");
    }
  }
  
  // acosh(cosh(r3d))
  {
    gmshfem::function::ScalarFunction< T_Scalar > test = gmshfem::function::acosh(gmshfem::function::cosh(gmshfem::function::r3d< T_Scalar >()));
      
    const T_Scalar check = gmshfem::post::evaluate(test, 0.9595, 0.6557, 0.0357);
    const T_Scalar expected = std::sqrt(0.9595 * 0.9595 + 0.6557 * 0.6557 + 0.0357 * 0.0357);
    gmshfem::msg::info << "\t- Check value = " << check << ", expected = " << expected << gmshfem::msg::endl;
    if(std::abs(check-expected) > 2. * gmshfem::scalar::Epsilon< T_Scalar >::value) {
      throw gmshfem::common::Exception("Error in " + std::string(__FILE__) + ", line " + std::to_string(__LINE__) + " -> function: hyperbolic< " + gmshfem::scalar::Name< T_Scalar >::name + " >");
    }
  }
}

template void hyperbolic< std::complex< double > >();
template void hyperbolic< std::complex< float > >();
template void hyperbolic< double >();
template void hyperbolic< float >();

template< class T_Scalar >
void logarithm()
{
  // log(x * y) - log(x) - 2 * log(y)
  {
    gmshfem::function::ScalarFunction< T_Scalar > test = gmshfem::function::log(gmshfem::function::x< T_Scalar >() * gmshfem::function::y< T_Scalar >()) - gmshfem::function::log(gmshfem::function::x< T_Scalar >()) - 2. * gmshfem::function::log(gmshfem::function::y< T_Scalar >());
    
    const T_Scalar check = gmshfem::post::evaluate(test, 0.8491, 0.9340, 0.6787);
    const T_Scalar expected = -std::log10(0.9340);
    gmshfem::msg::info << "\t- Check value = " << check << ", expected = " << expected << gmshfem::msg::endl;
    if(std::abs(check-expected) > 2. * gmshfem::scalar::Epsilon< T_Scalar >::value) {
      throw gmshfem::common::Exception("Error in " + std::string(__FILE__) + ", line " + std::to_string(__LINE__) + " -> function: logarithm< " + gmshfem::scalar::Name< T_Scalar >::name + " >");
    }
  }
  
  // ln(x / y) - ln(x) + 2 * ln(y)
  {
    gmshfem::function::ScalarFunction< T_Scalar > test = gmshfem::function::ln(gmshfem::function::x< T_Scalar >() / gmshfem::function::y< T_Scalar >()) - gmshfem::function::ln(gmshfem::function::x< T_Scalar >()) + 2. * gmshfem::function::ln(gmshfem::function::y< T_Scalar >());
    
    const T_Scalar check = gmshfem::post::evaluate(test, 0.7577, 0.7431, 0.3922);
    const T_Scalar expected = std::log(0.7431);
    gmshfem::msg::info << "\t- Check value = " << check << ", expected = " << expected << gmshfem::msg::endl;
    if(std::abs(check-expected) > 2. * gmshfem::scalar::Epsilon< T_Scalar >::value) {
      throw gmshfem::common::Exception("Error in " + std::string(__FILE__) + ", line " + std::to_string(__LINE__) + " -> function: logarithm< " + gmshfem::scalar::Name< T_Scalar >::name + " >");
    }
  }
  
  // ln(exp(r2d))
  {
    gmshfem::function::ScalarFunction< T_Scalar > test = gmshfem::function::ln(gmshfem::function::exp(gmshfem::function::r2d< T_Scalar >()));
    
    const T_Scalar check = gmshfem::post::evaluate(test, 0.6555, 0.1712, 0.7060);
    const T_Scalar expected = std::sqrt(0.6555 * 0.6555 + 0.1712 * 0.1712);
    gmshfem::msg::info << "\t- Check value = " << check << ", expected = " << expected << gmshfem::msg::endl;
    if(std::abs(check-expected) > 2. * gmshfem::scalar::Epsilon< T_Scalar >::value) {
      throw gmshfem::common::Exception("Error in " + std::string(__FILE__) + ", line " + std::to_string(__LINE__) + " -> function: logarithm< " + gmshfem::scalar::Name< T_Scalar >::name + " >");
    }
  }
}

template void logarithm< std::complex< double > >();
template void logarithm< std::complex< float > >();
template void logarithm< double >();
template void logarithm< float >();

template< class T_Scalar >
void root()
{
  // sqrt(x^2)
  {
    gmshfem::function::ScalarFunction< T_Scalar > test = gmshfem::function::sqrt(gmshfem::function::pow(gmshfem::function::x< T_Scalar >(), 2));
    
    const T_Scalar check = gmshfem::post::evaluate(test, 0.0318, 0.2769, 0.0462);
    const T_Scalar expected = 0.0318;
    gmshfem::msg::info << "\t- Check value = " << check << ", expected = " << expected << gmshfem::msg::endl;
    if(std::abs(check-expected) > 2. * gmshfem::scalar::Epsilon< T_Scalar >::value) {
      throw gmshfem::common::Exception("Error in " + std::string(__FILE__) + ", line " + std::to_string(__LINE__) + " -> function: root< " + gmshfem::scalar::Name< T_Scalar >::name + " >");
    }
  }
  
  // cbrt(z^3)
  {
    gmshfem::function::ScalarFunction< T_Scalar > test = gmshfem::function::cbrt(gmshfem::function::pow(gmshfem::function::z< T_Scalar >(), 3));
    
    const T_Scalar check = gmshfem::post::evaluate(test, 0.0971, 0.8235, 0.6948);
    const T_Scalar expected = 0.6948;
    gmshfem::msg::info << "\t- Check value = " << check << ", expected = " << expected << gmshfem::msg::endl;
    if(std::abs(check-expected) > 2. * gmshfem::scalar::Epsilon< T_Scalar >::value) {
      throw gmshfem::common::Exception("Error in " + std::string(__FILE__) + ", line " + std::to_string(__LINE__) + " -> function: root< " + gmshfem::scalar::Name< T_Scalar >::name + " >");
    }
  }
  
  // sqrt(y^4)
  {
    gmshfem::function::ScalarFunction< T_Scalar > test = gmshfem::function::sqrt(gmshfem::function::pow(gmshfem::function::y< T_Scalar >(), 4));
    
    const T_Scalar check = gmshfem::post::evaluate(test, 0.3171, 0.9502, 0.0344);
    const T_Scalar expected = 0.9502 * 0.9502;
    gmshfem::msg::info << "\t- Check value = " << check << ", expected = " << expected << gmshfem::msg::endl;
    if(std::abs(check-expected) > 2. * gmshfem::scalar::Epsilon< T_Scalar >::value) {
      throw gmshfem::common::Exception("Error in " + std::string(__FILE__) + ", line " + std::to_string(__LINE__) + " -> function: root< " + gmshfem::scalar::Name< T_Scalar >::name + " >");
    }
  }
}

template void root< std::complex< double > >();
template void root< std::complex< float > >();
template void root< double >();
template void root< float >();

template< class T_Scalar >
void bessel()
{
  // cylBesselJ
  {
    gmshfem::function::ScalarFunction< T_Scalar > r = gmshfem::function::r2d< T_Scalar >();
    gmshfem::function::ScalarFunction< T_Scalar > testCheck = gmshfem::function::cylBesselJ< T_Scalar >(3, r) + gmshfem::function::cylBesselJ< T_Scalar >(1, r);
    gmshfem::function::ScalarFunction< T_Scalar > testExpected = 2. * 2. * gmshfem::function::cylBesselJ< T_Scalar >(2, r) / r;
    
    const T_Scalar check = gmshfem::post::evaluate(testCheck, 0.4675, 0.9174, 0.1474);
    const T_Scalar expected = gmshfem::post::evaluate(testExpected, 0.4675, 0.9174, 0.1474);
    gmshfem::msg::info << gmshfem::msg::precision(18) << "\t- Check value = " << check << ", expected = " << expected << gmshfem::msg::endl;
    if(std::abs(check-expected) > 10. * gmshfem::scalar::Epsilon< T_Scalar >::value) {
      throw gmshfem::common::Exception("Error in " + std::string(__FILE__) + ", line " + std::to_string(__LINE__) + " -> function: bessel< " + gmshfem::scalar::Name< T_Scalar >::name + " >");
    }
  }
  
  // cylNeumann
  {
    gmshfem::function::ScalarFunction< T_Scalar > r = gmshfem::function::r2d< T_Scalar >();
    gmshfem::function::ScalarFunction< T_Scalar > testCheck = gmshfem::function::cylNeumann< T_Scalar >(3, r) + gmshfem::function::cylNeumann< T_Scalar >(1, r);
    gmshfem::function::ScalarFunction< T_Scalar > testExpected = 2. * 2. * gmshfem::function::cylNeumann< T_Scalar >(2, r) / r;
    
    const T_Scalar check = gmshfem::post::evaluate(testCheck, 0.4675, 0.9174, 0.1474);
    const T_Scalar expected = gmshfem::post::evaluate(testExpected, 0.4675, 0.9174, 0.1474);
    gmshfem::msg::info << gmshfem::msg::precision(18) << "\t- Check value = " << check << ", expected = " << expected << gmshfem::msg::endl;
    if(std::abs(check-expected) > 10. * gmshfem::scalar::Epsilon< T_Scalar >::value) {
      throw gmshfem::common::Exception("Error in " + std::string(__FILE__) + ", line " + std::to_string(__LINE__) + " -> function: bessel< " + gmshfem::scalar::Name< T_Scalar >::name + " >");
    }
  }
}

template void bessel< std::complex< double > >();
template void bessel< std::complex< float > >();
template void bessel< double >();
template void bessel< float >();

template< class T_Scalar >
void scalarFunction()
{
  // constructor and binary operations
  {
    gmshfem::function::ScalarFunction< T_Scalar > fct1 = 1.; // T_PScalar constructor
    gmshfem::function::ScalarFunction< T_Scalar > fct2 = T_Scalar(2.);  // T_Scalar constructor
    gmshfem::function::ScalarFunction< T_Scalar > fct3 = fct1;  // copy constructor
    gmshfem::function::ScalarFunction< T_Scalar > fct4 = std::move(fct2);  // move constructor
    
    gmshfem::function::ScalarFunction< T_Scalar > fct5 = 1.;
    gmshfem::function::ScalarFunction< T_Scalar > fct6 = 1.;
    fct5 = 3.; // T_PScalar =
    fct6 = T_Scalar(6.); // T_Scalar =
    gmshfem::function::ScalarFunction< T_Scalar > fct7 = 1.;
    gmshfem::function::ScalarFunction< T_Scalar > fct8 = 1.;
    fct7 = fct5; // copy =
    fct8 = std::move(fct6); // move =
    
    gmshfem::function::ScalarFunction< T_Scalar > test; // default constructor
    test = (1. - (2./(+fct4))) + (-((((3. * (1. + (fct3 * fct4) + 1. + (fct8 / fct7) -1. - fct4 - 1.)) + fct8) * 2.) - fct4) / 10.) / (fct5 / 3.);
    
    const gmshfem::function::ExecutionTreeWithDegree< T_Scalar, gmshfem::Degree::Degree0 > *tree = test.getTree();
    test.exportExecutionTree(gmshfem::domain::Domain(), "tree");
    
    const T_Scalar check = gmshfem::post::evaluate(test, 0.4387, 0.3816, 0.7655);
    const T_Scalar expected = -2.2;
    gmshfem::msg::info << "\t- Check value = " << check << ", expected = " << expected << gmshfem::msg::endl;
    if(std::abs(check-expected) > 2. * gmshfem::scalar::Epsilon< T_Scalar >::value || tree == nullptr) {
      throw gmshfem::common::Exception("Error in " + std::string(__FILE__) + ", line " + std::to_string(__LINE__) + " -> function: scalarFunction< " + gmshfem::scalar::Name< T_Scalar >::name + " >");
    }
  }
  
  // field
  {
    gmshfem::problem::Formulation< T_Scalar > formulation("scalarFunction");
    gmshfem::domain::Domain omega("omega");
    gmshfem::field::Field< T_Scalar, gmshfem::field::Form::Form0 > u0("u0", omega, gmshfem::field::FunctionSpaceTypeForm0::Lagrange);
    gmshfem::field::Field< T_Scalar, gmshfem::field::Form::Form0 > u1("u1", omega, gmshfem::field::FunctionSpaceTypeForm0::HierarchicalH1, 2);
    u0.addConstraint(omega, 1.);
    u1.addConstraint(omega, 2.);
    
    gmshfem::field::Field< T_Scalar, gmshfem::field::Form::Form3 > v("v", omega, gmshfem::field::FunctionSpaceTypeForm3::Constant);
    v.addConstraint(omega, 3.);
    
    formulation.integral(dof(u0), tf(u0), omega, "Gauss2");
    formulation.integral(dof(u1), tf(u1), omega, "Gauss2");
    formulation.integral(dof(v), tf(v), omega, "Gauss2");
    
    formulation.pre();
    
    gmshfem::function::ScalarFunction< T_Scalar > fct1 = u0;
    gmshfem::function::ScalarFunction< T_Scalar > fct2 = u1;
    gmshfem::function::ScalarFunction< T_Scalar > fct3 = v;
    
    gmshfem::function::ScalarFunction< T_Scalar > fct4 = 1.;
    gmshfem::function::ScalarFunction< T_Scalar > fct5 = 1.;
    gmshfem::function::ScalarFunction< T_Scalar > fct6 = 1.;
    
    fct4 = u0;
    fct5 = u1;
    fct6 = v;
    
    gmshfem::function::ScalarFunction< T_Scalar > test = (fct4 / fct1) + (fct5 / fct2) - (fct3 / fct6);
    test.exportExecutionTree(gmshfem::domain::Domain(), "tree");
    
    const T_Scalar check = gmshfem::post::integrate(test, omega, "Gauss2");
    const T_Scalar expected = 1.;
    gmshfem::msg::info << "\t- Check value = " << check << ", expected = " << expected << gmshfem::msg::endl;
    if(std::abs(check-expected) > 10. * gmshfem::scalar::Epsilon< T_Scalar >::value) {
      throw gmshfem::common::Exception("Error in " + std::string(__FILE__) + ", line " + std::to_string(__LINE__) + " -> function: scalarFunction< " + gmshfem::scalar::Name< T_Scalar >::name + " >");
    }
  }
}

template void scalarFunction< std::complex< double > >();
template void scalarFunction< std::complex< float > >();
template void scalarFunction< double >();
template void scalarFunction< float >();

template< class T_Scalar >
void vectorFunction()
{
  // constructor and binary operations
  {
    gmshfem::function::VectorFunction< T_Scalar > fct1 = typename gmshfem::MathObject< T_Scalar, gmshfem::Degree::Degree1 >::Object(1., 2., 3.); // T_Scalar constructor
    gmshfem::function::VectorFunction< T_Scalar > fct2 = fct1;  // copy constructor
    gmshfem::function::VectorFunction< T_Scalar > fct3 = std::move(fct1);  // move constructor
    
    gmshfem::function::VectorFunction< T_Scalar > fct4 = typename gmshfem::MathObject< T_Scalar, gmshfem::Degree::Degree1 >::Object(1., 0., 0.);
    fct4 = typename gmshfem::MathObject< T_Scalar, gmshfem::Degree::Degree1 >::Object(3., 2., 1.); // T_PScalar =
    gmshfem::function::VectorFunction< T_Scalar > fct5 = typename gmshfem::MathObject< T_Scalar, gmshfem::Degree::Degree1 >::Object(1., 0., 0.);
    gmshfem::function::VectorFunction< T_Scalar > fct6 = typename gmshfem::MathObject< T_Scalar, gmshfem::Degree::Degree1 >::Object(1., 0., 0.);
    fct5 = fct4; // copy =
    fct6 = std::move(fct4); // move =
    
    gmshfem::function::VectorFunction< T_Scalar > fct7 = fct2 + fct5; //(4, 4, 4)
    gmshfem::function::VectorFunction< T_Scalar > fct8 = typename gmshfem::MathObject< T_Scalar, gmshfem::Degree::Degree1 >::Object(2., 1., 0.) + fct3; //(3, 3, 3)
    gmshfem::function::VectorFunction< T_Scalar > fct9 = fct6 + typename gmshfem::MathObject< T_Scalar, gmshfem::Degree::Degree1 >::Object(0., 1., 2.); //(3, 3, 3)
    
    gmshfem::function::VectorFunction< T_Scalar > fct10 = fct7 - fct5; //(1, 2, 3)
    gmshfem::function::VectorFunction< T_Scalar > fct11 = typename gmshfem::MathObject< T_Scalar, gmshfem::Degree::Degree1 >::Object(2., 1., 0.) - fct8; //(-1, -2, -3)
    gmshfem::function::VectorFunction< T_Scalar > fct12 = fct9 - typename gmshfem::MathObject< T_Scalar, gmshfem::Degree::Degree1 >::Object(0., 1., 2.); //(3, 2, 1)
    
    gmshfem::function::ScalarFunction< T_Scalar > sca1 = 2.;
    gmshfem::function::VectorFunction< T_Scalar > fct13 = fct7 / sca1; //(2, 2, 2)
    gmshfem::function::VectorFunction< T_Scalar > fct13bis = typename gmshfem::MathObject< T_Scalar, gmshfem::Degree::Degree1 >::Object(4., 4., 4.) / sca1; //(2, 2, 2)
    gmshfem::function::VectorFunction< T_Scalar > fct14 = -fct11 / 1.; // (1, 2, 3)
    gmshfem::function::VectorFunction< T_Scalar > fct15 = +fct12 / 1.; // (3, 2, 1)
    
    gmshfem::function::VectorFunction< T_Scalar > fct16 = (fct10 % (fct13 + fct13bis)/2.) / 1.; //(-2, 4, -2)
    gmshfem::function::VectorFunction< T_Scalar > fct17 = typename gmshfem::MathObject< T_Scalar, gmshfem::Degree::Degree1 >::Object(2., 1., 0.) % fct14; // (3, -6, 3)
    gmshfem::function::VectorFunction< T_Scalar > fct18 = fct15 % typename gmshfem::MathObject< T_Scalar, gmshfem::Degree::Degree1 >::Object(2., 1., 0.); // (-1, 2, -1)
    
    gmshfem::function::ScalarFunction< T_Scalar > sca2 = (fct16 * fct17) + (typename gmshfem::MathObject< T_Scalar, gmshfem::Degree::Degree1 >::Object(1., 1., 1.) * fct18) + (fct18 * typename gmshfem::MathObject< T_Scalar, gmshfem::Degree::Degree1 >::Object(1., 2., 3.)); // -36
    
    gmshfem::function::TensorFunction< T_Scalar > ten1 = gmshfem::function::tensorDiag< T_Scalar >(1., 1., 1.);
    typename gmshfem::MathObject< T_Scalar, gmshfem::Degree::Degree2 >::Object m;
    m << 1., 0., 0., 0., 1., 0., 0., 0., 1.;
    
    gmshfem::function::VectorFunction< T_Scalar > test; // default constructor
    test = sca2 * fct16 + 1. * fct17 - sca2 * typename gmshfem::MathObject< T_Scalar, gmshfem::Degree::Degree1 >::Object(1., 1., 1.) +
           fct17 * sca2 - fct18 * 1. + typename gmshfem::MathObject< T_Scalar, gmshfem::Degree::Degree1 >::Object(1., 1., 1.) * sca2 +
           ten1 * fct16 + m * fct17 - ten1 * typename gmshfem::MathObject< T_Scalar, gmshfem::Degree::Degree1 >::Object(1., 1., 1.) +
           fct17 * ten1 - fct18 * m + typename gmshfem::MathObject< T_Scalar, gmshfem::Degree::Degree1 >::Object(1., 1., 1.) * ten1;
           // (111, -114, 111) + (-143, 178, -143) + (0, -3, 0) + (5, -7, 5) = (-27, 54, -27)
    
    const gmshfem::function::ExecutionTreeWithDegree< T_Scalar, gmshfem::Degree::Degree1 > *tree = test.getTree();
    test.exportExecutionTree(gmshfem::domain::Domain(), "tree");
    
    const typename gmshfem::MathObject< T_Scalar, gmshfem::Degree::Degree1 >::Object check = gmshfem::post::evaluate(test, 0.4456, 0.6463, 0.7094);
    const typename gmshfem::MathObject< T_Scalar, gmshfem::Degree::Degree1 >::Object expected = typename gmshfem::MathObject< T_Scalar, gmshfem::Degree::Degree1 >::Object(-27, 54, -27);
    gmshfem::msg::info << "\t- Check value = " << check.norm() << ", expected = " << expected.norm() << gmshfem::msg::endl;
    if((check-expected).norm() > 2. * gmshfem::scalar::Epsilon< T_Scalar >::value || tree == nullptr) {
      throw gmshfem::common::Exception("Error in " + std::string(__FILE__) + ", line " + std::to_string(__LINE__) + " -> function: vectorFunction< " + gmshfem::scalar::Name< T_Scalar >::name + " >");
    }
  }
  
  // field
  {
    gmshfem::problem::Formulation< T_Scalar > formulation("vectorFunction");
    gmshfem::domain::Domain omega("omega");
    gmshfem::field::Field< T_Scalar, gmshfem::field::Form::Form1 > u("u", omega, gmshfem::field::FunctionSpaceTypeForm1::HierarchicalHCurl, 2);
    u.addConstraint(omega, typename gmshfem::MathObject< T_Scalar, gmshfem::Degree::Degree1 >::Object(1., 3., 2.));
    
    formulation.integral(dof(u), tf(u), omega, "Gauss2");
    
    formulation.pre();
    
    gmshfem::function::VectorFunction< T_Scalar > fct1 = u;
    gmshfem::function::VectorFunction< T_Scalar > fct2 = typename gmshfem::MathObject< T_Scalar, gmshfem::Degree::Degree1 >::Object(1., 1., 1.);
    fct2 = u;
    
    gmshfem::function::VectorFunction< T_Scalar > test = (fct1 + fct2) / 2.;
        
    const typename gmshfem::MathObject< T_Scalar, gmshfem::Degree::Degree1 >::Object check = gmshfem::post::integrate( gmshfem::function::VectorFunction< T_Scalar >(typename gmshfem::MathObject< T_Scalar, gmshfem::Degree::Degree1 >::Object(1., 3., 2.)), omega, "Gauss2");
    const typename gmshfem::MathObject< T_Scalar, gmshfem::Degree::Degree1 >::Object expected = typename gmshfem::MathObject< T_Scalar, gmshfem::Degree::Degree1 >::Object(1., 3., 2.);
    gmshfem::msg::info << "\t- Check value = " << check.norm() << ", expected = " << expected.norm() << gmshfem::msg::endl;
    if((check-expected).norm() > 40. * gmshfem::scalar::Epsilon< T_Scalar >::value) {
      throw gmshfem::common::Exception("Error in " + std::string(__FILE__) + ", line " + std::to_string(__LINE__) + " -> function: vectorFunction< " + gmshfem::scalar::Name< T_Scalar >::name + " >");
    }
  }
  {
    gmshfem::problem::Formulation< T_Scalar > formulation("vectorFunction");
    gmshfem::domain::Domain dirichlet("dirichlet");
    gmshfem::field::Field< T_Scalar, gmshfem::field::Form::Form2 > u("u", dirichlet, gmshfem::field::FunctionSpaceTypeForm2::P_HierarchicalHCurl, 2);
    u.addConstraint(dirichlet, typename gmshfem::MathObject< T_Scalar, gmshfem::Degree::Degree1 >::Object(1., 2., 0.));
    
    formulation.integral(dof(u), tf(u), dirichlet, "Gauss2");
    
    formulation.pre();
    
    gmshfem::function::VectorFunction< T_Scalar > fct1 = u;
    gmshfem::function::VectorFunction< T_Scalar > fct2 = typename gmshfem::MathObject< T_Scalar, gmshfem::Degree::Degree1 >::Object(1., 1., 1.);
    fct2 = u;
    
    gmshfem::function::VectorFunction< T_Scalar > test = (fct1 + fct2) / 2.;
    test.exportExecutionTree(gmshfem::domain::Domain(), "tree");
    
    const typename gmshfem::MathObject< T_Scalar, gmshfem::Degree::Degree1 >::Object check = gmshfem::post::integrate(test, dirichlet, "Gauss2");
    const typename gmshfem::MathObject< T_Scalar, gmshfem::Degree::Degree1 >::Object expected = typename gmshfem::MathObject< T_Scalar, gmshfem::Degree::Degree1 >::Object(1., 2., 0.);
    gmshfem::msg::info << "\t- Check value = " << check.norm() << ", expected = " << expected.norm() << gmshfem::msg::endl;
    if((check-expected).norm() > 40. * gmshfem::scalar::Epsilon< T_Scalar >::value) {
      throw gmshfem::common::Exception("Error in " + std::string(__FILE__) + ", line " + std::to_string(__LINE__) + " -> function: vectorFunction< " + gmshfem::scalar::Name< T_Scalar >::name + " >");
    }
  }
  
  // comp
  {
    gmshfem::function::VectorFunction< T_Scalar > test = gmshfem::function::vector(gmshfem::function::x< T_Scalar >(), gmshfem::function::y< T_Scalar >(), gmshfem::function::z< T_Scalar >());
    test.exportExecutionTree(gmshfem::domain::Domain(), "tree");
    
    const T_Scalar check = gmshfem::post::evaluate(xComp(test) + yComp(test) - zComp(test), 0.1966, 0.2511, 0.6160);
    const T_Scalar expected = 0.1966 + 0.2511 - 0.6160;
    gmshfem::msg::info << "\t- Check value = " << check << ", expected = " << expected << gmshfem::msg::endl;
    if(std::abs(check-expected) > 2. * gmshfem::scalar::Epsilon< T_Scalar >::value) {
      throw gmshfem::common::Exception("Error in " + std::string(__FILE__) + ", line " + std::to_string(__LINE__) + " -> function: vectorFunction< " + gmshfem::scalar::Name< T_Scalar >::name + " >");
    }
  }
}

template void vectorFunction< std::complex< double > >();
template void vectorFunction< std::complex< float > >();
template void vectorFunction< double >();
template void vectorFunction< float >();

template< class T_Scalar >
void tensorFunction()
{
  // constructor and binary operations
  {
    gmshfem::domain::Domain omega("omega");
    typename gmshfem::MathObject< T_Scalar, gmshfem::Degree::Degree2 >::Object m;
    m << 1., 2., 3., -2., 3., 1., -3., 1., -2.;
    typename gmshfem::MathObject< T_Scalar, gmshfem::Degree::Degree2 >::Object n;
    n << 3., 1., 2., 0., 1., 3., -2., 0., -1.;
    gmshfem::function::TensorFunction< T_Scalar > fct1 = m; // T_PScalar constructor
    gmshfem::function::TensorFunction< T_Scalar > fct2 = fct1;  // copy constructor
    gmshfem::function::TensorFunction< T_Scalar > fct3 = std::move(fct1);  // move constructor
    
    gmshfem::function::TensorFunction< T_Scalar > fct4 = m;
    fct4 = n; // T_PScalar =
    gmshfem::function::TensorFunction< T_Scalar > fct5 = m;
    gmshfem::function::TensorFunction< T_Scalar > fct6 = m;
    fct5 = fct4; // copy =
    fct6 = std::move(fct4); // move =
    
    gmshfem::function::TensorFunction< T_Scalar > fct7 = fct2 + fct5; //(4, 3, 5; -2, 4, 4; -5, 1, -3)
    gmshfem::function::TensorFunction< T_Scalar > fct8 = n + fct3; //(4, 3, 5; -2, 4, 4; -5, 1, -3)
    gmshfem::function::TensorFunction< T_Scalar > fct9 = fct6 + m; //(4, 3, 5; -2, 4, 4; -5, 1, -3)
    
    gmshfem::function::TensorFunction< T_Scalar > fct10 = fct7 - fct5; //(3, 1, 2; 0, 1, 3; -2, 0, -1)
    gmshfem::function::TensorFunction< T_Scalar > fct11 = n - fct8; //(-1, -2, -3; 2, -3, -1; 3, -1, 2)
    gmshfem::function::TensorFunction< T_Scalar > fct12 = fct9 - m; //(3, 1, 2; 0, 1, 3; -2, 0, -1)

    gmshfem::function::ScalarFunction< T_Scalar > sca1 = 2.;
    gmshfem::function::TensorFunction< T_Scalar > fct13 = fct7 / sca1; //(2, 1.5, 2.5; -1, 2, 2; -2.5, 0.5, -1.5)
    gmshfem::function::TensorFunction< T_Scalar > fct13bis = m / sca1; //(0.5, 1, 1.5; -1, 1.5, 0.5; -1.5, 0.5, -1)
    gmshfem::function::TensorFunction< T_Scalar > fct14 = -fct11 / 1.; //(1, 2, 3; -2, 3, 1; -3, 1, -2)
    gmshfem::function::TensorFunction< T_Scalar > fct15 = +fct12 / 1.; //(3, 1, 2; 0, 1, 3; -2, 0, -1)

    gmshfem::function::TensorFunction< T_Scalar > test; // default constructor
    test = fct13 * fct15 + m * fct13bis - fct15 * n +
           sca1 * fct13 + fct14 * 1. - sca1 * n +
           fct14 * sca1 - 1. * fct15 + m * sca1;
           
    const gmshfem::function::ExecutionTreeWithDegree< T_Scalar, gmshfem::Degree::Degree2 > *tree = test.getTree();
    test.exportExecutionTree(gmshfem::domain::Domain(), "tree");
    
    const typename gmshfem::MathObject< T_Scalar, gmshfem::Degree::Degree2 >::Object check = gmshfem::post::evaluate(test, 0.6551, 0.1626, 0.1190);
    const typename gmshfem::MathObject< T_Scalar, gmshfem::Degree::Degree2 >::Object check2 = gmshfem::post::integrate(test, omega, "Gauss2");
    typename gmshfem::MathObject< T_Scalar, gmshfem::Degree::Degree2 >::Object expected;
    expected <<  -10., 15., 12.5, -18.5, 19., -0.5, -14., 3.5, -11.;
    gmshfem::msg::info << "\t- Check value = (" << check.norm() << " and " << check2.norm() << "), expected = " << expected.norm() << gmshfem::msg::endl;
    if((check-expected).norm() > 2. * gmshfem::scalar::Epsilon< T_Scalar >::value || (check2-expected).norm() > 200. * gmshfem::scalar::Epsilon< T_Scalar >::value || tree == nullptr) {
      throw gmshfem::common::Exception("Error in " + std::string(__FILE__) + ", line " + std::to_string(__LINE__) + " -> function: tensorFunction< " + gmshfem::scalar::Name< T_Scalar >::name + " >");
    }
  }
  
  // inv
  {
    gmshfem::function::TensorFunction< T_Scalar > test = gmshfem::function::tensor< T_Scalar >(gmshfem::function::x< T_Scalar >(), gmshfem::function::y< T_Scalar >(), gmshfem::function::z< T_Scalar >(), 0., 1., 0., 0., 0., 1.)  *gmshfem::function::tensorDiag(gmshfem::function::x< T_Scalar >(), gmshfem::function::y< T_Scalar >(), gmshfem::function::z< T_Scalar >());
    test.exportExecutionTree(gmshfem::domain::Domain(), "tree");
    
    const typename gmshfem::MathObject< T_Scalar, gmshfem::Degree::Degree2 >::Object check = gmshfem::post::evaluate(gmshfem::function::inv(test), 0.4733, 0.3517, 0.8308);
    typename gmshfem::MathObject< T_Scalar, gmshfem::Degree::Degree2 >::Object expected;
    expected << 1./(0.4733*0.4733), -0.3517/(0.4733*0.4733), -0.8308/(0.4733*0.4733), 0., 1./0.3517, 0., 0., 0., 1./0.8308;
    gmshfem::msg::info << "\t- Check value = " << check.norm() << ", expected = " << expected.norm() << gmshfem::msg::endl;
    if((check-expected).norm() > 10. * gmshfem::scalar::Epsilon< T_Scalar >::value) {
      throw gmshfem::common::Exception("Error in " + std::string(__FILE__) + ", line " + std::to_string(__LINE__) + " -> function: tensorFunction< " + gmshfem::scalar::Name< T_Scalar >::name + " >");
    }
  }
}

template void tensorFunction< std::complex< double > >();
template void tensorFunction< std::complex< float > >();
template void tensorFunction< double >();
template void tensorFunction< float >();

template< class T_Scalar >
void piecewiseFunction()
{
  // scalar
  {
    gmshfem::domain::Domain omega("omega");
    gmshfem::domain::Domain dirichlet("dirichlet");
    gmshfem::function::ScalarPiecewiseFunction< T_Scalar > test;
    test.addFunction(1., omega);
    test.addFunction(2., dirichlet);
    
    const T_Scalar check = gmshfem::post::integrate(+test, omega, "Gauss2");
    const T_Scalar expected = 1.;
    gmshfem::msg::info << "\t- Check value = " << check << ", expected = " << expected << gmshfem::msg::endl;
    if(std::abs(check-expected) > 10. * gmshfem::scalar::Epsilon< T_Scalar >::value) {
      throw gmshfem::common::Exception("Error in " + std::string(__FILE__) + ", line " + std::to_string(__LINE__) + " -> function: piecewiseFunction< " + gmshfem::scalar::Name< T_Scalar >::name + " >");
    }
  }
}

template void piecewiseFunction< std::complex< double > >();
template void piecewiseFunction< std::complex< float > >();
template void piecewiseFunction< double >();
template void piecewiseFunction< float >();

template< class T_Scalar >
void normalAndTangentFunction()
{
  // normal 2D
  {
    gmshfem::domain::Domain dirichlet("dirichlet");
    gmshfem::function::VectorFunction< T_Scalar > test = gmshfem::function::normal< T_Scalar >();
    
    const typename gmshfem::MathObject< T_Scalar, gmshfem::Degree::Degree1 >::Object check = gmshfem::post::integrate(test, dirichlet, "Gauss2");
    const typename gmshfem::MathObject< T_Scalar, gmshfem::Degree::Degree1 >::Object expected = typename gmshfem::MathObject< T_Scalar, gmshfem::Degree::Degree1 >::Object(0., 0., 1.);
    gmshfem::msg::info << "\t- Check value = " << check.norm() << ", expected = " << expected.norm() << gmshfem::msg::endl;
    if((check-expected).norm() > 40. * gmshfem::scalar::Epsilon< T_Scalar >::value) {
      throw gmshfem::common::Exception("Error in " + std::string(__FILE__) + ", line " + std::to_string(__LINE__) + " -> function: normalAndTangentFunction< " + gmshfem::scalar::Name< T_Scalar >::name + " >");
    }
  }
  
  // normal 1D
  {
    gmshfem::domain::Domain line("line");
    gmshfem::function::VectorFunction< T_Scalar > test = gmshfem::function::normal< T_Scalar >();
    
    const typename gmshfem::MathObject< T_Scalar, gmshfem::Degree::Degree1 >::Object check = gmshfem::post::integrate(test, line, "Gauss2");
    const typename gmshfem::MathObject< T_Scalar, gmshfem::Degree::Degree1 >::Object expected = typename gmshfem::MathObject< T_Scalar, gmshfem::Degree::Degree1 >::Object(0., -1., 0.);
    gmshfem::msg::info << "\t- Check value = " << check.norm() << ", expected = " << expected.norm() << gmshfem::msg::endl;
    if((check-expected).norm() > 40. * gmshfem::scalar::Epsilon< T_Scalar >::value) {
      throw gmshfem::common::Exception("Error in " + std::string(__FILE__) + ", line " + std::to_string(__LINE__) + " -> function: normalAndTangentFunction< " + gmshfem::scalar::Name< T_Scalar >::name + " >");
    }
  }
  
  // tangents 2D
  {
    gmshfem::domain::Domain dirichlet("dirichlet");
    gmshfem::function::VectorFunction< T_Scalar > test = gmshfem::function::tangent< T_Scalar >(0) % gmshfem::function::tangent< T_Scalar >(1);
    
    const typename gmshfem::MathObject< T_Scalar, gmshfem::Degree::Degree1 >::Object check = gmshfem::post::integrate(test, dirichlet, "Gauss2");
    const typename gmshfem::MathObject< T_Scalar, gmshfem::Degree::Degree1 >::Object expected = typename gmshfem::MathObject< T_Scalar, gmshfem::Degree::Degree1 >::Object(0., 0., 1.);
    gmshfem::msg::info << "\t- Check value = " << check.norm() << ", expected = " << expected.norm() << gmshfem::msg::endl;
    if((check-expected).norm() > 40. * gmshfem::scalar::Epsilon< T_Scalar >::value) {
      throw gmshfem::common::Exception("Error in " + std::string(__FILE__) + ", line " + std::to_string(__LINE__) + " -> function: normalAndTangentFunction< " + gmshfem::scalar::Name< T_Scalar >::name + " >");
    }
  }
  
  // tangents 1D
  {
    gmshfem::domain::Domain line("line");
    gmshfem::function::VectorFunction< T_Scalar > test = gmshfem::function::tangent< T_Scalar >(0);
    
    const typename gmshfem::MathObject< T_Scalar, gmshfem::Degree::Degree1 >::Object check = gmshfem::post::integrate(test, line, "Gauss2");
    const typename gmshfem::MathObject< T_Scalar, gmshfem::Degree::Degree1 >::Object expected = typename gmshfem::MathObject< T_Scalar, gmshfem::Degree::Degree1 >::Object(1., 0., 0.);
    gmshfem::msg::info << "\t- Check value = " << check.norm() << ", expected = " << expected.norm() << gmshfem::msg::endl;
    if((check-expected).norm() > 40. * gmshfem::scalar::Epsilon< T_Scalar >::value) {
      throw gmshfem::common::Exception("Error in " + std::string(__FILE__) + ", line " + std::to_string(__LINE__) + " -> function: normalAndTangentFunction< " + gmshfem::scalar::Name< T_Scalar >::name + " >");
    }
  }
}

template void normalAndTangentFunction< std::complex< double > >();
template void normalAndTangentFunction< std::complex< float > >();
template void normalAndTangentFunction< double >();
template void normalAndTangentFunction< float >();

template< class T_Scalar >
void interpolation()
{
  // linear type 1 scalar
  {
    std::vector< gmshfem::scalar::Precision< T_Scalar > > x(11);
    std::vector< T_Scalar > v(11);
    for(unsigned int i = 0; i <= 10; ++i) {
      x[i] = i/10.;
      v[i] = i/10.;
    }
    gmshfem::function::ScalarFunction< T_Scalar > test = gmshfem::function::linearInterpolation< T_Scalar >(std::move(x), std::move(v));
    
    const typename gmshfem::MathObject< T_Scalar, gmshfem::Degree::Degree0 >::Object check = gmshfem::post::evaluate(test, 0.2551, 0.5060, 0.6991);
    const typename gmshfem::MathObject< T_Scalar, gmshfem::Degree::Degree0 >::Object expected = 0.2551;
    gmshfem::msg::info << "\t- Check value = " << check << ", expected = " << expected << gmshfem::msg::endl;
    if(std::abs(check-expected) > 2. * gmshfem::scalar::Epsilon< T_Scalar >::value) {
      throw gmshfem::common::Exception("Error in " + std::string(__FILE__) + ", line " + std::to_string(__LINE__) + " -> function: interpolation< " + gmshfem::scalar::Name< T_Scalar >::name + " >");
    }
  }

  // linear type 1 vector
  {
    std::vector< gmshfem::scalar::Precision< T_Scalar > > x(11);
    std::vector< typename gmshfem::MathObject< T_Scalar, gmshfem::Degree::Degree1 >::Object > v(11);
    for(unsigned int i = 0; i <= 10; ++i) {
      x[i] = i/10.;
      v[i] = typename gmshfem::MathObject< T_Scalar, gmshfem::Degree::Degree1 >::Object(i/10., -1.*i/10., i/10.);
    }
    gmshfem::function::VectorFunction< T_Scalar > test = gmshfem::function::linearInterpolation< T_Scalar >(x, v);
    
    const typename gmshfem::MathObject< T_Scalar, gmshfem::Degree::Degree1 >::Object check = gmshfem::post::evaluate(test, 0.8407, 0.2543, 0.8143);
    const typename gmshfem::MathObject< T_Scalar, gmshfem::Degree::Degree1 >::Object expected = typename gmshfem::MathObject< T_Scalar, gmshfem::Degree::Degree1 >::Object(0.8407, -0.8407, 0.8407);
    gmshfem::msg::info << "\t- Check value = " << check.norm() << ", expected = " << expected.norm() << gmshfem::msg::endl;
    if((check-expected).norm() > 2. * gmshfem::scalar::Epsilon< T_Scalar >::value) {
      throw gmshfem::common::Exception("Error in " + std::string(__FILE__) + ", line " + std::to_string(__LINE__) + " -> function: interpolation< " + gmshfem::scalar::Name< T_Scalar >::name + " >");
    }
  }
  
  // linear type 1 tensor
  {
    std::vector< gmshfem::scalar::Precision< T_Scalar > > x(11);
    std::vector< typename gmshfem::MathObject< T_Scalar, gmshfem::Degree::Degree2 >::Object > v(11);
    for(unsigned int i = 0; i <= 10; ++i) {
      x[i] = i/10.;
      v[i] << i/10., -1.*i/10., i/10., -1.*i/10., i/10., i/10., i/10., i/10., -1.*i/10.;
    }
    gmshfem::function::TensorFunction< T_Scalar > test = gmshfem::function::linearInterpolation< T_Scalar >(x, v);
    
    const typename gmshfem::MathObject< T_Scalar, gmshfem::Degree::Degree2 >::Object check = gmshfem::post::evaluate(test, 0.2435, 0.9293, 0.3500);
    typename gmshfem::MathObject< T_Scalar, gmshfem::Degree::Degree2 >::Object expected;
    expected << 0.2435, -0.2435, 0.2435, -0.2435, 0.2435, 0.2435, 0.2435, 0.2435, -0.2435;
    gmshfem::msg::info << "\t- Check value = " << check.norm() << ", expected = " << expected.norm() << gmshfem::msg::endl;
    if((check-expected).norm() > 2. * gmshfem::scalar::Epsilon< T_Scalar >::value) {
      throw gmshfem::common::Exception("Error in " + std::string(__FILE__) + ", line " + std::to_string(__LINE__) + " -> function: interpolation< " + gmshfem::scalar::Name< T_Scalar >::name + " >");
    }
  }
  
  // linear type 2 scalar
  {
    std::vector< T_Scalar > v(11);
    for(unsigned int i = 0; i <= 10; ++i) {
      v[i] = i/10.;
    }
    gmshfem::function::ScalarFunction< T_Scalar > test = gmshfem::function::linearInterpolation< T_Scalar >(11, 0., 1., v);
    
    const typename gmshfem::MathObject< T_Scalar, gmshfem::Degree::Degree0 >::Object check = gmshfem::post::evaluate(test, 0.2551, 0.5060, 0.6991);
    const typename gmshfem::MathObject< T_Scalar, gmshfem::Degree::Degree0 >::Object expected = 0.2551;
    gmshfem::msg::info << "\t- Check value = " << check << ", expected = " << expected << gmshfem::msg::endl;
    if(std::abs(check-expected) > 2. * gmshfem::scalar::Epsilon< T_Scalar >::value) {
      throw gmshfem::common::Exception("Error in " + std::string(__FILE__) + ", line " + std::to_string(__LINE__) + " -> function: interpolation< " + gmshfem::scalar::Name< T_Scalar >::name + " >");
    }
  }

  // linear type 2 vector
  {
    std::vector< typename gmshfem::MathObject< T_Scalar, gmshfem::Degree::Degree1 >::Object > v(11);
    for(unsigned int i = 0; i <= 10; ++i) {
      v[i] = typename gmshfem::MathObject< T_Scalar, gmshfem::Degree::Degree1 >::Object(i/10., -1.*i/10., i/10.);
    }
    gmshfem::function::VectorFunction< T_Scalar > test = gmshfem::function::linearInterpolation< T_Scalar >(11, 0., 1., v);
    
    const typename gmshfem::MathObject< T_Scalar, gmshfem::Degree::Degree1 >::Object check = gmshfem::post::evaluate(test, 0.8407, 0.2543, 0.8143);
    const typename gmshfem::MathObject< T_Scalar, gmshfem::Degree::Degree1 >::Object expected = typename gmshfem::MathObject< T_Scalar, gmshfem::Degree::Degree1 >::Object(0.8407, -0.8407, 0.8407);
    gmshfem::msg::info << "\t- Check value = " << check.norm() << ", expected = " << expected.norm() << gmshfem::msg::endl;
    if((check-expected).norm() > 3. * gmshfem::scalar::Epsilon< T_Scalar >::value) {
      throw gmshfem::common::Exception("Error in " + std::string(__FILE__) + ", line " + std::to_string(__LINE__) + " -> function: interpolation< " + gmshfem::scalar::Name< T_Scalar >::name + " >");
    }
  }
  
  // linear type 2 tensor
  {
    std::vector< typename gmshfem::MathObject< T_Scalar, gmshfem::Degree::Degree2 >::Object > v(11);
    for(unsigned int i = 0; i <= 10; ++i) {
      v[i] << i/10., -1.*i/10., i/10., -1.*i/10., i/10., i/10., i/10., i/10., -1.*i/10.;
    }
    gmshfem::function::TensorFunction< T_Scalar > test = gmshfem::function::linearInterpolation< T_Scalar >(11, 0., 1., std::move(v));
    
    const typename gmshfem::MathObject< T_Scalar, gmshfem::Degree::Degree2 >::Object check = gmshfem::post::evaluate(test, 0.2435, 0.9293, 0.3500);
    typename gmshfem::MathObject< T_Scalar, gmshfem::Degree::Degree2 >::Object expected;
    expected << 0.2435, -0.2435, 0.2435, -0.2435, 0.2435, 0.2435, 0.2435, 0.2435, -0.2435;
    gmshfem::msg::info << "\t- Check value = " << check.norm() << ", expected = " << expected.norm() << gmshfem::msg::endl;
    if((check-expected).norm() > 2. * gmshfem::scalar::Epsilon< T_Scalar >::value) {
      throw gmshfem::common::Exception("Error in " + std::string(__FILE__) + ", line " + std::to_string(__LINE__) + " -> function: interpolation< " + gmshfem::scalar::Name< T_Scalar >::name + " >");
    }
  }
  
  // bilinear type 1 scalar
  {
    std::vector< gmshfem::scalar::Precision< T_Scalar > > x(11), y(11);
    std::vector< std::vector< T_Scalar > > v(11, std::vector< T_Scalar >(11));
    for(unsigned int i = 0; i <= 10; ++i) {
      x[i] = i/10.;
      for(unsigned int j = 0; j <= 10; ++j) {
        y[j] = j/10.;
        v[i][j] = i/10. + j/10.;
      }
    }
    gmshfem::function::ScalarFunction< T_Scalar > test = gmshfem::function::bilinearInterpolation< T_Scalar >(x, y, v);
    
    const typename gmshfem::MathObject< T_Scalar, gmshfem::Degree::Degree0 >::Object check = gmshfem::post::evaluate(test, 0.2551, 0.5060, 0.6991);
    const typename gmshfem::MathObject< T_Scalar, gmshfem::Degree::Degree0 >::Object expected = 0.2551 + 0.5060;
    gmshfem::msg::info << "\t- Check value = " << check << ", expected = " << expected << gmshfem::msg::endl;
    if(std::abs(check-expected) > 2. * gmshfem::scalar::Epsilon< T_Scalar >::value) {
      throw gmshfem::common::Exception("Error in " + std::string(__FILE__) + ", line " + std::to_string(__LINE__) + " -> function: interpolation< " + gmshfem::scalar::Name< T_Scalar >::name + " >");
    }
  }

  // bilinear type 1 vector
  {
    std::vector< gmshfem::scalar::Precision< T_Scalar > > x(11), y(11);
    std::vector< std::vector< typename gmshfem::MathObject< T_Scalar, gmshfem::Degree::Degree1 >::Object > > v(11, std::vector< typename gmshfem::MathObject< T_Scalar, gmshfem::Degree::Degree1 >::Object >(11));
    for(unsigned int i = 0; i <= 10; ++i) {
      x[i] = i/10.;
      for(unsigned int j = 0; j <= 10; ++j) {
        y[j] = j/10.;
        v[i][j] = typename gmshfem::MathObject< T_Scalar, gmshfem::Degree::Degree1 >::Object(i/10., -1.*i/10., i/10.) + typename gmshfem::MathObject< T_Scalar, gmshfem::Degree::Degree1 >::Object(j/10., -1.*j/10., j/10.);
      }
    }
    gmshfem::function::VectorFunction< T_Scalar > test = gmshfem::function::bilinearInterpolation< T_Scalar >(std::move(x), std::move(y), std::move(v));
    
    const typename gmshfem::MathObject< T_Scalar, gmshfem::Degree::Degree1 >::Object check = gmshfem::post::evaluate(test, 0.8407, 0.2543, 0.8143);
    const typename gmshfem::MathObject< T_Scalar, gmshfem::Degree::Degree1 >::Object expected = typename gmshfem::MathObject< T_Scalar, gmshfem::Degree::Degree1 >::Object(0.8407, -0.8407, 0.8407) + typename gmshfem::MathObject< T_Scalar, gmshfem::Degree::Degree1 >::Object(0.2543, -0.2543, 0.2543);
    gmshfem::msg::info << "\t- Check value = " << check.norm() << ", expected = " << expected.norm() << gmshfem::msg::endl;
    if((check-expected).norm() > 2. * gmshfem::scalar::Epsilon< T_Scalar >::value) {
      throw gmshfem::common::Exception("Error in " + std::string(__FILE__) + ", line " + std::to_string(__LINE__) + " -> function: interpolation< " + gmshfem::scalar::Name< T_Scalar >::name + " >");
    }
  }
  
  // bilinear type 1 tensor
  {
    std::vector< gmshfem::scalar::Precision< T_Scalar > > x(11), y(11);
    std::vector< std::vector< typename gmshfem::MathObject< T_Scalar, gmshfem::Degree::Degree2 >::Object > > v(11, std::vector< typename gmshfem::MathObject< T_Scalar, gmshfem::Degree::Degree2 >::Object >(11));
    for(unsigned int i = 0; i <= 10; ++i) {
      x[i] = i/10.;
      for(unsigned int j = 0; j <= 10; ++j) {
        y[j] = j/10.;
        v[i][j] << (i+j)/10., -1.*(i+j)/10., (i+j)/10., -1.*(i+j)/10., (i+j)/10., (i+j)/10., (i+j)/10., (i+j)/10., -1.*(i+j)/10.;
      }
    }
    gmshfem::function::TensorFunction< T_Scalar > test = gmshfem::function::bilinearInterpolation< T_Scalar >(x, y, v);
    
    const typename gmshfem::MathObject< T_Scalar, gmshfem::Degree::Degree2 >::Object check = gmshfem::post::evaluate(test, 0.2435, 0.9293, 0.3500);
    typename gmshfem::MathObject< T_Scalar, gmshfem::Degree::Degree2 >::Object expected;
    expected << 0.2435+0.9293, -0.2435-0.9293, 0.2435+0.9293, -0.2435-0.9293, 0.2435+0.9293, 0.2435+0.9293, 0.2435+0.9293, 0.2435+0.9293, -0.2435-0.9293;
    gmshfem::msg::info << "\t- Check value = " << check.norm() << ", expected = " << expected.norm() << gmshfem::msg::endl;
    if((check-expected).norm() > 2. * gmshfem::scalar::Epsilon< T_Scalar >::value) {
      throw gmshfem::common::Exception("Error in " + std::string(__FILE__) + ", line " + std::to_string(__LINE__) + " -> function: interpolation< " + gmshfem::scalar::Name< T_Scalar >::name + " >");
    }
  }
  
  // bilinear type 2 scalar
  {
    std::vector< std::vector< T_Scalar > > v(11, std::vector< T_Scalar >(11));
    for(unsigned int i = 0; i <= 10; ++i) {
      for(unsigned int j = 0; j <= 10; ++j) {
        v[i][j] = i/10. + j/10.;
      }
    }
    gmshfem::function::ScalarFunction< T_Scalar > test = gmshfem::function::bilinearInterpolation< T_Scalar >(11, 11, {0., 0.}, {1., 1.}, std::move(v));
    
    const typename gmshfem::MathObject< T_Scalar, gmshfem::Degree::Degree0 >::Object check = gmshfem::post::evaluate(test, 0.2551, 0.5060, 0.6991);
    const typename gmshfem::MathObject< T_Scalar, gmshfem::Degree::Degree0 >::Object expected = 0.2551 + 0.5060;
    gmshfem::msg::info << "\t- Check value = " << check << ", expected = " << expected << gmshfem::msg::endl;
    if(std::abs(check-expected) > 3. * gmshfem::scalar::Epsilon< T_Scalar >::value) {
      throw gmshfem::common::Exception("Error in " + std::string(__FILE__) + ", line " + std::to_string(__LINE__) + " -> function: interpolation< " + gmshfem::scalar::Name< T_Scalar >::name + " >");
    }
  }

  // bilinear type 2 vector
  {
    std::vector< std::vector< typename gmshfem::MathObject< T_Scalar, gmshfem::Degree::Degree1 >::Object > > v(11, std::vector< typename gmshfem::MathObject< T_Scalar, gmshfem::Degree::Degree1 >::Object >(11));
    for(unsigned int i = 0; i <= 10; ++i) {
      for(unsigned int j = 0; j <= 10; ++j) {
        v[i][j] = typename gmshfem::MathObject< T_Scalar, gmshfem::Degree::Degree1 >::Object(i/10., -1.*i/10., i/10.) + typename gmshfem::MathObject< T_Scalar, gmshfem::Degree::Degree1 >::Object(j/10., -1.*j/10., j/10.);
      }
    }
    gmshfem::function::VectorFunction< T_Scalar > test = gmshfem::function::bilinearInterpolation< T_Scalar >(11, 11, {0., 0.}, {1., 1.}, std::move(v));
    
    const typename gmshfem::MathObject< T_Scalar, gmshfem::Degree::Degree1 >::Object check = gmshfem::post::evaluate(test, 0.8407, 0.2543, 0.8143);
    const typename gmshfem::MathObject< T_Scalar, gmshfem::Degree::Degree1 >::Object expected = typename gmshfem::MathObject< T_Scalar, gmshfem::Degree::Degree1 >::Object(0.8407, -0.8407, 0.8407) + typename gmshfem::MathObject< T_Scalar, gmshfem::Degree::Degree1 >::Object(0.2543, -0.2543, 0.2543);
    gmshfem::msg::info << "\t- Check value = " << check.norm() << ", expected = " << expected.norm() << gmshfem::msg::endl;
    if((check-expected).norm() > 3. * gmshfem::scalar::Epsilon< T_Scalar >::value) {
      throw gmshfem::common::Exception("Error in " + std::string(__FILE__) + ", line " + std::to_string(__LINE__) + " -> function: interpolation< " + gmshfem::scalar::Name< T_Scalar >::name + " >");
    }
  }
  
  // bilinear type 2 tensor
  {
    std::vector< std::vector< typename gmshfem::MathObject< T_Scalar, gmshfem::Degree::Degree2 >::Object > > v(11, std::vector< typename gmshfem::MathObject< T_Scalar, gmshfem::Degree::Degree2 >::Object >(11));
    for(unsigned int i = 0; i <= 10; ++i) {
      for(unsigned int j = 0; j <= 10; ++j) {
        v[i][j] << (i+j)/10., -1.*(i+j)/10., (i+j)/10., -1.*(i+j)/10., (i+j)/10., (i+j)/10., (i+j)/10., (i+j)/10., -1.*(i+j)/10.;
      }
    }
    gmshfem::function::TensorFunction< T_Scalar > test = gmshfem::function::bilinearInterpolation< T_Scalar >(11, 11, {0., 0.}, {1., 1.}, v);
    
    const typename gmshfem::MathObject< T_Scalar, gmshfem::Degree::Degree2 >::Object check = gmshfem::post::evaluate(test, 0.2435, 0.9293, 0.3500);
    typename gmshfem::MathObject< T_Scalar, gmshfem::Degree::Degree2 >::Object expected;
    expected << 0.2435+0.9293, -0.2435-0.9293, 0.2435+0.9293, -0.2435-0.9293, 0.2435+0.9293, 0.2435+0.9293, 0.2435+0.9293, 0.2435+0.9293, -0.2435-0.9293;
    gmshfem::msg::info << "\t- Check value = " << check.norm() << ", expected = " << expected.norm() << gmshfem::msg::endl;
    if((check-expected).norm() > 9. * gmshfem::scalar::Epsilon< T_Scalar >::value) {
      throw gmshfem::common::Exception("Error in " + std::string(__FILE__) + ", line " + std::to_string(__LINE__) + " -> function: interpolation< " + gmshfem::scalar::Name< T_Scalar >::name + " >");
    }
  }
}

template void interpolation< std::complex< double > >();
template void interpolation< std::complex< float > >();
template void interpolation< double >();
template void interpolation< float >();

template< class T_Scalar >
void norm()
{
  // scalar
  {
    gmshfem::function::ScalarFunction< T_Scalar > test = gmshfem::function::norm(gmshfem::function::x< T_Scalar >());
    
    const T_Scalar check = gmshfem::post::evaluate(test, 0.5853, 0.5497, 0.9172);
    const T_Scalar expected = 0.5853;
    gmshfem::msg::info << "\t- Check value = " << check << ", expected = " << expected << gmshfem::msg::endl;
    if(std::abs(check-expected) > 1. * gmshfem::scalar::Epsilon< T_Scalar >::value) {
      throw gmshfem::common::Exception("Error in " + std::string(__FILE__) + ", line " + std::to_string(__LINE__) + " -> function: norm< " + gmshfem::scalar::Name< T_Scalar >::name + " >");
    }
  }
  
  // scalar (abs)
  {
    gmshfem::function::ScalarFunction< T_Scalar > test = gmshfem::function::abs(gmshfem::function::x< T_Scalar >());
    
    const T_Scalar check = gmshfem::post::evaluate(test, 0.2858, 0.7572, 0.7537);
    const T_Scalar expected = 0.2858;
    gmshfem::msg::info << "\t- Check value = " << check << ", expected = " << expected << gmshfem::msg::endl;
    if(std::abs(check-expected) > 1. * gmshfem::scalar::Epsilon< T_Scalar >::value) {
      throw gmshfem::common::Exception("Error in " + std::string(__FILE__) + ", line " + std::to_string(__LINE__) + " -> function: norm< " + gmshfem::scalar::Name< T_Scalar >::name + " >");
    }
  }
  
  // vector, 1
  {
    gmshfem::function::ScalarFunction< T_Scalar > test = gmshfem::function::norm(gmshfem::function::vector(gmshfem::function::x< T_Scalar >(), gmshfem::function::y< T_Scalar >(), gmshfem::function::z< T_Scalar >()), 1);
    
    const T_Scalar check = gmshfem::post::evaluate(test, 0.3804, 0.5678, 0.0759);
    const T_Scalar expected = 0.3804 + 0.5678 + 0.0759;
    gmshfem::msg::info << "\t- Check value = " << check << ", expected = " << expected << gmshfem::msg::endl;
    if(std::abs(check-expected) > 1. * gmshfem::scalar::Epsilon< T_Scalar >::value) {
      throw gmshfem::common::Exception("Error in " + std::string(__FILE__) + ", line " + std::to_string(__LINE__) + " -> function: norm< " + gmshfem::scalar::Name< T_Scalar >::name + " >");
    }
  }
  
  // vector, 2
  {
    gmshfem::function::ScalarFunction< T_Scalar > test = gmshfem::function::norm(gmshfem::function::vector(gmshfem::function::x< T_Scalar >(), gmshfem::function::y< T_Scalar >(), gmshfem::function::z< T_Scalar >()), 2);
    
    const T_Scalar check = gmshfem::post::evaluate(test, 0.3804, 0.5678, 0.0759);
    const T_Scalar expected = std::sqrt(0.3804*0.3804 + 0.5678*0.5678 + 0.0759*0.0759);
    gmshfem::msg::info << "\t- Check value = " << check << ", expected = " << expected << gmshfem::msg::endl;
    if(std::abs(check-expected) > 1. * gmshfem::scalar::Epsilon< T_Scalar >::value) {
      throw gmshfem::common::Exception("Error in " + std::string(__FILE__) + ", line " + std::to_string(__LINE__) + " -> function: norm< " + gmshfem::scalar::Name< T_Scalar >::name + " >");
    }
  }
  
  // vector, inf
  {
    gmshfem::function::ScalarFunction< T_Scalar > test = gmshfem::function::norm(gmshfem::function::vector(gmshfem::function::x< T_Scalar >(), gmshfem::function::y< T_Scalar >(), gmshfem::function::z< T_Scalar >()), -1);
    
    const T_Scalar check = gmshfem::post::evaluate(test, 0.3804, 0.5678, 0.0759);
    const T_Scalar expected = 0.5678;
    gmshfem::msg::info << "\t- Check value = " << check << ", expected = " << expected << gmshfem::msg::endl;
    if(std::abs(check-expected) > 1. * gmshfem::scalar::Epsilon< T_Scalar >::value) {
      throw gmshfem::common::Exception("Error in " + std::string(__FILE__) + ", line " + std::to_string(__LINE__) + " -> function: norm< " + gmshfem::scalar::Name< T_Scalar >::name + " >");
    }
  }
  
  // tensor, 1
  {
    gmshfem::function::ScalarFunction< T_Scalar > test = gmshfem::function::norm(gmshfem::function::tensorDiag(gmshfem::function::x< T_Scalar >(), gmshfem::function::y< T_Scalar >(), gmshfem::function::z< T_Scalar >()), 1);
    
    const T_Scalar check = gmshfem::post::evaluate(test, 0.0540, 0.5308, 0.7792);
    const T_Scalar expected = 0.7792;
    gmshfem::msg::info << "\t- Check value = " << check << ", expected = " << expected << gmshfem::msg::endl;
    if(std::abs(check-expected) > 1. * gmshfem::scalar::Epsilon< T_Scalar >::value) {
      throw gmshfem::common::Exception("Error in " + std::string(__FILE__) + ", line " + std::to_string(__LINE__) + " -> function: norm< " + gmshfem::scalar::Name< T_Scalar >::name + " >");
    }
  }
  
  // tensor, 2
  {
    gmshfem::function::ScalarFunction< T_Scalar > test = gmshfem::function::norm(gmshfem::function::tensorDiag(gmshfem::function::x< T_Scalar >(), gmshfem::function::y< T_Scalar >(), gmshfem::function::z< T_Scalar >()), 2);
    
    const T_Scalar check = gmshfem::post::evaluate(test, 0.0540, 0.5308, 0.7792);
    const T_Scalar expected = std::sqrt(0.0540*0.0540 + 0.5308*0.5308 + 0.7792*0.7792);
    gmshfem::msg::info << "\t- Check value = " << check << ", expected = " << expected << gmshfem::msg::endl;
    if(std::abs(check-expected) > 1. * gmshfem::scalar::Epsilon< T_Scalar >::value) {
      throw gmshfem::common::Exception("Error in " + std::string(__FILE__) + ", line " + std::to_string(__LINE__) + " -> function: norm< " + gmshfem::scalar::Name< T_Scalar >::name + " >");
    }
  }
  
  // tensor, inf
  {
    gmshfem::function::ScalarFunction< T_Scalar > test = gmshfem::function::norm(gmshfem::function::tensorDiag(gmshfem::function::x< T_Scalar >(), gmshfem::function::y< T_Scalar >(), gmshfem::function::z< T_Scalar >()), -1);
    
    const T_Scalar check = gmshfem::post::evaluate(test, 0.0540, 0.5308, 0.7792);
    const T_Scalar expected = 0.7792;
    gmshfem::msg::info << "\t- Check value = " << check << ", expected = " << expected << gmshfem::msg::endl;
    if(std::abs(check-expected) > 1. * gmshfem::scalar::Epsilon< T_Scalar >::value) {
      throw gmshfem::common::Exception("Error in " + std::string(__FILE__) + ", line " + std::to_string(__LINE__) + " -> function: norm< " + gmshfem::scalar::Name< T_Scalar >::name + " >");
    }
  }
}

template void norm< std::complex< double > >();
template void norm< std::complex< float > >();
template void norm< double >();
template void norm< float >();

template< class T_Scalar >
void complexFunction()
{
  // scalar: conj(x) * (real(x) + i imag(x))
  {
    const T_Scalar im(0., 1.);
    gmshfem::function::ScalarFunction< T_Scalar > x = gmshfem::function::x< T_Scalar >() + im * gmshfem::function::y< T_Scalar >();
    gmshfem::function::ScalarFunction< T_Scalar > test = gmshfem::function::conj(x) * (gmshfem::function::real(x) + im * gmshfem::function::imag(x));
    
    const T_Scalar check = gmshfem::post::evaluate(test, 0.9340, 0.1299, 0.5688);
    const T_Scalar expected = 0.9340*0.9340 + 0.1299*0.1299;
    gmshfem::msg::info << "\t- Check value = " << check << ", expected = " << expected << gmshfem::msg::endl;
    if(std::abs(check-expected) > 1. * gmshfem::scalar::Epsilon< T_Scalar >::value) {
      throw gmshfem::common::Exception("Error in " + std::string(__FILE__) + ", line " + std::to_string(__LINE__) + " -> function: complexFunction< " + gmshfem::scalar::Name< T_Scalar >::name + " >");
    }
  }
  
  // vector: conj(x)
  {
    const T_Scalar im(0., 1.);
    gmshfem::function::VectorFunction< T_Scalar > x = gmshfem::function::vector(gmshfem::function::x< T_Scalar >(), gmshfem::function::y< T_Scalar >(), gmshfem::function::z< T_Scalar >()) + im * gmshfem::function::vector(gmshfem::function::y< T_Scalar >(), gmshfem::function::z< T_Scalar >(), gmshfem::function::x< T_Scalar >());
    gmshfem::function::VectorFunction< T_Scalar > test = gmshfem::function::conj(x);
    
    const typename gmshfem::MathObject< T_Scalar, gmshfem::Degree::Degree1 >::Object check = gmshfem::post::evaluate(test, 0.5285, 0.1656, 0.6020);
    const typename gmshfem::MathObject< T_Scalar, gmshfem::Degree::Degree1 >::Object expected = typename gmshfem::MathObject< T_Scalar, gmshfem::Degree::Degree1 >::Object(T_Scalar(0.5285, -0.1656), T_Scalar(0.1656, -0.6020), T_Scalar(0.6020, -0.5285));
    gmshfem::msg::info << "\t- Check value = " << check.norm() << ", expected = " << expected.norm() << gmshfem::msg::endl;
    if((check-expected).norm() > 3. * gmshfem::scalar::Epsilon< T_Scalar >::value) {
      throw gmshfem::common::Exception("Error in " + std::string(__FILE__) + ", line " + std::to_string(__LINE__) + " -> function: complexFunction< " + gmshfem::scalar::Name< T_Scalar >::name + " >");
    }
  }
  
  // tensor: conj(x)
  {
    const T_Scalar im(0., 1.);
    gmshfem::function::TensorFunction< T_Scalar > x = gmshfem::function::tensorDiag(gmshfem::function::x< T_Scalar >(), gmshfem::function::y< T_Scalar >(), gmshfem::function::z< T_Scalar >()) + im * gmshfem::function::tensorDiag(gmshfem::function::y< T_Scalar >(), gmshfem::function::z< T_Scalar >(), gmshfem::function::x< T_Scalar >());
    gmshfem::function::TensorFunction< T_Scalar > test = gmshfem::function::conj(x);
    
    const typename gmshfem::MathObject< T_Scalar, gmshfem::Degree::Degree2 >::Object check = gmshfem::post::evaluate(test, 0.5285, 0.1656, 0.6020);
    typename gmshfem::MathObject< T_Scalar, gmshfem::Degree::Degree2 >::Object expected;
    expected << T_Scalar(0.5285, -0.1656), 0., 0., 0., T_Scalar(0.1656, -0.6020), 0., 0., 0., T_Scalar(0.6020, -0.5285);
    gmshfem::msg::info << "\t- Check value = " << check.norm() << ", expected = " << expected.norm() << gmshfem::msg::endl;
    if((check-expected).norm() > 3. * gmshfem::scalar::Epsilon< T_Scalar >::value) {
      throw gmshfem::common::Exception("Error in " + std::string(__FILE__) + ", line " + std::to_string(__LINE__) + " -> function: complexFunction< " + gmshfem::scalar::Name< T_Scalar >::name + " >");
    }
  }
  
  // realPart(x), imagPart(x)
  {
    const T_Scalar im(0., 1.);
    gmshfem::function::ScalarFunction< T_Scalar > x = T_Scalar(2.) + im;
    gmshfem::function::ScalarFunction< gmshfem::scalar::Precision< T_Scalar > > test1 = gmshfem::function::realPart(x);
    gmshfem::function::ScalarFunction< gmshfem::scalar::Precision< T_Scalar > > test2 = gmshfem::function::imagPart(x);
    
    const gmshfem::scalar::Precision< T_Scalar > check1 = gmshfem::post::evaluate(test1, 0.4752, 0.2759, 0.5983);
    const gmshfem::scalar::Precision< T_Scalar > check2 = gmshfem::post::evaluate(test2, 0.2858, 0.2957, 0.9573);
    const gmshfem::scalar::Precision< T_Scalar > expected1 = 2.;
    const gmshfem::scalar::Precision< T_Scalar > expected2 = 1.;
    gmshfem::msg::info << "\t- Check value = (" << check1 << ", " << check2 << "), expected = (" << expected1 << ", " << expected2 << ")" << gmshfem::msg::endl;
    if(std::abs(check1-expected1) > gmshfem::scalar::Epsilon< T_Scalar >::value || std::abs(check2-expected2) > gmshfem::scalar::Epsilon< T_Scalar >::value) {
      throw gmshfem::common::Exception("Error in " + std::string(__FILE__) + ", line " + std::to_string(__LINE__) + " -> function: complexFunction< " + gmshfem::scalar::Name< T_Scalar >::name + " >");
    }
  }
  
  // complex(x)
  {
    gmshfem::function::ScalarFunction< gmshfem::scalar::Precision< T_Scalar > > x = 2.;
    gmshfem::function::ScalarFunction< T_Scalar > test = gmshfem::function::complex(x);
    
    const T_Scalar check = gmshfem::post::evaluate(test, 0.2494, 0.5925, 0.1845);
    const T_Scalar expected = 2.;
    gmshfem::msg::info << "\t- Check value = " << check << ", expected = " << expected << gmshfem::msg::endl;
    if(std::abs(check-expected) > gmshfem::scalar::Epsilon< T_Scalar >::value) {
      throw gmshfem::common::Exception("Error in " + std::string(__FILE__) + ", line " + std::to_string(__LINE__) + " -> function: complexFunction< " + gmshfem::scalar::Name< T_Scalar >::name + " >");
    }
  }
}

template void complexFunction< std::complex< double > >();
template void complexFunction< std::complex< float > >();

template< class T_Scalar >
void coordinates()
{
  // polar coordinates
  {
    gmshfem::function::ScalarFunction< T_Scalar > r = gmshfem::function::r2d< T_Scalar >();
    gmshfem::function::ScalarFunction< T_Scalar > th = gmshfem::function::theta< T_Scalar >();
    
    gmshfem::function::ScalarFunction< T_Scalar > test1 = r * gmshfem::function::cos(th);
    gmshfem::function::ScalarFunction< T_Scalar > test2 = r * gmshfem::function::sin(th);
    
    const T_Scalar check1 = gmshfem::post::evaluate(test1, 0.7295, 0.2747, 0.1649);
    const T_Scalar check2 = gmshfem::post::evaluate(test2, 0.2747, 0.1649, 0.7295);
    const T_Scalar expected1 = 0.7295;
    const T_Scalar expected2 = 0.1649;
    gmshfem::msg::info << "\t- Check value = (" << check1 << ", " << check2 << ") expected = (" << expected1 << ", " << expected2 << ")" << gmshfem::msg::endl;
    if(std::abs(check1-expected1) > 2. * gmshfem::scalar::Epsilon< T_Scalar >::value || std::abs(check2-expected2) > 2. * gmshfem::scalar::Epsilon< T_Scalar >::value) {
      throw gmshfem::common::Exception("Error in " + std::string(__FILE__) + ", line " + std::to_string(__LINE__) + " -> function: coordinates< " + gmshfem::scalar::Name< T_Scalar >::name + " >");
    }
  }
  
  // polar coordinate
  {
    gmshfem::function::VectorFunction< T_Scalar > x = gmshfem::function::vector(gmshfem::function::x< T_Scalar >(), gmshfem::function::y< T_Scalar >(), gmshfem::function::z< T_Scalar >());
    gmshfem::function::VectorFunction< T_Scalar > test = gmshfem::function::vector< T_Scalar >(gmshfem::function::r2dComp(x), gmshfem::function::angularComp(x), T_Scalar(0.));
        
    const typename gmshfem::MathObject< T_Scalar, gmshfem::Degree::Degree1 >::Object check = gmshfem::post::evaluate(test, 0.7482, 0.4505, 0.0838);
    const typename gmshfem::MathObject< T_Scalar, gmshfem::Degree::Degree1 >::Object expected = typename gmshfem::MathObject< T_Scalar, gmshfem::Degree::Degree1 >::Object(std::sqrt(0.7482*0.7482 + 0.4505*0.4505), std::atan2(0.4505, 0.7482), 0.);
    gmshfem::msg::info << "\t- Check value = " << check.norm() << ", expected = " << expected.norm() << gmshfem::msg::endl;
    if((check-expected).norm() > 1. * gmshfem::scalar::Epsilon< T_Scalar >::value) {
      throw gmshfem::common::Exception("Error in " + std::string(__FILE__) + ", line " + std::to_string(__LINE__) + " -> function: coordinates< " + gmshfem::scalar::Name< T_Scalar >::name + " >");
    }
  }
  
  // domain coordinate
  {
    gmshfem::domain::Domain omega("omega");
    gmshfem::function::ScalarFunction< T_Scalar > test = gmshfem::function::x< T_Scalar >(omega) + gmshfem::function::y< T_Scalar >(omega) + gmshfem::function::z< T_Scalar >(omega);
        
    const T_Scalar check = gmshfem::post::evaluate(test, 0.2290, 0.9133, 0.1524);
    const T_Scalar expected = 0.2290 + 0.9133 + 0.1524;
    gmshfem::msg::info << "\t- Check value = " << check << ", expected = " << expected << gmshfem::msg::endl;
    if(std::abs(check-expected) > 1. * gmshfem::scalar::Epsilon< T_Scalar >::value) {
      throw gmshfem::common::Exception("Error in " + std::string(__FILE__) + ", line " + std::to_string(__LINE__) + " -> function: coordinates< " + gmshfem::scalar::Name< T_Scalar >::name + " >");
    }
  }
  
  // spherical coordinates
  {
    gmshfem::function::ScalarFunction< T_Scalar > r = gmshfem::function::r3d< T_Scalar >();
    gmshfem::function::ScalarFunction< T_Scalar > th = gmshfem::function::theta< T_Scalar >();
    gmshfem::function::ScalarFunction< T_Scalar > ph = gmshfem::function::phi< T_Scalar >();
    
    gmshfem::function::ScalarFunction< T_Scalar > test1 = r * gmshfem::function::cos(th) * gmshfem::function::sin(ph);
    gmshfem::function::ScalarFunction< T_Scalar > test2 = r * gmshfem::function::sin(th) * gmshfem::function::sin(ph);
    gmshfem::function::ScalarFunction< T_Scalar > test3 = r * gmshfem::function::cos(ph);
    
    const T_Scalar check1 = gmshfem::post::evaluate(test1, 0.2749, 0.1856, 0.9746);
    const T_Scalar check2 = gmshfem::post::evaluate(test2, 0.2653, 0.5632, 0.2759);
    const T_Scalar check3 = gmshfem::post::evaluate(test3, 0.4597, 0.2759, 0.5732);
    const T_Scalar expected1 = 0.2749;
    const T_Scalar expected2 = 0.5632;
    const T_Scalar expected3 = 0.5732;
    gmshfem::msg::info << "\t- Check value = (" << check1 << ", " << check2 << ", " << check3 << ") expected = (" << expected1 << ", " << expected2 << ", " << expected3 << ")" << gmshfem::msg::endl;
    if(std::abs(check1-expected1) > 2. * gmshfem::scalar::Epsilon< T_Scalar >::value || std::abs(check2-expected2) > 2. * gmshfem::scalar::Epsilon< T_Scalar >::value || std::abs(check3-expected3) > 2. * gmshfem::scalar::Epsilon< T_Scalar >::value) {
      throw gmshfem::common::Exception("Error in " + std::string(__FILE__) + ", line " + std::to_string(__LINE__) + " -> function: coordinates< " + gmshfem::scalar::Name< T_Scalar >::name + " >");
    }
  }
  
  // change of coordinates
  {
    gmshfem::function::ScalarFunction< T_Scalar > f = gmshfem::function::sin< T_Scalar >(gmshfem::function::x< T_Scalar >());
    
    gmshfem::function::ScalarFunction< T_Scalar > test = gmshfem::function::changeOfCoordinates< T_Scalar >(f, gmshfem::function::y< gmshfem::scalar::Precision< T_Scalar > >(), gmshfem::function::x< gmshfem::scalar::Precision< T_Scalar > >(), gmshfem::function::z< gmshfem::scalar::Precision< T_Scalar > >());
    
    const T_Scalar check = gmshfem::post::evaluate(test, 0.3758, 0.2756, 0.3656);
    const T_Scalar expected = 0.272124343161595; // sin(0.2756)
    gmshfem::msg::info << "\t- Check value = " << check << " expected = " << expected << gmshfem::msg::endl;
    if(std::abs(check-expected) > 2. * gmshfem::scalar::Epsilon< T_Scalar >::value) {
      throw gmshfem::common::Exception("Error in " + std::string(__FILE__) + ", line " + std::to_string(__LINE__) + " -> function: coordinates< " + gmshfem::scalar::Name< T_Scalar >::name + " >");
    }
  }
}

template void coordinates< std::complex< double > >();
template void coordinates< std::complex< float > >();
template void coordinates< double >();
template void coordinates< float >();

template< class T_Scalar >
void others()
{
  // Commutator, Heaviside
  {
    bool button = false;
    gmshfem::function::ScalarFunction< T_Scalar > x = gmshfem::function::commutator(gmshfem::function::real(gmshfem::function::x< T_Scalar >())-0.5+gmshfem::function::imag(gmshfem::function::x< T_Scalar >()), &button);
    gmshfem::function::ScalarFunction< T_Scalar > test = gmshfem::function::heaviside(gmshfem::function::conj(x));
    
    button = true;
    
    const T_Scalar check = gmshfem::post::evaluate(test, 0.7630, 0.6541, 0.6892);
    const T_Scalar expected = 1.;
    gmshfem::msg::info << "\t- Check value = " << check << ", expected = " << expected << gmshfem::msg::endl;
    if(std::abs(check-expected) > 1. * gmshfem::scalar::Epsilon< T_Scalar >::value) {
      throw gmshfem::common::Exception("Error in " + std::string(__FILE__) + ", line " + std::to_string(__LINE__) + " -> function: others< " + gmshfem::scalar::Name< T_Scalar >::name + " >");
    }
  }
}

template void others< std::complex< double > >();
template void others< std::complex< float > >();
template void others< double >();
template void others< float >();
