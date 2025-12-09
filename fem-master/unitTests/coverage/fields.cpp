// GmshFEM - Copyright (C) 2019-2021 Université de Liège
//
// See the LICENSE.txt file for license information. Please report all
// issues on https://gitlab.onelab.info/gmsh/fem/issues

#include <Formulation.h>
#include <Function.h>
#include <Post.h>

#include "fields.h"
#include "geo.h"


using gmshfem::equation::dof;
using gmshfem::equation::tf;
using gmshfem::equation::grad;

using namespace gmshfem::function;

void fields(gmshfem::common::GmshFem &gmshfem, unsigned int &numTest)
{
  gmshfem::msg::info << ++numTest << ") Fields" << gmshfem::msg::endl;

  Geo3D::tetrahedra();
  
  try {
    gmshfem::msg::info << "Trigonometric field functions:" << gmshfem::msg::endl;
    trigoField< std::complex< double > >();
    trigoField< std::complex< float > >();
    trigoField< double >();
    trigoField< float >();
  } catch(...) { throw; }

  try {
    gmshfem::msg::info << "Hyperbolic field functions:" << gmshfem::msg::endl;
    hyperbolicField< std::complex< double > >();
    hyperbolicField< std::complex< float > >();
    hyperbolicField< double >();
    hyperbolicField< float >();
  } catch(...) { throw; }

  try {
    gmshfem::msg::info << "Logarithm field functions:" << gmshfem::msg::endl;
    logarithmField< std::complex< double > >();
    logarithmField< std::complex< float > >();
    logarithmField< double >();
    logarithmField< float >();
  } catch(...) { throw; }

  try {
    gmshfem::msg::info << "Root field functions:" << gmshfem::msg::endl;
    rootField< std::complex< double > >();
    rootField< std::complex< float > >();
    rootField< double >();
    rootField< float >();
  } catch(...) { throw; }
  
  try {
    gmshfem::msg::info << "Norm field functions:" << gmshfem::msg::endl;
    normField< std::complex< double > >();
    normField< std::complex< float > >();
    normField< double >();
    normField< float >();
  } catch(...) { throw; }
  
  try {
    gmshfem::msg::info << "Complex field functions:" << gmshfem::msg::endl;
    complexFieldFunction< std::complex< double > >();
    complexFieldFunction< std::complex< float > >();
  } catch(...) { throw; }
  
  try {
    gmshfem::msg::info << "Derivative field functions:" << gmshfem::msg::endl;
    derivativeFieldFunction< std::complex< double > >();
    derivativeFieldFunction< std::complex< float > >();
    derivativeFieldFunction< double >();
    derivativeFieldFunction< float >();
  } catch(...) { throw; }
  
  removeGeo();
  
  Geo2D::triangle();
  
  try {
    gmshfem::msg::info << "Compound field:" << gmshfem::msg::endl;
    compoundField< std::complex< double > >();
    compoundField< std::complex< float > >();
    compoundField< double >();
    compoundField< float >();
  } catch(...) { throw; }

  removeGeo();
  
  return;
}

template< class T_Scalar >
void trigoField()
{
  gmshfem::domain::Domain omega("omega");
  gmshfem::field::Field< T_Scalar, gmshfem::field::Form::Form0 > fieldForm0("form0", omega, gmshfem::field::FunctionSpaceTypeForm0::Lagrange);
  fieldForm0.addConstraint(omega, 0.14827);
  gmshfem::field::Field< T_Scalar, gmshfem::field::Form::Form3 > fieldForm3("form3", omega, gmshfem::field::FunctionSpaceTypeForm3::Constant);
  fieldForm3.addConstraint(omega, 0.14827);
    
  gmshfem::problem::Formulation< T_Scalar > formulation("fake");

  formulation.integral(dof(fieldForm0), tf(fieldForm0), omega, "Gauss1");
  formulation.integral(dof(fieldForm3), tf(fieldForm3), omega, "Gauss1");

  formulation.pre();
    
  // cos^2(x) + sin^2(x)
  {
    gmshfem::function::ScalarFunction< T_Scalar > test = (gmshfem::function::pow(gmshfem::function::cos(fieldForm0), 2) + gmshfem::function::pow(gmshfem::function::sin(fieldForm3), 2) + gmshfem::function::pow(gmshfem::function::cos(fieldForm3), 2) + gmshfem::function::pow(gmshfem::function::sin(fieldForm0), 2) + fieldForm0 + fieldForm3 - fieldForm0 - fieldForm3) * fieldForm0 * fieldForm3 / fieldForm3 / fieldForm0;
    
    const T_Scalar check = gmshfem::post::evaluate(test, 0.3842, 0.3852, 0.5684);
    const T_Scalar expected = 2.;
    gmshfem::msg::info << "\t- Check value = " << check << ", expected = " << expected << gmshfem::msg::endl;
    if(std::abs(check-expected) > 2. * gmshfem::scalar::Epsilon< T_Scalar >::value) {
      throw gmshfem::common::Exception("Error in " + std::string(__FILE__) + ", line " + std::to_string(__LINE__) + " -> function: trigoField< " + gmshfem::scalar::Name< T_Scalar >::name + " >");
    }
  }
  
  // sin(y)/cos(y) - 2 * tan(y)
  {
    gmshfem::function::ScalarFunction< T_Scalar > test = gmshfem::function::sin(fieldForm0)/gmshfem::function::cos(fieldForm3) - gmshfem::function::tan(fieldForm3) - gmshfem::function::tan(fieldForm0);

    const T_Scalar check = gmshfem::post::evaluate(test, 0.2645, 0.4527, 0.1745);
    const T_Scalar expected = -std::tan(0.14827);
    gmshfem::msg::info << "\t- Check value = " << check << ", expected = " << expected << gmshfem::msg::endl;
    if(std::abs(check-expected) > 2. * gmshfem::scalar::Epsilon< T_Scalar >::value) {
      throw gmshfem::common::Exception("Error in " + std::string(__FILE__) + ", line " + std::to_string(__LINE__) + " -> function: trigoField< " + gmshfem::scalar::Name< T_Scalar >::name + " >");
    }
  }
}

template void trigoField< std::complex< double > >();
template void trigoField< std::complex< float > >();
template void trigoField< double >();
template void trigoField< float >();

template< class T_Scalar >
void hyperbolicField()
{
  gmshfem::domain::Domain omega("omega");
  gmshfem::field::Field< T_Scalar, gmshfem::field::Form::Form0 > fieldForm0("form0", omega, gmshfem::field::FunctionSpaceTypeForm0::Lagrange);
  fieldForm0.addConstraint(omega, 0.14827);
  gmshfem::field::Field< T_Scalar, gmshfem::field::Form::Form3 > fieldForm3("form3", omega, gmshfem::field::FunctionSpaceTypeForm3::Constant);
  fieldForm3.addConstraint(omega, 0.14827);
    
  gmshfem::problem::Formulation< T_Scalar > formulation("fake");

  formulation.integral(dof(fieldForm0), tf(fieldForm0), omega, "Gauss1");
  formulation.integral(dof(fieldForm3), tf(fieldForm3), omega, "Gauss1");

  formulation.pre();
    
  // cosh^2(x) - sinh^2(x)
  {
    gmshfem::function::ScalarFunction< T_Scalar > test = gmshfem::function::pow(gmshfem::function::cosh(fieldForm0), 2) - gmshfem::function::pow(gmshfem::function::sinh(fieldForm0), 2) + gmshfem::function::pow(gmshfem::function::cosh(fieldForm3), 2) - gmshfem::function::pow(gmshfem::function::sinh(fieldForm3), 2);
    
    const T_Scalar check = gmshfem::post::evaluate(test, 0.2895, 0.5835, 0.5927);
    const T_Scalar expected = 2.;
    gmshfem::msg::info << "\t- Check value = " << check << ", expected = " << expected << gmshfem::msg::endl;
    if(std::abs(check-expected) > 2. * gmshfem::scalar::Epsilon< T_Scalar >::value) {
      throw gmshfem::common::Exception("Error in " + std::string(__FILE__) + ", line " + std::to_string(__LINE__) + " -> function: hyperbolicField< " + gmshfem::scalar::Name< T_Scalar >::name + " >");
    }
  }
  
  // sinh(y)/cosh(y) - 2 * tanh(y)
  {
    gmshfem::function::ScalarFunction< T_Scalar > test = gmshfem::function::sinh(fieldForm0)/gmshfem::function::cosh(fieldForm3) - gmshfem::function::tanh(fieldForm3) - gmshfem::function::tanh(fieldForm0);
    
    const T_Scalar check = gmshfem::post::evaluate(test, 0.2756, 0.9653, 0.2658);
    const T_Scalar expected = -std::tanh(0.14827);
    gmshfem::msg::info << "\t- Check value = " << check << ", expected = " << expected << gmshfem::msg::endl;
    if(std::abs(check-expected) > 2. * gmshfem::scalar::Epsilon< T_Scalar >::value) {
      throw gmshfem::common::Exception("Error in " + std::string(__FILE__) + ", line " + std::to_string(__LINE__) + " -> function: hyperbolicField< " + gmshfem::scalar::Name< T_Scalar >::name + " >");
    }
  }
}

template void hyperbolicField< std::complex< double > >();
template void hyperbolicField< std::complex< float > >();
template void hyperbolicField< double >();
template void hyperbolicField< float >();

template< class T_Scalar >
void logarithmField()
{
  gmshfem::domain::Domain omega("omega");
  gmshfem::field::Field< T_Scalar, gmshfem::field::Form::Form0 > fieldForm0("form0", omega, gmshfem::field::FunctionSpaceTypeForm0::Lagrange);
  fieldForm0.addConstraint(omega, 0.25652);
  gmshfem::field::Field< T_Scalar, gmshfem::field::Form::Form3 > fieldForm3("form3", omega, gmshfem::field::FunctionSpaceTypeForm3::Constant);
  fieldForm3.addConstraint(omega, 4.3756);
    
  gmshfem::problem::Formulation< T_Scalar > formulation("fake");
  
  formulation.integral(dof(fieldForm0), tf(fieldForm0), omega, "Gauss1");
  formulation.integral(dof(fieldForm3), tf(fieldForm3), omega, "Gauss1");

  formulation.pre();
  
  // log(x * y) - log(x) - 2 * log(y)
  {
    gmshfem::function::ScalarFunction< T_Scalar > test = gmshfem::function::log(fieldForm0 * fieldForm3) - gmshfem::function::log(fieldForm0) - 2. * gmshfem::function::log(fieldForm3);
    
    const T_Scalar check = gmshfem::post::evaluate(test, 0.2755, 0.5735, 0.9474);
    const T_Scalar expected = -std::log10(4.3756);
    gmshfem::msg::info << "\t- Check value = " << check << ", expected = " << expected << gmshfem::msg::endl;
    if(std::abs(check-expected) > 10. * gmshfem::scalar::Epsilon< T_Scalar >::value) {
      throw gmshfem::common::Exception("Error in " + std::string(__FILE__) + ", line " + std::to_string(__LINE__) + " -> function: logarithmField< " + gmshfem::scalar::Name< T_Scalar >::name + " >");
    }
  }
  
  // ln(x / y) - ln(x) + 2 * ln(y)
  {
    gmshfem::function::ScalarFunction< T_Scalar > test = gmshfem::function::ln(fieldForm0 / fieldForm3) - gmshfem::function::ln(fieldForm0) + 2. * gmshfem::function::ln(fieldForm3);
    
    const T_Scalar check = gmshfem::post::evaluate(test, 0.3756, 0.3628, 0.5936);
    const T_Scalar expected = std::log(4.3756);
    gmshfem::msg::info << "\t- Check value = " << check << ", expected = " << expected << gmshfem::msg::endl;
    if(std::abs(check-expected) > 10. * gmshfem::scalar::Epsilon< T_Scalar >::value) {
      throw gmshfem::common::Exception("Error in " + std::string(__FILE__) + ", line " + std::to_string(__LINE__) + " -> function: logarithmField< " + gmshfem::scalar::Name< T_Scalar >::name + " >");
    }
  }
  
  // ln(exp(x))
  {
    gmshfem::function::ScalarFunction< T_Scalar > test = gmshfem::function::ln(gmshfem::function::exp(fieldForm0)) + gmshfem::function::ln(gmshfem::function::exp(fieldForm3));
    
    const T_Scalar check = gmshfem::post::evaluate(test, 0.3856, 0.8624, 0.1746);
    const T_Scalar expected = 0.25652 + 4.3756;
    gmshfem::msg::info << "\t- Check value = " << check << ", expected = " << expected << gmshfem::msg::endl;
    if(std::abs(check-expected) > 10. * gmshfem::scalar::Epsilon< T_Scalar >::value) {
      throw gmshfem::common::Exception("Error in " + std::string(__FILE__) + ", line " + std::to_string(__LINE__) + " -> function: logarithmField< " + gmshfem::scalar::Name< T_Scalar >::name + " >");
    }
  }
}

template void logarithmField< std::complex< double > >();
template void logarithmField< std::complex< float > >();
template void logarithmField< double >();
template void logarithmField< float >();

template< class T_Scalar >
void rootField()
{
  gmshfem::domain::Domain omega("omega");
  gmshfem::field::Field< T_Scalar, gmshfem::field::Form::Form0 > fieldForm0("form0", omega, gmshfem::field::FunctionSpaceTypeForm0::Lagrange);
  fieldForm0.addConstraint(omega, 36.375);
  gmshfem::field::Field< T_Scalar, gmshfem::field::Form::Form3 > fieldForm3("form3", omega, gmshfem::field::FunctionSpaceTypeForm3::Constant);
  fieldForm3.addConstraint(omega, 4.2568);
    
  gmshfem::problem::Formulation< T_Scalar > formulation("fake");
  
  formulation.integral(dof(fieldForm0), tf(fieldForm0), omega, "Gauss1");
  formulation.integral(dof(fieldForm3), tf(fieldForm3), omega, "Gauss1");

  formulation.pre();
  
  // x^2 + sqrt(x)
  {
    gmshfem::function::ScalarFunction< T_Scalar > test = gmshfem::function::pow(fieldForm0, 2) + gmshfem::function::sqrt(fieldForm0) - gmshfem::function::pow(fieldForm3, 2) - gmshfem::function::sqrt(fieldForm3);
    
    const T_Scalar check = gmshfem::post::evaluate(test, 0.3853, 0.0572, 0.2756);
    const T_Scalar expected = 36.375 * 36.375 + std::sqrt(36.375) - 4.2568 * 4.2568 - std::sqrt(4.2568);
    gmshfem::msg::info << "\t- Check value = " << check << ", expected = " << expected << gmshfem::msg::endl;
    if(std::abs(check-expected) > 1.e5 * gmshfem::scalar::Epsilon< T_Scalar >::value) {
      throw gmshfem::common::Exception("Error in " + std::string(__FILE__) + ", line " + std::to_string(__LINE__) + " -> function: rootField< " + gmshfem::scalar::Name< T_Scalar >::name + " >");
    }
  }
  
  // cbrt(x)
  {
    gmshfem::function::ScalarFunction< T_Scalar > test = gmshfem::function::pow(gmshfem::function::cbrt(fieldForm0), 3) - gmshfem::function::pow(gmshfem::function::cbrt(fieldForm3), 3);
    
    const T_Scalar check = gmshfem::post::evaluate(test, 0.0971, 0.8235, 0.6948);
    const T_Scalar expected = 36.375 - 4.2568;
    gmshfem::msg::info << "\t- Check value = " << check << ", expected = " << expected << gmshfem::msg::endl;
    if(std::abs(check-expected) > 1.e5 * gmshfem::scalar::Epsilon< T_Scalar >::value) {
      throw gmshfem::common::Exception("Error in " + std::string(__FILE__) + ", line " + std::to_string(__LINE__) + " -> function: rootField< " + gmshfem::scalar::Name< T_Scalar >::name + " >");
    }
  }
}

template void rootField< std::complex< double > >();
template void rootField< std::complex< float > >();
template void rootField< double >();
template void rootField< float >();

template< class T_Scalar >
void normField()
{
  gmshfem::domain::Domain omega("omega");
  gmshfem::domain::Domain dirichlet("dirichlet");
  gmshfem::field::Field< T_Scalar, gmshfem::field::Form::Form0 > fieldForm0("form0", omega, gmshfem::field::FunctionSpaceTypeForm0::Lagrange);
  fieldForm0.addConstraint(omega, 0.2756);
  gmshfem::field::Field< T_Scalar, gmshfem::field::Form::Form1 > fieldForm1("form1", omega, gmshfem::field::FunctionSpaceTypeForm1::HierarchicalHCurl, 1);
  fieldForm1.addConstraint(omega, typename gmshfem::MathObject< T_Scalar, gmshfem::Degree::Degree1 >::Object(0.3824, 3.5735, 2.55));
  gmshfem::field::Field< T_Scalar, gmshfem::field::Form::Form2 > fieldForm2("form2", dirichlet, gmshfem::field::FunctionSpaceTypeForm2::P_HierarchicalHCurl, 1);
  fieldForm2.addConstraint(dirichlet, typename gmshfem::MathObject< T_Scalar, gmshfem::Degree::Degree1 >::Object(0.3824, 3.5735, 0.));
  gmshfem::field::Field< T_Scalar, gmshfem::field::Form::Form3 > fieldForm3("form3", omega, gmshfem::field::FunctionSpaceTypeForm3::Constant);
  fieldForm3.addConstraint(omega, 47.265);
    
  gmshfem::problem::Formulation< T_Scalar > formulation("fake");
  
  formulation.integral(dof(fieldForm0), tf(fieldForm0), omega, "Gauss1");
  formulation.integral(dof(fieldForm1), tf(fieldForm1), omega, "Gauss1");
  formulation.integral(dof(fieldForm2), tf(fieldForm2), dirichlet, "Gauss1");
  formulation.integral(dof(fieldForm3), tf(fieldForm3), omega, "Gauss1");

  formulation.pre();
  
  // scalar
  {
    gmshfem::function::ScalarFunction< T_Scalar > test = gmshfem::function::norm(fieldForm0) + gmshfem::function::norm(fieldForm3);
    
    const T_Scalar check = gmshfem::post::evaluate(test, 0.5853, 0.5497, 0.9172);
    const T_Scalar expected = 0.2756 + 47.265;
    gmshfem::msg::info << "\t- Check value = " << check << ", expected = " << expected << gmshfem::msg::endl;
    if(std::abs(check-expected) > 50. * gmshfem::scalar::Epsilon< T_Scalar >::value) {
      throw gmshfem::common::Exception("Error in " + std::string(__FILE__) + ", line " + std::to_string(__LINE__) + " -> function: normField< " + gmshfem::scalar::Name< T_Scalar >::name + " >");
    }
  }
  
  // scalar (abs)
  {
    gmshfem::function::ScalarFunction< T_Scalar > test = gmshfem::function::abs(fieldForm0) - gmshfem::function::abs(fieldForm3);
    
    const T_Scalar check = gmshfem::post::evaluate(test, 0.2858, 0.7572, 0.7537);
    const T_Scalar expected = 0.2756 - 47.265;
    gmshfem::msg::info << "\t- Check value = " << check << ", expected = " << expected << gmshfem::msg::endl;
    if(std::abs(check-expected) > 50. * gmshfem::scalar::Epsilon< T_Scalar >::value) {
      throw gmshfem::common::Exception("Error in " + std::string(__FILE__) + ", line " + std::to_string(__LINE__) + " -> function: normField< " + gmshfem::scalar::Name< T_Scalar >::name + " >");
    }
  }
  
  // vector, 1
  {
    gmshfem::function::ScalarFunction< T_Scalar > test1 = gmshfem::function::norm(fieldForm1, 1);
    gmshfem::function::ScalarFunction< T_Scalar > test2 = gmshfem::function::norm(fieldForm2, 1);
    
    const T_Scalar check1 = gmshfem::post::evaluate(test1, 0.3804, 0.5678, 0.);
    const T_Scalar check2 = gmshfem::post::evaluate(test2, 0.3804, 0.5678, 0.);
    const T_Scalar expected1 = 0.3824 + 3.5735 + 2.55;
    const T_Scalar expected2 = 0.3824 + 3.5735;
    gmshfem::msg::info << "\t- Check value = (" << check1 << ", " << check2 << "), expected = (" << expected1 << ", " << expected2 << ")" << gmshfem::msg::endl;
    if(std::abs(check1-expected1) > 50. * gmshfem::scalar::Epsilon< T_Scalar >::value || std::abs(check2-expected2) > 50. * gmshfem::scalar::Epsilon< T_Scalar >::value) {
      throw gmshfem::common::Exception("Error in " + std::string(__FILE__) + ", line " + std::to_string(__LINE__) + " -> function: normField< " + gmshfem::scalar::Name< T_Scalar >::name + " >");
    }
  }
}

template void normField< std::complex< double > >();
template void normField< std::complex< float > >();
template void normField< double >();
template void normField< float >();

template< class T_Scalar >
void complexFieldFunction()
{
  gmshfem::domain::Domain omega("omega");
  gmshfem::field::Field< T_Scalar, gmshfem::field::Form::Form0 > fieldForm0("form0", omega, gmshfem::field::FunctionSpaceTypeForm0::Lagrange);
  fieldForm0.addConstraint(omega, T_Scalar(9.483, 0.3985));
  gmshfem::field::Field< T_Scalar, gmshfem::field::Form::Form3 > fieldForm3("form3", omega, gmshfem::field::FunctionSpaceTypeForm3::Constant);
  fieldForm3.addConstraint(omega, T_Scalar(0.347, 25.257));
    
  gmshfem::problem::Formulation< T_Scalar > formulation("fake");
  
  formulation.integral(dof(fieldForm0), tf(fieldForm0), omega, "Gauss1");
  formulation.integral(dof(fieldForm3), tf(fieldForm3), omega, "Gauss1");

  formulation.pre();
  
  // scalar: conj(x)
  {
    gmshfem::function::ScalarFunction< T_Scalar > test = gmshfem::function::conj(fieldForm0) * gmshfem::function::conj(fieldForm3);
    
    const T_Scalar check = gmshfem::post::evaluate(test, 0.9340, 0.1299, 0.5688);
    const T_Scalar expected = T_Scalar(9.483 * 0.347 - 0.3985 * 25.257, - (9.483 * 25.257 + 0.3985 * 0.347));
    gmshfem::msg::info << "\t- Check value = " << check << ", expected = " << expected << gmshfem::msg::endl;
    if(std::abs(check-expected) > 500. * gmshfem::scalar::Epsilon< T_Scalar >::value) {
      throw gmshfem::common::Exception("Error in " + std::string(__FILE__) + ", line " + std::to_string(__LINE__) + " -> function: complexFieldFunction< " + gmshfem::scalar::Name< T_Scalar >::name + " >");
    }
  }
  
  // vector: real(x)
  {
    gmshfem::function::ScalarFunction< T_Scalar > test = gmshfem::function::real(fieldForm0) + gmshfem::function::real(fieldForm3);
    
    const T_Scalar check = gmshfem::post::evaluate(test, 0.9340, 0.1299, 0.5688);
    const T_Scalar expected = 9.483 + 0.347;
    gmshfem::msg::info << "\t- Check value = " << check << ", expected = " << expected << gmshfem::msg::endl;
    if(std::abs(check-expected) > 50. * gmshfem::scalar::Epsilon< T_Scalar >::value) {
      throw gmshfem::common::Exception("Error in " + std::string(__FILE__) + ", line " + std::to_string(__LINE__) + " -> function: complexFieldFunction< " + gmshfem::scalar::Name< T_Scalar >::name + " >");
    }
  }
  
  // tensor: imag(x)
  {
    gmshfem::function::ScalarFunction< T_Scalar > test = gmshfem::function::imag(fieldForm0) - gmshfem::function::imag(fieldForm3);
    
    const T_Scalar check = gmshfem::post::evaluate(test, 0.9340, 0.1299, 0.5688);
    const T_Scalar expected = 0.3985 - 25.257;
    gmshfem::msg::info << "\t- Check value = " << check << ", expected = " << expected << gmshfem::msg::endl;
    if(std::abs(check-expected) > 50. * gmshfem::scalar::Epsilon< T_Scalar >::value) {
      throw gmshfem::common::Exception("Error in " + std::string(__FILE__) + ", line " + std::to_string(__LINE__) + " -> function: complexFieldFunction< " + gmshfem::scalar::Name< T_Scalar >::name + " >");
    }
  }
}

template void complexFieldFunction< std::complex< double > >();
template void complexFieldFunction< std::complex< float > >();

template< class T_Scalar >
void derivativeFieldFunction()
{
  gmshfem::domain::Domain omega("omega");
  gmshfem::domain::Domain dirichlet("dirichlet");
  gmshfem::field::Field< T_Scalar, gmshfem::field::Form::Form0 > fieldForm0("form0", omega, gmshfem::field::FunctionSpaceTypeForm0::Lagrange);
  fieldForm0.addConstraint(omega, 0.14827);
  gmshfem::field::Field< T_Scalar, gmshfem::field::Form::Form1 > fieldForm1("form1", omega, gmshfem::field::FunctionSpaceTypeForm1::HierarchicalHCurl, 1);
  fieldForm1.addConstraint(omega, typename gmshfem::MathObject< T_Scalar, gmshfem::Degree::Degree1 >::Object(0.3824, 3.5735, 2.55));
  gmshfem::field::Field< T_Scalar, gmshfem::field::Form::Form2 > fieldForm2("form2", dirichlet, gmshfem::field::FunctionSpaceTypeForm2::P_HierarchicalHCurl, 1);
  fieldForm2.addConstraint(dirichlet, typename gmshfem::MathObject< T_Scalar, gmshfem::Degree::Degree1 >::Object(0.3824, 3.5735, 0.));
    
  gmshfem::problem::Formulation< T_Scalar > formulation("fake");

  formulation.integral(dof(fieldForm0), tf(fieldForm0), omega, "Gauss1");
  formulation.integral(dof(fieldForm1), tf(fieldForm1), omega, "Gauss1");
  formulation.integral(dof(fieldForm2), tf(fieldForm2), dirichlet, "Gauss1");

  formulation.pre();
    
  // d(form0) + d(form1)
  {
    gmshfem::function::ScalarFunction< T_Scalar > test = gmshfem::function::norm(gmshfem::function::d(fieldForm0) + gmshfem::function::d(fieldForm1)) + d(fieldForm2);
    
    const T_Scalar check = gmshfem::post::evaluate(test, 0.3842, 0.3852, 0.5684);
    const T_Scalar expected = 0.;
    gmshfem::msg::info << "\t- Check value = " << check << ", expected = " << expected << gmshfem::msg::endl;
    if(std::abs(check-expected) > 1.e4 * gmshfem::scalar::Epsilon< T_Scalar >::value) {
      throw gmshfem::common::Exception("Error in " + std::string(__FILE__) + ", line " + std::to_string(__LINE__) + " -> function: derivativeFieldFunction< " + gmshfem::scalar::Name< T_Scalar >::name + " >");
    }
  }
}

template void derivativeFieldFunction< std::complex< double > >();
template void derivativeFieldFunction< std::complex< float > >();
template void derivativeFieldFunction< double >();
template void derivativeFieldFunction< float >();

template< class T_Scalar >
void compoundField()
{
  gmshfem::domain::Domain omega("omega");
  gmshfem::domain::Domain domega = omega.getBoundary();

  gmshfem::field::CompoundField< T_Scalar, gmshfem::field::Form::Form0, 2 > v1("v1", omega, gmshfem::field::FunctionSpaceTypeForm0::Lagrange);
  gmshfem::field::CompoundField< T_Scalar, gmshfem::field::Form::Form0, 2 > v2("v2", omega, gmshfem::field::FunctionSpaceTypeForm0::Lagrange);
  gmshfem::field::CompoundField< T_Scalar, gmshfem::field::Form::Form0, 2 > v3("v3", omega, gmshfem::field::FunctionSpaceTypeForm0::Lagrange);
  gmshfem::field::CompoundField< T_Scalar, gmshfem::field::Form::Form0, 2 > v4("v4", omega, gmshfem::field::FunctionSpaceTypeForm0::Lagrange);
  
  gmshfem::function::VectorFunction< T_Scalar > nullv = gmshfem::function::vector< T_Scalar >(0.,0.,0.);
  v1.addConstraint(domega, nullv);
  v2.addConstraint(domega, nullv);
  v3.addConstraint(domega, nullv);
  v4.addConstraint(domega, nullv);

  gmshfem::function::ScalarFunction< T_Scalar > X = gmshfem::function::x< T_Scalar >() * (gmshfem::function::x< T_Scalar >() - 1.);
  gmshfem::function::ScalarFunction< T_Scalar > Y = gmshfem::function::y< T_Scalar >() * (gmshfem::function::y< T_Scalar >() - 1.);
  gmshfem::function::ScalarFunction< T_Scalar > F = X*Y;
  gmshfem::function::TensorFunction< T_Scalar > L = gmshfem::function::tensor< T_Scalar >
  (
    F,2.*F,0.,
    3.*F, 4.*F, 0.,
    0., 0., 0.
  );
  gmshfem::function::ScalarFunction< T_Scalar > dX = 2. * gmshfem::function::x< T_Scalar >() - 1.;
  gmshfem::function::ScalarFunction< T_Scalar > dY = 2. * gmshfem::function::y< T_Scalar >() - 1.;
  gmshfem::function::VectorFunction< T_Scalar > divL = gmshfem::function::vector< T_Scalar > ( dX*Y+3.*X*dY, 2.*dX*Y+4.*X*dY, 0. );

  // Version 1: reference
  gmshfem::problem::Formulation< T_Scalar > formulation1("formulation1");

  formulation1.integral(grad(dof(v1)), grad(tf(v1)), omega, "Gauss6");
  formulation1.integral(divL, tf(v1), omega, "Gauss6");

  formulation1.pre();
  formulation1.assemble();
  formulation1.solve();

  // Version 2: test linear term with a grad(tf())
  gmshfem::problem::Formulation< T_Scalar > formulation2("formulation2");

  formulation2.integral(grad(dof(v2)), grad(tf(v2)), omega, "Gauss6");
  formulation2.integral(-L, grad(tf(v2)), omega, "Gauss6");

  formulation2.pre();
  formulation2.assemble();
  formulation2.solve();

  // Version 3: test a double dot product (C * L)
  TensorFunction< T_Scalar > null = tensor< T_Scalar >
  (
    0., 0., 0.,
    0., 0., 0.,
    0., 0., 0.
  );
  TensorFunction< T_Scalar > C_00 = tensor< T_Scalar >
  (
    1.,0.,0.,
    0., 1., 0.,
    0., 0., 0.
  );

  TensorFunction< T_Scalar, 4 > C = tensor< T_Scalar, 4 >
  (
    C_00, null, null,
    null, C_00, null,
    null, null, null
  );
  gmshfem::problem::Formulation< T_Scalar > formulation3("formulation3");

  formulation3.integral( grad(dof(v3)), grad(tf(v3)), omega, "Gauss6");
  formulation3.integral( C * L, grad(tf(v3)), omega, "Gauss6");

  formulation3.pre();
  formulation3.assemble();
  formulation3.solve();

  // Version 4: Comparaison with the analytic double dot product CL
  gmshfem::function::VectorFunction< T_Scalar > ex = gmshfem::function::vector< T_Scalar >(1., 0., 0.);
  gmshfem::function::VectorFunction< T_Scalar > ey = gmshfem::function::vector< T_Scalar >(0., 1., 0.);
  gmshfem::function::TensorFunction< T_Scalar > CL = gmshfem::function::tensor< T_Scalar >
  (
    ex*L*ex+ey*L*ey, 0., 0.,
    0., ex*L*ex+ey*L*ey, 0.,
    0., 0., 0.
  );
    
  gmshfem::problem::Formulation< T_Scalar > formulation4("formulation4");

  formulation4.integral( grad(dof(v4)), grad(tf(v4)), omega, "Gauss6");
  formulation4.integral(CL, grad(tf(v4)), omega, "Gauss6");

  formulation4.pre();
  formulation4.assemble();
  formulation4.solve();
  
  {
    gmshfem::function::VectorFunction< T_Scalar > test = v1 - v2;
    
    const T_Scalar check = gmshfem::post::integrate(norm(test), omega, "Gauss6");
    const T_Scalar expected = 0.;
    gmshfem::msg::info << "\t- Check value = " << check << ", expected = " << expected << gmshfem::msg::endl;
    if(std::abs(check-expected) > gmshfem::scalar::Epsilon< T_Scalar >::value) {
      throw gmshfem::common::Exception("Error in " + std::string(__FILE__) + ", line " + std::to_string(__LINE__) + " -> function: compoundField< " + gmshfem::scalar::Name< T_Scalar >::name + " >");
    }
  }
  
  {
    gmshfem::function::VectorFunction< T_Scalar > test = v3 - v4;
    
    const T_Scalar check = gmshfem::post::integrate(norm(test), omega, "Gauss6");
    const T_Scalar expected = 0.;
    gmshfem::msg::info << "\t- Check value = " << check << ", expected = " << expected << gmshfem::msg::endl;
    if(std::abs(check-expected) > gmshfem::scalar::Epsilon< T_Scalar >::value) {
      throw gmshfem::common::Exception("Error in " + std::string(__FILE__) + ", line " + std::to_string(__LINE__) + " -> function: compoundField< " + gmshfem::scalar::Name< T_Scalar >::name + " >");
    }
  }
}

template void compoundField< std::complex< double > >();
template void compoundField< std::complex< float > >();
template void compoundField< double >();
template void compoundField< float >();
