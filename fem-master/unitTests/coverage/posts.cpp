// GmshFEM - Copyright (C) 2019-2021 Université de Liège
//
// See the LICENSE.txt file for license information. Please report all
// issues on https://gitlab.onelab.info/gmsh/fem/issues

#include <Formulation.h>
#include <Function.h>
#include <Post.h>

#include "posts.h"
#include "geo.h"


using gmshfem::equation::dof;
using gmshfem::equation::dt_dof;
using gmshfem::equation::dt2_dof;
using gmshfem::equation::tf;
using namespace gmshfem::function;

void posts(gmshfem::common::GmshFem &gmshfem, unsigned int &numTest)
{
  gmshfem::msg::info << ++numTest << ") Posts" << gmshfem::msg::endl;

  Geo2D::triangle();

  try {
    gmshfem::msg::info << "Point evaluation and field copy:" << gmshfem::msg::endl;
    pointEvaluationAndFieldCopy< std::complex< double > >();
    pointEvaluationAndFieldCopy< std::complex< float > >();
    pointEvaluationAndFieldCopy< double >();
    pointEvaluationAndFieldCopy< float >();
  } catch(...) { throw; }
  
  try {
    gmshfem::msg::info << "Algebra:" << gmshfem::msg::endl;
    algebra< std::complex< double > >();
    algebra< std::complex< float > >();
    algebra< double >();
    algebra< float >();
  } catch(...) { throw; }
  
  removeGeo();
  
  Geo3D::tetrahedra();
  
  try {
    gmshfem::msg::info << "Integrate:" << gmshfem::msg::endl;
    integrate< std::complex< double > >();
    integrate< std::complex< float > >();
    integrate< double >();
    integrate< float >();
  } catch(...) { throw; }
  
  try {
    gmshfem::msg::info << "Save:" << gmshfem::msg::endl;
    save< std::complex< double > >();
    save< std::complex< float > >();
    save< double >();
    save< float >();
  } catch(...) { throw; }
  
  
  removeGeo();
  
  try {
    gmshfem::msg::info << "Basis functions draw:" << gmshfem::msg::endl;
    basisFunctionsDraw< std::complex< double > >();
    basisFunctionsDraw< std::complex< float > >();
    basisFunctionsDraw< double >();
    basisFunctionsDraw< float >();
  } catch(...) { throw; }
  
  return;
}

template< class T_Scalar >
void integrate()
{
  {
    gmshfem::domain::Domain omega("omega");
        
    gmshfem::problem::Formulation< T_Scalar > formulation("integrate");

    gmshfem::field::Field< T_Scalar, gmshfem::field::Form::Form0 > u("u", omega, gmshfem::field::FunctionSpaceTypeForm0::Lagrange);
    
    formulation.integral(dof(u), tf(u), omega, "Gauss2");
    formulation.integral(-1., tf(u), omega, "Gauss2");
    
    formulation.pre();
    formulation.assemble();
    formulation.solve();
        
    const T_Scalar check1 = gmshfem::post::integrate(u, omega, "Gauss2");
    gmshfem::post::Circle circle("circle", 0.5, 0.5, 0., 0.2, 0., 0., 1., 100);
    const T_Scalar check2 = gmshfem::post::integrate(2. * u, circle, "Gauss2");
    gmshfem::post::Line line("line", 0.2, 0.2, 0., 0.8, 0.8, 0., 50);
    const T_Scalar check3 = gmshfem::post::integrate(u+1., line, "Gauss2");
    gmshfem::post::Plane plane("plane", 0.2, 0.2, 0., 0.6, 0., 0., 0., 0.6, 0., 40, 40);
    const T_Scalar check4 = gmshfem::post::integrate(u/2., plane, "Gauss2");
    
    const double pi = 3.14159265358979323846264338327950288419716939937510;
    
    const T_Scalar expected1 = 1.;
    const T_Scalar expected2 = 2. * 2. * 100. * 0.2 * std::sin(pi/100.);
    const T_Scalar expected3 = 2. * (0.8-0.2) * std::sqrt(2.);
    const T_Scalar expected4 = 0.5 * (0.8-0.2) * (0.8-0.2);
    gmshfem::msg::info << "\t- Check value = (" << check1 << ", " << check2 << ", " << check3 << ", " << check4 << "), expected = (" << expected1 << ", " << expected2 << ", " << expected3 << ", " << expected4 << ")" << gmshfem::msg::endl;
    if(std::abs(check1-expected1) > 10. * gmshfem::scalar::Epsilon< T_Scalar >::value || std::abs(check2-expected2) > 10. * gmshfem::scalar::Epsilon< T_Scalar >::value || std::abs(check3-expected3) > 10. * gmshfem::scalar::Epsilon< T_Scalar >::value) {
      throw gmshfem::common::Exception("Error in " + std::string(__FILE__) + ", line " + std::to_string(__LINE__) + " -> function: integrate< " + gmshfem::scalar::Name< T_Scalar >::name + " >");
    }
  }
}

template void integrate< std::complex< double > >();
template void integrate< std::complex< float > >();
template void integrate< double >();
template void integrate< float >();

template< class T_Scalar >
void save()
{
  {
    gmshfem::domain::Domain omega("omega");
    gmshfem::function::TensorFunction< T_Scalar > fct = gmshfem::function::tensor< T_Scalar >(1., 2., 3.,  4., 5., 6.,  7., 8., 9.);
        
    gmshfem::post::save(fct, omega, "fct");
    gmshfem::post::Circle circle("circle", 0.5, 0.5, 0., 0.2, 0., 0., 1., 100);
    gmshfem::post::save(2. * fct, circle, "fct");
    gmshfem::post::Line line("line", 0.2, 0.2, 0., 0.8, 0.8, 0., 50);
    gmshfem::post::save(fct + tensorDiag< T_Scalar >(1., 2., 3.), line, "fct");
    gmshfem::post::Plane plane("plane", 0.2, 0.2, 0., 0.6, 0., 0., 0., 0.6, 0., 40, 40);
    gmshfem::post::save(fct/2., plane, "fct");
    gmshfem::post::Disk disk("disk", 0.5, 0.5, 0., 0.2, 0., 0., 1., 100);
    gmshfem::post::save(fct * 3., disk, "fct");
    gmshfem::post::Sphere sphere("sphere", 0.5, 0.5, 0.5, 0.2, 100);
    gmshfem::post::save(fct / 3., sphere, "fct");

    // Unable to check ...
  }
}

template void save< std::complex< double > >();
template void save< std::complex< float > >();
template void save< double >();
template void save< float >();

template< class T_Scalar >
void basisFunctionsDraw()
{
  bool derivative[2] = {false, true};
  std::string element[8] = {"point", "line", "triangle", "quadrangle", "tetrahedron", "hexahedron", "prism", "pyramid"};
  std::string form[4] = {"form0", "form1", "form2", "form3"};
  std::string basisFunction[4] = {"isoParametric", "hierarchical", "p_isoParametric", "p_hierarchical"};
  
  gmshfem::domain::Domain empty;
  
  for(unsigned int i = 0; i < 2; ++i) {
    for(unsigned int j = 0; j < 8; ++j) {
      for(unsigned int k = 0; k < 4; ++k) {
        for(unsigned int l = 0; l < 4; ++l) {
          if(form[k] == "form0") {
            if(basisFunction[l] == "isoParametric") {
              gmshfem::field::Field< T_Scalar, gmshfem::field::Form::Form0 > field("field", empty, gmshfem::functionSpaceH1::Lagrange, 1);
              field.getFunctionSpace()->draw(derivative[i], element[j], 0, 0.5);
            }
            else if(basisFunction[l] == "hierarchical" && j != 7) {
              gmshfem::field::Field< T_Scalar, gmshfem::field::Form::Form0 > field("field", empty, gmshfem::functionSpaceH1::HierarchicalH1, 1);
              field.getFunctionSpace()->draw(derivative[i], element[j], 0, 0.5);
            }
            else {
              continue;
            }
          }
          else if(form[k] == "form1") {
            if(basisFunction[l] == "hierarchical" && j >= 1 && j != 7) {
              gmshfem::field::Field< T_Scalar, gmshfem::field::Form::Form1 > field("field", empty, gmshfem::functionSpaceHCurl::HierarchicalHCurl, 1);
              field.getFunctionSpace()->draw(derivative[i], element[j], 0, 0.5);
            }
            else if(basisFunction[l] == "p_isoParametric" && j < 4) {
              gmshfem::field::Field< T_Scalar, gmshfem::field::Form::Form1 > field("field", empty, gmshfem::functionSpaceHCurl::P_Lagrange, 1);
              field.getFunctionSpace()->draw(derivative[i], element[j], 0, 0.5);
            }
            else if(basisFunction[l] == "p_hierarchical" && j < 4 && j != 7) {
              gmshfem::field::Field< T_Scalar, gmshfem::field::Form::Form1 > field("field", empty, gmshfem::functionSpaceHCurl::P_HierarchicalH1, 1);
              field.getFunctionSpace()->draw(derivative[i], element[j], 0, 0.5);
            }
            else {
              continue;
            }
          }
          else if(form[k] == "form2") {
            if(basisFunction[l] == "p_hierarchical" && j >= 2 && j < 4 && j != 7) {
              gmshfem::field::Field< T_Scalar, gmshfem::field::Form::Form2 > field("field", empty, gmshfem::functionSpaceHDiv::P_HierarchicalHCurl, 1);
              field.getFunctionSpace()->draw(derivative[i], element[j], 0, 0.5);
            }
            else {
              continue;
            }
          }
          else if(form[k] == "form3") {
            if(i != 1) {
              if(basisFunction[l] == "p_hierarchical" && j >= 2 && j < 4) {
                gmshfem::field::Field< T_Scalar, gmshfem::field::Form::Form3 > field("field", empty, gmshfem::functionSpaceH0::P_Constant, 1);
                field.getFunctionSpace()->draw(derivative[i], element[j], 0, 0.5);
              }
              else if(basisFunction[l] == "hierarchical" && j >= 4) {
                gmshfem::field::Field< T_Scalar, gmshfem::field::Form::Form3 > field("field", empty, gmshfem::functionSpaceH0::Constant, 1);
                field.getFunctionSpace()->draw(derivative[i], element[j], 0, 0.5);
              }
              else {
                continue;
              }
            }
          }
          else {
            throw gmshfem::common::Exception("Error in " + std::string(__FILE__) + ", line " + std::to_string(__LINE__) + " -> function: basisFunctionsDraw< " + gmshfem::scalar::Name< T_Scalar >::name + " >");
          }
        }
      }
    }
  }
}

template void basisFunctionsDraw< std::complex< double > >();
template void basisFunctionsDraw< std::complex< float > >();
template void basisFunctionsDraw< double >();
template void basisFunctionsDraw< float >();

template< class T_Scalar >
void pointEvaluationAndFieldCopy()
{
  {
    gmshfem::domain::Domain omega("omega");
    gmshfem::problem::Formulation< T_Scalar > formulation("pointEvaluationAndFieldCopy");
    
    gmshfem::field::Field< T_Scalar, gmshfem::field::Form::Form0 > u("u", omega, gmshfem::field::FunctionSpaceTypeForm0::HierarchicalH1, 1);
    formulation.integral(dof(u), tf(u), omega, "Gauss4");
    formulation.integral(- gmshfem::function::x< T_Scalar >(), tf(u), omega, "Gauss4");
    
    formulation.pre();
    formulation.assemble();
    formulation.solve();
    
    gmshfem::field::Field< T_Scalar, gmshfem::field::Form::Form0 > v("v", omega, gmshfem::field::FunctionSpaceTypeForm0::HierarchicalH1, 1);
    v = u;
    const T_Scalar check = gmshfem::post::evaluate(v, 0.8147, 0.9058, 0.);
    const T_Scalar expected = 0.8147;
    gmshfem::msg::info << "\t- Check value = " << check << ", expected = " << expected << gmshfem::msg::endl;
    if(std::abs(check-expected) > 10. * gmshfem::scalar::Epsilon< T_Scalar >::value) {
      throw gmshfem::common::Exception("Error in " + std::string(__FILE__) + ", line " + std::to_string(__LINE__) + " -> function: pointEvaluation< " + gmshfem::scalar::Name< T_Scalar >::name + " >");
    }
  }
  
  {
    gmshfem::domain::Domain omega("omega");
    gmshfem::problem::Formulation< T_Scalar > formulation("pointEvaluationAndFieldCopy");
    
    gmshfem::field::Field< T_Scalar, gmshfem::field::Form::Form0 > u("u", omega, gmshfem::field::FunctionSpaceTypeForm0::HierarchicalH1, 1);
    formulation.integral(dof(u), tf(u), omega, "Gauss4");
    formulation.integral(- gmshfem::function::x< T_Scalar >(), tf(u), omega, "Gauss4");
    
    formulation.pre();
    formulation.assemble();
    formulation.solve();
    
    gmshfem::field::Field< T_Scalar, gmshfem::field::Form::Form0 > v("v", omega, gmshfem::field::FunctionSpaceTypeForm0::HierarchicalH1, 1);
    v = u;
    const typename gmshfem::MathObject< T_Scalar, gmshfem::Degree::Degree1 >::Object check = gmshfem::post::evaluate(grad(v), 0.8147, 0.9058, 0.);
    const typename gmshfem::MathObject< T_Scalar, gmshfem::Degree::Degree1 >::Object expected = typename gmshfem::MathObject< T_Scalar, gmshfem::Degree::Degree1 >::Object(1., 0., 0.);
    gmshfem::msg::info << "\t- Check value = " << check.norm() << ", expected = " << expected.norm() << gmshfem::msg::endl;
    if((check-expected).norm() > 100. * gmshfem::scalar::Epsilon< T_Scalar >::value) {
      throw gmshfem::common::Exception("Error in " + std::string(__FILE__) + ", line " + std::to_string(__LINE__) + " -> function: pointEvaluation< " + gmshfem::scalar::Name< T_Scalar >::name + " >");
    }
  }
  
  {
    gmshfem::domain::Domain omega("omega");
    gmshfem::problem::Formulation< T_Scalar > formulation("pointEvaluationAndFieldCopy");
    
    gmshfem::field::Field< T_Scalar, gmshfem::field::Form::Form1 > u("u", omega, gmshfem::field::FunctionSpaceTypeForm1::HierarchicalHCurl, 1);
    formulation.integral(dof(u), tf(u), omega, "Gauss4");
    formulation.integral(- gmshfem::function::vector< T_Scalar >(gmshfem::function::x< T_Scalar >(), gmshfem::function::y< T_Scalar >(), 0.), tf(u), omega, "Gauss4");
    
    formulation.pre();
    formulation.assemble();
    formulation.solve();
    
    gmshfem::field::Field< T_Scalar, gmshfem::field::Form::Form1 > v("v", omega, gmshfem::field::FunctionSpaceTypeForm1::HierarchicalHCurl, 1);
    v = u;
    const typename gmshfem::MathObject< T_Scalar, gmshfem::Degree::Degree1 >::Object check = gmshfem::post::evaluate(v, 0.9134, 0.6324, 0.);
    const typename gmshfem::MathObject< T_Scalar, gmshfem::Degree::Degree1 >::Object expected = typename gmshfem::MathObject< T_Scalar, gmshfem::Degree::Degree1 >::Object(0.9134, 0.6324, 0.);
    gmshfem::msg::info << "\t- Check value = " << check.norm() << ", expected = " << expected.norm() << gmshfem::msg::endl;
    if((check-expected).norm() > 10. * gmshfem::scalar::Epsilon< T_Scalar >::value) {
      throw gmshfem::common::Exception("Error in " + std::string(__FILE__) + ", line " + std::to_string(__LINE__) + " -> function: pointEvaluation< " + gmshfem::scalar::Name< T_Scalar >::name + " >");
    }
  }
  
  {
    gmshfem::domain::Domain omega("omega");
    gmshfem::problem::Formulation< T_Scalar > formulation("pointEvaluationAndFieldCopy");
    
    gmshfem::field::Field< T_Scalar, gmshfem::field::Form::Form1 > u("u", omega, gmshfem::field::FunctionSpaceTypeForm1::HierarchicalHCurl, 2);
    formulation.integral(dof(u), tf(u), omega, "Gauss4");
    formulation.integral(- gmshfem::function::vector< T_Scalar >(gmshfem::function::x< T_Scalar >(), gmshfem::function::y< T_Scalar >(), 0.), tf(u), omega, "Gauss4");
    
    formulation.pre();
    formulation.assemble();
    formulation.solve();
    
    gmshfem::field::Field< T_Scalar, gmshfem::field::Form::Form1 > v("v", omega, gmshfem::field::FunctionSpaceTypeForm1::HierarchicalHCurl, 2);
    v = u;
    const typename gmshfem::MathObject< T_Scalar, gmshfem::Degree::Degree1 >::Object check = gmshfem::post::evaluate(curl(v), 0.9134, 0.6324, 0.);
    const typename gmshfem::MathObject< T_Scalar, gmshfem::Degree::Degree1 >::Object expected = typename gmshfem::MathObject< T_Scalar, gmshfem::Degree::Degree1 >::Object(0., 0., 0.);
    gmshfem::msg::info << "\t- Check value = " << check.norm() << ", expected = " << expected.norm() << gmshfem::msg::endl;
    if((check-expected).norm() > 1000. * gmshfem::scalar::Epsilon< T_Scalar >::value) {
      throw gmshfem::common::Exception("Error in " + std::string(__FILE__) + ", line " + std::to_string(__LINE__) + " -> function: pointEvaluation< " + gmshfem::scalar::Name< T_Scalar >::name + " >");
    }
  }
  
  {
    gmshfem::domain::Domain omega("omega");
    gmshfem::problem::Formulation< T_Scalar > formulation("pointEvaluationAndFieldCopy");
    
    gmshfem::field::Field< T_Scalar, gmshfem::field::Form::Form2 > u("u", omega, gmshfem::field::FunctionSpaceTypeForm2::P_HierarchicalHCurl, 1);
    formulation.integral(dof(u), tf(u), omega, "Gauss4");
    formulation.integral(- gmshfem::function::vector< T_Scalar >(gmshfem::function::x< T_Scalar >(), gmshfem::function::y< T_Scalar >(), 0.), tf(u), omega, "Gauss4");
    
    formulation.pre();
    formulation.assemble();
    formulation.solve();
    
    gmshfem::field::Field< T_Scalar, gmshfem::field::Form::Form2 > v("v", omega, gmshfem::field::FunctionSpaceTypeForm2::P_HierarchicalHCurl, 1);
    v = u;
    const typename gmshfem::MathObject< T_Scalar, gmshfem::Degree::Degree1 >::Object check = gmshfem::post::evaluate(v, 0.2785, 0.5469, 0.);
    const typename gmshfem::MathObject< T_Scalar, gmshfem::Degree::Degree1 >::Object expected = typename gmshfem::MathObject< T_Scalar, gmshfem::Degree::Degree1 >::Object(0.2785, 0.5469, 0.);
    gmshfem::msg::info << "\t- Check value = " << check.norm() << ", expected = " << expected.norm() << gmshfem::msg::endl;
    if((check-expected).norm() > 10. * gmshfem::scalar::Epsilon< T_Scalar >::value) {
      throw gmshfem::common::Exception("Error in " + std::string(__FILE__) + ", line " + std::to_string(__LINE__) + " -> function: pointEvaluation< " + gmshfem::scalar::Name< T_Scalar >::name + " >");
    }
  }
  
  {
    gmshfem::domain::Domain omega("omega");
    gmshfem::problem::Formulation< T_Scalar > formulation("pointEvaluationAndFieldCopy");
    
    gmshfem::field::Field< T_Scalar, gmshfem::field::Form::Form2 > u("u", omega, gmshfem::field::FunctionSpaceTypeForm2::P_HierarchicalHCurl, 1);
    formulation.integral(dof(u), tf(u), omega, "Gauss4");
    formulation.integral(- gmshfem::function::vector< T_Scalar >(gmshfem::function::x< T_Scalar >(), gmshfem::function::y< T_Scalar >(), 0.), tf(u), omega, "Gauss4");
    
    formulation.pre();
    formulation.assemble();
    formulation.solve();
    
    gmshfem::field::Field< T_Scalar, gmshfem::field::Form::Form2 > v("v", omega, gmshfem::field::FunctionSpaceTypeForm2::P_HierarchicalHCurl, 1);
    v = u;
    const T_Scalar check = gmshfem::post::evaluate(div(v), 0.2785, 0.5469, 0.);
    const T_Scalar expected = 2.;
    gmshfem::msg::info << "\t- Check value = " << check << ", expected = " << expected << gmshfem::msg::endl;
    if(std::abs(check-expected) > 500. * gmshfem::scalar::Epsilon< T_Scalar >::value) {
      throw gmshfem::common::Exception("Error in " + std::string(__FILE__) + ", line " + std::to_string(__LINE__) + " -> function: pointEvaluation< " + gmshfem::scalar::Name< T_Scalar >::name + " >");
    }
  }
  
  {
    gmshfem::domain::Domain omega("omega");
    gmshfem::problem::Formulation< T_Scalar > formulation("pointEvaluationAndFieldCopy");
    
    gmshfem::field::Field< T_Scalar, gmshfem::field::Form::Form3 > u("u", omega, gmshfem::field::FunctionSpaceTypeForm3::P_Constant, 1);
    formulation.integral(dof(u), tf(u), omega, "Gauss4");
    formulation.integral(- 1., tf(u), omega, "Gauss4");
    
    formulation.pre();
    formulation.assemble();
    formulation.solve();
    
    gmshfem::field::Field< T_Scalar, gmshfem::field::Form::Form3 > v("v", omega, gmshfem::field::FunctionSpaceTypeForm3::P_Constant, 1);
    v = u;
    const T_Scalar check = gmshfem::post::evaluate(v, 0.5634, 0.1735, 0.);
    const T_Scalar expected = 1.;
    gmshfem::msg::info << "\t- Check value = " << check << ", expected = " << expected << gmshfem::msg::endl;
    if(std::abs(check-expected) > 10. * gmshfem::scalar::Epsilon< T_Scalar >::value) {
      throw gmshfem::common::Exception("Error in " + std::string(__FILE__) + ", line " + std::to_string(__LINE__) + " -> function: pointEvaluation< " + gmshfem::scalar::Name< T_Scalar >::name + " >");
    }
  }
}

template void pointEvaluationAndFieldCopy< std::complex< double > >();
template void pointEvaluationAndFieldCopy< std::complex< float > >();
template void pointEvaluationAndFieldCopy< double >();
template void pointEvaluationAndFieldCopy< float >();

template< class T_Scalar >
void algebra()
{
  {
    gmshfem::domain::Domain omega("omega");
    gmshfem::domain::Domain dirichlet("dirichlet");
    gmshfem::problem::Formulation< T_Scalar > formulation("algebra");
    
    gmshfem::field::Field< T_Scalar, gmshfem::field::Form::Form0 > u("u", omega, gmshfem::field::FunctionSpaceTypeForm0::Lagrange);
    u.addConstraint(dirichlet, 0.);
    formulation.integral(-grad(dof(u)), grad(tf(u)), omega, "Gauss4");
    formulation.integral(1., tf(u), omega, "Gauss4");
    
    formulation.pre();
    formulation.assemble();
    formulation.solve();
    
    gmshfem::algebra::MatrixCRSFast< T_Scalar > matCRSFast;
    formulation.getLHS(matCRSFast);
    gmshfem::algebra::MatrixCRS< T_Scalar > matCRS;
    matCRS = matCRSFast;
    gmshfem::algebra::MatrixCCS< T_Scalar > matCCS;
    matCCS = matCRSFast;
    
    gmshfem::algebra::MatrixCRSFast< T_Scalar > matCRSFast2(matCRSFast);
    
    gmshfem::algebra::MatrixCRS< T_Scalar > matCRS2(matCRS);
    formulation.getLHS(matCRS2);
    matCCS = matCRS2;

    gmshfem::algebra::MatrixCCS< T_Scalar > matCCS2(matCCS);
    formulation.getLHS(matCCS2);
    matCRS = matCCS2;
    
    gmshfem::algebra::MatrixCRSFast< T_Scalar > matCRSFast3(std::move(matCRSFast2));
    gmshfem::algebra::MatrixCRS< T_Scalar > matCRS3(std::move(matCRS2));
    gmshfem::algebra::MatrixCCS< T_Scalar > matCCS3(std::move(matCCS2));
    
    matCRSFast2 = matCRSFast3;
    matCRS2 = matCRS3;
    matCCS2 = matCCS3;
    
    matCRSFast3 = std::move(matCRSFast2);
    matCRS3 = std::move(matCRS2);
    matCCS3 = std::move(matCCS2);
        
    {
      const bool check1 = matCRSFast3.isSymmetric();
      const bool expected1 = true;
      const bool check2 = matCRSFast3.isHermitian();
      const bool expected2 = true;
      gmshfem::msg::info << "\t- Check value = (" << check1 << ", " << check2 << "), expected = (" << expected1 << ", " << expected2 << ")" << gmshfem::msg::endl;
      if(check1 != expected1 || check2 != expected2) {
        throw gmshfem::common::Exception("Error in " + std::string(__FILE__) + ", line " + std::to_string(__LINE__) + " -> function: algebra< " + gmshfem::scalar::Name< T_Scalar >::name + " >");
      }
    }
    {
      const bool check1 = matCRS3.isSymmetric();
      const bool expected1 = true;
      const bool check2 = matCRS3.isHermitian();
      const bool expected2 = true;
      gmshfem::msg::info << "\t- Check value = (" << check1 << ", " << check2 << "), expected = (" << expected1 << ", " << expected2 << ")" << gmshfem::msg::endl;
      if(check1 != expected1 || check2 != expected2) {
        throw gmshfem::common::Exception("Error in " + std::string(__FILE__) + ", line " + std::to_string(__LINE__) + " -> function: algebra< " + gmshfem::scalar::Name< T_Scalar >::name + " >");
      }
    }
    {
      const bool check1 = matCCS3.isSymmetric();
      const bool expected1 = true;
      const bool check2 = matCCS3.isHermitian();
      const bool expected2 = true;
      gmshfem::msg::info << "\t- Check value = (" << check1 << ", " << check2 << "), expected = (" << expected1 << ", " << expected2 << ")" << gmshfem::msg::endl;
      if(check1 != expected1 || check2 != expected2) {
        throw gmshfem::common::Exception("Error in " + std::string(__FILE__) + ", line " + std::to_string(__LINE__) + " -> function: algebra< " + gmshfem::scalar::Name< T_Scalar >::name + " >");
      }
    }
    
    {
      const unsigned int check1 = matCRS3.numberOfNonZero();
      const unsigned int check2 = matCCS3.numberOfNonZero();
      const unsigned int expected = matCRSFast3.numberOfNonZero();
      gmshfem::msg::info << "\t- Check value = (" << check1 << ", " << check2 << "), expected = (" << expected << ", " << expected << ")" << gmshfem::msg::endl;
      if(check1 != expected || check2 != expected) {
        throw gmshfem::common::Exception("Error in " + std::string(__FILE__) + ", line " + std::to_string(__LINE__) + " -> function: algebra< " + gmshfem::scalar::Name< T_Scalar >::name + " >");
      }
    }
    
    {
      const double check1 = matCRS3.sparsity();
      const double check2 = matCCS3.sparsity();
      const double expected = matCRSFast3.sparsity();
      gmshfem::msg::info << "\t- Check value = (" << check1 << ", " << check2 << "), expected = (" << expected << ", " << expected << ")" << gmshfem::msg::endl;
      if(check1 != expected || check2 != expected) {
        throw gmshfem::common::Exception("Error in " + std::string(__FILE__) + ", line " + std::to_string(__LINE__) + " -> function: algebra< " + gmshfem::scalar::Name< T_Scalar >::name + " >");
      }
    }
    
    matCRSFast3.save("matCRSFast");
    matCRSFast3.saveSpyPlot("plot_matCRSFast");
    
    matCRS3.save("matCRS");
    matCRS3.saveSpyPlot("plot_matCRS");
    
    matCCS3.save("matCCS");
    matCCS3.saveSpyPlot("plot_matCCS");
    
    gmshfem::algebra::Vector< T_Scalar > b;
    formulation.getRHS(b);
    
    gmshfem::algebra::Vector< T_Scalar > x;
    u.getUnknownVector(x);
    
    {
      const gmshfem::scalar::Precision< T_Scalar > check = gmshfem::algebra::residual(b, matCRS3, x);
      const gmshfem::scalar::Precision< T_Scalar > expected = formulation.getResidual(); // return b - Ax;
      gmshfem::msg::info << "\t- Check value = " << check << ", expected = " << expected << gmshfem::msg::endl;
      if(std::abs(check-expected) > 2. * gmshfem::scalar::Epsilon< T_Scalar >::value) {
        throw gmshfem::common::Exception("Error in " + std::string(__FILE__) + ", line " + std::to_string(__LINE__) + " -> function: algebra< " + gmshfem::scalar::Name< T_Scalar >::name + " >");
      }
    }
    
    {
      gmshfem::algebra::Vector< T_Scalar > b2;
      gmshfem::algebra::linear(b2, matCRSFast3, x);
      T_Scalar check = 0.;
      for(unsigned int i = 0; i < b2.size(); ++i) {
        check += std::abs((b2[i] + b[i]) * (b2[i] + b[i]));
      }
      check = std::sqrt(check);
      const T_Scalar expected = 0.;
      gmshfem::msg::info << "\t- Check value = " << check << ", expected = " << expected << gmshfem::msg::endl;
      if(std::abs(check-expected) > 50. * gmshfem::scalar::Epsilon< T_Scalar >::value) {
        throw gmshfem::common::Exception("Error in " + std::string(__FILE__) + ", line " + std::to_string(__LINE__) + " -> function: algebra< " + gmshfem::scalar::Name< T_Scalar >::name + " >");
      }
    }
    
    {
      gmshfem::algebra::Vector< T_Scalar > b3;
      gmshfem::algebra::linear(b3, matCRS3, x);
      T_Scalar check = 0.;
      for(unsigned int i = 0; i < b3.size(); ++i) {
        check += std::abs((b3[i] + b[i]) * (b3[i] + b[i]));
      }
      check = std::sqrt(check);
      const T_Scalar expected = 0.;
      gmshfem::msg::info << "\t- Check value = " << check << ", expected = " << expected << gmshfem::msg::endl;
      if(std::abs(check-expected) > 50. * gmshfem::scalar::Epsilon< T_Scalar >::value) {
        throw gmshfem::common::Exception("Error in " + std::string(__FILE__) + ", line " + std::to_string(__LINE__) + " -> function: algebra< " + gmshfem::scalar::Name< T_Scalar >::name + " >");
      }
    }
  }
}

template void algebra< std::complex< double > >();
template void algebra< std::complex< float > >();
template void algebra< double >();
template void algebra< float >();
