// GmshFEM - Copyright (C) 2019-2022, A. Royer, E. Béchet, C. Geuzaine, Université de Liège
//
// See the LICENSE.txt file for license information. Please report all
// issues on https://gitlab.onelab.info/gmsh/fem/issues

#include "CompoundField.h"

#include "Exception.h"
#include "instantiate.h"

namespace gmshfem::field
{


  template< class T_Scalar, field::Form T_Form, unsigned int T_NumFields >
  bool CompoundField< T_Scalar, T_Form, T_NumFields >::_checkDomain() const
  {
    if constexpr(T_Form == field::Form::Form0) {
      return true;
    }
    else if constexpr(T_Form == field::Form::Form1) {
      bool error = false;

      // Check dimension of domain
      if(this->_domain.minDim() < 1) {
        if(this->_functionSpace) {
          if(this->_functionSpace->isFormP()) {
            // Ok
          }
          else {
            msg::warning << "CompoundField '" << this->_name << "': " + std::string(NameOfForm< T_Form >::value) << " is not defined on points" << msg::endl;
            error = true;
          }
        }
      }

      // Check that FormP is not defined on 3D mesh
      if(this->_functionSpace) {
        if(this->_functionSpace->isFormP()) {
          if(this->_domain.maxDim() == 3) {
            msg::warning << "CompoundField '" << this->_name << "': Perpendicular forms are meaningless on 3D domains" << msg::endl;
            error = true;
          }
        }
      }

      return !error;
    }
    else if constexpr(T_Form == field::Form::Form2) {
      bool error = false;

      if(this->_domain.minDim() < 2) {
        if(this->_functionSpace) {
          if(this->_functionSpace->isFormP()) {
            if(this->_domain.minDim() < 1) {
              msg::warning << "CompoundField '" << this->_name << "': " << std::string(NameOfForm< T_Form >::value) << " is not defined on points" << msg::endl;
              error = true;
            }
          }
          else {
            msg::warning << "CompoundField '" << this->_name << "': " << std::string(NameOfForm< T_Form >::value) << " is not defined on lines, neither on points" << msg::endl;
            error = true;
          }
        }
      }

      // Check that FormP is not defined on 3D mesh
      if(this->_functionSpace) {
        if(this->_functionSpace->isFormP()) {
          if(this->_domain.maxDim() == 3) {
            msg::warning << "CompoundField '" << this->_name << "': Parallel forms are meaningless on 3D domains" << msg::endl;
            error = true;
          }
        }
      }

      return !error;
    }
    else if constexpr(T_Form == field::Form::Form3) {
      bool error = false;

      // Check dimension of domain
      if(this->_domain.minDim() < 3) {
        if(this->_functionSpace) {
          if(this->_functionSpace->isFormP()) {
            // Ok
          }
          else {
            msg::warning << "CompoundField '" << this->_name << "': " + std::string(NameOfForm< T_Form >::value) << " is not defined on surfaces, lines, or points" << msg::endl;
            error = true;
          }
        }
      }

      // Check that FormP is not defined on 3D mesh
      if(this->_functionSpace) {
        if(this->_functionSpace->isFormP()) {
          if(this->_domain.maxDim() == 3) {
            msg::warning << "CompoundField '" << this->_name << "': Perpendicular forms are meaningless on 3D domains" << msg::endl;
            error = true;
          }
        }
      }

      return !error;
    }
  }

  template< class T_Scalar, field::Form T_Form, unsigned int T_NumFields >
  CompoundField< T_Scalar, T_Form, T_NumFields >::CompoundField() :
    FieldInterface< T_Scalar >(), _constraints(), _constraintsOnComponent(), _myself(*this)
  {
    for(unsigned int i = 1; i < T_NumFields; ++i) {
      FieldInterface< T_Scalar >::s_incrementNbrOfFields();
    }
  }

  template< class T_Scalar, field::Form T_Form, unsigned int T_NumFields >
  CompoundField< T_Scalar, T_Form, T_NumFields >::CompoundField(const std::string &name, const domain::Domain &domain, const field::FunctionSpaceOfForm< T_Form > &type, const unsigned int order, const std::string &model) :
    FieldInterface< T_Scalar >(name, domain, model), _constraints(), _constraintsOnComponent(), _myself(*this)
  {
    for(unsigned int i = 1; i < T_NumFields; ++i) {
      FieldInterface< T_Scalar >::s_incrementNbrOfFields();
    }
    
    if constexpr(T_Form == field::Form::Form0) {
      switch(type) {
      case field::FunctionSpaceTypeForm0::Lagrange:
        this->_functionSpace = new field::FunctionSpaceLagrange< scalar::Precision< T_Scalar > >(order);
        break;
      case field::FunctionSpaceTypeForm0::HierarchicalH1:
        this->_functionSpace = new field::FunctionSpaceHierarchicalH1< scalar::Precision< T_Scalar > >(order);
        break;
      default:
        throw common::Exception("Unknown function space for " + std::string(NameOfForm< T_Form >::value) + " field");
        break;
      }
    }
    else if constexpr(T_Form == field::Form::Form1) {
      switch(type) {
      case field::FunctionSpaceTypeForm1::HierarchicalHCurl:
        this->_functionSpace = new field::FunctionSpaceHierarchicalHCurl< scalar::Precision< T_Scalar > >(order);
        break;
      case field::FunctionSpaceTypeForm1::P_Lagrange:
        this->_functionSpace = new field::FunctionSpacePerpendicular< scalar::Precision< T_Scalar > >(field::FunctionSpaceLagrange< scalar::Precision< T_Scalar > >(order));
        break;
      case field::FunctionSpaceTypeForm1::P_HierarchicalH1:
        this->_functionSpace = new field::FunctionSpacePerpendicular< scalar::Precision< T_Scalar > >(field::FunctionSpaceHierarchicalH1< scalar::Precision< T_Scalar > >(order));
        break;
      case field::FunctionSpaceTypeForm1::D_Lagrange:
        this->_functionSpace = new field::FunctionSpaceGradH1< scalar::Precision< T_Scalar > >(field::FunctionSpaceLagrange< scalar::Precision< T_Scalar > >(order));
        break;
      case field::FunctionSpaceTypeForm1::D_HierarchicalH1:
        this->_functionSpace = new field::FunctionSpaceGradH1< scalar::Precision< T_Scalar > >(field::FunctionSpaceHierarchicalH1< scalar::Precision< T_Scalar > >(order));
        break;
      default:
        throw common::Exception("Unknown function space for " + std::string(NameOfForm< T_Form >::value) + " field");
        break;
      }
    }
    else if constexpr(T_Form == field::Form::Form2) {
      switch(type) {
      case field::FunctionSpaceTypeForm2::P_HierarchicalHCurl:
        this->_functionSpace = new field::FunctionSpaceParallel< scalar::Precision< T_Scalar > >(field::FunctionSpaceHierarchicalHCurl< scalar::Precision< T_Scalar > >(order));
        break;
      case field::FunctionSpaceTypeForm2::D_HierarchicalHCurl:
        this->_functionSpace = new field::FunctionSpaceCurlHCurl< scalar::Precision< T_Scalar > >(field::FunctionSpaceHierarchicalHCurl< scalar::Precision< T_Scalar > >(order));
        break;
      case field::FunctionSpaceTypeForm2::DP_Lagrange:
        this->_functionSpace = new field::FunctionSpaceCurlHCurl< scalar::Precision< T_Scalar > >(field::FunctionSpacePerpendicular< scalar::Precision< T_Scalar > >(field::FunctionSpaceLagrange< scalar::Precision< T_Scalar > >(order)));
        break;
      case field::FunctionSpaceTypeForm2::DP_HierarchicalH1:
        this->_functionSpace = new field::FunctionSpaceCurlHCurl< scalar::Precision< T_Scalar > >(field::FunctionSpacePerpendicular< scalar::Precision< T_Scalar > >(field::FunctionSpaceHierarchicalH1< scalar::Precision< T_Scalar > >(order)));
        break;
      default:
        throw common::Exception("Unknown function space for " + std::string(NameOfForm< T_Form >::value) + " field");
        break;
      }
    }
    else if constexpr(T_Form == field::Form::Form3) {
      switch(type) {
      case field::FunctionSpaceTypeForm3::Constant:
        this->_functionSpace = new field::FunctionSpaceConstant< scalar::Precision< T_Scalar > >(false);
        break;
      case field::FunctionSpaceTypeForm3::P_Constant:
        this->_functionSpace = new field::FunctionSpaceConstant< scalar::Precision< T_Scalar > >(true);
        break;
      case field::FunctionSpaceTypeForm3::DP_HierarchicalHCurl:
        this->_functionSpace = new field::FunctionSpaceDivHDiv< scalar::Precision< T_Scalar > >(field::FunctionSpaceParallel< scalar::Precision< T_Scalar > >(field::FunctionSpaceHierarchicalHCurl< scalar::Precision< T_Scalar > >(order)));
        break;
      default:
        throw common::Exception("Unknown function space for " + std::string(NameOfForm< T_Form >::value) + " field");
        break;
      }
    }
    
    _checkDomain();
  }

  template< class T_Scalar, field::Form T_Form, unsigned int T_NumFields >
  CompoundField< T_Scalar, T_Form, T_NumFields >::CompoundField(const CompoundField< T_Scalar, T_Form, T_NumFields > &other) :
    FieldInterface< T_Scalar >(other), _constraints(other._constraints), _constraintsOnComponent(), _myself(*this)
  {
    for(unsigned int i = 1; i < T_NumFields; ++i) {
      FieldInterface< T_Scalar >::s_incrementNbrOfFields();
    }
  }

  template< class T_Scalar, field::Form T_Form, unsigned int T_NumFields >
  CompoundField< T_Scalar, T_Form, T_NumFields >::~CompoundField()
  {
  }

  template< class T_Scalar, field::Form T_Form, unsigned int T_NumFields >
  void CompoundField< T_Scalar, T_Form, T_NumFields >::clear()
  {
    FieldInterface< T_Scalar >::clear();
    _constraints.clear();
    for(unsigned int i = 0; i < T_NumFields; ++i) {
      _constraintsOnComponent[i].clear();
    }
  }

  template< class T_Scalar, field::Form T_Form, unsigned int T_NumFields >
  field::Form CompoundField< T_Scalar, T_Form, T_NumFields >::form() const
  {
    return T_Form;
  }

  template< class T_Scalar, field::Form T_Form, unsigned int T_NumFields >
  unsigned int CompoundField< T_Scalar, T_Form, T_NumFields >::multiplicity() const
  {
    return T_NumFields;
  }

  template< class T_Scalar, field::Form T_Form, unsigned int T_NumFields >
  field::FunctionSpace< scalar::Precision< T_Scalar >, T_Form > *CompoundField< T_Scalar, T_Form, T_NumFields >::getFunctionSpace() const
  {
    return static_cast< field::FunctionSpace< scalar::Precision< T_Scalar >, T_Form > * >(this->_functionSpace);
  }

  template< class T_Scalar, field::Form T_Form, unsigned int T_NumFields >
  CompoundField< T_Scalar, T_Form, T_NumFields > &CompoundField< T_Scalar, T_Form, T_NumFields >::operator=(const CompoundField< T_Scalar, T_Form, T_NumFields > &other)
  {
    this->_copy(other);

    _constraints.clear();
    for(auto i = 0ULL; i < other._constraints.size(); ++i) {
      Constraint< T_Scalar, DegreeOfCompoundForm< T_Form >::value > constraint(other._constraints[i]._domain, other._constraints[i]._function);
      _constraints.push_back(constraint);
    }
    for(unsigned int i = 0; i < T_NumFields; ++i) {
      _constraintsOnComponent[i].clear();
      for(auto j = 0ULL; j < other._constraintsOnComponent[i].size(); ++j) {
        Constraint< T_Scalar, DegreeOfForm< T_Form >::value > constraint(other._constraintsOnComponent[i][j]._domain, other._constraintsOnComponent[i][j]._function);
        _constraintsOnComponent[i].push_back(constraint);
      }
    }
    _checkDomain();

    return *this;
  }

  template< class T_Scalar, field::Form T_Form, unsigned int T_NumFields >
  void CompoundField< T_Scalar, T_Form, T_NumFields >::addConstraint(const domain::Domain &domain, const function::Function< T_Scalar, DegreeOfCompoundForm< T_Form >::value > &function)
  {
    Constraint< T_Scalar, DegreeOfCompoundForm< T_Form >::value > constraint(domain, function);
    _constraints.push_back(constraint);
  }

  template< class T_Scalar, field::Form T_Form, unsigned int T_NumFields >
  std::vector< Constraint< T_Scalar, DegreeOfCompoundForm< T_Form >::value > > CompoundField< T_Scalar, T_Form, T_NumFields >::getConstraints() const
  {
    return _constraints;
  }
  
  template< class T_Scalar, field::Form T_Form, unsigned int T_NumFields >
  template< unsigned int T_Component, class >
  void CompoundField< T_Scalar, T_Form, T_NumFields >::addConstraintOnComponent(const domain::Domain &domain, const function::Function< T_Scalar, DegreeOfForm< T_Form >::value > &function)
  {
    Constraint< T_Scalar, DegreeOfForm< T_Form >::value > constraint(domain, function);
    _constraintsOnComponent[T_Component].push_back(constraint);
  }
  
  template< class T_Scalar, field::Form T_Form, unsigned int T_NumFields >
  template< unsigned int T_Component, class >
  std::vector< Constraint< T_Scalar, DegreeOfForm< T_Form >::value > > CompoundField< T_Scalar, T_Form, T_NumFields >::getConstraintsOnComponent() const
  {
    return _constraintsOnComponent[T_Component];
  }
  
  template< class T_Scalar, field::Form T_Form, unsigned int T_NumFields >
  void CompoundField< T_Scalar, T_Form, T_NumFields >::removeConstraints()
  {
    _constraints.clear();
    for(std::size_t i = 0; i < T_NumFields; ++i) {
      _constraintsOnComponent[i].clear();
    }
  }

  template< class T_Scalar, field::Form T_Form, unsigned int T_NumFields >
  const function::Function< T_Scalar, DegreeOfCompoundForm< T_Form >::value > &CompoundField< T_Scalar, T_Form, T_NumFields >::getEvaluableFunction() const
  {
    return _myself;
  }

  INSTANTIATE_CLASS_3(CompoundField, 4, 4, 2, TEMPLATE_ARGS(std::complex< double >, std::complex< float >, double, float), TEMPLATE_ARGS(field::Form::Form0, field::Form::Form1, field::Form::Form2, field::Form::Form3), TEMPLATE_ARGS(2, 3))
  
  template void CompoundField< std::complex< double >, field::Form::Form0, 2 >::addConstraintOnComponent< 0 >(const domain::Domain &domain, const function::Function< std::complex< double >, DegreeOfForm< field::Form::Form0 >::value > &function);
  template void CompoundField< std::complex< double >, field::Form::Form0, 2 >::addConstraintOnComponent< 1 >(const domain::Domain &domain, const function::Function< std::complex< double >, DegreeOfForm< field::Form::Form0 >::value > &function);
  template void CompoundField< std::complex< double >, field::Form::Form0, 3 >::addConstraintOnComponent< 0 >(const domain::Domain &domain, const function::Function< std::complex< double >, DegreeOfForm< field::Form::Form0 >::value > &function);
  template void CompoundField< std::complex< double >, field::Form::Form0, 3 >::addConstraintOnComponent< 1 >(const domain::Domain &domain, const function::Function< std::complex< double >, DegreeOfForm< field::Form::Form0 >::value > &function);
  template void CompoundField< std::complex< double >, field::Form::Form0, 3 >::addConstraintOnComponent< 2 >(const domain::Domain &domain, const function::Function< std::complex< double >, DegreeOfForm< field::Form::Form0 >::value > &function);
  
  template void CompoundField< std::complex< double >, field::Form::Form1, 2 >::addConstraintOnComponent< 0 >(const domain::Domain &domain, const function::Function< std::complex< double >, DegreeOfForm< field::Form::Form1 >::value > &function);
  template void CompoundField< std::complex< double >, field::Form::Form1, 2 >::addConstraintOnComponent< 1 >(const domain::Domain &domain, const function::Function< std::complex< double >, DegreeOfForm< field::Form::Form1 >::value > &function);
  template void CompoundField< std::complex< double >, field::Form::Form1, 3 >::addConstraintOnComponent< 0 >(const domain::Domain &domain, const function::Function< std::complex< double >, DegreeOfForm< field::Form::Form1 >::value > &function);
  template void CompoundField< std::complex< double >, field::Form::Form1, 3 >::addConstraintOnComponent< 1 >(const domain::Domain &domain, const function::Function< std::complex< double >, DegreeOfForm< field::Form::Form1 >::value > &function);
  template void CompoundField< std::complex< double >, field::Form::Form1, 3 >::addConstraintOnComponent< 2 >(const domain::Domain &domain, const function::Function< std::complex< double >, DegreeOfForm< field::Form::Form1 >::value > &function);
  
  template void CompoundField< std::complex< double >, field::Form::Form2, 2 >::addConstraintOnComponent< 0 >(const domain::Domain &domain, const function::Function< std::complex< double >, DegreeOfForm< field::Form::Form2 >::value > &function);
  template void CompoundField< std::complex< double >, field::Form::Form2, 2 >::addConstraintOnComponent< 1 >(const domain::Domain &domain, const function::Function< std::complex< double >, DegreeOfForm< field::Form::Form2 >::value > &function);
  template void CompoundField< std::complex< double >, field::Form::Form2, 3 >::addConstraintOnComponent< 0 >(const domain::Domain &domain, const function::Function< std::complex< double >, DegreeOfForm< field::Form::Form2 >::value > &function);
  template void CompoundField< std::complex< double >, field::Form::Form2, 3 >::addConstraintOnComponent< 1 >(const domain::Domain &domain, const function::Function< std::complex< double >, DegreeOfForm< field::Form::Form2 >::value > &function);
  template void CompoundField< std::complex< double >, field::Form::Form2, 3 >::addConstraintOnComponent< 2 >(const domain::Domain &domain, const function::Function< std::complex< double >, DegreeOfForm< field::Form::Form2 >::value > &function);
  
  template void CompoundField< std::complex< double >, field::Form::Form3, 2 >::addConstraintOnComponent< 0 >(const domain::Domain &domain, const function::Function< std::complex< double >, DegreeOfForm< field::Form::Form3 >::value > &function);
  template void CompoundField< std::complex< double >, field::Form::Form3, 2 >::addConstraintOnComponent< 1 >(const domain::Domain &domain, const function::Function< std::complex< double >, DegreeOfForm< field::Form::Form3 >::value > &function);
  template void CompoundField< std::complex< double >, field::Form::Form3, 3 >::addConstraintOnComponent< 0 >(const domain::Domain &domain, const function::Function< std::complex< double >, DegreeOfForm< field::Form::Form3 >::value > &function);
  template void CompoundField< std::complex< double >, field::Form::Form3, 3 >::addConstraintOnComponent< 1 >(const domain::Domain &domain, const function::Function< std::complex< double >, DegreeOfForm< field::Form::Form3 >::value > &function);
  template void CompoundField< std::complex< double >, field::Form::Form3, 3 >::addConstraintOnComponent< 2 >(const domain::Domain &domain, const function::Function< std::complex< double >, DegreeOfForm< field::Form::Form3 >::value > &function);
  
  
  template void CompoundField< std::complex< float >, field::Form::Form0, 2 >::addConstraintOnComponent< 0 >(const domain::Domain &domain, const function::Function< std::complex< float >, DegreeOfForm< field::Form::Form0 >::value > &function);
  template void CompoundField< std::complex< float >, field::Form::Form0, 2 >::addConstraintOnComponent< 1 >(const domain::Domain &domain, const function::Function< std::complex< float >, DegreeOfForm< field::Form::Form0 >::value > &function);
  template void CompoundField< std::complex< float >, field::Form::Form0, 3 >::addConstraintOnComponent< 0 >(const domain::Domain &domain, const function::Function< std::complex< float >, DegreeOfForm< field::Form::Form0 >::value > &function);
  template void CompoundField< std::complex< float >, field::Form::Form0, 3 >::addConstraintOnComponent< 1 >(const domain::Domain &domain, const function::Function< std::complex< float >, DegreeOfForm< field::Form::Form0 >::value > &function);
  template void CompoundField< std::complex< float >, field::Form::Form0, 3 >::addConstraintOnComponent< 2 >(const domain::Domain &domain, const function::Function< std::complex< float >, DegreeOfForm< field::Form::Form0 >::value > &function);
  
  template void CompoundField< std::complex< float >, field::Form::Form1, 2 >::addConstraintOnComponent< 0 >(const domain::Domain &domain, const function::Function< std::complex< float >, DegreeOfForm< field::Form::Form1 >::value > &function);
  template void CompoundField< std::complex< float >, field::Form::Form1, 2 >::addConstraintOnComponent< 1 >(const domain::Domain &domain, const function::Function< std::complex< float >, DegreeOfForm< field::Form::Form1 >::value > &function);
  template void CompoundField< std::complex< float >, field::Form::Form1, 3 >::addConstraintOnComponent< 0 >(const domain::Domain &domain, const function::Function< std::complex< float >, DegreeOfForm< field::Form::Form1 >::value > &function);
  template void CompoundField< std::complex< float >, field::Form::Form1, 3 >::addConstraintOnComponent< 1 >(const domain::Domain &domain, const function::Function< std::complex< float >, DegreeOfForm< field::Form::Form1 >::value > &function);
  template void CompoundField< std::complex< float >, field::Form::Form1, 3 >::addConstraintOnComponent< 2 >(const domain::Domain &domain, const function::Function< std::complex< float >, DegreeOfForm< field::Form::Form1 >::value > &function);
  
  template void CompoundField< std::complex< float >, field::Form::Form2, 2 >::addConstraintOnComponent< 0 >(const domain::Domain &domain, const function::Function< std::complex< float >, DegreeOfForm< field::Form::Form2 >::value > &function);
  template void CompoundField< std::complex< float >, field::Form::Form2, 2 >::addConstraintOnComponent< 1 >(const domain::Domain &domain, const function::Function< std::complex< float >, DegreeOfForm< field::Form::Form2 >::value > &function);
  template void CompoundField< std::complex< float >, field::Form::Form2, 3 >::addConstraintOnComponent< 0 >(const domain::Domain &domain, const function::Function< std::complex< float >, DegreeOfForm< field::Form::Form2 >::value > &function);
  template void CompoundField< std::complex< float >, field::Form::Form2, 3 >::addConstraintOnComponent< 1 >(const domain::Domain &domain, const function::Function< std::complex< float >, DegreeOfForm< field::Form::Form2 >::value > &function);
  template void CompoundField< std::complex< float >, field::Form::Form2, 3 >::addConstraintOnComponent< 2 >(const domain::Domain &domain, const function::Function< std::complex< float >, DegreeOfForm< field::Form::Form2 >::value > &function);
  
  template void CompoundField< std::complex< float >, field::Form::Form3, 2 >::addConstraintOnComponent< 0 >(const domain::Domain &domain, const function::Function< std::complex< float >, DegreeOfForm< field::Form::Form3 >::value > &function);
  template void CompoundField< std::complex< float >, field::Form::Form3, 2 >::addConstraintOnComponent< 1 >(const domain::Domain &domain, const function::Function< std::complex< float >, DegreeOfForm< field::Form::Form3 >::value > &function);
  template void CompoundField< std::complex< float >, field::Form::Form3, 3 >::addConstraintOnComponent< 0 >(const domain::Domain &domain, const function::Function< std::complex< float >, DegreeOfForm< field::Form::Form3 >::value > &function);
  template void CompoundField< std::complex< float >, field::Form::Form3, 3 >::addConstraintOnComponent< 1 >(const domain::Domain &domain, const function::Function< std::complex< float >, DegreeOfForm< field::Form::Form3 >::value > &function);
  template void CompoundField< std::complex< float >, field::Form::Form3, 3 >::addConstraintOnComponent< 2 >(const domain::Domain &domain, const function::Function< std::complex< float >, DegreeOfForm< field::Form::Form3 >::value > &function);
  
  
  template void CompoundField< double, field::Form::Form0, 2 >::addConstraintOnComponent< 0 >(const domain::Domain &domain, const function::Function< double, DegreeOfForm< field::Form::Form0 >::value > &function);
  template void CompoundField< double, field::Form::Form0, 2 >::addConstraintOnComponent< 1 >(const domain::Domain &domain, const function::Function< double, DegreeOfForm< field::Form::Form0 >::value > &function);
  template void CompoundField< double, field::Form::Form0, 3 >::addConstraintOnComponent< 0 >(const domain::Domain &domain, const function::Function< double, DegreeOfForm< field::Form::Form0 >::value > &function);
  template void CompoundField< double, field::Form::Form0, 3 >::addConstraintOnComponent< 1 >(const domain::Domain &domain, const function::Function< double, DegreeOfForm< field::Form::Form0 >::value > &function);
  template void CompoundField< double, field::Form::Form0, 3 >::addConstraintOnComponent< 2 >(const domain::Domain &domain, const function::Function< double, DegreeOfForm< field::Form::Form0 >::value > &function);
  
  template void CompoundField< double, field::Form::Form1, 2 >::addConstraintOnComponent< 0 >(const domain::Domain &domain, const function::Function< double, DegreeOfForm< field::Form::Form1 >::value > &function);
  template void CompoundField< double, field::Form::Form1, 2 >::addConstraintOnComponent< 1 >(const domain::Domain &domain, const function::Function< double, DegreeOfForm< field::Form::Form1 >::value > &function);
  template void CompoundField< double, field::Form::Form1, 3 >::addConstraintOnComponent< 0 >(const domain::Domain &domain, const function::Function< double, DegreeOfForm< field::Form::Form1 >::value > &function);
  template void CompoundField< double, field::Form::Form1, 3 >::addConstraintOnComponent< 1 >(const domain::Domain &domain, const function::Function< double, DegreeOfForm< field::Form::Form1 >::value > &function);
  template void CompoundField< double, field::Form::Form1, 3 >::addConstraintOnComponent< 2 >(const domain::Domain &domain, const function::Function< double, DegreeOfForm< field::Form::Form1 >::value > &function);
  
  template void CompoundField< double, field::Form::Form2, 2 >::addConstraintOnComponent< 0 >(const domain::Domain &domain, const function::Function< double, DegreeOfForm< field::Form::Form2 >::value > &function);
  template void CompoundField< double, field::Form::Form2, 2 >::addConstraintOnComponent< 1 >(const domain::Domain &domain, const function::Function< double, DegreeOfForm< field::Form::Form2 >::value > &function);
  template void CompoundField< double, field::Form::Form2, 3 >::addConstraintOnComponent< 0 >(const domain::Domain &domain, const function::Function< double, DegreeOfForm< field::Form::Form2 >::value > &function);
  template void CompoundField< double, field::Form::Form2, 3 >::addConstraintOnComponent< 1 >(const domain::Domain &domain, const function::Function< double, DegreeOfForm< field::Form::Form2 >::value > &function);
  template void CompoundField< double, field::Form::Form2, 3 >::addConstraintOnComponent< 2 >(const domain::Domain &domain, const function::Function< double, DegreeOfForm< field::Form::Form2 >::value > &function);
  
  template void CompoundField< double, field::Form::Form3, 2 >::addConstraintOnComponent< 0 >(const domain::Domain &domain, const function::Function< double, DegreeOfForm< field::Form::Form3 >::value > &function);
  template void CompoundField< double, field::Form::Form3, 2 >::addConstraintOnComponent< 1 >(const domain::Domain &domain, const function::Function< double, DegreeOfForm< field::Form::Form3 >::value > &function);
  template void CompoundField< double, field::Form::Form3, 3 >::addConstraintOnComponent< 0 >(const domain::Domain &domain, const function::Function< double, DegreeOfForm< field::Form::Form3 >::value > &function);
  template void CompoundField< double, field::Form::Form3, 3 >::addConstraintOnComponent< 1 >(const domain::Domain &domain, const function::Function< double, DegreeOfForm< field::Form::Form3 >::value > &function);
  template void CompoundField< double, field::Form::Form3, 3 >::addConstraintOnComponent< 2 >(const domain::Domain &domain, const function::Function< double, DegreeOfForm< field::Form::Form3 >::value > &function);
  
  
  template void CompoundField< float, field::Form::Form0, 2 >::addConstraintOnComponent< 0 >(const domain::Domain &domain, const function::Function< float, DegreeOfForm< field::Form::Form0 >::value > &function);
  template void CompoundField< float, field::Form::Form0, 2 >::addConstraintOnComponent< 1 >(const domain::Domain &domain, const function::Function< float, DegreeOfForm< field::Form::Form0 >::value > &function);
  template void CompoundField< float, field::Form::Form0, 3 >::addConstraintOnComponent< 0 >(const domain::Domain &domain, const function::Function< float, DegreeOfForm< field::Form::Form0 >::value > &function);
  template void CompoundField< float, field::Form::Form0, 3 >::addConstraintOnComponent< 1 >(const domain::Domain &domain, const function::Function< float, DegreeOfForm< field::Form::Form0 >::value > &function);
  template void CompoundField< float, field::Form::Form0, 3 >::addConstraintOnComponent< 2 >(const domain::Domain &domain, const function::Function< float, DegreeOfForm< field::Form::Form0 >::value > &function);
  
  template void CompoundField< float, field::Form::Form1, 2 >::addConstraintOnComponent< 0 >(const domain::Domain &domain, const function::Function< float, DegreeOfForm< field::Form::Form1 >::value > &function);
  template void CompoundField< float, field::Form::Form1, 2 >::addConstraintOnComponent< 1 >(const domain::Domain &domain, const function::Function< float, DegreeOfForm< field::Form::Form1 >::value > &function);
  template void CompoundField< float, field::Form::Form1, 3 >::addConstraintOnComponent< 0 >(const domain::Domain &domain, const function::Function< float, DegreeOfForm< field::Form::Form1 >::value > &function);
  template void CompoundField< float, field::Form::Form1, 3 >::addConstraintOnComponent< 1 >(const domain::Domain &domain, const function::Function< float, DegreeOfForm< field::Form::Form1 >::value > &function);
  template void CompoundField< float, field::Form::Form1, 3 >::addConstraintOnComponent< 2 >(const domain::Domain &domain, const function::Function< float, DegreeOfForm< field::Form::Form1 >::value > &function);
  
  template void CompoundField< float, field::Form::Form2, 2 >::addConstraintOnComponent< 0 >(const domain::Domain &domain, const function::Function< float, DegreeOfForm< field::Form::Form2 >::value > &function);
  template void CompoundField< float, field::Form::Form2, 2 >::addConstraintOnComponent< 1 >(const domain::Domain &domain, const function::Function< float, DegreeOfForm< field::Form::Form2 >::value > &function);
  template void CompoundField< float, field::Form::Form2, 3 >::addConstraintOnComponent< 0 >(const domain::Domain &domain, const function::Function< float, DegreeOfForm< field::Form::Form2 >::value > &function);
  template void CompoundField< float, field::Form::Form2, 3 >::addConstraintOnComponent< 1 >(const domain::Domain &domain, const function::Function< float, DegreeOfForm< field::Form::Form2 >::value > &function);
  template void CompoundField< float, field::Form::Form2, 3 >::addConstraintOnComponent< 2 >(const domain::Domain &domain, const function::Function< float, DegreeOfForm< field::Form::Form2 >::value > &function);
  
  template void CompoundField< float, field::Form::Form3, 2 >::addConstraintOnComponent< 0 >(const domain::Domain &domain, const function::Function< float, DegreeOfForm< field::Form::Form3 >::value > &function);
  template void CompoundField< float, field::Form::Form3, 2 >::addConstraintOnComponent< 1 >(const domain::Domain &domain, const function::Function< float, DegreeOfForm< field::Form::Form3 >::value > &function);
  template void CompoundField< float, field::Form::Form3, 3 >::addConstraintOnComponent< 0 >(const domain::Domain &domain, const function::Function< float, DegreeOfForm< field::Form::Form3 >::value > &function);
  template void CompoundField< float, field::Form::Form3, 3 >::addConstraintOnComponent< 1 >(const domain::Domain &domain, const function::Function< float, DegreeOfForm< field::Form::Form3 >::value > &function);
  template void CompoundField< float, field::Form::Form3, 3 >::addConstraintOnComponent< 2 >(const domain::Domain &domain, const function::Function< float, DegreeOfForm< field::Form::Form3 >::value > &function);
  
  
  template std::vector< Constraint< std::complex< double >, DegreeOfForm< field::Form::Form0 >::value > > CompoundField< std::complex< double >, field::Form::Form0, 2 >::getConstraintsOnComponent< 0 >() const;
  template std::vector< Constraint< std::complex< double >, DegreeOfForm< field::Form::Form0 >::value > > CompoundField< std::complex< double >, field::Form::Form0, 2 >::getConstraintsOnComponent< 1 >() const;
  template std::vector< Constraint< std::complex< double >, DegreeOfForm< field::Form::Form0 >::value > > CompoundField< std::complex< double >, field::Form::Form0, 3 >::getConstraintsOnComponent< 0 >() const;
  template std::vector< Constraint< std::complex< double >, DegreeOfForm< field::Form::Form0 >::value > > CompoundField< std::complex< double >, field::Form::Form0, 3 >::getConstraintsOnComponent< 1 >() const;
  template std::vector< Constraint< std::complex< double >, DegreeOfForm< field::Form::Form0 >::value > > CompoundField< std::complex< double >, field::Form::Form0, 3 >::getConstraintsOnComponent< 2 >() const;
  
  template std::vector< Constraint< std::complex< double >, DegreeOfForm< field::Form::Form1 >::value > > CompoundField< std::complex< double >, field::Form::Form1, 2 >::getConstraintsOnComponent< 0 >() const;
  template std::vector< Constraint< std::complex< double >, DegreeOfForm< field::Form::Form1 >::value > > CompoundField< std::complex< double >, field::Form::Form1, 2 >::getConstraintsOnComponent< 1 >() const;
  template std::vector< Constraint< std::complex< double >, DegreeOfForm< field::Form::Form1 >::value > > CompoundField< std::complex< double >, field::Form::Form1, 3 >::getConstraintsOnComponent< 0 >() const;
  template std::vector< Constraint< std::complex< double >, DegreeOfForm< field::Form::Form1 >::value > > CompoundField< std::complex< double >, field::Form::Form1, 3 >::getConstraintsOnComponent< 1 >() const;
  template std::vector< Constraint< std::complex< double >, DegreeOfForm< field::Form::Form1 >::value > > CompoundField< std::complex< double >, field::Form::Form1, 3 >::getConstraintsOnComponent< 2 >() const;
  
  template std::vector< Constraint< std::complex< double >, DegreeOfForm< field::Form::Form2 >::value > > CompoundField< std::complex< double >, field::Form::Form2, 2 >::getConstraintsOnComponent< 0 >() const;
  template std::vector< Constraint< std::complex< double >, DegreeOfForm< field::Form::Form2 >::value > > CompoundField< std::complex< double >, field::Form::Form2, 2 >::getConstraintsOnComponent< 1 >() const;
  template std::vector< Constraint< std::complex< double >, DegreeOfForm< field::Form::Form2 >::value > > CompoundField< std::complex< double >, field::Form::Form2, 3 >::getConstraintsOnComponent< 0 >() const;
  template std::vector< Constraint< std::complex< double >, DegreeOfForm< field::Form::Form2 >::value > > CompoundField< std::complex< double >, field::Form::Form2, 3 >::getConstraintsOnComponent< 1 >() const;
  template std::vector< Constraint< std::complex< double >, DegreeOfForm< field::Form::Form2 >::value > > CompoundField< std::complex< double >, field::Form::Form2, 3 >::getConstraintsOnComponent< 2 >() const;
  
  template std::vector< Constraint< std::complex< double >, DegreeOfForm< field::Form::Form3 >::value > > CompoundField< std::complex< double >, field::Form::Form3, 2 >::getConstraintsOnComponent< 0 >() const;
  template std::vector< Constraint< std::complex< double >, DegreeOfForm< field::Form::Form3 >::value > > CompoundField< std::complex< double >, field::Form::Form3, 2 >::getConstraintsOnComponent< 1 >() const;
  template std::vector< Constraint< std::complex< double >, DegreeOfForm< field::Form::Form3 >::value > > CompoundField< std::complex< double >, field::Form::Form3, 3 >::getConstraintsOnComponent< 0 >() const;
  template std::vector< Constraint< std::complex< double >, DegreeOfForm< field::Form::Form3 >::value > > CompoundField< std::complex< double >, field::Form::Form3, 3 >::getConstraintsOnComponent< 1 >() const;
  template std::vector< Constraint< std::complex< double >, DegreeOfForm< field::Form::Form3 >::value > > CompoundField< std::complex< double >, field::Form::Form3, 3 >::getConstraintsOnComponent< 2 >() const;
  
  
  template std::vector< Constraint< std::complex< float >, DegreeOfForm< field::Form::Form0 >::value > > CompoundField< std::complex< float >, field::Form::Form0, 2 >::getConstraintsOnComponent< 0 >() const;
  template std::vector< Constraint< std::complex< float >, DegreeOfForm< field::Form::Form0 >::value > > CompoundField< std::complex< float >, field::Form::Form0, 2 >::getConstraintsOnComponent< 1 >() const;
  template std::vector< Constraint< std::complex< float >, DegreeOfForm< field::Form::Form0 >::value > > CompoundField< std::complex< float >, field::Form::Form0, 3 >::getConstraintsOnComponent< 0 >() const;
  template std::vector< Constraint< std::complex< float >, DegreeOfForm< field::Form::Form0 >::value > > CompoundField< std::complex< float >, field::Form::Form0, 3 >::getConstraintsOnComponent< 1 >() const;
  template std::vector< Constraint< std::complex< float >, DegreeOfForm< field::Form::Form0 >::value > > CompoundField< std::complex< float >, field::Form::Form0, 3 >::getConstraintsOnComponent< 2 >() const;
  
  template std::vector< Constraint< std::complex< float >, DegreeOfForm< field::Form::Form1 >::value > > CompoundField< std::complex< float >, field::Form::Form1, 2 >::getConstraintsOnComponent< 0 >() const;
  template std::vector< Constraint< std::complex< float >, DegreeOfForm< field::Form::Form1 >::value > > CompoundField< std::complex< float >, field::Form::Form1, 2 >::getConstraintsOnComponent< 1 >() const;
  template std::vector< Constraint< std::complex< float >, DegreeOfForm< field::Form::Form1 >::value > > CompoundField< std::complex< float >, field::Form::Form1, 3 >::getConstraintsOnComponent< 0 >() const;
  template std::vector< Constraint< std::complex< float >, DegreeOfForm< field::Form::Form1 >::value > > CompoundField< std::complex< float >, field::Form::Form1, 3 >::getConstraintsOnComponent< 1 >() const;
  template std::vector< Constraint< std::complex< float >, DegreeOfForm< field::Form::Form1 >::value > > CompoundField< std::complex< float >, field::Form::Form1, 3 >::getConstraintsOnComponent< 2 >() const;
  
  template std::vector< Constraint< std::complex< float >, DegreeOfForm< field::Form::Form2 >::value > > CompoundField< std::complex< float >, field::Form::Form2, 2 >::getConstraintsOnComponent< 0 >() const;
  template std::vector< Constraint< std::complex< float >, DegreeOfForm< field::Form::Form2 >::value > > CompoundField< std::complex< float >, field::Form::Form2, 2 >::getConstraintsOnComponent< 1 >() const;
  template std::vector< Constraint< std::complex< float >, DegreeOfForm< field::Form::Form2 >::value > > CompoundField< std::complex< float >, field::Form::Form2, 3 >::getConstraintsOnComponent< 0 >() const;
  template std::vector< Constraint< std::complex< float >, DegreeOfForm< field::Form::Form2 >::value > > CompoundField< std::complex< float >, field::Form::Form2, 3 >::getConstraintsOnComponent< 1 >() const;
  template std::vector< Constraint< std::complex< float >, DegreeOfForm< field::Form::Form2 >::value > > CompoundField< std::complex< float >, field::Form::Form2, 3 >::getConstraintsOnComponent< 2 >() const;
  
  template std::vector< Constraint< std::complex< float >, DegreeOfForm< field::Form::Form3 >::value > > CompoundField< std::complex< float >, field::Form::Form3, 2 >::getConstraintsOnComponent< 0 >() const;
  template std::vector< Constraint< std::complex< float >, DegreeOfForm< field::Form::Form3 >::value > > CompoundField< std::complex< float >, field::Form::Form3, 2 >::getConstraintsOnComponent< 1 >() const;
  template std::vector< Constraint< std::complex< float >, DegreeOfForm< field::Form::Form3 >::value > > CompoundField< std::complex< float >, field::Form::Form3, 3 >::getConstraintsOnComponent< 0 >() const;
  template std::vector< Constraint< std::complex< float >, DegreeOfForm< field::Form::Form3 >::value > > CompoundField< std::complex< float >, field::Form::Form3, 3 >::getConstraintsOnComponent< 1 >() const;
  template std::vector< Constraint< std::complex< float >, DegreeOfForm< field::Form::Form3 >::value > > CompoundField< std::complex< float >, field::Form::Form3, 3 >::getConstraintsOnComponent< 2 >() const;
  
  
  template std::vector< Constraint< double, DegreeOfForm< field::Form::Form0 >::value > > CompoundField< double, field::Form::Form0, 2 >::getConstraintsOnComponent< 0 >() const;
  template std::vector< Constraint< double, DegreeOfForm< field::Form::Form0 >::value > > CompoundField< double, field::Form::Form0, 2 >::getConstraintsOnComponent< 1 >() const;
  template std::vector< Constraint< double, DegreeOfForm< field::Form::Form0 >::value > > CompoundField< double, field::Form::Form0, 3 >::getConstraintsOnComponent< 0 >() const;
  template std::vector< Constraint< double, DegreeOfForm< field::Form::Form0 >::value > > CompoundField< double, field::Form::Form0, 3 >::getConstraintsOnComponent< 1 >() const;
  template std::vector< Constraint< double, DegreeOfForm< field::Form::Form0 >::value > > CompoundField< double, field::Form::Form0, 3 >::getConstraintsOnComponent< 2 >() const;
  
  template std::vector< Constraint< double, DegreeOfForm< field::Form::Form1 >::value > > CompoundField< double, field::Form::Form1, 2 >::getConstraintsOnComponent< 0 >() const;
  template std::vector< Constraint< double, DegreeOfForm< field::Form::Form1 >::value > > CompoundField< double, field::Form::Form1, 2 >::getConstraintsOnComponent< 1 >() const;
  template std::vector< Constraint< double, DegreeOfForm< field::Form::Form1 >::value > > CompoundField< double, field::Form::Form1, 3 >::getConstraintsOnComponent< 0 >() const;
  template std::vector< Constraint< double, DegreeOfForm< field::Form::Form1 >::value > > CompoundField< double, field::Form::Form1, 3 >::getConstraintsOnComponent< 1 >() const;
  template std::vector< Constraint< double, DegreeOfForm< field::Form::Form1 >::value > > CompoundField< double, field::Form::Form1, 3 >::getConstraintsOnComponent< 2 >() const;
  
  template std::vector< Constraint< double, DegreeOfForm< field::Form::Form2 >::value > > CompoundField< double, field::Form::Form2, 2 >::getConstraintsOnComponent< 0 >() const;
  template std::vector< Constraint< double, DegreeOfForm< field::Form::Form2 >::value > > CompoundField< double, field::Form::Form2, 2 >::getConstraintsOnComponent< 1 >() const;
  template std::vector< Constraint< double, DegreeOfForm< field::Form::Form2 >::value > > CompoundField< double, field::Form::Form2, 3 >::getConstraintsOnComponent< 0 >() const;
  template std::vector< Constraint< double, DegreeOfForm< field::Form::Form2 >::value > > CompoundField< double, field::Form::Form2, 3 >::getConstraintsOnComponent< 1 >() const;
  template std::vector< Constraint< double, DegreeOfForm< field::Form::Form2 >::value > > CompoundField< double, field::Form::Form2, 3 >::getConstraintsOnComponent< 2 >() const;
  
  template std::vector< Constraint< double, DegreeOfForm< field::Form::Form3 >::value > > CompoundField< double, field::Form::Form3, 2 >::getConstraintsOnComponent< 0 >() const;
  template std::vector< Constraint< double, DegreeOfForm< field::Form::Form3 >::value > > CompoundField< double, field::Form::Form3, 2 >::getConstraintsOnComponent< 1 >() const;
  template std::vector< Constraint< double, DegreeOfForm< field::Form::Form3 >::value > > CompoundField< double, field::Form::Form3, 3 >::getConstraintsOnComponent< 0 >() const;
  template std::vector< Constraint< double, DegreeOfForm< field::Form::Form3 >::value > > CompoundField< double, field::Form::Form3, 3 >::getConstraintsOnComponent< 1 >() const;
  template std::vector< Constraint< double, DegreeOfForm< field::Form::Form3 >::value > > CompoundField< double, field::Form::Form3, 3 >::getConstraintsOnComponent< 2 >() const;
  
  
  template std::vector< Constraint< float, DegreeOfForm< field::Form::Form0 >::value > > CompoundField< float, field::Form::Form0, 2 >::getConstraintsOnComponent< 0 >() const;
  template std::vector< Constraint< float, DegreeOfForm< field::Form::Form0 >::value > > CompoundField< float, field::Form::Form0, 2 >::getConstraintsOnComponent< 1 >() const;
  template std::vector< Constraint< float, DegreeOfForm< field::Form::Form0 >::value > > CompoundField< float, field::Form::Form0, 3 >::getConstraintsOnComponent< 0 >() const;
  template std::vector< Constraint< float, DegreeOfForm< field::Form::Form0 >::value > > CompoundField< float, field::Form::Form0, 3 >::getConstraintsOnComponent< 1 >() const;
  template std::vector< Constraint< float, DegreeOfForm< field::Form::Form0 >::value > > CompoundField< float, field::Form::Form0, 3 >::getConstraintsOnComponent< 2 >() const;
  
  template std::vector< Constraint< float, DegreeOfForm< field::Form::Form1 >::value > > CompoundField< float, field::Form::Form1, 2 >::getConstraintsOnComponent< 0 >() const;
  template std::vector< Constraint< float, DegreeOfForm< field::Form::Form1 >::value > > CompoundField< float, field::Form::Form1, 2 >::getConstraintsOnComponent< 1 >() const;
  template std::vector< Constraint< float, DegreeOfForm< field::Form::Form1 >::value > > CompoundField< float, field::Form::Form1, 3 >::getConstraintsOnComponent< 0 >() const;
  template std::vector< Constraint< float, DegreeOfForm< field::Form::Form1 >::value > > CompoundField< float, field::Form::Form1, 3 >::getConstraintsOnComponent< 1 >() const;
  template std::vector< Constraint< float, DegreeOfForm< field::Form::Form1 >::value > > CompoundField< float, field::Form::Form1, 3 >::getConstraintsOnComponent< 2 >() const;
  
  template std::vector< Constraint< float, DegreeOfForm< field::Form::Form2 >::value > > CompoundField< float, field::Form::Form2, 2 >::getConstraintsOnComponent< 0 >() const;
  template std::vector< Constraint< float, DegreeOfForm< field::Form::Form2 >::value > > CompoundField< float, field::Form::Form2, 2 >::getConstraintsOnComponent< 1 >() const;
  template std::vector< Constraint< float, DegreeOfForm< field::Form::Form2 >::value > > CompoundField< float, field::Form::Form2, 3 >::getConstraintsOnComponent< 0 >() const;
  template std::vector< Constraint< float, DegreeOfForm< field::Form::Form2 >::value > > CompoundField< float, field::Form::Form2, 3 >::getConstraintsOnComponent< 1 >() const;
  template std::vector< Constraint< float, DegreeOfForm< field::Form::Form2 >::value > > CompoundField< float, field::Form::Form2, 3 >::getConstraintsOnComponent< 2 >() const;
  
  template std::vector< Constraint< float, DegreeOfForm< field::Form::Form3 >::value > > CompoundField< float, field::Form::Form3, 2 >::getConstraintsOnComponent< 0 >() const;
  template std::vector< Constraint< float, DegreeOfForm< field::Form::Form3 >::value > > CompoundField< float, field::Form::Form3, 2 >::getConstraintsOnComponent< 1 >() const;
  template std::vector< Constraint< float, DegreeOfForm< field::Form::Form3 >::value > > CompoundField< float, field::Form::Form3, 3 >::getConstraintsOnComponent< 0 >() const;
  template std::vector< Constraint< float, DegreeOfForm< field::Form::Form3 >::value > > CompoundField< float, field::Form::Form3, 3 >::getConstraintsOnComponent< 1 >() const;
  template std::vector< Constraint< float, DegreeOfForm< field::Form::Form3 >::value > > CompoundField< float, field::Form::Form3, 3 >::getConstraintsOnComponent< 2 >() const;


} // namespace gmshfem::field
