// GmshFEM - Copyright (C) 2019-2022, A. Royer, E. Béchet, C. Geuzaine, Université de Liège
//
// See the LICENSE.txt file for license information. Please report all
// issues on https://gitlab.onelab.info/gmsh/fem/issues

#include "Field.h"

#include "Exception.h"
#include "instantiate.h"

namespace gmshfem::field
{


  template< class T_Scalar, field::Form T_Form >
  bool Field< T_Scalar, T_Form >::_checkDomain() const
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
            msg::warning << "Field '" << this->_name << "': " + std::string(NameOfForm< T_Form >::value) << " is not defined on points" << msg::endl;
            error = true;
          }
        }
      }

      // Check that FormP is not defined on 3D mesh
      if(this->_functionSpace) {
        if(this->_functionSpace->isFormP()) {
          if(this->_domain.maxDim() == 3) {
            msg::warning << "Field '" << this->_name << "': Perpendicular forms are meaningless on 3D domains" << msg::endl;
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
              msg::warning << "Field '" << this->_name << "': " << std::string(NameOfForm< T_Form >::value) << " is not defined on points" << msg::endl;
              error = true;
            }
          }
          else {
            msg::warning << "Field '" << this->_name << "': " << std::string(NameOfForm< T_Form >::value) << " is not defined on lines, neither on points" << msg::endl;
            error = true;
          }
        }
      }

      // Check that FormP is not defined on 3D mesh
      if(this->_functionSpace) {
        if(this->_functionSpace->isFormP()) {
          if(this->_domain.maxDim() == 3) {
            msg::warning << "Field '" << this->_name << "': Parallel forms are meaningless on 3D domains" << msg::endl;
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
            msg::warning << "Field '" << this->_name << "': " + std::string(NameOfForm< T_Form >::value) << " is not defined on surfaces, lines, or points" << msg::endl;
            error = true;
          }
        }
      }

      // Check that FormP is not defined on 3D mesh
      if(this->_functionSpace) {
        if(this->_functionSpace->isFormP()) {
          if(this->_domain.maxDim() == 3) {
            msg::warning << "Field '" << this->_name << "': Perpendicular forms are meaningless on 3D domains" << msg::endl;
            error = true;
          }
        }
      }

      return !error;
    }
  }

  template< class T_Scalar, field::Form T_Form >
  Field< T_Scalar, T_Form >::Field() :
    FieldInterface< T_Scalar >(), _constraints(), _myself(*this)
  {
  }

  template< class T_Scalar, field::Form T_Form >
  Field< T_Scalar, T_Form >::Field(const std::string &name, const domain::Domain &domain, const field::FunctionSpaceOfForm< T_Form > &type, const unsigned int order, const std::string &model) :
    FieldInterface< T_Scalar >(name, domain, model), _constraints(), _myself(*this)
  {
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

  template< class T_Scalar, field::Form T_Form >
  Field< T_Scalar, T_Form >::Field(const std::string &name, const domain::Domain &domain, const FunctionSpace< scalar::Precision< T_Scalar >, T_Form > &functionSpace, const std::string &model) :
    FieldInterface< T_Scalar >(name, domain, model), _constraints(), _myself(*this)
  {
    this->_functionSpace = functionSpace.copy();
    _checkDomain();
  }

  template< class T_Scalar, field::Form T_Form >
  Field< T_Scalar, T_Form >::Field(const Field< T_Scalar, T_Form > &other) :
    FieldInterface< T_Scalar >(other), _constraints(other._constraints), _myself(*this)
  {
    _checkDomain();
  }

  template< class T_Scalar, field::Form T_Form >
  Field< T_Scalar, T_Form >::~Field()
  {
  }

  template< class T_Scalar, field::Form T_Form >
  void Field< T_Scalar, T_Form >::clear()
  {
    FieldInterface< T_Scalar >::clear();
    _constraints.clear();
  }

  template< class T_Scalar, field::Form T_Form >
  field::Form Field< T_Scalar, T_Form >::form() const
  {
    return T_Form;
  }

  template< class T_Scalar, field::Form T_Form >
  unsigned int Field< T_Scalar, T_Form >::multiplicity() const
  {
    return 1;
  }

  template< class T_Scalar, field::Form T_Form >
  field::FunctionSpace< scalar::Precision< T_Scalar >, T_Form > *Field< T_Scalar, T_Form >::getFunctionSpace() const
  {
    return static_cast< field::FunctionSpace< scalar::Precision< T_Scalar >, T_Form > * >(this->_functionSpace);
  }

  template< class T_Scalar, field::Form T_Form >
  Field< T_Scalar, T_Form > &Field< T_Scalar, T_Form >::operator=(const Field< T_Scalar, T_Form > &other)
  {
    this->_copy(other);

    _constraints.clear();
    for(auto i = 0ULL; i < other._constraints.size(); ++i) {
      Constraint< T_Scalar, DegreeOfForm< T_Form >::value > constraint(other._constraints[i]._domain, other._constraints[i]._function);
      _constraints.push_back(constraint);
    }
    _checkDomain();

    return *this;
  }

  template< class T_Scalar, field::Form T_Form >
  void Field< T_Scalar, T_Form >::addConstraint(const domain::Domain &domain, const function::Function< T_Scalar, DegreeOfForm< T_Form >::value > &function)
  {
    Constraint< T_Scalar, DegreeOfForm< T_Form >::value > constraint(domain, function);
    _constraints.push_back(constraint);
  }

  template< class T_Scalar, field::Form T_Form >
  std::vector< Constraint< T_Scalar, DegreeOfForm< T_Form >::value > > Field< T_Scalar, T_Form >::getConstraints() const
  {
    return _constraints;
  }
  
  template< class T_Scalar, field::Form T_Form >
  void Field< T_Scalar, T_Form >::removeConstraints()
  {
    _constraints.clear();
  }

  template< class T_Scalar, field::Form T_Form >
  const function::Function< T_Scalar, field::DegreeOfForm< T_Form >::value > &Field< T_Scalar, T_Form >::getEvaluableFunction() const
  {
    return _myself;
  }

  INSTANTIATE_CLASS_2(Field, 4, 4, TEMPLATE_ARGS(std::complex< double >, std::complex< float >, double, float), TEMPLATE_ARGS(field::Form::Form0, field::Form::Form1, field::Form::Form2, field::Form::Form3))


} // namespace gmshfem::field
