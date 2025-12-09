// GmshFEM - Copyright (C) 2019-2022, A. Royer, E. Béchet, C. Geuzaine, Université de Liège
//
// See the LICENSE.txt file for license information. Please report all
// issues on https://gitlab.onelab.info/gmsh/fem/issues

#include "Formulation.h"

#include "Assembler.h"
#include "ElementBucket.h"
#include "Exception.h"
#include "FieldEvaluator.h"
#include "FunctionSpaceBucket.h"
#include "Hilbert.h"
#include "IndiceBucket.h"
#include "MatrixFactory.h"
#include "Memory.h"
#include "Message.h"
#include "OmpInterface.h"
#include "Options.h"
#include "Solver.h"
#include "Term.h"
#include "VectorFactory.h"
#include "instantiate.h"
#include "scalar.h"

#include <algorithm>
#include <cstdint>
#include <gmsh.h>
#include <iostream>
#include <map>
#include <unordered_map>

namespace gmshfem::problem
{


  template< class T_Scalar >
  Formulation< T_Scalar >::Formulation(const std::string &name, const std::string &matrixOptions) :
    _name(name), _terms(), _dofs(), _A(nullptr), _multiB(1), _solver(nullptr), _unknownFields(), _entities{}, _elementTypes(), _frequency(0.), _attributes()
  {
    initSystem(matrixOptions);
  }

  template< class T_Scalar >
  Formulation< T_Scalar >::~Formulation()
  {
    removeTerms();
    removeSystem();
    for(auto it = _attributes.begin(); it != _attributes.end(); ++it) {
      std::free(it->second);
    }
  }

  template< class T_Scalar >
  Formulation< T_Scalar >::Formulation(Formulation< T_Scalar > &&other) :
    _name(other._name), _terms(std::move(other._terms)), _dofs(std::move(other._dofs)), _A(other._A), _multiB(std::move(other._multiB)), _solver(other._solver), _unknownFields(std::move(other._unknownFields)), _entities{}, _elementTypes(), _frequency(other._frequency), _attributes(other._attributes)
  {
    other._terms.clear();
    other._dofs.clear();
    other._A = nullptr;
    other._multiB.clear();
    other._solver = nullptr;
    other._unknownFields.clear();
  }

  template< class T_Scalar >
  void Formulation< T_Scalar >::initSystem(const std::string &matrixOptions)
  {
    try {
      _A = new system::MatrixFactory< T_Scalar >(matrixOptions);
    }
    catch(const std::exception &exc) {
      msg::error << "Unable to allocate the matrix" << msg::endl;
      msg::error << "Origin of error: " << exc.what() << msg::endl;
      _A = nullptr;
      removeSystem();
      throw;
    }

    try {
      _multiB.resize(numRHS());
    }
    catch(const std::exception &exc) {
      msg::error << "Unable to allocate the right hand side vectors" << msg::endl;
      msg::error << "Origin of error: " << exc.what() << msg::endl;
      removeSystem();
      throw;
    }


    try {
      _solver = new system::Solver< T_Scalar >(_A, _multiB);
    }
    catch(const std::exception &exc) {
      msg::error << "Unable to allocate the solver" << msg::endl;
      msg::error << "Origin of error: " << exc.what() << msg::endl;
      _solver = nullptr;
      removeSystem();
      throw;
    }
  }

  template< class T_Scalar >
  void Formulation< T_Scalar >::removeTerms()
  {
    for(auto i = 0ULL; i < _terms.size(); ++i) {
      delete _terms[i];
    }
    _terms.clear();
  }

  template< class T_Scalar >
  void Formulation< T_Scalar >::removeSystem()
  {
    if(_A) {
      delete _A;
      _A = nullptr;
    }
    if(_solver) {
      delete _solver;
      _solver = nullptr;
    }

    _multiB.clear();
  }

  template< class T_Scalar >
  void Formulation< T_Scalar >::setSystemToZero()
  {
    _A->setToZero();
    for (auto& bi: _multiB)
      bi.setToZero();
  }

  template< class T_Scalar >
  void Formulation< T_Scalar >::setRHSToZero()
  {
    for (auto& bi: _multiB)
      bi.setToZero();
  }

  template< class T_Scalar >
  std::string Formulation< T_Scalar >::name() const
  {
    return _name;
  }

  template< class T_Scalar >
  void Formulation< T_Scalar >::_addField(field::FieldInterface< T_Scalar > *field)
  {
    for(auto i = 0ULL; i < _unknownFields.size(); ++i) {
      if(field::fieldIsStillValid(_unknownFields[i].first)) {
        if(field->tag() == _unknownFields[i].second->tag()) {
          return;
        }
      }
    }
    _unknownFields.push_back(std::make_pair(field->tag(), field));
  }

  template< class T_Scalar >
  typename std::vector< term::Term< T_Scalar > * >::iterator Formulation< T_Scalar >::begin()
  {
    return _terms.begin();
  }

  template< class T_Scalar >
  typename std::vector< term::Term< T_Scalar > * >::iterator Formulation< T_Scalar >::end()
  {
    return _terms.end();
  }

  template< class T_Scalar >
  template< class T_Object >
  void Formulation< T_Scalar >::setAttribute(const std::string &name, const T_Object &attribute)
  {
    auto it = _attributes.find(name);
    if(it == _attributes.end()) {
      T_Object *ptr = (T_Object *)std::malloc(sizeof(T_Object));
      *ptr = attribute;
      _attributes.insert(std::make_pair(name, ptr));
      return;
    }
    *static_cast< T_Object * >(it->second) = attribute;
  }

  INSTANTIATE_CLASS_FCT(void, , Formulation, setAttribute, 0, 4, class, TEMPLATE_ARGS(std::complex< double >, std::complex< float >, double, float), 9, class, TEMPLATE_ARGS(bool, int, unsigned int, long long, unsigned long long, std::complex< double >, scalar::Precision< std::complex< double > >, char, std::string), TEMPLATE_PARAMS(const std::string &, const TEMPLATE_PARAM_1 &))


  template< class T_Scalar >
  template< class T_Object >
  void Formulation< T_Scalar >::getAttribute(const std::string &name, T_Object &attribute) const
  {
    auto it = _attributes.find(name);
    if(it == _attributes.end()) {
      msg::error << "There is no attribute '" + name + "' defined in formulation '" + _name + "'" << msg::endl;
      return;
    }
    attribute = *static_cast< T_Object * >(it->second);
  }

  INSTANTIATE_CLASS_FCT_CONST(void, , Formulation, getAttribute, 0, 4, class, TEMPLATE_ARGS(std::complex< double >, std::complex< float >, double, float), 9, class, TEMPLATE_ARGS(bool, int, unsigned int, long long, unsigned long long, std::complex< double >, scalar::Precision< std::complex< double > >, char, std::string), TEMPLATE_PARAMS(const std::string &, TEMPLATE_PARAM_1 &))


  template< class T_Scalar >
  unsigned Formulation< T_Scalar >::numRHS() const
  {
    return _numRHS;
  }

  template< class T_Scalar >
  void Formulation< T_Scalar >::setCurrentRHS(unsigned idx)
  {
    _currentRHS = idx;
    _numRHS = std::max(idx + 1, _numRHS);
  }

  // Bilinear terms
  template< class T_Scalar >
  unsigned int Formulation< T_Scalar >::integral(const equation::EquationInterface< T_Scalar, Degree::Degree0, field::Form::Form0 > &equationLhs, const equation::EquationInterface< T_Scalar, Degree::Degree0, field::Form::Form0 > &equationRhs, const domain::GeometricObject &domain, const std::string &integrationType, const term::ProductType productType)
  {
    term::BilinearTermInterface< T_Scalar > *newTerm = new term::BilinearTerm< T_Scalar, field::Form::Form0, field::Form::Form0 >(domain, integrationType, equationLhs, equationRhs, productType);
    _terms.push_back(newTerm);
    _addField(equationLhs.getField());
    return newTerm->tag();
  }

  template< class T_Scalar >
  unsigned int Formulation< T_Scalar >::integral(const equation::EquationInterface< T_Scalar, Degree::Degree0, field::Form::Form0 > &equationLhs, const equation::EquationInterface< T_Scalar, Degree::Degree0, field::Form::Form1 > &equationRhs, const domain::GeometricObject &domain, const std::string &integrationType, const term::ProductType productType)
  {
    term::BilinearTermInterface< T_Scalar > *newTerm = new term::BilinearTerm< T_Scalar, field::Form::Form0, field::Form::Form1 >(domain, integrationType, equationLhs, equationRhs, productType);
    _terms.push_back(newTerm);
    _addField(equationLhs.getField());
    return newTerm->tag();
  }

  template< class T_Scalar >
  unsigned int Formulation< T_Scalar >::integral(const equation::EquationInterface< T_Scalar, Degree::Degree0, field::Form::Form0 > &equationLhs, const equation::EquationInterface< T_Scalar, Degree::Degree0, field::Form::Form2 > &equationRhs, const domain::GeometricObject &domain, const std::string &integrationType, const term::ProductType productType)
  {
    term::BilinearTermInterface< T_Scalar > *newTerm = new term::BilinearTerm< T_Scalar, field::Form::Form0, field::Form::Form2 >(domain, integrationType, equationLhs, equationRhs, productType);
    _terms.push_back(newTerm);
    _addField(equationLhs.getField());
    return newTerm->tag();
  }

  template< class T_Scalar >
  unsigned int Formulation< T_Scalar >::integral(const equation::EquationInterface< T_Scalar, Degree::Degree0, field::Form::Form0 > &equationLhs, const equation::EquationInterface< T_Scalar, Degree::Degree0, field::Form::Form3 > &equationRhs, const domain::GeometricObject &domain, const std::string &integrationType, const term::ProductType productType)
  {
    term::BilinearTermInterface< T_Scalar > *newTerm = new term::BilinearTerm< T_Scalar, field::Form::Form0, field::Form::Form3 >(domain, integrationType, equationLhs, equationRhs, productType);
    _terms.push_back(newTerm);
    _addField(equationLhs.getField());
    return newTerm->tag();
  }

  template< class T_Scalar >
  unsigned int Formulation< T_Scalar >::integral(const equation::EquationInterface< T_Scalar, Degree::Degree0, field::Form::Form1 > &equationLhs, const equation::EquationInterface< T_Scalar, Degree::Degree0, field::Form::Form0 > &equationRhs, const domain::GeometricObject &domain, const std::string &integrationType, const term::ProductType productType)
  {
    term::BilinearTermInterface< T_Scalar > *newTerm = new term::BilinearTerm< T_Scalar, field::Form::Form1, field::Form::Form0 >(domain, integrationType, equationLhs, equationRhs, productType);
    _terms.push_back(newTerm);
    _addField(equationLhs.getField());
    return newTerm->tag();
  }

  template< class T_Scalar >
  unsigned int Formulation< T_Scalar >::integral(const equation::EquationInterface< T_Scalar, Degree::Degree0, field::Form::Form1 > &equationLhs, const equation::EquationInterface< T_Scalar, Degree::Degree0, field::Form::Form1 > &equationRhs, const domain::GeometricObject &domain, const std::string &integrationType, const term::ProductType productType)
  {
    term::BilinearTermInterface< T_Scalar > *newTerm = new term::BilinearTerm< T_Scalar, field::Form::Form1, field::Form::Form1 >(domain, integrationType, equationLhs, equationRhs, productType);
    _terms.push_back(newTerm);
    _addField(equationLhs.getField());
    return newTerm->tag();
  }

  template< class T_Scalar >
  unsigned int Formulation< T_Scalar >::integral(const equation::EquationInterface< T_Scalar, Degree::Degree0, field::Form::Form1 > &equationLhs, const equation::EquationInterface< T_Scalar, Degree::Degree0, field::Form::Form2 > &equationRhs, const domain::GeometricObject &domain, const std::string &integrationType, const term::ProductType productType)
  {
    term::BilinearTermInterface< T_Scalar > *newTerm = new term::BilinearTerm< T_Scalar, field::Form::Form1, field::Form::Form2 >(domain, integrationType, equationLhs, equationRhs, productType);
    _terms.push_back(newTerm);
    _addField(equationLhs.getField());
    return newTerm->tag();
  }

  template< class T_Scalar >
  unsigned int Formulation< T_Scalar >::integral(const equation::EquationInterface< T_Scalar, Degree::Degree0, field::Form::Form1 > &equationLhs, const equation::EquationInterface< T_Scalar, Degree::Degree0, field::Form::Form3 > &equationRhs, const domain::GeometricObject &domain, const std::string &integrationType, const term::ProductType productType)
  {
    term::BilinearTermInterface< T_Scalar > *newTerm = new term::BilinearTerm< T_Scalar, field::Form::Form1, field::Form::Form3 >(domain, integrationType, equationLhs, equationRhs, productType);
    _terms.push_back(newTerm);
    _addField(equationLhs.getField());
    return newTerm->tag();
  }

  template< class T_Scalar >
  unsigned int Formulation< T_Scalar >::integral(const equation::EquationInterface< T_Scalar, Degree::Degree0, field::Form::Form2 > &equationLhs, const equation::EquationInterface< T_Scalar, Degree::Degree0, field::Form::Form0 > &equationRhs, const domain::GeometricObject &domain, const std::string &integrationType, const term::ProductType productType)
  {
    term::BilinearTermInterface< T_Scalar > *newTerm = new term::BilinearTerm< T_Scalar, field::Form::Form2, field::Form::Form0 >(domain, integrationType, equationLhs, equationRhs, productType);
    _terms.push_back(newTerm);
    _addField(equationLhs.getField());
    return newTerm->tag();
  }

  template< class T_Scalar >
  unsigned int Formulation< T_Scalar >::integral(const equation::EquationInterface< T_Scalar, Degree::Degree0, field::Form::Form2 > &equationLhs, const equation::EquationInterface< T_Scalar, Degree::Degree0, field::Form::Form1 > &equationRhs, const domain::GeometricObject &domain, const std::string &integrationType, const term::ProductType productType)
  {
    term::BilinearTermInterface< T_Scalar > *newTerm = new term::BilinearTerm< T_Scalar, field::Form::Form2, field::Form::Form1 >(domain, integrationType, equationLhs, equationRhs, productType);
    _terms.push_back(newTerm);
    _addField(equationLhs.getField());
    return newTerm->tag();
  }

  template< class T_Scalar >
  unsigned int Formulation< T_Scalar >::integral(const equation::EquationInterface< T_Scalar, Degree::Degree0, field::Form::Form2 > &equationLhs, const equation::EquationInterface< T_Scalar, Degree::Degree0, field::Form::Form2 > &equationRhs, const domain::GeometricObject &domain, const std::string &integrationType, const term::ProductType productType)
  {
    term::BilinearTermInterface< T_Scalar > *newTerm = new term::BilinearTerm< T_Scalar, field::Form::Form2, field::Form::Form2 >(domain, integrationType, equationLhs, equationRhs, productType);
    _terms.push_back(newTerm);
    _addField(equationLhs.getField());
    return newTerm->tag();
  }

  template< class T_Scalar >
  unsigned int Formulation< T_Scalar >::integral(const equation::EquationInterface< T_Scalar, Degree::Degree0, field::Form::Form2 > &equationLhs, const equation::EquationInterface< T_Scalar, Degree::Degree0, field::Form::Form3 > &equationRhs, const domain::GeometricObject &domain, const std::string &integrationType, const term::ProductType productType)
  {
    term::BilinearTermInterface< T_Scalar > *newTerm = new term::BilinearTerm< T_Scalar, field::Form::Form2, field::Form::Form3 >(domain, integrationType, equationLhs, equationRhs, productType);
    _terms.push_back(newTerm);
    _addField(equationLhs.getField());
    return newTerm->tag();
  }

  template< class T_Scalar >
  unsigned int Formulation< T_Scalar >::integral(const equation::EquationInterface< T_Scalar, Degree::Degree0, field::Form::Form3 > &equationLhs, const equation::EquationInterface< T_Scalar, Degree::Degree0, field::Form::Form0 > &equationRhs, const domain::GeometricObject &domain, const std::string &integrationType, const term::ProductType productType)
  {
    term::BilinearTermInterface< T_Scalar > *newTerm = new term::BilinearTerm< T_Scalar, field::Form::Form3, field::Form::Form0 >(domain, integrationType, equationLhs, equationRhs, productType);
    _terms.push_back(newTerm);
    _addField(equationLhs.getField());
    return newTerm->tag();
  }

  template< class T_Scalar >
  unsigned int Formulation< T_Scalar >::integral(const equation::EquationInterface< T_Scalar, Degree::Degree0, field::Form::Form3 > &equationLhs, const equation::EquationInterface< T_Scalar, Degree::Degree0, field::Form::Form1 > &equationRhs, const domain::GeometricObject &domain, const std::string &integrationType, const term::ProductType productType)
  {
    term::BilinearTermInterface< T_Scalar > *newTerm = new term::BilinearTerm< T_Scalar, field::Form::Form3, field::Form::Form1 >(domain, integrationType, equationLhs, equationRhs, productType);
    _terms.push_back(newTerm);
    _addField(equationLhs.getField());
    return newTerm->tag();
  }

  template< class T_Scalar >
  unsigned int Formulation< T_Scalar >::integral(const equation::EquationInterface< T_Scalar, Degree::Degree0, field::Form::Form3 > &equationLhs, const equation::EquationInterface< T_Scalar, Degree::Degree0, field::Form::Form2 > &equationRhs, const domain::GeometricObject &domain, const std::string &integrationType, const term::ProductType productType)
  {
    term::BilinearTermInterface< T_Scalar > *newTerm = new term::BilinearTerm< T_Scalar, field::Form::Form3, field::Form::Form2 >(domain, integrationType, equationLhs, equationRhs, productType);
    _terms.push_back(newTerm);
    _addField(equationLhs.getField());
    return newTerm->tag();
  }

  template< class T_Scalar >
  unsigned int Formulation< T_Scalar >::integral(const equation::EquationInterface< T_Scalar, Degree::Degree0, field::Form::Form3 > &equationLhs, const equation::EquationInterface< T_Scalar, Degree::Degree0, field::Form::Form3 > &equationRhs, const domain::GeometricObject &domain, const std::string &integrationType, const term::ProductType productType)
  {
    term::BilinearTermInterface< T_Scalar > *newTerm = new term::BilinearTerm< T_Scalar, field::Form::Form3, field::Form::Form3 >(domain, integrationType, equationLhs, equationRhs, productType);
    _terms.push_back(newTerm);
    _addField(equationLhs.getField());
    return newTerm->tag();
  }

  template< class T_Scalar >
  unsigned int Formulation< T_Scalar >::integral(const equation::EquationInterface< T_Scalar, Degree::Degree1, field::Form::Form0 > &equationLhs, const equation::EquationInterface< T_Scalar, Degree::Degree1, field::Form::Form0 > &equationRhs, const domain::GeometricObject &domain, const std::string &integrationType, const term::ProductType productType)
  {
    term::BilinearTermInterface< T_Scalar > *newTerm = new term::BilinearTerm< T_Scalar, field::Form::Form0, field::Form::Form0 >(domain, integrationType, equationLhs, equationRhs, productType);
    _terms.push_back(newTerm);
    _addField(equationLhs.getField());
    return newTerm->tag();
  }

  template< class T_Scalar >
  unsigned int Formulation< T_Scalar >::integral(const equation::EquationInterface< T_Scalar, Degree::Degree1, field::Form::Form0 > &equationLhs, const equation::EquationInterface< T_Scalar, Degree::Degree1, field::Form::Form1 > &equationRhs, const domain::GeometricObject &domain, const std::string &integrationType, const term::ProductType productType)
  {
    term::BilinearTermInterface< T_Scalar > *newTerm = new term::BilinearTerm< T_Scalar, field::Form::Form0, field::Form::Form1 >(domain, integrationType, equationLhs, equationRhs, productType);
    _terms.push_back(newTerm);
    _addField(equationLhs.getField());
    return newTerm->tag();
  }

  template< class T_Scalar >
  unsigned int Formulation< T_Scalar >::integral(const equation::EquationInterface< T_Scalar, Degree::Degree1, field::Form::Form0 > &equationLhs, const equation::EquationInterface< T_Scalar, Degree::Degree1, field::Form::Form2 > &equationRhs, const domain::GeometricObject &domain, const std::string &integrationType, const term::ProductType productType)
  {
    term::BilinearTermInterface< T_Scalar > *newTerm = new term::BilinearTerm< T_Scalar, field::Form::Form0, field::Form::Form2 >(domain, integrationType, equationLhs, equationRhs, productType);
    _terms.push_back(newTerm);
    _addField(equationLhs.getField());
    return newTerm->tag();
  }

  template< class T_Scalar >
  unsigned int Formulation< T_Scalar >::integral(const equation::EquationInterface< T_Scalar, Degree::Degree1, field::Form::Form0 > &equationLhs, const equation::EquationInterface< T_Scalar, Degree::Degree1, field::Form::Form3 > &equationRhs, const domain::GeometricObject &domain, const std::string &integrationType, const term::ProductType productType)
  {
    term::BilinearTermInterface< T_Scalar > *newTerm = new term::BilinearTerm< T_Scalar, field::Form::Form0, field::Form::Form3 >(domain, integrationType, equationLhs, equationRhs, productType);
    _terms.push_back(newTerm);
    _addField(equationLhs.getField());
    return newTerm->tag();
  }

  template< class T_Scalar >
  unsigned int Formulation< T_Scalar >::integral(const equation::EquationInterface< T_Scalar, Degree::Degree1, field::Form::Form1 > &equationLhs, const equation::EquationInterface< T_Scalar, Degree::Degree1, field::Form::Form0 > &equationRhs, const domain::GeometricObject &domain, const std::string &integrationType, const term::ProductType productType)
  {
    term::BilinearTermInterface< T_Scalar > *newTerm = new term::BilinearTerm< T_Scalar, field::Form::Form1, field::Form::Form0 >(domain, integrationType, equationLhs, equationRhs, productType);
    _terms.push_back(newTerm);
    _addField(equationLhs.getField());
    return newTerm->tag();
  }

  template< class T_Scalar >
  unsigned int Formulation< T_Scalar >::integral(const equation::EquationInterface< T_Scalar, Degree::Degree1, field::Form::Form1 > &equationLhs, const equation::EquationInterface< T_Scalar, Degree::Degree1, field::Form::Form1 > &equationRhs, const domain::GeometricObject &domain, const std::string &integrationType, const term::ProductType productType)
  {
    term::BilinearTermInterface< T_Scalar > *newTerm = new term::BilinearTerm< T_Scalar, field::Form::Form1, field::Form::Form1 >(domain, integrationType, equationLhs, equationRhs, productType);
    _terms.push_back(newTerm);
    _addField(equationLhs.getField());
    return newTerm->tag();
  }

  template< class T_Scalar >
  unsigned int Formulation< T_Scalar >::integral(const equation::EquationInterface< T_Scalar, Degree::Degree1, field::Form::Form1 > &equationLhs, const equation::EquationInterface< T_Scalar, Degree::Degree1, field::Form::Form2 > &equationRhs, const domain::GeometricObject &domain, const std::string &integrationType, const term::ProductType productType)
  {
    term::BilinearTermInterface< T_Scalar > *newTerm = new term::BilinearTerm< T_Scalar, field::Form::Form1, field::Form::Form2 >(domain, integrationType, equationLhs, equationRhs, productType);
    _terms.push_back(newTerm);
    _addField(equationLhs.getField());
    return newTerm->tag();
  }

  template< class T_Scalar >
  unsigned int Formulation< T_Scalar >::integral(const equation::EquationInterface< T_Scalar, Degree::Degree1, field::Form::Form1 > &equationLhs, const equation::EquationInterface< T_Scalar, Degree::Degree1, field::Form::Form3 > &equationRhs, const domain::GeometricObject &domain, const std::string &integrationType, const term::ProductType productType)
  {
    term::BilinearTermInterface< T_Scalar > *newTerm = new term::BilinearTerm< T_Scalar, field::Form::Form1, field::Form::Form3 >(domain, integrationType, equationLhs, equationRhs, productType);
    _terms.push_back(newTerm);
    _addField(equationLhs.getField());
    return newTerm->tag();
  }

  template< class T_Scalar >
  unsigned int Formulation< T_Scalar >::integral(const equation::EquationInterface< T_Scalar, Degree::Degree1, field::Form::Form2 > &equationLhs, const equation::EquationInterface< T_Scalar, Degree::Degree1, field::Form::Form0 > &equationRhs, const domain::GeometricObject &domain, const std::string &integrationType, const term::ProductType productType)
  {
    term::BilinearTermInterface< T_Scalar > *newTerm = new term::BilinearTerm< T_Scalar, field::Form::Form2, field::Form::Form0 >(domain, integrationType, equationLhs, equationRhs, productType);
    _terms.push_back(newTerm);
    _addField(equationLhs.getField());
    return newTerm->tag();
  }

  template< class T_Scalar >
  unsigned int Formulation< T_Scalar >::integral(const equation::EquationInterface< T_Scalar, Degree::Degree1, field::Form::Form2 > &equationLhs, const equation::EquationInterface< T_Scalar, Degree::Degree1, field::Form::Form1 > &equationRhs, const domain::GeometricObject &domain, const std::string &integrationType, const term::ProductType productType)
  {
    term::BilinearTermInterface< T_Scalar > *newTerm = new term::BilinearTerm< T_Scalar, field::Form::Form2, field::Form::Form1 >(domain, integrationType, equationLhs, equationRhs, productType);
    _terms.push_back(newTerm);
    _addField(equationLhs.getField());
    return newTerm->tag();
  }

  template< class T_Scalar >
  unsigned int Formulation< T_Scalar >::integral(const equation::EquationInterface< T_Scalar, Degree::Degree1, field::Form::Form2 > &equationLhs, const equation::EquationInterface< T_Scalar, Degree::Degree1, field::Form::Form2 > &equationRhs, const domain::GeometricObject &domain, const std::string &integrationType, const term::ProductType productType)
  {
    term::BilinearTermInterface< T_Scalar > *newTerm = new term::BilinearTerm< T_Scalar, field::Form::Form2, field::Form::Form2 >(domain, integrationType, equationLhs, equationRhs, productType);
    _terms.push_back(newTerm);
    _addField(equationLhs.getField());
    return newTerm->tag();
  }

  template< class T_Scalar >
  unsigned int Formulation< T_Scalar >::integral(const equation::EquationInterface< T_Scalar, Degree::Degree1, field::Form::Form2 > &equationLhs, const equation::EquationInterface< T_Scalar, Degree::Degree1, field::Form::Form3 > &equationRhs, const domain::GeometricObject &domain, const std::string &integrationType, const term::ProductType productType)
  {
    term::BilinearTermInterface< T_Scalar > *newTerm = new term::BilinearTerm< T_Scalar, field::Form::Form2, field::Form::Form3 >(domain, integrationType, equationLhs, equationRhs, productType);
    _terms.push_back(newTerm);
    _addField(equationLhs.getField());
    return newTerm->tag();
  }

  template< class T_Scalar >
  unsigned int Formulation< T_Scalar >::integral(const equation::EquationInterface< T_Scalar, Degree::Degree1, field::Form::Form3 > &equationLhs, const equation::EquationInterface< T_Scalar, Degree::Degree1, field::Form::Form0 > &equationRhs, const domain::GeometricObject &domain, const std::string &integrationType, const term::ProductType productType)
  {
    term::BilinearTermInterface< T_Scalar > *newTerm = new term::BilinearTerm< T_Scalar, field::Form::Form3, field::Form::Form0 >(domain, integrationType, equationLhs, equationRhs, productType);
    _terms.push_back(newTerm);
    _addField(equationLhs.getField());
    return newTerm->tag();
  }

  template< class T_Scalar >
  unsigned int Formulation< T_Scalar >::integral(const equation::EquationInterface< T_Scalar, Degree::Degree1, field::Form::Form3 > &equationLhs, const equation::EquationInterface< T_Scalar, Degree::Degree1, field::Form::Form1 > &equationRhs, const domain::GeometricObject &domain, const std::string &integrationType, const term::ProductType productType)
  {
    term::BilinearTermInterface< T_Scalar > *newTerm = new term::BilinearTerm< T_Scalar, field::Form::Form3, field::Form::Form1 >(domain, integrationType, equationLhs, equationRhs, productType);
    _terms.push_back(newTerm);
    _addField(equationLhs.getField());
    return newTerm->tag();
  }

  template< class T_Scalar >
  unsigned int Formulation< T_Scalar >::integral(const equation::EquationInterface< T_Scalar, Degree::Degree1, field::Form::Form3 > &equationLhs, const equation::EquationInterface< T_Scalar, Degree::Degree1, field::Form::Form2 > &equationRhs, const domain::GeometricObject &domain, const std::string &integrationType, const term::ProductType productType)
  {
    term::BilinearTermInterface< T_Scalar > *newTerm = new term::BilinearTerm< T_Scalar, field::Form::Form3, field::Form::Form2 >(domain, integrationType, equationLhs, equationRhs, productType);
    _terms.push_back(newTerm);
    _addField(equationLhs.getField());
    return newTerm->tag();
  }

  template< class T_Scalar >
  unsigned int Formulation< T_Scalar >::integral(const equation::EquationInterface< T_Scalar, Degree::Degree1, field::Form::Form3 > &equationLhs, const equation::EquationInterface< T_Scalar, Degree::Degree1, field::Form::Form3 > &equationRhs, const domain::GeometricObject &domain, const std::string &integrationType, const term::ProductType productType)
  {
    term::BilinearTermInterface< T_Scalar > *newTerm = new term::BilinearTerm< T_Scalar, field::Form::Form3, field::Form::Form3 >(domain, integrationType, equationLhs, equationRhs, productType);
    _terms.push_back(newTerm);
    _addField(equationLhs.getField());
    return newTerm->tag();
  }

  template< class T_Scalar >
  unsigned int Formulation< T_Scalar >::integral(const equation::EquationInterface< T_Scalar, Degree::Degree2, field::Form::Form0 > &equationLhs, const equation::EquationInterface< T_Scalar, Degree::Degree2, field::Form::Form0 > &equationRhs, const domain::GeometricObject &domain, const std::string &integrationType, const term::ProductType productType)
  {
    term::BilinearTermInterface< T_Scalar > *newTerm = new term::BilinearTerm< T_Scalar, field::Form::Form0, field::Form::Form0 >(domain, integrationType, equationLhs, equationRhs, productType);
    _terms.push_back(newTerm);
    _addField(equationLhs.getField());
    return newTerm->tag();
  }

  template< class T_Scalar >
  unsigned int Formulation< T_Scalar >::integral(const equation::EquationInterface< T_Scalar, Degree::Degree2, field::Form::Form0 > &equationLhs, const equation::EquationInterface< T_Scalar, Degree::Degree2, field::Form::Form1 > &equationRhs, const domain::GeometricObject &domain, const std::string &integrationType, const term::ProductType productType)
  {
    term::BilinearTermInterface< T_Scalar > *newTerm = new term::BilinearTerm< T_Scalar, field::Form::Form0, field::Form::Form1 >(domain, integrationType, equationLhs, equationRhs, productType);
    _terms.push_back(newTerm);
    _addField(equationLhs.getField());
    return newTerm->tag();
  }

  template< class T_Scalar >
  unsigned int Formulation< T_Scalar >::integral(const equation::EquationInterface< T_Scalar, Degree::Degree2, field::Form::Form0 > &equationLhs, const equation::EquationInterface< T_Scalar, Degree::Degree2, field::Form::Form2 > &equationRhs, const domain::GeometricObject &domain, const std::string &integrationType, const term::ProductType productType)
  {
    term::BilinearTermInterface< T_Scalar > *newTerm = new term::BilinearTerm< T_Scalar, field::Form::Form0, field::Form::Form2 >(domain, integrationType, equationLhs, equationRhs, productType);
    _terms.push_back(newTerm);
    _addField(equationLhs.getField());
    return newTerm->tag();
  }

  template< class T_Scalar >
  unsigned int Formulation< T_Scalar >::integral(const equation::EquationInterface< T_Scalar, Degree::Degree2, field::Form::Form0 > &equationLhs, const equation::EquationInterface< T_Scalar, Degree::Degree2, field::Form::Form3 > &equationRhs, const domain::GeometricObject &domain, const std::string &integrationType, const term::ProductType productType)
  {
    term::BilinearTermInterface< T_Scalar > *newTerm = new term::BilinearTerm< T_Scalar, field::Form::Form0, field::Form::Form3 >(domain, integrationType, equationLhs, equationRhs, productType);
    _terms.push_back(newTerm);
    _addField(equationLhs.getField());
    return newTerm->tag();
  }

  template< class T_Scalar >
  unsigned int Formulation< T_Scalar >::integral(const equation::EquationInterface< T_Scalar, Degree::Degree2, field::Form::Form1 > &equationLhs, const equation::EquationInterface< T_Scalar, Degree::Degree2, field::Form::Form0 > &equationRhs, const domain::GeometricObject &domain, const std::string &integrationType, const term::ProductType productType)
  {
    term::BilinearTermInterface< T_Scalar > *newTerm = new term::BilinearTerm< T_Scalar, field::Form::Form1, field::Form::Form0 >(domain, integrationType, equationLhs, equationRhs, productType);
    _terms.push_back(newTerm);
    _addField(equationLhs.getField());
    return newTerm->tag();
  }

  template< class T_Scalar >
  unsigned int Formulation< T_Scalar >::integral(const equation::EquationInterface< T_Scalar, Degree::Degree2, field::Form::Form1 > &equationLhs, const equation::EquationInterface< T_Scalar, Degree::Degree2, field::Form::Form1 > &equationRhs, const domain::GeometricObject &domain, const std::string &integrationType, const term::ProductType productType)
  {
    term::BilinearTermInterface< T_Scalar > *newTerm = new term::BilinearTerm< T_Scalar, field::Form::Form1, field::Form::Form1 >(domain, integrationType, equationLhs, equationRhs, productType);
    _terms.push_back(newTerm);
    _addField(equationLhs.getField());
    return newTerm->tag();
  }

  template< class T_Scalar >
  unsigned int Formulation< T_Scalar >::integral(const equation::EquationInterface< T_Scalar, Degree::Degree2, field::Form::Form1 > &equationLhs, const equation::EquationInterface< T_Scalar, Degree::Degree2, field::Form::Form2 > &equationRhs, const domain::GeometricObject &domain, const std::string &integrationType, const term::ProductType productType)
  {
    term::BilinearTermInterface< T_Scalar > *newTerm = new term::BilinearTerm< T_Scalar, field::Form::Form1, field::Form::Form2 >(domain, integrationType, equationLhs, equationRhs, productType);
    _terms.push_back(newTerm);
    _addField(equationLhs.getField());
    return newTerm->tag();
  }

  template< class T_Scalar >
  unsigned int Formulation< T_Scalar >::integral(const equation::EquationInterface< T_Scalar, Degree::Degree2, field::Form::Form1 > &equationLhs, const equation::EquationInterface< T_Scalar, Degree::Degree2, field::Form::Form3 > &equationRhs, const domain::GeometricObject &domain, const std::string &integrationType, const term::ProductType productType)
  {
    term::BilinearTermInterface< T_Scalar > *newTerm = new term::BilinearTerm< T_Scalar, field::Form::Form1, field::Form::Form3 >(domain, integrationType, equationLhs, equationRhs, productType);
    _terms.push_back(newTerm);
    _addField(equationLhs.getField());
    return newTerm->tag();
  }

  template< class T_Scalar >
  unsigned int Formulation< T_Scalar >::integral(const equation::EquationInterface< T_Scalar, Degree::Degree2, field::Form::Form2 > &equationLhs, const equation::EquationInterface< T_Scalar, Degree::Degree2, field::Form::Form0 > &equationRhs, const domain::GeometricObject &domain, const std::string &integrationType, const term::ProductType productType)
  {
    term::BilinearTermInterface< T_Scalar > *newTerm = new term::BilinearTerm< T_Scalar, field::Form::Form2, field::Form::Form0 >(domain, integrationType, equationLhs, equationRhs, productType);
    _terms.push_back(newTerm);
    _addField(equationLhs.getField());
    return newTerm->tag();
  }

  template< class T_Scalar >
  unsigned int Formulation< T_Scalar >::integral(const equation::EquationInterface< T_Scalar, Degree::Degree2, field::Form::Form2 > &equationLhs, const equation::EquationInterface< T_Scalar, Degree::Degree2, field::Form::Form1 > &equationRhs, const domain::GeometricObject &domain, const std::string &integrationType, const term::ProductType productType)
  {
    term::BilinearTermInterface< T_Scalar > *newTerm = new term::BilinearTerm< T_Scalar, field::Form::Form2, field::Form::Form1 >(domain, integrationType, equationLhs, equationRhs, productType);
    _terms.push_back(newTerm);
    _addField(equationLhs.getField());
    return newTerm->tag();
  }

  template< class T_Scalar >
  unsigned int Formulation< T_Scalar >::integral(const equation::EquationInterface< T_Scalar, Degree::Degree2, field::Form::Form2 > &equationLhs, const equation::EquationInterface< T_Scalar, Degree::Degree2, field::Form::Form2 > &equationRhs, const domain::GeometricObject &domain, const std::string &integrationType, const term::ProductType productType)
  {
    term::BilinearTermInterface< T_Scalar > *newTerm = new term::BilinearTerm< T_Scalar, field::Form::Form2, field::Form::Form2 >(domain, integrationType, equationLhs, equationRhs, productType);
    _terms.push_back(newTerm);
    _addField(equationLhs.getField());
    return newTerm->tag();
  }

  template< class T_Scalar >
  unsigned int Formulation< T_Scalar >::integral(const equation::EquationInterface< T_Scalar, Degree::Degree2, field::Form::Form2 > &equationLhs, const equation::EquationInterface< T_Scalar, Degree::Degree2, field::Form::Form3 > &equationRhs, const domain::GeometricObject &domain, const std::string &integrationType, const term::ProductType productType)
  {
    term::BilinearTermInterface< T_Scalar > *newTerm = new term::BilinearTerm< T_Scalar, field::Form::Form2, field::Form::Form3 >(domain, integrationType, equationLhs, equationRhs, productType);
    _terms.push_back(newTerm);
    _addField(equationLhs.getField());
    return newTerm->tag();
  }

  template< class T_Scalar >
  unsigned int Formulation< T_Scalar >::integral(const equation::EquationInterface< T_Scalar, Degree::Degree2, field::Form::Form3 > &equationLhs, const equation::EquationInterface< T_Scalar, Degree::Degree2, field::Form::Form0 > &equationRhs, const domain::GeometricObject &domain, const std::string &integrationType, const term::ProductType productType)
  {
    term::BilinearTermInterface< T_Scalar > *newTerm = new term::BilinearTerm< T_Scalar, field::Form::Form3, field::Form::Form0 >(domain, integrationType, equationLhs, equationRhs, productType);
    _terms.push_back(newTerm);
    _addField(equationLhs.getField());
    return newTerm->tag();
  }

  template< class T_Scalar >
  unsigned int Formulation< T_Scalar >::integral(const equation::EquationInterface< T_Scalar, Degree::Degree2, field::Form::Form3 > &equationLhs, const equation::EquationInterface< T_Scalar, Degree::Degree2, field::Form::Form1 > &equationRhs, const domain::GeometricObject &domain, const std::string &integrationType, const term::ProductType productType)
  {
    term::BilinearTermInterface< T_Scalar > *newTerm = new term::BilinearTerm< T_Scalar, field::Form::Form3, field::Form::Form1 >(domain, integrationType, equationLhs, equationRhs, productType);
    _terms.push_back(newTerm);
    _addField(equationLhs.getField());
    return newTerm->tag();
  }

  template< class T_Scalar >
  unsigned int Formulation< T_Scalar >::integral(const equation::EquationInterface< T_Scalar, Degree::Degree2, field::Form::Form3 > &equationLhs, const equation::EquationInterface< T_Scalar, Degree::Degree2, field::Form::Form2 > &equationRhs, const domain::GeometricObject &domain, const std::string &integrationType, const term::ProductType productType)
  {
    term::BilinearTermInterface< T_Scalar > *newTerm = new term::BilinearTerm< T_Scalar, field::Form::Form3, field::Form::Form2 >(domain, integrationType, equationLhs, equationRhs, productType);
    _terms.push_back(newTerm);
    _addField(equationLhs.getField());
    return newTerm->tag();
  }

  template< class T_Scalar >
  unsigned int Formulation< T_Scalar >::integral(const equation::EquationInterface< T_Scalar, Degree::Degree2, field::Form::Form3 > &equationLhs, const equation::EquationInterface< T_Scalar, Degree::Degree2, field::Form::Form3 > &equationRhs, const domain::GeometricObject &domain, const std::string &integrationType, const term::ProductType productType)
  {
    term::BilinearTermInterface< T_Scalar > *newTerm = new term::BilinearTerm< T_Scalar, field::Form::Form3, field::Form::Form3 >(domain, integrationType, equationLhs, equationRhs, productType);
    _terms.push_back(newTerm);
    _addField(equationLhs.getField());
    return newTerm->tag();
  }

  // Linear terms
  template< class T_Scalar >
  unsigned int Formulation< T_Scalar >::integral(const function::Function< T_Scalar, Degree::Degree0 > &functionLhs, const equation::EquationInterface< T_Scalar, Degree::Degree0, field::Form::Form0 > &equationRhs, const domain::GeometricObject &domain, const std::string &integrationType, const term::ProductType productType)
  {
    term::LinearTermInterface< T_Scalar > *newTerm = new term::LinearTerm< T_Scalar, Degree::Degree0, field::Form::Form0 >(domain, integrationType, functionLhs, equationRhs, productType, _currentRHS);
    _terms.push_back(newTerm);
    return newTerm->tag();
  }

  template< class T_Scalar >
  unsigned int Formulation< T_Scalar >::integral(const function::Function< T_Scalar, Degree::Degree0 > &functionLhs, const equation::EquationInterface< T_Scalar, Degree::Degree0, field::Form::Form1 > &equationRhs, const domain::GeometricObject &domain, const std::string &integrationType, const term::ProductType productType)
  {
    term::LinearTermInterface< T_Scalar > *newTerm = new term::LinearTerm< T_Scalar, Degree::Degree0, field::Form::Form1 >(domain, integrationType, functionLhs, equationRhs, productType, _currentRHS);
    _terms.push_back(newTerm);
    return newTerm->tag();
  }

  template< class T_Scalar >
  unsigned int Formulation< T_Scalar >::integral(const function::Function< T_Scalar, Degree::Degree0 > &functionLhs, const equation::EquationInterface< T_Scalar, Degree::Degree0, field::Form::Form2 > &equationRhs, const domain::GeometricObject &domain, const std::string &integrationType, const term::ProductType productType)
  {
    term::LinearTermInterface< T_Scalar > *newTerm = new term::LinearTerm< T_Scalar, Degree::Degree0, field::Form::Form2 >(domain, integrationType, functionLhs, equationRhs, productType, _currentRHS);
    _terms.push_back(newTerm);
    return newTerm->tag();
  }

  template< class T_Scalar >
  unsigned int Formulation< T_Scalar >::integral(const function::Function< T_Scalar, Degree::Degree0 > &functionLhs, const equation::EquationInterface< T_Scalar, Degree::Degree0, field::Form::Form3 > &equationRhs, const domain::GeometricObject &domain, const std::string &integrationType, const term::ProductType productType)
  {
    term::LinearTermInterface< T_Scalar > *newTerm = new term::LinearTerm< T_Scalar, Degree::Degree0, field::Form::Form3 >(domain, integrationType, functionLhs, equationRhs, productType, _currentRHS);
    _terms.push_back(newTerm);
    return newTerm->tag();
  }

  template< class T_Scalar >
  unsigned int Formulation< T_Scalar >::integral(const function::Function< T_Scalar, Degree::Degree1 > &functionLhs, const equation::EquationInterface< T_Scalar, Degree::Degree1, field::Form::Form0 > &equationRhs, const domain::GeometricObject &domain, const std::string &integrationType, const term::ProductType productType)
  {
    term::LinearTermInterface< T_Scalar > *newTerm = new term::LinearTerm< T_Scalar, Degree::Degree1, field::Form::Form0 >(domain, integrationType, functionLhs, equationRhs, productType, _currentRHS);
    _terms.push_back(newTerm);
    return newTerm->tag();
  }

  template< class T_Scalar >
  unsigned int Formulation< T_Scalar >::integral(const function::Function< T_Scalar, Degree::Degree1 > &functionLhs, const equation::EquationInterface< T_Scalar, Degree::Degree1, field::Form::Form1 > &equationRhs, const domain::GeometricObject &domain, const std::string &integrationType, const term::ProductType productType)
  {
    term::LinearTermInterface< T_Scalar > *newTerm = new term::LinearTerm< T_Scalar, Degree::Degree1, field::Form::Form1 >(domain, integrationType, functionLhs, equationRhs, productType, _currentRHS);
    _terms.push_back(newTerm);
    return newTerm->tag();
  }

  template< class T_Scalar >
  unsigned int Formulation< T_Scalar >::integral(const function::Function< T_Scalar, Degree::Degree1 > &functionLhs, const equation::EquationInterface< T_Scalar, Degree::Degree1, field::Form::Form2 > &equationRhs, const domain::GeometricObject &domain, const std::string &integrationType, const term::ProductType productType)
  {
    term::LinearTermInterface< T_Scalar > *newTerm = new term::LinearTerm< T_Scalar, Degree::Degree1, field::Form::Form2 >(domain, integrationType, functionLhs, equationRhs, productType, _currentRHS);
    _terms.push_back(newTerm);
    return newTerm->tag();
  }

  template< class T_Scalar >
  unsigned int Formulation< T_Scalar >::integral(const function::Function< T_Scalar, Degree::Degree1 > &functionLhs, const equation::EquationInterface< T_Scalar, Degree::Degree1, field::Form::Form3 > &equationRhs, const domain::GeometricObject &domain, const std::string &integrationType, const term::ProductType productType)
  {
    term::LinearTermInterface< T_Scalar > *newTerm = new term::LinearTerm< T_Scalar, Degree::Degree1, field::Form::Form3 >(domain, integrationType, functionLhs, equationRhs, productType, _currentRHS);
    _terms.push_back(newTerm);
    return newTerm->tag();
  }

  template< class T_Scalar >
  unsigned int Formulation< T_Scalar >::integral(const function::Function< T_Scalar, Degree::Degree2 > &functionLhs, const equation::EquationInterface< T_Scalar, Degree::Degree2, field::Form::Form0 > &equationRhs, const domain::GeometricObject &domain, const std::string &integrationType, const term::ProductType productType)
  {
    term::LinearTermInterface< T_Scalar > *newTerm = new term::LinearTerm< T_Scalar, Degree::Degree2, field::Form::Form0 >(domain, integrationType, functionLhs, equationRhs, productType, _currentRHS);
    _terms.push_back(newTerm);
    return newTerm->tag();
  }

  template< class T_Scalar >
  unsigned int Formulation< T_Scalar >::integral(const function::Function< T_Scalar, Degree::Degree2 > &functionLhs, const equation::EquationInterface< T_Scalar, Degree::Degree2, field::Form::Form1 > &equationRhs, const domain::GeometricObject &domain, const std::string &integrationType, const term::ProductType productType)
  {
    term::LinearTermInterface< T_Scalar > *newTerm = new term::LinearTerm< T_Scalar, Degree::Degree2, field::Form::Form1 >(domain, integrationType, functionLhs, equationRhs, productType, _currentRHS);
    _terms.push_back(newTerm);
    return newTerm->tag();
  }

  template< class T_Scalar >
  unsigned int Formulation< T_Scalar >::integral(const function::Function< T_Scalar, Degree::Degree2 > &functionLhs, const equation::EquationInterface< T_Scalar, Degree::Degree2, field::Form::Form2 > &equationRhs, const domain::GeometricObject &domain, const std::string &integrationType, const term::ProductType productType)
  {
    term::LinearTermInterface< T_Scalar > *newTerm = new term::LinearTerm< T_Scalar, Degree::Degree2, field::Form::Form2 >(domain, integrationType, functionLhs, equationRhs, productType, _currentRHS);
    _terms.push_back(newTerm);
    return newTerm->tag();
  }

  template< class T_Scalar >
  unsigned int Formulation< T_Scalar >::integral(const function::Function< T_Scalar, Degree::Degree2 > &functionLhs, const equation::EquationInterface< T_Scalar, Degree::Degree2, field::Form::Form3 > &equationRhs, const domain::GeometricObject &domain, const std::string &integrationType, const term::ProductType productType)
  {
    term::LinearTermInterface< T_Scalar > *newTerm = new term::LinearTerm< T_Scalar, Degree::Degree2, field::Form::Form3 >(domain, integrationType, functionLhs, equationRhs, productType, _currentRHS);
    _terms.push_back(newTerm);
    return newTerm->tag();
  }

  // Global quantity terms
  template< class T_Scalar >
  void Formulation< T_Scalar >::globalTerm(field::GlobalQuantity< T_Scalar > &globalQuantity, const field::FixedComponent &fixComponent, const T_Scalar &fixedValue)
  {
    field::FieldInterface< T_Scalar > *associatedDualField = globalQuantity.getAssociatedDualField();
    field::FieldInterface< T_Scalar > *associatedPrimalField = globalQuantity.getAssociatedPrimalField();
    if(associatedDualField == nullptr || associatedPrimalField == nullptr) {
      msg::warning << "Global quantity '" << globalQuantity.name() << "' is not assigned to a field" << msg::endl;
      msg::warning << "The corresponding global term is ignored." << msg::endl;
      return;
    }

    globalQuantity._setFixedComponent(fixComponent);
    if(fixComponent == field::FixedComponent::Dual) {
      globalQuantity._setDualValue(fixedValue);
    }
    else if(fixComponent == field::FixedComponent::Primal) {
      globalQuantity._setPrimalValue(fixedValue);
    }
    else if(fixComponent == field::FixedComponent::None) {
      // circuit equations
    }
  }

  template< class T_Scalar >
  std::vector< term::Term< T_Scalar > * > Formulation< T_Scalar >::_getTermOnEntity(const std::pair< int, int > &entity) const
  {
    std::vector< term::Term< T_Scalar > * > terms;
    for(auto i = 0ULL; i < _terms.size(); ++i) {
      if(_terms[i]->domain().have(entity)) terms.push_back(_terms[i]);
    }

    return terms;
  }

  template< class T_Scalar >
  unsigned long long Formulation< T_Scalar >::getTotalNumberOfDof() const
  {
    return _dofs.nbrDofs();
  }

  template< class T_Scalar >
  unsigned long long Formulation< T_Scalar >::getNumberOfUnknownDof() const
  {
    return _dofs.nbrUnknownDofs();
  }

  template< class T_Scalar >
  unsigned long long Formulation< T_Scalar >::getNumberOfFixedDof() const
  {
    return _dofs.nbrFixedDofs();
  }

  template< class T_Scalar >
  void Formulation< T_Scalar >::setAngularFrequency(const scalar::Precision< T_Scalar > &frequency)
  {
    _frequency = frequency;
  }

  template< class T_Scalar >
  scalar::Precision< T_Scalar > Formulation< T_Scalar >::getAngularFrequency() const
  {
    return _frequency;
  }

  template< class T_Scalar >
  field::FieldInterface< T_Scalar > *Formulation< T_Scalar >::getField(const std::string &name) const
  {
    for(auto i = 0ULL; i < _unknownFields.size(); ++i) {
      if(field::fieldIsStillValid(_unknownFields[i].first)) {
        if(_unknownFields[i].second->name() == name) {
          return _unknownFields[i].second;
        }
      }
    }
    return nullptr;
  }

  template< class T_Scalar >
  common::Memory Formulation< T_Scalar >::getEstimatedFactorizationMemoryUsage() const
  {
    try {
      return _solver->getEstimatedFactorizationMemoryUsage();
    }
    catch(const std::exception &exc) {
      msg::error << "Unable to compute the estimated factoriation memory usage" << msg::endl;
      return common::Memory();
    }
  }

  template< class T_Scalar >
  unsigned long long Formulation< T_Scalar >::getNumberOfNonZeros() const
  {
    return _A->numberOfNonZeros();
  }

  template< class T_Scalar >
  static bool s_checkFieldsValidity(const std::vector< std::pair< unsigned int, field::FieldInterface< T_Scalar > * > > &unknownFields, const std::string &name)
  {
    std::string currentModel;
    gmsh::model::getCurrent(currentModel);

    for(auto i = 0ULL; i < unknownFields.size(); ++i) {
      if(!field::fieldIsStillValid(unknownFields[i].first)) {
        msg::error << "Field '" << unknownFields[i].first << "' is not valid" << msg::endl;
        msg::error << "Perhaps you use a destroyed field in the formulation '" << name << "'" << msg::endl;

        return false;
      }

      if(unknownFields[i].second->model() != currentModel) {
        msg::error << "Field '" << unknownFields[i].first << "' is not defined in the same model as the current one ('" << currentModel << "')" << msg::endl;
        return false;
      }
    }
    return true;
  }

  template< class T_Scalar >
  static void s_setGlobalValuePattern(const std::vector< std::pair< unsigned int, field::FieldInterface< T_Scalar > * > > &unknownFields, system::VectorFactory< T_Scalar > *b, system::MatrixFactory< T_Scalar > *A)
  {
    for(auto i = 0ULL; i < unknownFields.size(); ++i) {
      for(auto it = unknownFields[i].second->firstGlobalQuantity(); it != unknownFields[i].second->lastGlobalQuantity(); ++it) {
        if(it->second->fixedComponent() == field::FixedComponent::Primal) {
          A->addPatternGlobalDof(it->second->getNumDualDof());
        }
      }
    }
  }

  template< class T_Scalar >
  static void s_setLinkedIndices(const unsigned long long nbrLinkedDofs, const std::vector< std::pair< unsigned int, field::FieldInterface< T_Scalar > * > > &unknownFields, system::MatrixFactory< T_Scalar > *A, std::vector<system::VectorFactory< T_Scalar >>& b)
  {
    std::vector< unsigned long long > indicesLC(nbrLinkedDofs);
    if(indicesLC.size() != 0) {
      for(auto i = 0ULL; i < unknownFields.size(); ++i) {
        unknownFields[i].second->getLinkedIndices(indicesLC);
      }
    }
    std::vector< unsigned long long > indicesLCcopy = indicesLC;
    A->setIndicesLC(indicesLC);
    for (auto& bi: b)
      {
        auto copy = indicesLCcopy;
        bi.setIndicesLC(copy);
      }
  }

  template< class T_Scalar >
  common::Timer Formulation< T_Scalar >::pre(const DofsSort::Algorithm algo)
  {
    common::Timer time;
    time.tick();

    msg::info << "Pre-processing " << _name << "..." << msg::endl;

    if(_A == nullptr || _solver == nullptr) {
      throw common::Exception("This system is not initialized: did you forgot to call 'Formulation::initSystem()'?");
    }

    if(!s_checkFieldsValidity(_unknownFields, _name)) {
      msg::info << "Pre-processing aborted" << msg::endl;
      time.tock();
      return time;
    }

    // clean _dofs and fields before pre-processing the formulation
    _dofs.clear();
    for(auto i = 0ULL; i < _unknownFields.size(); ++i) {
      _unknownFields[i].second->clearDofs();
      _unknownFields[i].second->invalidateOrderedDofCache();
    }



    common::Timer prepro;
    prepro.tick();

    try {
      _dofs.built(_terms, _unknownFields);
    }
    catch(const std::exception &exc) {
      msg::error << "Unable to create the dof dictionary" << msg::endl;
      msg::error << "Origin of error: " << exc.what() << msg::endl;

      _dofs.clear();
      for(auto i = 0ULL; i < _unknownFields.size(); ++i) {
        _unknownFields[i].second->clearDofs();
      }

      msg::info << "Pre-processing aborted" << msg::endl;
      time.tock();
      return time;
    }
    prepro.tock();

    double bubbleRatio = _dofs.nbrUnknownDofs() == 0 ? 0 : 100. * double(_dofs.nbrBubbleDofs()) / _dofs.nbrUnknownDofs();
    msg::info << _dofs.nbrDofs() << " dofs created in " << prepro << "s:" << msg::endl;
    msg::info << " - " << _dofs.nbrUnknownDofs() << " unknown dofs" << msg::endl;
    msg::info << "  * " << _dofs.nbrBubbleDofs() << " bubble unknown dofs " << msg::fill(40, '.') << " " << msg::precision(3) << bubbleRatio << "%" << msg::endl;
    if(_dofs.nbrUnknownGlobalDofs()) {
      msg::info << "  * " << _dofs.nbrUnknownGlobalDofs() << " unknown global dofs" << msg::endl;
    }
    msg::info << " - " << _dofs.nbrFixedDofs() << " fixed dofs" << msg::endl;
    if(_dofs.nbrFixedGlobalDofs()) {
      msg::info << "  * " << _dofs.nbrFixedGlobalDofs() << " fixed global dofs" << msg::endl;
    }
    if(_dofs.nbrLinkedDofs()) {
      double bubbleLinkedRatio = _dofs.nbrLinkedDofs() == 0 ? 0 : 100. * double(_dofs.nbrBubbleLinkedDofs()) / _dofs.nbrLinkedDofs();
      msg::info << " - " << _dofs.nbrLinkedDofs() << " linked dofs" << msg::endl;
      msg::info << "  * " << _dofs.nbrBubbleLinkedDofs() << " bubble linked dofs " << msg::fill(40, '.') << " " << msg::precision(3) << bubbleLinkedRatio << "%" << msg::endl;
    }

    if(common::Options::instance()->memory) {
      common::Memory memory;
      for(auto i = 0ULL; i < _unknownFields.size(); ++i) {
        memory += _unknownFields[i].second->memory();
      }
      msg::info << "Memory footprint of fields: " << memory << msg::endl;
    }

    for(auto i = 0ULL; i < _unknownFields.size(); ++i) {
      field::FieldInterface< T_Scalar > *primal = nullptr, *dual = nullptr;
      for(auto it = _unknownFields[i].second->firstGlobalQuantity(); it != _unknownFields[i].second->lastGlobalQuantity(); ++it) {
        if(it->second->fixedComponent() == field::FixedComponent::Primal) {
          primal = it->second->getAssociatedPrimalField();
          dual = it->second->getAssociatedDualField();
        }
      }
      if(primal != nullptr && dual != nullptr) {
        std::vector< term::Term< T_Scalar > * > newTerm;
        for(auto itTerm = _terms.begin(); itTerm != _terms.end(); ++itTerm) {
          term::Term< T_Scalar > *dualTerm = (*itTerm)->switchTestFunctionField(primal, dual);
          if(dualTerm) {
            newTerm.push_back(dualTerm);
          }
        }
        for(auto itTerm = newTerm.begin(); itTerm != newTerm.end(); ++itTerm) {
          _terms.push_back(*itTerm);
        }
      }
    }

    try {
      const int dim = gmsh::model::getDimension();
      bool useBubble = (algo == DofsSort::Algorithm::Default) && (dim > 1);
      _A->init(_dofs.nbrUnknownDofs(), useBubble ? _dofs.nbrBubbleDofs() : 0);
      // b initialization deferred to assembly


      if(algo == DofsSort::Algorithm::Default || algo == DofsSort::Algorithm::Hilbert) {
        _dofs.reorderWithHilbert(true);
      }
      else if(algo != DofsSort::Algorithm::RCM && algo != DofsSort::Algorithm::None) {
        msg::warning << "Unknown Dof sort algorithm, default algorithm is used instead" << msg::endl;
      }

      // Linked indices deferred to assembly
      s_setLinkedIndices(_dofs.nbrLinkedDofs(), _unknownFields, _A, _multiB);

      _buildPattern(common::Options::instance()->elementsSortAlgorithm);
      // Used to require "b", but actually doesn't use it
      s_setGlobalValuePattern<T_Scalar>(_unknownFields, nullptr, _A);
      _A->finalizePattern();

      if(algo == DofsSort::Algorithm::RCM) {
        _dofs.reorderWithRCM(_A);
      }
    }
    catch(const std::exception &exc) {
      msg::error << "Unable to initalize the system" << msg::endl;
      msg::info << "Pre-processing aborted" << msg::endl;
      time.tock();
      return time;
    }

    for(auto dim = 0; dim <= 3; dim++) {
      gmsh::model::getEntities(_entities[dim], dim);
      gmsh::model::mesh::getElementTypes(_elementTypes[dim], dim);
    }

    time.tock();
    msg::info << "Done pre-processing in " << time << "s" << msg::endl;
    if(common::Options::instance()->memory) {
      auto memory = _A->memory();
      for (auto& bi: _multiB)
        memory += bi.memory();
      msg::info << "Memory footprint of system: " << memory << msg::endl;
    }

    return time;
  }

  template< class T_Scalar >
  static void s_setDirichletConditions(const unsigned long long nbrFixedDofs, const std::vector< std::pair< unsigned int, field::FieldInterface< T_Scalar > * > > &unknownFields, std::vector<system::VectorFactory< T_Scalar >>& b)
  {
    std::vector< T_Scalar > valuesDC(nbrFixedDofs);
    if(valuesDC.size() != 0) {
      for(auto i = 0ULL; i < unknownFields.size(); ++i) {
        unknownFields[i].second->getFixedValues(valuesDC);
      }
    }
    for(auto &bi : b) {
      std::vector< T_Scalar > copy = valuesDC;
      bi.setValuesDC(copy);
    }
  }

  template< class T_Scalar >
  static void s_setLinkedConditions(const unsigned long long nbrLinkedDofs, const std::vector< std::pair< unsigned int, field::FieldInterface< T_Scalar > * > > &unknownFields, system::MatrixFactory< T_Scalar > *A, std::vector<system::VectorFactory< T_Scalar >>& b)
  {
    std::vector< T_Scalar > valuesLC(nbrLinkedDofs);
    if(valuesLC.size() != 0) {
      for(auto i = 0ULL; i < unknownFields.size(); ++i) {
        unknownFields[i].second->getLinkedValues(valuesLC);
      }
    }
    std::vector< T_Scalar > valuesLCcopy = valuesLC;
    A->setValuesLC(valuesLCcopy);
    for (auto& bi: b)
    {
      valuesLCcopy = valuesLC;
      bi.setValuesLC(valuesLCcopy);
    }
  }

  template< class T_Scalar >
  static void s_setGlobalValueConditions(const std::vector< std::pair< unsigned int, field::FieldInterface< T_Scalar > * > > &unknownFields, std::vector<system::VectorFactory< T_Scalar >> &b, system::MatrixFactory< T_Scalar > *A)
  {
    for(auto i = 0ULL; i < unknownFields.size(); ++i) {
      if(field::fieldIsStillValid(unknownFields[i].first)) {
        for(auto it = unknownFields[i].second->firstGlobalQuantity(); it != unknownFields[i].second->lastGlobalQuantity(); ++it) {
          if(field::globalValueIsStillValid(it->first)) {
            if(it->second->fixedComponent() == field::FixedComponent::Dual) {
              for (auto &bi: b)
                bi.addValue(it->second->getNumPrimalDof(), -it->second->getDualValue());
            }
            else if(it->second->fixedComponent() == field::FixedComponent::Primal) {
              A->addValue(it->second->getNumDualDof(), it->second->getNumDualDof(), T_Scalar(-1.));
              unknownFields[i].second->assignValueTo(it->second->getNumPrimalDof(), dofs::Type::Fixed, it->second->getPrimalValue());
            }
          }
        }
      }
    }
  }

  template< class T_Scalar >
  void s_setModule(const bool separate, const std::string &name, const scalar::Precision< T_Scalar > &frequency, const std::vector< term::Term< T_Scalar > * > &terms, system::MatrixFactory< T_Scalar > *A)
  {
    bool haveStiffness = false;
    bool haveDamping = false;
    bool haveMass = false;
    for(auto i = 0ULL; i < terms.size(); ++i) {
      if(terms[i]->isBilinear()) {
        if(static_cast< term::BilinearTermInterface< T_Scalar > * >(terms[i])->unknownFieldType() == equation::UnknownFieldType::NotDt) {
          haveStiffness = true;
        }
        else if(static_cast< term::BilinearTermInterface< T_Scalar > * >(terms[i])->unknownFieldType() == equation::UnknownFieldType::Dt) {
          haveDamping = true;
        }
        else if(static_cast< term::BilinearTermInterface< T_Scalar > * >(terms[i])->unknownFieldType() == equation::UnknownFieldType::DtDt) {
          haveMass = true;
        }
      }
    }

    system::MatrixModule< T_Scalar > *matrixModule = nullptr;
    if(separate) {
      if(haveMass && haveDamping && haveStiffness) { // Module M-C-K
        matrixModule = new system::MCKModule< T_Scalar >(A);
      }
      else if(haveMass && !haveDamping && haveStiffness) { // Module M-K
        matrixModule = new system::MKModule< T_Scalar >(A);
      }
      else if(!haveMass && haveDamping && haveStiffness) { // Module C-K
        matrixModule = new system::CKModule< T_Scalar >(A);
      }
      else {
        msg::warning << "There are no mass or damping terms in formulation'" << name << "'" << msg::endl;
        matrixModule = new system::AModule< T_Scalar >(A);
      }
    }
    else {
      if(!haveMass && !haveDamping) {
        matrixModule = new system::AModule< T_Scalar >(A);
      }
      else {
        if constexpr(scalar::IsComplex< T_Scalar >::value) {
          if(frequency == 0.) {
            msg::warning << "The frequency is set to zero" << msg::endl;
          }
          matrixModule = new system::AFrequencyModule< T_Scalar >(A);
          static_cast< system::AFrequencyModule< T_Scalar > * >(matrixModule)->setFrequency(frequency);
        }
        else {
          throw common::Exception("Only a complex formulation can be assembled into a 'Ax = b' system when mass or damping terms are defined");
        }
      }
    }
    A->setModule(matrixModule);
  }

  template< class T_Scalar >
  common::Timer Formulation< T_Scalar >::assemble(const bool separate, const ElementsSort::Algorithm algo)
  {
    common::Timer time;
    time.tick();

    msg::info << "Assembling " << _name << "..." << msg::endl;

    if(numRHS() != _multiB.size()) {
      _multiB.clear();
      _multiB.resize(numRHS());
    }
    for(auto &bi : _multiB) {
      bi.init(_dofs.nbrUnknownDofs());
    }
    s_setLinkedIndices(_dofs.nbrLinkedDofs(), _unknownFields, _A, _multiB);


    if(_A->getModule() == nullptr) {
      s_setModule(separate, _name, _frequency, _terms, _A);
    }

    if(!s_checkFieldsValidity(_unknownFields, _name)) {
      msg::info << "Assembling aborted" << msg::endl;
      time.tock();
      return time;
    }

    if(_dofs.nbrUnknownDofs() == 0) {
      msg::warning << "There are no unknown dofs to assemble. Did you forget to pre-process the formulation?" << msg::endl;
      time.tock();
      msg::info << "Done assembling in " << time << "s" << msg::endl;
      return time;
    }

    s_setGlobalValueConditions(_unknownFields, _multiB, _A);
    s_setDirichletConditions(_dofs.nbrFixedDofs(), _unknownFields, _multiB);
    s_setLinkedConditions(_dofs.nbrLinkedDofs(), _unknownFields, _A, _multiB);

    try {
      msg::info << " - 3D entities" << msg::endl;
      _assembleDim(3, algo);
      msg::info << " - 2D entities" << msg::endl;
      _assembleDim(2, algo);
      msg::info << " - 1D entities" << msg::endl;
      _assembleDim(1, algo);
      msg::info << " - 0D entities" << msg::endl;
      _assembleDim(0, algo);
    }
    catch(const std::exception &exce) {
      msg::error << "Unable to assemble the system" << msg::endl;
      msg::error << exce.what() << msg::endl;
      _A->setToZero();
      for (auto& bi: _multiB)
        bi.setToZero();
      msg::info << "Assembling aborted" << msg::endl;
      time.tock();
      return time;
    }

    time.tock();
    msg::info << "Done assembling in " << time << "s" << msg::endl;

    return time;
  }

  template< class T_Scalar >
  common::Timer Formulation< T_Scalar >::solveAll(const bool reusePreconditioner)
  {
    msg::info << "Solving " << _name << " with " << _numRHS << " right-hand sides..." << msg::endl;

    common::Timer time;
    time.tick();

    _solutions.resize(_numRHS);

    try {
      _solver->solveAll(_solutions, reusePreconditioner);
    }
    catch(const std::exception &exc) {
      msg::error << "Unable to solve the system" << msg::endl;
      msg::error << exc.what() << msg::endl;
      msg::info << "Solving aborted" << msg::endl;
      time.tock();
      return time;
    }

    if(_solutions.back().size() != _dofs.nbrUnknownDofs()) {
      msg::error << "An error occurred while solving the system" << msg::endl;
      msg::info << "Solving aborted" << msg::endl;
      time.tock();
      return time;
    }

    //setSolutionIntoFields(values);
    //rawSolution = std::move(values); // Move values to rawSolution; values not needed anymore: can be safely moved

    time.tock();
    msg::info << "Done solving in " << time << "s" << msg::endl;
    return time;
  }

  template< class T_Scalar >
  void Formulation< T_Scalar >::loadSolution(unsigned idx)
  {
    if (idx >= numRHS()) {
      throw common::Exception("Tried to load a solution whose index (" + std::to_string(idx) + ") is higher than the number of fields (" + std::to_string(numRHS()) + ").");
    }
    if (idx >= _solutions.size()) {
      throw common::Exception("Tried to load a solution whose index (" + std::to_string(idx) + ") is higher than the number of computed solutions (" + std::to_string(numRHS()) + ").");
    }
    setSolutionIntoFields(_solutions.at(idx));
  }

  static void s_hilbertSort(const int elementType, const unsigned int tag)
  {
    const int problemDim = gmsh::model::getDimension();
    double xmin, ymin, zmin, xmax, ymax, zmax;
    std::vector< double > barycenters;

    gmsh::model::mesh::preallocateBarycenters(elementType, barycenters, tag);

#pragma omp parallel num_threads(omp::getMaxThreads())
    {
      const unsigned int numThreads = omp::getNumThreads();
      const unsigned int myThreadID = omp::getThreadNum();

      gmsh::model::mesh::getBarycenters(elementType, tag, true, true, barycenters, myThreadID, numThreads);
    }

    const unsigned long long nbrElements = barycenters.size() / 3;

    xmin = barycenters[0];
    xmax = barycenters[0];
    ymin = barycenters[1];
    ymax = barycenters[1];
    zmin = barycenters[2];
    zmax = barycenters[2];

    for(auto i = 0ULL; i < nbrElements; ++i) {
      xmin = (barycenters[3 * i + 0] < xmin ? barycenters[3 * i + 0] : xmin);
      xmax = (barycenters[3 * i + 0] > xmax ? barycenters[3 * i + 0] : xmax);
      ymin = (barycenters[3 * i + 1] < ymin ? barycenters[3 * i + 1] : ymin);
      ymax = (barycenters[3 * i + 1] > ymax ? barycenters[3 * i + 1] : ymax);
      zmin = (barycenters[3 * i + 2] < zmin ? barycenters[3 * i + 2] : zmin);
      zmax = (barycenters[3 * i + 2] > zmax ? barycenters[3 * i + 2] : zmax);
    }

    if(problemDim == 2) {
      std::vector< reorder::SortedEntity< 2 > > sortedEntities;
      sortedEntities.resize(nbrElements);

#pragma omp parallel num_threads(omp::getMaxThreads())
      {
        reorder::SortedEntity< 2 > se;
#pragma omp for
        for(auto i = 0ULL; i < nbrElements; ++i) {
          se.f[0] = barycenters[3 * i + 0];
          se.f[1] = barycenters[3 * i + 1];
          se.ptrA = reinterpret_cast< void * >(i);
          sortedEntities[i] = se;
        }
      }

      float min[2] = {static_cast< float >(xmin - 0.000001 * std::abs(xmin)),
                      static_cast< float >(ymin - 0.000001 * std::abs(ymin))};
      float max[2] = {static_cast< float >(xmax + 0.000001 * std::abs(xmax)),
                      static_cast< float >(ymax + 0.000001 * std::abs(ymax))};
      reorder::Hilbert< 2, 4 > hilbert;
      hilbert.apply(sortedEntities, min, max);

      barycenters.clear();
      barycenters.shrink_to_fit();
      std::vector< std::size_t > degree(nbrElements);
#pragma omp parallel for num_threads(omp::getMaxThreads())
      for(auto i = 0ULL; i < nbrElements; ++i) {
        degree[i] = reinterpret_cast< std::size_t >(sortedEntities[i].ptrB);
      }

      gmsh::model::mesh::reorderElements(elementType, tag, degree);
    }
    else if(problemDim == 3) {
      std::vector< reorder::SortedEntity< 3 > > sortedEntities;
      sortedEntities.resize(nbrElements);

#pragma omp parallel num_threads(omp::getMaxThreads())
      {
        reorder::SortedEntity< 3 > se;
#pragma omp for
        for(auto i = 0ULL; i < nbrElements; ++i) {
          se.f[0] = barycenters[3 * i + 0];
          se.f[1] = barycenters[3 * i + 1];
          se.f[2] = barycenters[3 * i + 2];
          se.ptrA = reinterpret_cast< void * >(i);
          sortedEntities[i] = se;
        }
      }

      float min[3] = {static_cast< float >(xmin - 0.000001 * std::abs(xmin)),
                      static_cast< float >(ymin - 0.000001 * std::abs(ymin)),
                      static_cast< float >(zmin - 0.000001 * std::abs(zmin))};
      float max[3] = {static_cast< float >(xmax + 0.000001 * std::abs(xmax)),
                      static_cast< float >(ymax + 0.000001 * std::abs(ymax)),
                      static_cast< float >(zmax + 0.000001 * std::abs(zmax))};
      reorder::Hilbert< 3, 8 > hilbert;
      hilbert.apply(sortedEntities, min, max);

      barycenters.clear();
      barycenters.shrink_to_fit();
      std::vector< std::size_t > degree(nbrElements);
#pragma omp parallel for num_threads(omp::getMaxThreads())
      for(auto i = 0ULL; i < nbrElements; ++i) {
        degree[i] = reinterpret_cast< std::size_t >(sortedEntities[i].ptrB);
      }

      gmsh::model::mesh::reorderElements(elementType, tag, degree);
    }
  }

  template< class T_Scalar >
  common::Timer Formulation< T_Scalar >::_buildPattern(const ElementsSort::Algorithm algo)
  {
    common::Timer time;
    time.tick();

    gmsh::vectorpair entities;
    gmsh::model::getEntities(entities);

    for(auto i = 0ULL; i < entities.size(); ++i) {
      std::vector< term::Term< T_Scalar > * > terms = _getTermOnEntity(entities[i]);
      if(terms.size() == 0) continue;

      term::Pattern< T_Scalar > pattern;
      for(auto j = 0ULL; j < terms.size(); ++j) {
        if(terms[j]->isBilinear()) {
          pattern.addBilinearTerm(static_cast< term::BilinearTermInterface< T_Scalar > * >(terms[j]));
        }
        else if(terms[j]->isLinear()) {
          pattern.addLinearTerm(static_cast< term::LinearTermInterface< T_Scalar > * >(terms[j]));
        }
      }
      pattern.sort(entities[i]);

      std::vector< int > elementTypes;
      gmsh::model::mesh::getElementTypes(elementTypes, entities[i].first, entities[i].second);

      for(auto typeIndex = 0ULL; typeIndex < elementTypes.size(); ++typeIndex) {
        if(entities[i].first > 1) {
          switch(algo) {
          case ElementsSort::Algorithm::Hilbert:
            s_hilbertSort(elementTypes[typeIndex], entities[i].second);
            break;
          default:
            break;
          }
        }
        for (auto& bi: _multiB)
          pattern.build(&bi, _A, elementTypes[typeIndex], entities[i]);

      }
    }

    time.tock();
    return time;
  }

  template< class T_Scalar >
  common::Timer Formulation< T_Scalar >::_assembleDim(const int dimToAssemble, const ElementsSort::Algorithm algo)
  {
    common::Timer time;
    time.tick();

    ElementBucket< scalar::Precision< T_Scalar > > elements;
    IndiceBucket indices;
    FunctionSpaceBucket< scalar::Precision< T_Scalar > > functionSpaces;

    std::vector< int > elementTypesOfEntity;

    for(auto elemType : _elementTypes[dimToAssemble]) {
      for(auto entity : _entities[dimToAssemble]) {

        // Check this entity has some elements of this type
        gmsh::model::mesh::getElementTypes(elementTypesOfEntity, dimToAssemble, entity.second);
        if (std::count(elementTypesOfEntity.begin(), elementTypesOfEntity.end(), elemType) == 0) continue;


        std::vector< term::Term< T_Scalar > * > terms = _getTermOnEntity(entity);
        if(terms.size() == 0) continue;

        term::Assembler< T_Scalar > assembler;
        for(auto j = 0ULL; j < terms.size(); ++j) {
          if(terms[j]->isActivated()) {
            if(terms[j]->isBilinear()) {
              assembler.addBilinearTerm(static_cast< term::BilinearTermInterface< T_Scalar > * >(terms[j]));
            }
            else if(terms[j]->isLinear()) {
              assembler.addLinearTerm(static_cast< term::LinearTermInterface< T_Scalar > * >(terms[j]));
            }
          }
        }
        assembler.sort();

        if(entity.first > 1) {
          switch(algo) {
          case ElementsSort::Algorithm::Hilbert:
            s_hilbertSort(elemType, entity.second);
            for(auto i = 0ULL; i < _unknownFields.size(); ++i) {
              _unknownFields[i].second->invalidateOrderedDofCache();
            }
            break;
          default: break;
          }
        }

        assembler.assemblyInitialization(indices, elements, functionSpaces, elemType, entity);

        if(common::Options::instance()->memory) {
          std::string elementTypesName;
          int dim, order, numNodes, numPrimaryNodes;
          std::vector< double > nodeCoord;
          gmsh::model::mesh::getElementProperties(elemType, elementTypesName, dim, order, numNodes, nodeCoord, numPrimaryNodes);
          msg::info << "On '" << elementTypesName << "' of entity (" << entity.first << ", " << entity.second << "):" << msg::endl;
          msg::info << " - Memory footprint of geometric quantities: " << elements.memory() << msg::endl;
          msg::info << " - Memory footprint of dof indices: " << indices.memory() << msg::endl;
          msg::info << " - Memory footprint of basis functions: " << functionSpaces.memory() << msg::endl;
          common::Memory memoryTerm;
          for(auto j = 0ULL; j < terms.size(); ++j) {
            memoryTerm += terms[j]->memory();
          }
          msg::info << " - Memory footprint of terms: " << memoryTerm << msg::endl;
        }

        assembler.assemble(elements, indices, _multiB, _A, elemType, entity);

        elements.clear();
        indices.clear();
        functionSpaces.clear();
      }
    }

    time.tock();
    return time;
  }

  template< class T_Scalar >
  void Formulation< T_Scalar >::setSolutionIntoFields(const std::vector< T_Scalar > &solution)
  {
    for(auto i = 0ULL; i < _unknownFields.size(); ++i) {
      if(field::fieldIsStillValid(_unknownFields[i].first)) {
        _unknownFields[i].second->assignValues(solution);
      }
      for(auto it = _unknownFields[i].second->firstGlobalQuantity(); it != _unknownFields[i].second->lastGlobalQuantity(); ++it) {
        if(field::globalValueIsStillValid(it->first)) {
          if(it->second->isActivated()) {
            if(it->second->fixedComponent() == field::FixedComponent::Dual) {
              it->second->_setPrimalValue(solution[it->second->getNumPrimalDof()]);
            }
            else if(it->second->fixedComponent() == field::FixedComponent::Primal) {
              it->second->_setDualValue(solution[it->second->getNumDualDof()]);
            }
          }
        }
      }
    }
  }

  template< class T_Scalar >
  common::Timer Formulation< T_Scalar >::solve(const bool reusePreconditioner)
  {
    algebra::Vector< T_Scalar > dummy;
    return solve(reusePreconditioner, dummy);
  }

  template< class T_Scalar >
  common::Timer Formulation< T_Scalar >::solve(const bool reusePreconditioner,
                                               algebra::Vector< T_Scalar > &rawSolution)
  {
    msg::info << "Solving " << _name << "..." << msg::endl;

    common::Timer time;
    time.tick();

    std::vector< T_Scalar > values;

    try {
      _solver->solve(values, reusePreconditioner);
    }
    catch(const std::exception &exc) {
      msg::error << "Unable to solve the system" << msg::endl;
      msg::error << exc.what() << msg::endl;
      msg::info << "Solving aborted" << msg::endl;
      time.tock();
      return time;
    }

    if(values.size() != _dofs.nbrUnknownDofs()) {
      msg::error << "An error occurred while solving the system" << msg::endl;
      msg::info << "Solving aborted" << msg::endl;
      time.tock();
      return time;
    }

    setSolutionIntoFields(values);
    rawSolution = std::move(values); // Move values to rawSolution; values not needed anymore: can be safely moved

    time.tock();
    msg::info << "Done solving in " << time << "s" << msg::endl;
    return time;
  }

  template< class T_Scalar >
  common::Timer Formulation< T_Scalar >::eigensolve(algebra::Vector< scalar::ComplexPrecision< T_Scalar > > &eigenvalues, const bool computeEigenvectors, const unsigned long long numberOfEigenvalues, const scalar::ComplexPrecision< T_Scalar > target)
  {
    msg::info << "Solving " << _name << " for eigenvalues... " << msg::endl;

    common::Timer time;
    time.tick();

    try {
      std::vector< scalar::ComplexPrecision< T_Scalar > > stdEigenvalues;
      std::vector< std::vector< scalar::ComplexPrecision< T_Scalar > > > eigenvectors;
      _solver->eigensolve(stdEigenvalues, eigenvectors, computeEigenvectors, numberOfEigenvalues, target);
      eigenvalues = std::move(stdEigenvalues);

      msg::info << "Found " << eigenvalues.size() << (computeEigenvectors ? " eigenpairs." : " eigenvalues.") << msg::endl;

      if(computeEigenvectors) {
        for(auto i = 0ULL; i < _unknownFields.size(); ++i) {
          if(field::fieldIsStillValid(_unknownFields[i].first)) {
            field::FieldModule< T_Scalar > *module = _unknownFields[i].second->getModule("Eigenpair");
            if(module == nullptr) {
              module = new field::EigenpairModule< T_Scalar >(_unknownFields[i].second);
              _unknownFields[i].second->setModule(module);
            }
            else {
              msg::debug << "A eigenpair module already exists in field '" << _unknownFields[i].second->name() << "', it is then overridden" << msg::endl;
              module->clear();
            }

            for(auto j = 0ULL; j < eigenvalues.size(); ++j) {
              static_cast< field::EigenpairModule< T_Scalar > * >(module)->assignEigenpair(eigenvalues[j], eigenvectors[j]);
            }
          }
        }
      }
    }
    catch(const std::exception &exc) {
      eigenvalues.clear();
      msg::error << "Unable to find" << (computeEigenvectors ? " eigenpairs" : " eigenvalues") << msg::endl;
      msg::error << exc.what() << msg::endl;
      msg::info << "Eigensolve aborted" << msg::endl;
      time.tock();
      return time;
    }

    time.tock();
    msg::info << "Done solving for eigenvalues in " << time << "s" << msg::endl;
    return time;
  }

  template< class T_Scalar >
  scalar::Precision< T_Scalar > Formulation< T_Scalar >::getResidual() const
  {
    algebra::MatrixCRSFast< T_Scalar > A;
    algebra::Vector< T_Scalar > x;
    algebra::Vector< T_Scalar > b;

    getLHS(A);
    getRHS(b);

    std::vector< T_Scalar > values(_dofs.nbrUnknownDofs());
    for(auto i = 0ULL; i < _unknownFields.size(); ++i) {
      if(field::fieldIsStillValid(_unknownFields[i].first)) {
        _unknownFields[i].second->getUnknownValues(values);
      }
    }

    x = std::move(values);

    return algebra::residual(b, A, x);
  }

  template< class T_Scalar >
  void Formulation< T_Scalar >::getLHS(algebra::Matrix< T_Scalar > &matrix) const
  {
    if(_A->getModule() != nullptr) {
      if(_A->getModule()->name() != "A" && _A->getModule()->name() != "AFrequency") {
        msg::error << "'getLHS' cannot be used with system assembled with separated matrices: use 'getMass', 'getDamping' or 'getStiffness'" << msg::endl;
        return;
      }
      matrix.extract(_A);
    }
  }

  template< class T_Scalar >
  void Formulation< T_Scalar >::getMass(algebra::Matrix< T_Scalar > &matrix) const
  {
    if(_A->getModule() != nullptr) {
      if(_A->getModule()->name() == "A" || _A->getModule()->name() == "AFrequency" || _A->getModule()->name() == "CK") {
        msg::error << "There is no mass matrix associated to system of type '" << _A->getModule()->name() << "'" << msg::endl;
        return;
      }
      _A->getModule()->activate('M');
      matrix.extract(_A);
    }
  }

  template< class T_Scalar >
  void Formulation< T_Scalar >::getDamping(algebra::Matrix< T_Scalar > &matrix) const
  {
    if(_A->getModule() != nullptr) {
      if(_A->getModule()->name() == "A" || _A->getModule()->name() == "AFrequency" || _A->getModule()->name() == "MK") {
        msg::error << "There is no damping matrix associated to system of type '" << _A->getModule()->name() << "'" << msg::endl;
        return;
      }
      _A->getModule()->activate('C');
      matrix.extract(_A);
    }
  }

  template< class T_Scalar >
  void Formulation< T_Scalar >::getStiffness(algebra::Matrix< T_Scalar > &matrix) const
  {
    if(_A->getModule() != nullptr) {
      if(_A->getModule()->name() == "A" || _A->getModule()->name() == "AFrequency") {
        msg::error << "There is no stiffness matrix associated to system of type '" << _A->getModule()->name() << "'" << msg::endl;
        return;
      }
      _A->getModule()->activate('K');
      matrix.extract(_A);
    }
  }

  template< class T_Scalar >
  static void s_getLHSBlock(algebra::Matrix< T_Scalar > &matrix, const std::vector< unsigned long long > &dofTags, const std::vector< unsigned long long > &tfTags, const system::MatrixFactory< T_Scalar > *A, const int tag)
  {
    algebra::MatrixCRSFast< T_Scalar > globalMatrix;
    globalMatrix.extract(A);

    std::vector< unsigned long long > ai;
    ai.reserve(dofTags.size());
    std::vector< unsigned long long > aj;
    aj.reserve(tfTags.size());
    std::vector< T_Scalar > a;
    a.reserve(tfTags.size());

    ai.push_back(0);
    unsigned long long currentLine = 0;
    for(auto i = 0ULL; i < globalMatrix.size(0); ++i) {
      if(i == dofTags[currentLine]) {
        ai.push_back(ai.back());
        auto it = tfTags.begin();
        for(auto jj = globalMatrix.ai()[i]; jj < globalMatrix.ai()[i + 1]; ++jj) {
          auto it2 = std::find(it, tfTags.end(), globalMatrix.aj()[jj]);
          if(it2 != tfTags.end()) {
            it = it2;
            aj.push_back(it2 - tfTags.begin());
            a.push_back(globalMatrix.a()[jj]);
            ai[ai.size() - 1]++;
          }
        }
        currentLine++;
      }
    }

    algebra::MatrixCRSFast< T_Scalar > blockMatrix(dofTags.size(), tfTags.size(), ai, aj, a);
    matrix = blockMatrix;
  }

  template< class T_Scalar >
  static void s_getLHSIndices(std::vector< unsigned long long > &dofTags, std::vector< unsigned long long > &tfTags, const field::FieldInterface< T_Scalar > &dof, const field::FieldInterface< T_Scalar > &tf)
  {
    dof.getUnknownIndices(dofTags);
    tf.getUnknownIndices(tfTags);
    std::sort(dofTags.begin(), dofTags.end());
    std::sort(tfTags.begin(), tfTags.end());
  }

  template< class T_Scalar >
  static void s_getLHSIndices(std::vector< unsigned long long > &dofTags, std::vector< unsigned long long > &tfTags, const std::vector< const field::FieldInterface< T_Scalar > * > &dofs, const std::vector< const field::FieldInterface< T_Scalar > * > &tfs)
  {
    for(auto i = 0ULL; i < dofs.size(); ++i) {
      std::vector< unsigned long long > localDofTags;
      dofs[i]->getUnknownIndices(localDofTags);
      dofTags.insert(dofTags.end(), localDofTags.begin(), localDofTags.end());
    }
    for(auto i = 0ULL; i < tfs.size(); ++i) {
      std::vector< unsigned long long > localTfTags;
      tfs[i]->getUnknownIndices(localTfTags);
      tfTags.insert(tfTags.end(), localTfTags.begin(), localTfTags.end());
    }
    std::sort(dofTags.begin(), dofTags.end());
    std::sort(tfTags.begin(), tfTags.end());
  }

  template< class T_Scalar >
  void Formulation< T_Scalar >::getLHSBlock(algebra::Matrix< T_Scalar > &matrix, const field::FieldInterface< T_Scalar > &dof, const field::FieldInterface< T_Scalar > &tf) const
  {
    if(_A->getModule()->name() != "A" && _A->getModule()->name() != "AFrequency") {
      msg::error << "'s_getLHSBlock' cannot be used with system assembled with separated matrices: use 'getMass', 'getDamping' or 'getStiffness'" << msg::endl;
      return;
    }
    if(matrix.format() == algebra::MatrixFormat::CRSFast) {
      msg::error << "Cannot use CRSFast format to get a block a matrix" << msg::endl;
      return;
    }

    std::vector< unsigned long long > dofTags, tfTags;
    s_getLHSIndices(dofTags, tfTags, dof, tf);
    s_getLHSBlock(matrix, dofTags, tfTags, _A, 0);
  }

  template< class T_Scalar >
  void Formulation< T_Scalar >::getLHSBlock(algebra::Matrix< T_Scalar > &matrix, const std::vector< const field::FieldInterface< T_Scalar > * > &dofs, const std::vector< const field::FieldInterface< T_Scalar > * > &tfs) const
  {
    if(_A->getModule()->name() != "A" && _A->getModule()->name() != "AFrequency") {
      msg::error << "'s_getLHSBlock' cannot be used with system assembled with separated matrices: use 'getMass', 'getDamping' or 'getStiffness'" << msg::endl;
      return;
    }
    if(matrix.format() == algebra::MatrixFormat::CRSFast) {
      msg::error << "Cannot use CRSFast format to get a block a matrix" << msg::endl;
      return;
    }

    std::vector< unsigned long long > dofTags, tfTags;
    s_getLHSIndices(dofTags, tfTags, dofs, tfs);
    s_getLHSBlock(matrix, dofTags, tfTags, _A, 0);
  }

  template< class T_Scalar >
  void Formulation< T_Scalar >::getMassBlock(algebra::Matrix< T_Scalar > &matrix, const field::FieldInterface< T_Scalar > &dof, const field::FieldInterface< T_Scalar > &tf) const
  {
    if(_A->getModule()->name() == "A" || _A->getModule()->name() == "AFrequency" || _A->getModule()->name() == "CK") {
      msg::error << "There is no mass matrix associated to system of type '" << _A->getModule()->name() << "'" << msg::endl;
      return;
    }
    if(matrix.format() == algebra::MatrixFormat::CRSFast) {
      msg::error << "Cannot use CRSFast format to get a block a matrix" << msg::endl;
      return;
    }

    std::vector< unsigned long long > dofTags, tfTags;
    s_getLHSIndices(dofTags, tfTags, dof, tf);
    _A->getModule()->activate('M');
    s_getLHSBlock(matrix, dofTags, tfTags, _A, 0);
  }

  template< class T_Scalar >
  void Formulation< T_Scalar >::getMassBlock(algebra::Matrix< T_Scalar > &matrix, const std::vector< const field::FieldInterface< T_Scalar > * > &dofs, const std::vector< const field::FieldInterface< T_Scalar > * > &tfs) const
  {
    if(_A->getModule()->name() != "A" || _A->getModule()->name() != "AFrequency" || _A->getModule()->name() == "CK") {
      msg::error << "There is no mass matrix associated to system of type '" << _A->getModule()->name() << "'" << msg::endl;
      return;
    }
    if(matrix.format() == algebra::MatrixFormat::CRSFast) {
      msg::error << "Cannot use CRSFast format to get a block a matrix" << msg::endl;
      return;
    }

    std::vector< unsigned long long > dofTags, tfTags;
    s_getLHSIndices(dofTags, tfTags, dofs, tfs);
    _A->getModule()->activate('M');
    s_getLHSBlock(matrix, dofTags, tfTags, _A, 0);
  }

  template< class T_Scalar >
  void Formulation< T_Scalar >::getDampingBlock(algebra::Matrix< T_Scalar > &matrix, const field::FieldInterface< T_Scalar > &dof, const field::FieldInterface< T_Scalar > &tf) const
  {
    if(_A->getModule()->name() == "A" || _A->getModule()->name() == "AFrequency" || _A->getModule()->name() == "MK") {
      msg::error << "There is no damping matrix associated to system of type '" << _A->getModule()->name() << "'" << msg::endl;
      return;
    }
    if(matrix.format() == algebra::MatrixFormat::CRSFast) {
      msg::error << "Cannot use CRSFast format to get a block a matrix" << msg::endl;
      return;
    }

    std::vector< unsigned long long > dofTags, tfTags;
    s_getLHSIndices(dofTags, tfTags, dof, tf);
    _A->getModule()->activate('C');
    s_getLHSBlock(matrix, dofTags, tfTags, _A, 0);
  }

  template< class T_Scalar >
  void Formulation< T_Scalar >::getDampingBlock(algebra::Matrix< T_Scalar > &matrix, const std::vector< const field::FieldInterface< T_Scalar > * > &dofs, const std::vector< const field::FieldInterface< T_Scalar > * > &tfs) const
  {
    if(_A->getModule()->name() != "A" || _A->getModule()->name() != "AFrequency" || _A->getModule()->name() == "MK") {
      msg::error << "There is no damping matrix associated to system of type '" << _A->getModule()->name() << "'" << msg::endl;
      return;
    }
    if(matrix.format() == algebra::MatrixFormat::CRSFast) {
      msg::error << "Cannot use CRSFast format to get a block a matrix" << msg::endl;
      return;
    }

    std::vector< unsigned long long > dofTags, tfTags;
    s_getLHSIndices(dofTags, tfTags, dofs, tfs);
    _A->getModule()->activate('C');
    s_getLHSBlock(matrix, dofTags, tfTags, _A, 0);
  }

  template< class T_Scalar >
  void Formulation< T_Scalar >::getStiffnessBlock(algebra::Matrix< T_Scalar > &matrix, const field::FieldInterface< T_Scalar > &dof, const field::FieldInterface< T_Scalar > &tf) const
  {
    if(_A->getModule()->name() == "A" || _A->getModule()->name() == "AFrequency") {
      msg::error << "There is no stiffness matrix associated to system of type '" << _A->getModule()->name() << "'" << msg::endl;
      return;
    }
    if(matrix.format() == algebra::MatrixFormat::CRSFast) {
      msg::error << "Cannot use CRSFast format to get a block a matrix" << msg::endl;
      return;
    }

    std::vector< unsigned long long > dofTags, tfTags;
    s_getLHSIndices(dofTags, tfTags, dof, tf);
    _A->getModule()->activate('K');
    s_getLHSBlock(matrix, dofTags, tfTags, _A, 0);
  }

  template< class T_Scalar >
  void Formulation< T_Scalar >::getStiffnessBlock(algebra::Matrix< T_Scalar > &matrix, const std::vector< const field::FieldInterface< T_Scalar > * > &dofs, const std::vector< const field::FieldInterface< T_Scalar > * > &tfs) const
  {
    if(_A->getModule()->name() != "A" || _A->getModule()->name() != "AFrequency") {
      msg::error << "There is no stiffness matrix associated to system of type '" << _A->getModule()->name() << "'" << msg::endl;
      return;
    }
    if(matrix.format() == algebra::MatrixFormat::CRSFast) {
      msg::error << "Cannot use CRSFast format to get a block a matrix" << msg::endl;
      return;
    }

    std::vector< unsigned long long > dofTags, tfTags;
    s_getLHSIndices(dofTags, tfTags, dofs, tfs);
    _A->getModule()->activate('K');
    s_getLHSBlock(matrix, dofTags, tfTags, _A, 0);
  }

  template< class T_Scalar >
  void Formulation< T_Scalar >::getRHS(algebra::Vector< T_Scalar > &vector) const
  {
    vector.extract(&_multiB.at(0));
  }

  template< class T_Scalar >
  static void s_getRHSBlock(algebra::Vector< T_Scalar > &vector, const std::vector< unsigned long long > &tfTags, const system::VectorFactory< T_Scalar > *b)
  {
    std::vector< T_Scalar > vec(tfTags.size());
    algebra::Vector< T_Scalar > globalVector;
    globalVector.extract(b);
    for(auto i = 0ULL; i < tfTags.size(); ++i) {
      vec[i] = globalVector[tfTags[i]];
    }
    vector = std::move(vec);
  }

  template< class T_Scalar >
  void Formulation< T_Scalar >::getRHSBlock(algebra::Vector< T_Scalar > &vector, const field::FieldInterface< T_Scalar > &tf) const
  {
    std::vector< unsigned long long > tfTags;
    tf.getUnknownIndices(tfTags);
    std::sort(tfTags.begin(), tfTags.end());
    s_getRHSBlock(vector, tfTags, &_multiB.at(0));
  }

  template< class T_Scalar >
  void Formulation< T_Scalar >::getRHSBlock(algebra::Vector< T_Scalar > &vector, const std::vector< const field::FieldInterface< T_Scalar > * > &tfs) const
  {
    std::vector< unsigned long long > tfTags;
    for(auto i = 0ULL; i < tfs.size(); ++i) {
      std::vector< unsigned long long > localTfTags;
      tfs[i]->getUnknownIndices(localTfTags);
      tfTags.insert(tfTags.end(), localTfTags.begin(), localTfTags.end());
    }
    std::sort(tfTags.begin(), tfTags.end());
    s_getRHSBlock(vector, tfTags, &_multiB.at(0));
  }

  template< class T_Scalar >
  void Formulation< T_Scalar >::setRHS(const algebra::Vector< T_Scalar > &vector)
  {
    _multiB.at(0).setToZero();
    for(auto i = 0ULL; i < vector.size(); ++i) {
      _multiB.at(0).addValue(i, vector[i]);
    }
  }

  template< class T_Scalar >
  static void s_setRHSBlock(const algebra::Vector< T_Scalar > &vector, const std::vector< unsigned long long > &tfTags, system::VectorFactory< T_Scalar > *b)
  {
    b->setToZero(tfTags);
    for(auto i = 0ULL; i < tfTags.size(); ++i) {
      b->addValue(tfTags[i], vector[i]);
    }
  }

  template< class T_Scalar >
  void Formulation< T_Scalar >::setRHSBlock(const algebra::Vector< T_Scalar > &vector, const field::FieldInterface< T_Scalar > &tf)
  {
    std::vector< unsigned long long > tfTags;
    tf.getUnknownIndices(tfTags);
    std::sort(tfTags.begin(), tfTags.end());
    s_setRHSBlock(vector, tfTags, &_multiB.at(0));
  }

  template< class T_Scalar >
  void Formulation< T_Scalar >::setRHSBlock(const algebra::Vector< T_Scalar > &vector, const std::vector< const field::FieldInterface< T_Scalar > * > &tfs)
  {
    std::vector< unsigned long long > tfTags;
    for(auto i = 0ULL; i < tfs.size(); ++i) {
      std::vector< unsigned long long > localTfTags;
      tfs[i]->getUnknownIndices(localTfTags);
      tfTags.insert(tfTags.end(), localTfTags.begin(), localTfTags.end());
    }
    std::sort(tfTags.begin(), tfTags.end());
    s_setRHSBlock(vector, tfTags, &_multiB.at(0));
  }

  INSTANTIATE_CLASS(Formulation, 4, TEMPLATE_ARGS(std::complex< double >, std::complex< float >, double, float))

} // namespace gmshfem::problem
