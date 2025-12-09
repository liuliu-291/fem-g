// GmshFEM - Copyright (C) 2019-2022, A. Royer, E. Béchet, C. Geuzaine, Université de Liège
//
// See the LICENSE.txt file for license information. Please report all
// issues on https://gitlab.onelab.info/gmsh/fem/issues

#include "FieldModule.h"

#include "Dof.h"
#include "Exception.h"
#include "FieldInterface.h"
#include "instantiate.h"

namespace gmshfem::field
{


  //
  // class FieldModule
  //

  template< class T_Scalar >
  FieldModule< T_Scalar >::FieldModule(FieldInterface< T_Scalar > *parentField) :
    _parentField(parentField)
  {
  }

  template< class T_Scalar >
  FieldModule< T_Scalar >::FieldModule(const FieldModule< T_Scalar > &other) :
    _parentField(other._parentField)
  {
  }

  template< class T_Scalar >
  FieldModule< T_Scalar >::~FieldModule()
  {
  }

  INSTANTIATE_CLASS(FieldModule, 4, TEMPLATE_ARGS(std::complex< double >, std::complex< float >, double, float))


  //
  // class EigenpairModule : public FieldModule
  //

  template< class T_Scalar >
  EigenpairModule< T_Scalar >::EigenpairModule(FieldInterface< T_Scalar > *parentField) :
    FieldModule< T_Scalar >(parentField), _mapping(), _eigenpairs()
  {
  }

  template< class T_Scalar >
  EigenpairModule< T_Scalar >::EigenpairModule(const EigenpairModule< T_Scalar > &other) :
    FieldModule< T_Scalar >(other), _mapping(other._mapping), _eigenpairs(other._eigenpairs)
  {
  }

  template< class T_Scalar >
  EigenpairModule< T_Scalar >::~EigenpairModule()
  {
  }

  template< class T_Scalar >
  std::string EigenpairModule< T_Scalar >::name() const
  {
    return "Eigenpair";
  }

  template< class T_Scalar >
  common::Memory EigenpairModule< T_Scalar >::memory() const
  {
    common::Memory mem = _eigenpairs.size() + sizeof(std::pair< scalar::ComplexPrecision< T_Scalar >, std::vector< scalar::ComplexPrecision< T_Scalar > > >);
    for(auto i = 0ULL; i < _eigenpairs.size(); ++i) {
      mem += _eigenpairs[i].second.size() * sizeof(scalar::ComplexPrecision< T_Scalar >);
    }
    mem += sizeof(_mapping.size());
    unsigned int n = _mapping.bucket_count();
    double m = _mapping.max_load_factor();
    if(m > 1.0) {
      mem += n * m * sizeof(std::pair< unsigned int, unsigned int >);
    }
    else {
      mem += n * sizeof(std::pair< unsigned int, unsigned int >);
    }

    return mem;
  }

  template< class T_Scalar >
  void EigenpairModule< T_Scalar >::clear()
  {
    _eigenpairs.clear();
    _eigenpairs.shrink_to_fit();
    _mapping.clear();
  }

  template< class T_Scalar >
  EigenpairModule< T_Scalar > *EigenpairModule< T_Scalar >::copy() const
  {
    return new EigenpairModule< T_Scalar >(*this);
  }

  template< class T_Scalar >
  void EigenpairModule< T_Scalar >::assignEigenpair(const scalar::ComplexPrecision< T_Scalar > &eigenvalue, const std::vector< scalar::ComplexPrecision< T_Scalar > > &eigenvector)
  {
    std::vector< scalar::ComplexPrecision< T_Scalar > > myEigenvector;
    myEigenvector.reserve(this->_parentField->numberOfUnknownDofs());
    if(_mapping.size() == 0) {
      _mapping.reserve(this->_parentField->numberOfUnknownDofs());
      for(auto it = this->_parentField->begin(); it != this->_parentField->end(); ++it) {
        if(it->first->type() == dofs::Type::Unknown) {
          _mapping.insert(std::make_pair(it->first->numDof(), myEigenvector.size()));
          myEigenvector.push_back(eigenvector[it->first->numDof() - 1]);
        }
      }
    }
    else {
      for(auto it = this->_parentField->begin(); it != this->_parentField->end(); ++it) {
        if(it->first->type() == dofs::Type::Unknown) {
          myEigenvector.push_back(eigenvector[it->first->numDof() - 1]);
        }
      }
    }
    _eigenpairs.push_back(std::make_pair(eigenvalue, myEigenvector));
  }

  template< class T_Scalar >
  typename std::vector< std::pair< scalar::ComplexPrecision< T_Scalar >, std::vector< scalar::ComplexPrecision< T_Scalar > > > >::const_iterator EigenpairModule< T_Scalar >::begin() const
  {
    return _eigenpairs.begin();
  }

  template< class T_Scalar >
  typename std::vector< std::pair< scalar::ComplexPrecision< T_Scalar >, std::vector< scalar::ComplexPrecision< T_Scalar > > > >::const_iterator EigenpairModule< T_Scalar >::end() const
  {
    return _eigenpairs.end();
  }

  template< class T_Scalar >
  void EigenpairModule< T_Scalar >::getValues(const unsigned int eigenTag, const std::vector< int > &typeKeys, const std::vector< unsigned long long > &entityKeys, std::vector< scalar::ComplexPrecision< T_Scalar > > &values, const unsigned int begin, const unsigned int end) const
  {
    if(eigenTag > _eigenpairs.size()) {
      throw common::Exception("Trying to access the eigenpair number " + std::to_string(eigenTag) + " while there are only " + std::to_string(_eigenpairs.size()) + " stored in field '" + this->_parentField->name() + "'");
    }
    for(auto i = begin; i < end; ++i) {
      const dofs::Dof *dof = this->_parentField->searchDof(typeKeys[i], entityKeys[i]);
      if(dof != nullptr) {
        if(dof->type() == dofs::Type::Unknown) {
          values[i] = _eigenpairs[eigenTag].second[_mapping.at(dof->numDof())];
        }
        else {
          values[i] = 0.;
        }
      }
      else {
        values[i] = 0.;
      }
    }
  }

  INSTANTIATE_CLASS(EigenpairModule, 4, TEMPLATE_ARGS(std::complex< double >, std::complex< float >, double, float))


} // namespace gmshfem::field
