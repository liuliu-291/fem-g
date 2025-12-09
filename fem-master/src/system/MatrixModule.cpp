// GmshFEM - Copyright (C) 2019-2022, A. Royer, E. Béchet, C. Geuzaine, Université de Liège
//
// See the LICENSE.txt file for license information. Please report all
// issues on https://gitlab.onelab.info/gmsh/fem/issues

#include "MatrixModule.h"

#include "Exception.h"
#include "MatrixFactory.h"
#include "instantiate.h"

namespace gmshfem::system
{


  //
  // class MatrixModule
  //

  template< class T_Scalar >
  MatrixModule< T_Scalar >::MatrixModule(MatrixFactory< T_Scalar > *parentMatrix) :
    _parentMatrix(parentMatrix)
  {
  }

  template< class T_Scalar >
  MatrixModule< T_Scalar >::~MatrixModule()
  {
  }

  INSTANTIATE_CLASS(MatrixModule, 4, TEMPLATE_ARGS(std::complex< double >, std::complex< float >, double, float))


  //
  // class AModule : public MatrixModule
  //

  template< class T_Scalar >
  AModule< T_Scalar >::AModule(MatrixFactory< T_Scalar > *parentMatrix) :
    MatrixModule< T_Scalar >(parentMatrix), _a(parentMatrix->numberOfNonZeros(), 0.)
  {
  }

  template< class T_Scalar >
  AModule< T_Scalar >::~AModule()
  {
  }

  template< class T_Scalar >
  void AModule< T_Scalar >::setToZero()
  {
    for(auto i = 0ULL; i < _a.size(); ++i) {
      _a[i] = 0.;
    }
  }

  template< class T_Scalar >
  std::string AModule< T_Scalar >::name() const
  {
    return "A";
  }

  template< class T_Scalar >
  common::Memory AModule< T_Scalar >::memory() const
  {
    common::Memory mem = sizeof(this->_parentMatrix) + sizeof(_a) + _a.size() * sizeof(T_Scalar);
    return mem;
  }

  template< class T_Scalar >
  const T_Scalar &AModule< T_Scalar >::operator[](const unsigned int index) const
  {
    return _a[index];
  }

  template< class T_Scalar >
  T_Scalar &AModule< T_Scalar >::operator[](const unsigned int index)
  {
    return _a[index];
  }

  template< class T_Scalar >
  void AModule< T_Scalar >::activate(const char &param) const
  {
  }

  template< class T_Scalar >
  const T_Scalar *AModule< T_Scalar >::getMatrix() const
  {
    return _a.data();
  }

  INSTANTIATE_CLASS(AModule, 4, TEMPLATE_ARGS(std::complex< double >, std::complex< float >, double, float))


  //
  // class AFrequencyModule : public MatrixModule
  //

  template< class T_CPScalar >
  AFrequencyModule< T_CPScalar >::AFrequencyModule(MatrixFactory< T_CPScalar > *parentMatrix) :
    MatrixModule< T_CPScalar >(parentMatrix), _a(parentMatrix->numberOfNonZeros(), 0.), _frequency(0.)
  {
  }

  template< class T_CPScalar >
  AFrequencyModule< T_CPScalar >::AFrequencyModule(MatrixFactory< T_CPScalar > *parentMatrix, std::vector< T_CPScalar > &matrix) :
    MatrixModule< T_CPScalar >(parentMatrix), _a(std::move(matrix)), _frequency(0.)
  {
  }

  template< class T_CPScalar >
  AFrequencyModule< T_CPScalar >::~AFrequencyModule()
  {
  }

  template< class T_CPScalar >
  void AFrequencyModule< T_CPScalar >::setToZero()
  {
    for(auto i = 0ULL; i < _a.size(); ++i) {
      _a[i] = 0.;
    }
  }

  template< class T_CPScalar >
  std::string AFrequencyModule< T_CPScalar >::name() const
  {
    return "AFrequency";
  }

  template< class T_CPScalar >
  common::Memory AFrequencyModule< T_CPScalar >::memory() const
  {
    common::Memory mem = sizeof(this->_parentMatrix) + sizeof(_a) + _a.size() * sizeof(T_CPScalar);
    return mem;
  }

  template< class T_CPScalar >
  const T_CPScalar &AFrequencyModule< T_CPScalar >::operator[](const unsigned int index) const
  {
    return _a[index];
  }

  template< class T_CPScalar >
  T_CPScalar &AFrequencyModule< T_CPScalar >::operator[](const unsigned int index)
  {
    return _a[index];
  }

  template< class T_CPScalar >
  void AFrequencyModule< T_CPScalar >::activate(const char &param) const
  {
  }

  template< class T_CPScalar >
  void AFrequencyModule< T_CPScalar >::setFrequency(const scalar::Precision< T_CPScalar > &frequency)
  {
    _frequency = frequency;
  }

  template< class T_CPScalar >
  scalar::Precision< T_CPScalar > AFrequencyModule< T_CPScalar >::getFrequency() const
  {
    return _frequency;
  }

  template< class T_CPScalar >
  const T_CPScalar *AFrequencyModule< T_CPScalar >::getMatrix() const
  {
    return _a.data();
  }

  INSTANTIATE_CLASS(AFrequencyModule, 4, TEMPLATE_ARGS(std::complex< double >, std::complex< float >, double, float))


  //
  // class MKModule : public MatrixModule
  //

  template< class T_Scalar >
  MKModule< T_Scalar >::MKModule(MatrixFactory< T_Scalar > *parentMatrix) :
    MatrixModule< T_Scalar >(parentMatrix), _a(2, std::vector< T_Scalar >(parentMatrix->numberOfNonZeros(), 0.)), _activeComponents(0)
  {
  }

  template< class T_Scalar >
  MKModule< T_Scalar >::~MKModule()
  {
  }

  template< class T_Scalar >
  void MKModule< T_Scalar >::setToZero()
  {
    for(auto i = 0; i < 2; ++i) {
      for(auto j = 0ULL; j < _a[i].size(); ++j) {
        _a[i][j] = 0.;
      }
    }
  }

  template< class T_Scalar >
  std::string MKModule< T_Scalar >::name() const
  {
    return "MK";
  }

  template< class T_Scalar >
  common::Memory MKModule< T_Scalar >::memory() const
  {
    common::Memory mem = sizeof(this->_parentMatrix) + sizeof(_a) + _a.size() * _a[0].size() * sizeof(T_Scalar);
    return mem;
  }

  template< class T_Scalar >
  const T_Scalar &MKModule< T_Scalar >::operator[](const unsigned int index) const
  {
    return _a[_activeComponents][index];
  }

  template< class T_Scalar >
  T_Scalar &MKModule< T_Scalar >::operator[](const unsigned int index)
  {
    return _a[_activeComponents][index];
  }

  template< class T_Scalar >
  void MKModule< T_Scalar >::activate(const char &param) const
  {
    if(param == 'M') {
      _activeComponents = 0;
    }
    else if(param == 'K') {
      _activeComponents = 1;
    }
    else {
      throw common::Exception("Parameter '" + std::string(1, param) + "' has no meaning for the MKModule");
    }
  }

  template< class T_Scalar >
  const T_Scalar *MKModule< T_Scalar >::getMatrix() const
  {
    return _a[_activeComponents].data();
  }

  INSTANTIATE_CLASS(MKModule, 4, TEMPLATE_ARGS(std::complex< double >, std::complex< float >, double, float))


  //
  // class CKModule : public MatrixModule
  //

  template< class T_Scalar >
  CKModule< T_Scalar >::CKModule(MatrixFactory< T_Scalar > *parentMatrix) :
    MatrixModule< T_Scalar >(parentMatrix), _a(2, std::vector< T_Scalar >(parentMatrix->numberOfNonZeros(), 0.)), _activeComponents(0)
  {
  }

  template< class T_Scalar >
  CKModule< T_Scalar >::~CKModule()
  {
  }

  template< class T_Scalar >
  void CKModule< T_Scalar >::setToZero()
  {
    for(auto i = 0; i < 2; ++i) {
      for(auto j = 0ULL; j < _a[i].size(); ++j) {
        _a[i][j] = 0.;
      }
    }
  }

  template< class T_Scalar >
  std::string CKModule< T_Scalar >::name() const
  {
    return "CK";
  }

  template< class T_Scalar >
  common::Memory CKModule< T_Scalar >::memory() const
  {
    common::Memory mem = sizeof(this->_parentMatrix) + sizeof(_a) + _a.size() * _a[0].size() * sizeof(T_Scalar);
    return mem;
  }

  template< class T_Scalar >
  const T_Scalar &CKModule< T_Scalar >::operator[](const unsigned int index) const
  {
    return _a[_activeComponents][index];
  }

  template< class T_Scalar >
  T_Scalar &CKModule< T_Scalar >::operator[](const unsigned int index)
  {
    return _a[_activeComponents][index];
  }

  template< class T_Scalar >
  void CKModule< T_Scalar >::activate(const char &param) const
  {
    if(param == 'C') {
      _activeComponents = 0;
    }
    else if(param == 'K') {
      _activeComponents = 1;
    }
    else {
      throw common::Exception("Parameter '" + std::string(1, param) + "' has no meaning for the CKModule");
    }
  }

  template< class T_Scalar >
  const T_Scalar *CKModule< T_Scalar >::getMatrix() const
  {
    return _a[_activeComponents].data();
  }

  INSTANTIATE_CLASS(CKModule, 4, TEMPLATE_ARGS(std::complex< double >, std::complex< float >, double, float))


  //
  // class MCKModule : public MatrixModule
  //

  template< class T_Scalar >
  MCKModule< T_Scalar >::MCKModule(MatrixFactory< T_Scalar > *parentMatrix) :
    MatrixModule< T_Scalar >(parentMatrix), _a(3, std::vector< T_Scalar >(parentMatrix->numberOfNonZeros(), 0.)), _activeComponents(0)
  {
  }

  template< class T_Scalar >
  MCKModule< T_Scalar >::~MCKModule()
  {
  }

  template< class T_Scalar >
  void MCKModule< T_Scalar >::setToZero()
  {
    for(auto i = 0; i < 3; ++i) {
      for(auto j = 0ULL; j < _a[i].size(); ++j) {
        _a[i][j] = 0.;
      }
    }
  }

  template< class T_Scalar >
  std::string MCKModule< T_Scalar >::name() const
  {
    return "MCK";
  }

  template< class T_Scalar >
  common::Memory MCKModule< T_Scalar >::memory() const
  {
    common::Memory mem = sizeof(this->_parentMatrix) + sizeof(_a) + _a.size() * _a[0].size() * sizeof(T_Scalar);
    return mem;
  }

  template< class T_Scalar >
  const T_Scalar &MCKModule< T_Scalar >::operator[](const unsigned int index) const
  {
    return _a[_activeComponents][index];
  }

  template< class T_Scalar >
  T_Scalar &MCKModule< T_Scalar >::operator[](const unsigned int index)
  {
    return _a[_activeComponents][index];
  }

  template< class T_Scalar >
  void MCKModule< T_Scalar >::activate(const char &param) const
  {
    if(param == 'M') {
      _activeComponents = 0;
    }
    else if(param == 'C') {
      _activeComponents = 1;
    }
    else if(param == 'K') {
      _activeComponents = 2;
    }
    else {
      throw common::Exception("Parameter '" + std::string(1, param) + "' has no meaning for the CKModule");
    }
  }

  template< class T_Scalar >
  const T_Scalar *MCKModule< T_Scalar >::getMatrix() const
  {
    return _a[_activeComponents].data();
  }

  INSTANTIATE_CLASS(MCKModule, 4, TEMPLATE_ARGS(std::complex< double >, std::complex< float >, double, float))


} // namespace gmshfem::system
