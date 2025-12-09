// GmshFEM - Copyright (C) 2019-2022, A. Royer, E. Béchet, C. Geuzaine, Université de Liège
//
// See the LICENSE.txt file for license information. Please report all
// issues on https://gitlab.onelab.info/gmsh/fem/issues

#include "VectorFactory.h"

#include "CSVio.h"
#include "Dof.h"
#include "Exception.h"
#include "IndiceBucket.h"
#include "Message.h"
#include "OmpInterface.h"
#include "PetscInterface.h"
#include "gmshfemDefines.h"
#include "instantiate.h"

#include <complex>

namespace gmshfem::system
{


  template< class T_Scalar >
  VectorFactory< T_Scalar >::VectorFactory() :
    _size(0), _vec(), _valuesDC(), _valuesLC(), _indicesLC(), _shouldBeDestroyedWithPetsc()
  {
  }

  template< class T_Scalar >
  VectorFactory< T_Scalar >::VectorFactory(const unsigned long long size) :
    _size(size), _vec(size, 0.), _valuesDC(), _valuesLC(), _indicesLC(), _shouldBeDestroyedWithPetsc()
  {
  }

  template< class T_Scalar >
  VectorFactory< T_Scalar >::~VectorFactory()
  {
    removeVector();
  }

  template< class T_Scalar >
  void VectorFactory< T_Scalar >::init(const unsigned long long size)
  {
    _size = size;
    removeVector();
    _vec.resize(_size, 0.);
  }

  template< class T_Scalar >
  void VectorFactory< T_Scalar >::removeVector()
  {
    _vec.clear();
    _vec.shrink_to_fit();
    _valuesDC.clear();
    _valuesDC.shrink_to_fit();
    _valuesLC.clear();
    _valuesLC.shrink_to_fit();
    _indicesLC.clear();
    _indicesLC.shrink_to_fit();
    removePetscData();
  }

  template< class T_Scalar >
  void VectorFactory< T_Scalar >::removePetscData()
  {
    for(auto i = 0ULL; i < _shouldBeDestroyedWithPetsc.size(); ++i) {
      std::free(_shouldBeDestroyedWithPetsc[i]);
    }
    _shouldBeDestroyedWithPetsc.clear();
  }

  template< class T_Scalar >
  void VectorFactory< T_Scalar >::setToZero()
  {
    for(auto i = 0ULL; i < _vec.size(); ++i) {
      _vec[i] = 0.;
    }
    _valuesDC.clear();
    _valuesLC.clear();
    _indicesLC.clear();
  }

  template< class T_Scalar >
  void VectorFactory< T_Scalar >::setToZero(const std::vector< unsigned long long > &indices)
  {
    for(auto i = 0ULL; i < indices.size(); ++i) {
      _vec[indices[i]] = 0.;
    }
  }

  template< class T_Scalar >
  void VectorFactory< T_Scalar >::setValuesDC(std::vector< T_Scalar > &values)
  {
    _valuesDC = std::move(values);
  }

  template< class T_Scalar >
  void VectorFactory< T_Scalar >::setValuesLC(std::vector< T_Scalar > &values)
  {
    _valuesLC = std::move(values);
  }

  template< class T_Scalar >
  void VectorFactory< T_Scalar >::setIndicesLC(std::vector< unsigned long long > &indices)
  {
    _indicesLC = std::move(indices);
  }

  template< class T_Scalar >
  void VectorFactory< T_Scalar >::addValue(const unsigned long long n, const T_Scalar &value)
  {
    _vec[n] += value;
  }

  template< class T_Scalar >
  struct HackComplex {
    T_Scalar r, i;
  };

  template< class T_Scalar >
  void VectorFactory< T_Scalar >::addValues(const unsigned long long n, const std::pair< unsigned long long, int > *const index, const T_Scalar *const values)
  {
    for(auto i = 0ULL; i < n; ++i) {
      if(index[i].second == dofs::AssembleType::Unknown || index[i].second == dofs::AssembleType::UnknownBubble || index[i].second == dofs::AssembleType::UnknownGlobal) {
        HackComplex< scalar::Precision< T_Scalar > > *hack = reinterpret_cast< HackComplex< scalar::Precision< T_Scalar > > * >(&_vec[index[i].first]);
        const scalar::Precision< T_Scalar > re = std::real(values[i]);
        const scalar::Precision< T_Scalar > im = std::imag(values[i]);
#pragma omp atomic update
        hack->r += re;
#pragma omp atomic update
        hack->i += im;
      }
      else if(index[i].second == dofs::AssembleType::Linked || index[i].second == dofs::AssembleType::LinkedBubble) {
        HackComplex< scalar::Precision< T_Scalar > > *hack = reinterpret_cast< HackComplex< scalar::Precision< T_Scalar > > * >(&_vec[_indicesLC[index[i].first]]);
        const T_Scalar value = std::conj(_valuesLC[index[i].first]) * values[i];
        const scalar::Precision< T_Scalar > re = std::real(value);
        const scalar::Precision< T_Scalar > im = std::imag(value);
#pragma omp atomic update
        hack->r += re;
#pragma omp atomic update
        hack->i += im;
      }
    }
  }

  template<>
  void VectorFactory< double >::addValues(const unsigned long long n, const std::pair< unsigned long long, int > *const index, const double *const values)
  {
    for(auto i = 0ULL; i < n; ++i) {
      if(index[i].second == dofs::AssembleType::Unknown || index[i].second == dofs::AssembleType::UnknownBubble || index[i].second == dofs::AssembleType::UnknownGlobal) {
#pragma omp atomic update
        _vec[index[i].first] += values[i];
      }
      else if(index[i].second == dofs::AssembleType::Linked || index[i].second == dofs::AssembleType::LinkedBubble) {
#pragma omp atomic update
        _vec[_indicesLC[index[i].first]] += _valuesLC[index[i].first] * values[i];
      }
    }
  }

  template<>
  void VectorFactory< float >::addValues(const unsigned long long n, const std::pair< unsigned long long, int > *const index, const float *const values)
  {
    for(auto i = 0ULL; i < n; ++i) {
      if(index[i].second == dofs::AssembleType::Unknown || index[i].second == dofs::AssembleType::UnknownBubble || index[i].second == dofs::AssembleType::UnknownGlobal) {
#pragma omp atomic update
        _vec[index[i].first] += values[i];
      }
      else if(index[i].second == dofs::AssembleType::Linked || index[i].second == dofs::AssembleType::LinkedBubble) {
#pragma omp atomic update
        _vec[_indicesLC[index[i].first]] += _valuesLC[index[i].first] * values[i];
      }
    }
  }

  template< class T_Scalar >
  void VectorFactory< T_Scalar >::addValuesDC(const unsigned long long m, const unsigned long long n, const std::pair< unsigned long long, int > *const indexI, const std::pair< unsigned long long, int > *const indexJ, T_Scalar *const values)
  {
    for(auto i = 0ULL; i < m; ++i) {
      if(indexI[i].second == dofs::AssembleType::Unknown || indexI[i].second == dofs::AssembleType::UnknownBubble || indexI[i].second == dofs::AssembleType::UnknownGlobal) {
        for(auto j = 0ULL; j < n; ++j) {
          if(indexJ[j].second == dofs::AssembleType::Fixed) {
            const T_Scalar value = values[j * n + i] * _valuesDC[indexJ[j].first];
            HackComplex< scalar::Precision< T_Scalar > > *hack = reinterpret_cast< HackComplex< scalar::Precision< T_Scalar > > * >(&_vec[indexI[i].first]);
            const scalar::Precision< scalar::Precision< T_Scalar > > re = std::real(value);
            const scalar::Precision< scalar::Precision< T_Scalar > > im = std::imag(value);
#pragma omp atomic update
            hack->r += re;
#pragma omp atomic update
            hack->i += im;
          }
        }
      }
      else if(indexI[i].second == dofs::AssembleType::Linked || indexI[i].second == dofs::AssembleType::LinkedBubble) {
        for(auto j = 0ULL; j < n; ++j) {
          if(indexJ[j].second == dofs::AssembleType::Fixed) {
            const T_Scalar value = std::conj(_valuesLC[indexI[i].first]) * values[j * n + i] * _valuesDC[indexJ[j].first];
            HackComplex< scalar::Precision< T_Scalar > > *hack = reinterpret_cast< HackComplex< scalar::Precision< T_Scalar > > * >(&_vec[_indicesLC[indexI[i].first]]);
            const scalar::Precision< scalar::Precision< T_Scalar > > re = std::real(value);
            const scalar::Precision< scalar::Precision< T_Scalar > > im = std::imag(value);
#pragma omp atomic update
            hack->r += re;
#pragma omp atomic update
            hack->i += im;
          }
        }
      }
    }
  }

  template<>
  void VectorFactory< double >::addValuesDC(const unsigned long long m, const unsigned long long n, const std::pair< unsigned long long, int > *const indexI, const std::pair< unsigned long long, int > *const indexJ, double *const values)
  {
    for(auto i = 0ULL; i < m; ++i) {
      if(indexI[i].second == dofs::AssembleType::Unknown || indexI[i].second == dofs::AssembleType::UnknownBubble || indexI[i].second == dofs::AssembleType::UnknownGlobal) {
        for(auto j = 0ULL; j < n; ++j) {
          if(indexJ[j].second == dofs::AssembleType::Fixed) {
            const double value = values[j * n + i] * _valuesDC[indexJ[j].first];
#pragma omp atomic update
            _vec[indexI[i].first] += value;
          }
        }
      }
      else if(indexI[i].second == dofs::AssembleType::Linked || indexI[i].second == dofs::AssembleType::LinkedBubble) {
        for(auto j = 0ULL; j < n; ++j) {
          if(indexJ[j].second == dofs::AssembleType::Fixed) {
            const double value = values[j * n + i] * _valuesDC[indexJ[j].first];
#pragma omp atomic update
            _vec[_indicesLC[indexI[i].first]] += _valuesLC[indexI[i].first] * value;
          }
        }
      }
    }
  }

  template<>
  void VectorFactory< float >::addValuesDC(const unsigned long long m, const unsigned long long n, const std::pair< unsigned long long, int > *const indexI, const std::pair< unsigned long long, int > *const indexJ, float *const values)
  {
    for(auto i = 0ULL; i < m; ++i) {
      if(indexI[i].second == dofs::AssembleType::Unknown || indexI[i].second == dofs::AssembleType::UnknownBubble || indexI[i].second == dofs::AssembleType::UnknownGlobal) {
        for(auto j = 0ULL; j < n; ++j) {
          if(indexJ[j].second == dofs::AssembleType::Fixed) {
            const float value = values[j * n + i] * _valuesDC[indexJ[j].first];
#pragma omp atomic update
            _vec[indexI[i].first] += value;
          }
        }
      }
      else if(indexI[i].second == dofs::AssembleType::Linked || indexI[i].second == dofs::AssembleType::LinkedBubble) {
        for(auto j = 0ULL; j < n; ++j) {
          if(indexJ[j].second == dofs::AssembleType::Fixed) {
            const double value = values[j * n + i] * _valuesDC[indexJ[j].first];
#pragma omp atomic update
            _vec[_indicesLC[indexI[i].first]] += _valuesLC[indexI[i].first] * value;
          }
        }
      }
    }
  }

#ifdef HAVE_PETSC
  template< class T_Scalar >
  Vec VectorFactory< T_Scalar >::getPetsc() const
  {
    Vec vecPetsc;
    PetscInterface< T_Scalar, PetscScalar > interface;
    interface(_size, &_vec[0], &vecPetsc, _shouldBeDestroyedWithPetsc);
    return vecPetsc;
  }
#endif

  template< class T_Scalar >
  common::Memory VectorFactory< T_Scalar >::memory() const
  {
    return common::Memory((_vec.size() + _valuesDC.size()) * sizeof(T_Scalar));
  }

  template< class T_Scalar >
  unsigned long long VectorFactory< T_Scalar >::size() const
  {
    return _size;
  }

  INSTANTIATE_CLASS(VectorFactory, 4, TEMPLATE_ARGS(std::complex< double >, std::complex< float >, double, float))


} // namespace gmshfem::system
