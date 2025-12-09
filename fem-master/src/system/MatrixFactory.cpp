// GmshFEM - Copyright (C) 2019-2022, A. Royer, E. Béchet, C. Geuzaine, Université de Liège
//
// See the LICENSE.txt file for license information. Please report all
// issues on https://gitlab.onelab.info/gmsh/fem/issues

#include "MatrixFactory.h"

#include "CSVio.h"
#include "Dof.h"
#include "Exception.h"
#include "KahanSum.h"
#include "Message.h"
#include "PPMio.h"
#include "PetscInterface.h"
#include "gmshfemDefines.h"
#include "instantiate.h"

namespace gmshfem::system
{


  // *****************************
  // MatrixFactory
  // *****************************

  // private
  template< class T_Scalar >
  void MatrixFactory< T_Scalar >::_reorder(const std::vector< unsigned long long > &degree)
  {
    unsigned long long *newAJ = new unsigned long long[_ai[_size]];
    unsigned long long *newAI = new unsigned long long[_size + 1];
#pragma omp parallel num_threads(omp::getMaxThreads())
    {
#pragma omp for
      for(auto i = 0ULL; i < _size; ++i) {
        newAI[degree[i] - 1] = _ai[i + 1] - _ai[i];
      }
#pragma omp single
      {
        for(auto i = 1ULL; i < _size; ++i) {
          newAI[i] += newAI[i - 1];
        }
        for(auto i = _size; i > 0; --i) {
          newAI[i] = newAI[i - 1];
        }
        newAI[0] = 0;
      }
#pragma omp for
      for(auto i = 0ULL; i < _size; ++i) {
        for(auto j = _ai[i]; j < _ai[i + 1]; ++j) {
          newAJ[newAI[degree[i] - 1] + j - _ai[i]] = degree[_aj[j]] - 1;
        }
      }
#pragma omp single
      {
        std::swap(_aj, newAJ);
        std::swap(_ai, newAI);
        delete[] newAJ;
        delete[] newAI;
      }
#pragma omp for
      for(auto i = 0ULL; i < _size; ++i) {
        std::sort(&_aj[_ai[i]], &_aj[_ai[i + 1]]);
      }
    }
  }

  template< class T_Scalar >
  MatrixFactory< T_Scalar >::MatrixFactory(const std::string &options) :
    _size(0), _options(options), _ai(nullptr), _aj(nullptr), _module(nullptr), _valuesLC(), _indicesLC(), _pattern(nullptr), _patternMemory(nullptr), _havePattern(false), _locks(), _shouldBeDestroyedWithPetsc()
  {
  }

  template< class T_Scalar >
  MatrixFactory< T_Scalar >::MatrixFactory(const unsigned long long size, const std::string &options, const unsigned long long nbrBubble) :
    _size(size), _options(options), _ai(nullptr), _aj(nullptr), _module(nullptr), _valuesLC(), _indicesLC(), _pattern(nullptr), _patternMemory(nullptr), _havePattern(false), _locks(), _shouldBeDestroyedWithPetsc()
  {
    init(size, nbrBubble);
  }

  template< class T_Scalar >
  MatrixFactory< T_Scalar >::~MatrixFactory()
  {
    removeMatrix();
  }

  template< class T_Scalar >
  void MatrixFactory< T_Scalar >::setModule(MatrixModule< T_Scalar > *module)
  {
    _module = module;
  }

  template< class T_Scalar >
  const MatrixModule< T_Scalar > *MatrixFactory< T_Scalar >::getModule() const
  {
    return _module;
  }

  template< class T_Scalar >
  void MatrixFactory< T_Scalar >::init(const unsigned long long size, const unsigned long long nbrBubble)
  {
    removeMatrix();
    _size = size;
    _locks.clear();
    _locks.resize(_size - nbrBubble);
    _patternMemory = new std::vector< unsigned long long >[size];
    _pattern = new std::vector< unsigned long long > *[size];
    for(auto i = 0ULL; i < size; ++i) {
      _pattern[i] = &_patternMemory[i];
    }
    _havePattern = false;
  }

  template< class T_Scalar >
  void MatrixFactory< T_Scalar >::removeMatrix()
  {
    if(_ai) {
      delete[] _ai;
      _ai = nullptr;
    }

    if(_aj) {
      delete[] _aj;
      _aj = nullptr;
    }

    if(_module) {
      delete _module;
      _module = nullptr;
    }

    removePetscData();

    _options.clean();
    _size = 0;
    _locks.clear();

    if(_patternMemory) {
      delete[] _patternMemory;
      _patternMemory = nullptr;
    }

    if(_pattern) {
      delete[] _pattern;
      _pattern = nullptr;
    }
    _havePattern = false;
    _valuesLC.clear();
    _valuesLC.shrink_to_fit();
    _indicesLC.clear();
    _indicesLC.shrink_to_fit();
  }

  template< class T_Scalar >
  void MatrixFactory< T_Scalar >::removePetscData()
  {
    for(auto i = 0ULL; i < _shouldBeDestroyedWithPetsc.size(); ++i) {
      std::free(_shouldBeDestroyedWithPetsc[i]);
    }
    _shouldBeDestroyedWithPetsc.clear();
  }

  template< class T_Scalar >
  void MatrixFactory< T_Scalar >::setPatternToZero()
  {
    setToZero();
    if(_havePattern) {
      for(auto i = 0ULL; i < _ai[_size]; ++i) {
        _aj[i] = 0;
      }
      for(auto i = 0ULL; i < _size; ++i) {
        _ai[i] = 0;
      }
    }

    _havePattern = false;
  }

  template< class T_Scalar >
  void MatrixFactory< T_Scalar >::setToZero()
  {
    if(_module) {
      _module->setToZero();
    }
  }

  template< class T_Scalar >
  void MatrixFactory< T_Scalar >::finalizePattern()
  {
    _ai = new unsigned long long[_size + 1];

    _ai[0] = 0;
    for(auto i = 1ULL; i <= _size; ++i) {
      _ai[i] = _pattern[i - 1]->size() + _ai[i - 1];
    }

    _aj = new unsigned long long[_ai[_size]];

#pragma omp parallel for num_threads(omp::getMaxThreads())
    for(auto i = 0ULL; i < _size; ++i) {
      unsigned long long *ptr = &(*_pattern[i])[0];
      unsigned long long *ptrA = &_aj[_ai[i]];
      for(auto j = 0ULL; j < _ai[i + 1] - _ai[i]; ++j) {
        ptrA[j] = ptr[j];
      }
    }
    delete[] _pattern;
    _pattern = nullptr;
    delete[] _patternMemory;
    _patternMemory = nullptr;

    _havePattern = true;
  }

  template< class T_Scalar >
  void MatrixFactory< T_Scalar >::setValuesLC(std::vector< T_Scalar > &values)
  {
    _valuesLC = std::move(values);
  }

  template< class T_Scalar >
  void MatrixFactory< T_Scalar >::setIndicesLC(std::vector< unsigned long long > &indices)
  {
    _indicesLC = std::move(indices);
  }

  template< class T_Scalar >
  struct HackComplex {
    T_Scalar r, i;
  };

  template< class T_Scalar >
  bool MatrixFactory< T_Scalar >::addValues(const unsigned long long m, const unsigned long long n, const std::pair< unsigned long long, int > *const indexI, const std::pair< unsigned long long, int > *const indexJ, T_Scalar *const values)
  {
    bool haveDC = false;
    const bool sym = _options.symmetric() || _options.hermitian();

    for(auto i = 0ULL; i < m; ++i) {
      unsigned long long iIndex = 0;
      if(indexI[i].second == dofs::AssembleType::Unknown || indexI[i].second == dofs::AssembleType::UnknownBubble || indexI[i].second == dofs::AssembleType::UnknownGlobal) {
        iIndex = indexI[i].first;
      }
      else if(indexI[i].second == dofs::AssembleType::Linked || indexI[i].second == dofs::AssembleType::LinkedBubble) {
        iIndex = _indicesLC[indexI[i].first];
      }
      else if(indexI[i].second == dofs::AssembleType::Fixed) {
        haveDC = true;
        continue;
      }

      if(indexI[i].second == dofs::AssembleType::Unknown || indexI[i].second == dofs::AssembleType::UnknownGlobal) {
        for(auto j = 0ULL; j < n; ++j) {
          unsigned long long jIndex = 0;
          if(indexJ[j].second == dofs::AssembleType::Unknown || indexJ[j].second == dofs::AssembleType::UnknownBubble || indexJ[j].second == dofs::AssembleType::UnknownGlobal) {
            jIndex = indexJ[j].first;
          }
          else if(indexJ[j].second == dofs::AssembleType::Linked || indexJ[j].second == dofs::AssembleType::LinkedBubble) {
            jIndex = _indicesLC[indexJ[j].first];
          }
          else if(indexJ[j].second == dofs::AssembleType::Fixed) {
            haveDC = true;
            continue;
          }

          if(sym && jIndex < iIndex) continue; // upper triangular matrix

          const unsigned long long *place = std::lower_bound(&_aj[_ai[iIndex]], &_aj[_ai[iIndex + 1]], jIndex);
          const unsigned long long pos = place - _aj;
          if(indexJ[j].second == dofs::AssembleType::Unknown || indexJ[j].second == dofs::AssembleType::UnknownGlobal) {
            if constexpr(scalar::IsComplex< T_Scalar >::value) {
              HackComplex< scalar::Precision< T_Scalar > > *hack = reinterpret_cast< HackComplex< scalar::Precision< T_Scalar > > * >(&(*_module)[pos]);
              const scalar::Precision< T_Scalar > re = std::real(values[i * n + j]);
              const scalar::Precision< T_Scalar > im = std::imag(values[i * n + j]);
#pragma omp atomic update
              hack->r += re;
#pragma omp atomic update
              hack->i += im;
            }
            else {
              T_Scalar &update = (*_module)[pos];
#pragma omp atomic update
              update += values[i * n + j];
            }
          }
          else if(indexJ[j].second == dofs::AssembleType::Linked) {
            if constexpr(scalar::IsComplex< T_Scalar >::value) {
              HackComplex< scalar::Precision< T_Scalar > > *hack = reinterpret_cast< HackComplex< scalar::Precision< T_Scalar > > * >(&(*_module)[pos]);
              const T_Scalar value = _valuesLC[indexJ[j].first] * values[i * n + j];
              const scalar::Precision< T_Scalar > re = std::real(value);
              const scalar::Precision< T_Scalar > im = std::imag(value);
#pragma omp atomic update
              hack->r += re;
#pragma omp atomic update
              hack->i += im;
            }
            else {
              T_Scalar &update = (*_module)[pos];
              const T_Scalar value = _valuesLC[indexJ[j].first] * values[i * n + j];
#pragma omp atomic update
              update += value;
            }
          }
          else if(indexJ[j].second == dofs::AssembleType::UnknownBubble) {
            (*_module)[pos] += values[i * n + j];
          }
          else if(indexJ[j].second == dofs::AssembleType::LinkedBubble) {
            (*_module)[pos] += _valuesLC[indexJ[j].first] * values[i * n + j];
          }
        }
      }
      else if(indexI[i].second == dofs::AssembleType::Linked) {
        for(auto j = 0ULL; j < n; ++j) {
          unsigned long long jIndex = 0;
          if(indexJ[j].second == dofs::AssembleType::Unknown || indexJ[j].second == dofs::AssembleType::UnknownBubble || indexJ[j].second == dofs::AssembleType::UnknownGlobal) {
            jIndex = indexJ[j].first;
          }
          else if(indexJ[j].second == dofs::AssembleType::Linked || indexJ[j].second == dofs::AssembleType::LinkedBubble) {
            jIndex = _indicesLC[indexJ[j].first];
          }
          else if(indexJ[j].second == dofs::AssembleType::Fixed) {
            haveDC = true;
            continue;
          }

          if(sym && jIndex < iIndex) continue; // upper triangular matrix

          const unsigned long long *place = std::lower_bound(&_aj[_ai[iIndex]], &_aj[_ai[iIndex + 1]], jIndex);
          const unsigned long long pos = place - _aj;
          if(indexJ[j].second == dofs::AssembleType::Unknown || indexJ[j].second == dofs::AssembleType::UnknownGlobal) {
            if constexpr(scalar::IsComplex< T_Scalar >::value) {
              HackComplex< scalar::Precision< T_Scalar > > *hack = reinterpret_cast< HackComplex< scalar::Precision< T_Scalar > > * >(&(*_module)[pos]);
              const T_Scalar value = std::conj(_valuesLC[indexI[i].first]) * values[i * n + j];
              const scalar::Precision< T_Scalar > re = std::real(value);
              const scalar::Precision< T_Scalar > im = std::imag(value);
#pragma omp atomic update
              hack->r += re;
#pragma omp atomic update
              hack->i += im;
            }
            else {
              T_Scalar &update = (*_module)[pos];
              const T_Scalar value = _valuesLC[indexI[i].first] * values[i * n + j];
#pragma omp atomic update
              update += value;
            }
          }
          else if(indexJ[j].second == dofs::AssembleType::Linked) {
            if constexpr(scalar::IsComplex< T_Scalar >::value) {
              HackComplex< scalar::Precision< T_Scalar > > *hack = reinterpret_cast< HackComplex< scalar::Precision< T_Scalar > > * >(&(*_module)[pos]);
              const T_Scalar value = std::conj(_valuesLC[indexI[i].first]) * _valuesLC[indexJ[j].first] * values[i * n + j];
              const scalar::Precision< T_Scalar > re = std::real(value);
              const scalar::Precision< T_Scalar > im = std::imag(value);
#pragma omp atomic update
              hack->r += re;
#pragma omp atomic update
              hack->i += im;
            }
            else {
              T_Scalar &update = (*_module)[pos];
              const T_Scalar value = _valuesLC[indexI[i].first] * _valuesLC[indexJ[j].first] * values[i * n + j];
#pragma omp atomic update
              update += value;
            }
          }
          else if(indexJ[j].second == dofs::AssembleType::UnknownBubble) {
            if constexpr(scalar::IsComplex< T_Scalar >::value) {
              (*_module)[pos] += std::conj(_valuesLC[indexI[i].first]) * values[i * n + j];
            }
            else {
              (*_module)[pos] += _valuesLC[indexI[i].first] * values[i * n + j];
            }
          }
          else if(indexJ[j].second == dofs::AssembleType::LinkedBubble) {
            if constexpr(scalar::IsComplex< T_Scalar >::value) {
              (*_module)[pos] += std::conj(_valuesLC[indexI[i].first]) * _valuesLC[indexJ[j].first] * values[i * n + j];
            }
            else {
              (*_module)[pos] += _valuesLC[indexI[i].first] * _valuesLC[indexJ[j].first] * values[i * n + j];
            }
          }
        }
      }
      else if(indexI[i].second == dofs::AssembleType::UnknownBubble) {
        for(auto j = 0ULL; j < n; ++j) {
          unsigned long long jIndex = 0;
          if(indexJ[j].second == dofs::AssembleType::Unknown || indexJ[j].second == dofs::AssembleType::UnknownBubble || indexJ[j].second == dofs::AssembleType::UnknownGlobal) {
            jIndex = indexJ[j].first;
          }
          else if(indexJ[j].second == dofs::AssembleType::Linked || indexJ[j].second == dofs::AssembleType::LinkedBubble) {
            jIndex = _indicesLC[indexJ[j].first];
          }
          else if(indexJ[j].second == dofs::AssembleType::Fixed) {
            haveDC = true;
            continue;
          }

          if(sym && jIndex < iIndex) continue; // upper triangular matrix

          const unsigned long long *place = std::lower_bound(&_aj[_ai[iIndex]], &_aj[_ai[iIndex + 1]], jIndex);
          const unsigned long long pos = place - _aj;
          if(indexJ[j].second == dofs::AssembleType::Unknown || indexJ[j].second == dofs::AssembleType::UnknownBubble || indexJ[j].second == dofs::AssembleType::UnknownGlobal) {
            (*_module)[pos] += values[i * n + j];
          }
          else if(indexJ[j].second == dofs::AssembleType::Linked || indexJ[j].second == dofs::AssembleType::LinkedBubble) {
            (*_module)[pos] += _valuesLC[indexJ[j].first] * values[i * n + j];
          }
        }
      }
      else if(indexI[i].second == dofs::AssembleType::LinkedBubble) {
        for(auto j = 0ULL; j < n; ++j) {
          unsigned long long jIndex = 0;
          if(indexJ[j].second == dofs::AssembleType::Unknown || indexJ[j].second == dofs::AssembleType::UnknownBubble || indexJ[j].second == dofs::AssembleType::UnknownGlobal) {
            jIndex = indexJ[j].first;
          }
          else if(indexJ[j].second == dofs::AssembleType::Linked || indexJ[j].second == dofs::AssembleType::LinkedBubble) {
            jIndex = _indicesLC[indexJ[j].first];
          }
          else if(indexJ[j].second == dofs::AssembleType::Fixed) {
            haveDC = true;
            continue;
          }

          if(sym && jIndex < iIndex) continue; // upper triangular matrix

          const unsigned long long *place = std::lower_bound(&_aj[_ai[iIndex]], &_aj[_ai[iIndex + 1]], jIndex);
          const unsigned long long pos = place - _aj;
          if(indexJ[j].second == dofs::AssembleType::Unknown || indexJ[j].second == dofs::AssembleType::UnknownBubble || indexJ[j].second == dofs::AssembleType::UnknownGlobal) {
            if constexpr(scalar::IsComplex< T_Scalar >::value) {
              (*_module)[pos] += std::conj(_valuesLC[indexI[i].first]) * values[i * n + j];
            }
            else {
              (*_module)[pos] += _valuesLC[indexI[i].first] * values[i * n + j];
            }
          }
          else if(indexJ[j].second == dofs::AssembleType::Linked || indexJ[j].second == dofs::AssembleType::LinkedBubble) {
            if constexpr(scalar::IsComplex< T_Scalar >::value) {
              (*_module)[pos] += std::conj(_valuesLC[indexI[i].first]) * _valuesLC[indexJ[j].first] * values[i * n + j];
            }
            else {
              (*_module)[pos] += _valuesLC[indexI[i].first] * _valuesLC[indexJ[j].first] * values[i * n + j];
            }
          }
        }
      }
    }

    return haveDC;
  }

  template< class T_Scalar >
  void MatrixFactory< T_Scalar >::addValue(const unsigned long long i, const unsigned long long j, const T_Scalar &values)
  {
    unsigned long long *place = std::lower_bound(&_aj[_ai[i]], &_aj[_ai[i + 1]], j);
    const unsigned long long pos = place - _aj;
    (*_module)[pos] += values;
  }

  template< class T_Scalar >
  void MatrixFactory< T_Scalar >::addPattern(const unsigned long long m, const unsigned long long n, const std::pair< unsigned long long, int > *const indexI, const std::pair< unsigned long long, int > *const indexJ)
  {
    const bool sym = _options.symmetric() || _options.hermitian();
    std::vector< unsigned long long > *bubble = nullptr;

    for(auto i = 0ULL; i < m; ++i) {
      unsigned long long iIndex = 0;
      if(indexI[i].second == dofs::AssembleType::Unknown || indexI[i].second == dofs::AssembleType::UnknownBubble || indexI[i].second == dofs::AssembleType::UnknownGlobal) {
        iIndex = indexI[i].first;
      }
      else if(indexI[i].second == dofs::AssembleType::Linked || indexI[i].second == dofs::AssembleType::LinkedBubble) {
        iIndex = _indicesLC[indexI[i].first];
      }
      else if(indexI[i].second == dofs::AssembleType::Fixed) {
        continue;
      }

      if(indexI[i].second == dofs::AssembleType::Unknown || indexI[i].second == dofs::AssembleType::UnknownGlobal || indexI[i].second == dofs::AssembleType::Linked) {
        _locks[iIndex].lock();
        _pattern[iIndex]->reserve(_pattern[iIndex]->size() + n);
        for(auto j = 0ULL; j < n; ++j) {
          unsigned long long jIndex = 0;
          if(indexJ[j].second == dofs::AssembleType::Unknown || indexJ[j].second == dofs::AssembleType::UnknownBubble || indexJ[j].second == dofs::AssembleType::UnknownGlobal) {
            jIndex = indexJ[j].first;
          }
          else if(indexJ[j].second == dofs::AssembleType::Linked || indexJ[j].second == dofs::AssembleType::LinkedBubble) {
            jIndex = _indicesLC[indexJ[j].first];
          }
          else if(indexJ[j].second == dofs::AssembleType::Fixed) {
            continue;
          }

          if(sym && jIndex < iIndex) continue; // upper triangular matrix

          const auto place = std::lower_bound(_pattern[iIndex]->begin(), _pattern[iIndex]->end(), jIndex);
          if(place == _pattern[iIndex]->end()) {
            _pattern[iIndex]->push_back(jIndex);
          }
          else {
            if(*place != jIndex) {
              _pattern[iIndex]->emplace(place, jIndex);
            }
          }
        }
        _locks[iIndex].unlock();
      }
      else if(indexI[i].second == dofs::AssembleType::UnknownBubble || indexI[i].second == dofs::AssembleType::LinkedBubble) {
        if(bubble == nullptr) {
          bubble = _pattern[iIndex];
          _pattern[iIndex]->reserve(_pattern[iIndex]->size() + n);
          for(auto j = 0ULL; j < n; ++j) {
            unsigned long long jIndex = 0;
            if(indexJ[j].second == dofs::AssembleType::Unknown || indexJ[j].second == dofs::AssembleType::UnknownBubble || indexJ[j].second == dofs::AssembleType::UnknownGlobal) {
              jIndex = indexJ[j].first;
            }
            else if(indexJ[j].second == dofs::AssembleType::Linked || indexJ[j].second == dofs::AssembleType::LinkedBubble) {
              jIndex = _indicesLC[indexJ[j].first];
            }
            else if(indexJ[j].second == dofs::AssembleType::Fixed) {
              continue;
            }

            if(sym && jIndex < iIndex) continue; // upper triangular matrix

            const auto place = std::lower_bound(_pattern[iIndex]->begin(), _pattern[iIndex]->end(), jIndex);
            if(place == _pattern[iIndex]->end()) {
              _pattern[iIndex]->push_back(jIndex);
            }
            else {
              if(*place != jIndex) {
                _pattern[iIndex]->emplace(place, jIndex);
              }
            }
          }
        }
        else {
          _pattern[iIndex] = bubble;
        }
      }
    }
  }

  template< class T_Scalar >
  void MatrixFactory< T_Scalar >::addPatternGlobalDof(const unsigned long long i)
  {
    std::pair< unsigned long long, int > indexI = std::make_pair(i, 4);
    addPattern(1, 1, &indexI, &indexI);
  }

  template< class T_Scalar >
  bool MatrixFactory< T_Scalar >::havePattern() const
  {
    return _havePattern;
  }

#ifdef HAVE_PETSC
  template< class T_Scalar >
  Mat MatrixFactory< T_Scalar >::getPetsc() const
  {
    PetscInterface< T_Scalar, PetscScalar > interface;
    Mat matPetsc;
    interface(_size, _size, _ai, _aj, &(*_module)[0], &matPetsc, _shouldBeDestroyedWithPetsc, _options.symmetric() || _options.hermitian());

    if(_options.symmetric()) {
      MatSetOption(matPetsc, MAT_SYMMETRIC, PETSC_TRUE);
      MatSetOption(matPetsc, MAT_SYMMETRY_ETERNAL, PETSC_TRUE);
    }
    if(_options.hermitian()) {
      MatSetOption(matPetsc, MAT_HERMITIAN, PETSC_TRUE);
      MatSetOption(matPetsc, MAT_SYMMETRY_ETERNAL, PETSC_TRUE);
    }
    return matPetsc;
  }
#endif

  template< class T_Scalar >
  common::Memory MatrixFactory< T_Scalar >::memory() const
  {
    common::Memory mem = (_ai[_size] + _size + 1) * sizeof(unsigned long long);
    if(_module) {
      mem += _module->memory();
    }
    return mem;
  }

  template< class T_Scalar >
  unsigned long long MatrixFactory< T_Scalar >::size() const
  {
    return _size;
  }

  template< class T_Scalar >
  unsigned long long MatrixFactory< T_Scalar >::numberOfNonZeros() const
  {
    if(!_havePattern) {
      return 0;
    }
    return _ai[_size];
  }

  template< class T_Scalar >
  system::MatrixOptions MatrixFactory< T_Scalar >::getOptions() const
  {
    return _options;
  }

  INSTANTIATE_CLASS(MatrixFactory, 4, TEMPLATE_ARGS(std::complex< double >, std::complex< float >, double, float))


} // namespace gmshfem::system
