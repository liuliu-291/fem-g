// GmshFEM - Copyright (C) 2019-2022, A. Royer, E. Béchet, C. Geuzaine, Université de Liège
//
// See the LICENSE.txt file for license information. Please report all
// issues on https://gitlab.onelab.info/gmsh/fem/issues

#include "MatrixCRS.h"

#include "KahanSum.h"
#include "MatrixFactory.h"
#include "Message.h"
#include "OmpInterface.h"
#include "PetscInterface.h"
#include "instantiate.h"
#include "io.h"

namespace gmshfem::algebra
{


  template< class T_Scalar >
  void MatrixCRS< T_Scalar >::_buildPetsc()
  {
#ifdef HAVE_PETSC
    system::PetscInterface< T_Scalar, PetscScalar > interface;
    interface(this->_size[0], this->_size[1], &_ai[0], &_aj[0], &_a[0], &this->_matPetsc, this->_shouldBeDestroyedWithPetsc, false);
#endif
  }

  template< class T_Scalar >
  MatrixCRS< T_Scalar >::MatrixCRS() :
    Matrix< T_Scalar >(), _ai(), _aj(), _a()
  {
  }

  template< class T_Scalar >
  MatrixCRS< T_Scalar >::MatrixCRS(const MatrixCRS &other) :
    Matrix< T_Scalar >(other), _ai(other._ai), _aj(other._aj), _a(other._a)
  {
  }

  template< class T_Scalar >
  MatrixCRS< T_Scalar >::MatrixCRS(MatrixCRS &&other) :
    Matrix< T_Scalar >(other), _ai(std::move(other._ai)), _aj(std::move(other._aj)), _a(std::move(other._a))
  {
  }

  template< class T_Scalar >
  MatrixCRS< T_Scalar >::MatrixCRS(std::vector< unsigned long long > &&ai, std::vector< unsigned long long > &&aj, std::vector< T_Scalar > &&a) :
    Matrix< T_Scalar >(ai.size() - 1, *std::max_element(aj.begin(), aj.end()) + 1), _ai(std::move(ai)), _aj(std::move(aj)), _a(std::move(a))
  {
  }

  template< class T_Scalar >
  MatrixCRS< T_Scalar >::~MatrixCRS()
  {
  }

  template< class T_Scalar >
  Eigen::SparseMatrix< T_Scalar, Eigen::RowMajor, long long > MatrixCRS< T_Scalar >::getEigen()
  {
    return Eigen::Map< Eigen::SparseMatrix< T_Scalar, Eigen::RowMajor, long long > >(this->_size[0], this->_size[1], numberOfNonZero(), reinterpret_cast< long long * >(_ai.data()), reinterpret_cast< long long * >(_aj.data()), _a.data());
  }

  template< class T_Scalar >
  MatrixCRS< T_Scalar > &MatrixCRS< T_Scalar >::operator=(const MatrixCRS &other)
  {
    this->copySize(other);
    _ai = other._ai;
    _aj = other._aj;
    _a = other._a;
    return *this;
  }

  template< class T_Scalar >
  MatrixCRS< T_Scalar > &MatrixCRS< T_Scalar >::operator=(MatrixCRS &&other)
  {
    this->copySize(other);
    _ai = std::move(other._ai);
    _aj = std::move(other._aj);
    _a = std::move(other._a);
    return *this;
  }

  template< class T_Scalar >
  MatrixCRS< T_Scalar > &MatrixCRS< T_Scalar >::operator=(const Matrix< T_Scalar > &other)
  {
    clear();
    switch(other.format()) {
    case MatrixFormat::CRSFast: {
      const MatrixCRSFast< T_Scalar > &mat = static_cast< const MatrixCRSFast< T_Scalar > & >(other);
      this->copySize(other);
      _ai = std::vector< unsigned long long >(mat.ai(), mat.ai() + this->_size[0] + 1);
      _aj = std::vector< unsigned long long >(mat.aj(), mat.aj() + _ai.back());
      _a = std::vector< T_Scalar >(mat.a(), mat.a() + _ai.back());
    } break;
    case MatrixFormat::CCS: {
      const MatrixCCS< T_Scalar > &mat = static_cast< const MatrixCCS< T_Scalar > & >(other);
      this->copySize(other);

      _ai.resize(this->_size[0] + 1, 0);
#pragma omp parallel for num_threads(omp::getMaxThreads())
      for(auto i = 0ULL; i < mat.aj()[this->_size[1]]; ++i) {
#pragma omp atomic update
        _ai[mat.ai()[i]]++;
      }

      unsigned long long tmp = _ai[0];
      _ai[0] = 0;
      for(auto i = 1ULL; i < _ai.size(); ++i) {
        const unsigned long long tmp2 = _ai[i];
        _ai[i] = _ai[i - 1] + tmp;
        tmp = tmp2;
      }

      _aj.resize(_ai[this->_size[0]], 0.);
      _a.resize(_ai[this->_size[0]], 0.);

      for(auto j = 0ULL; j < this->_size[1]; ++j) {
        for(auto i = mat.aj()[j]; i < mat.aj()[j + 1]; ++i) {
          _aj[_ai[mat.ai()[i]]] = j;
          _a[_ai[mat.ai()[i]]] = mat.a()[i];
          _ai[mat.ai()[i]]++;
        }
      }

      for(auto i = _ai.size() - 1; i > 0; --i) {
        _ai[i] = _ai[i - 1];
      }
      _ai[0] = 0;
    } break;
    default:
      common::Exception("Conversion of matrix not implement.");
      break;
    }
    return *this;
  }

  template< class T_Scalar >
  void MatrixCRS< T_Scalar >::extract(const system::MatrixFactory< T_Scalar > *matrix)
  {
    if(matrix->getModule() == nullptr) {
      return;
    }
    this->_size[0] = this->_size[1] = matrix->size();
    _ai = std::vector< unsigned long long >(matrix->_ai, matrix->_ai + this->_size[0] + 1);
    _aj = std::vector< unsigned long long >(matrix->_aj, matrix->_aj + _ai.back());
    _a = std::vector< T_Scalar >(matrix->getModule()->getMatrix(), matrix->getModule()->getMatrix() + _ai.back());
  }

  template< class T_Scalar >
  void MatrixCRS< T_Scalar >::clear()
  {
    this->_size[0] = this->_size[1] = 0;
    _ai.clear();
    _aj.clear();
    _a.clear();
  }

  template< class T_Scalar >
  MatrixFormat MatrixCRS< T_Scalar >::format() const
  {
    return MatrixFormat::CRS;
  }

  template< class T_Scalar >
  unsigned long long MatrixCRS< T_Scalar >::numberOfNonZero() const
  {
    return _ai.back();
  }

  template< class T_Scalar >
  bool MatrixCRS< T_Scalar >::isSymmetric(const scalar::Precision< T_Scalar > tolerance) const
  {
    if(!this->isSquare()) return false;

    for(auto i = 0ULL; i < this->_size[0]; ++i) {
      for(auto j = _ai[i]; j < _ai[i + 1]; ++j) {
        const unsigned long long *place = std::lower_bound(&_aj[_ai[_aj[j]]], &_aj[_ai[_aj[j] + 1]], i);
        const std::size_t pos = place - &_aj[0];
        if(std::abs(_a[pos] - _a[j]) > tolerance * std::abs(_a[pos])) {
          msg::error << _a[pos] << " " << _a[j] << msg::endl;
          return false;
        }
      }
    }

    return true;
  }

  template< class T_Scalar >
  bool MatrixCRS< T_Scalar >::isHermitian(const scalar::Precision< T_Scalar > tolerance) const
  {
    if(!this->isSquare()) return false;

    if(scalar::IsComplex< T_Scalar >::value) {
      for(auto i = 0ULL; i < this->_size[0]; ++i) {
        for(auto j = _ai[i]; j < _ai[i + 1]; ++j) {
          const unsigned long long *place = std::lower_bound(&_aj[_ai[_aj[j]]], &_aj[_ai[_aj[j] + 1]], i);
          const std::size_t pos = place - &_aj[0];
          if(std::abs(_a[pos] - std::conj(_a[j])) > tolerance * std::abs(_a[pos])) {
            return false;
          }
        }
      }
      return true;
    }
    else {
      return isSymmetric(tolerance);
    }

    return false;
  }

  template< class T_Scalar >
  void MatrixCRS< T_Scalar >::save(const std::string &path) const
  {
    common::CSVio file(path);

    if(!file.isOpen()) {
      throw common::Exception("Cannot open file " + path + ".csv");
    }

    file << csv::precision(scalar::PrecisionDigits< T_Scalar >::value) << csv::scientific;
    for(auto i = 0ULL; i < this->_size[0]; ++i) {
      for(auto j = _ai[i]; j < _ai[i + 1]; ++j) {
        file << i << _aj[j] << _a[j] << csv::endl;
      }
    }

    file.close();
  }

  template< class T_Scalar >
  void MatrixCRS< T_Scalar >::saveSpyPlot(const std::string &path, const unsigned int pointSize, const common::Color &zero, const common::Color &noZero) const
  {
    MatrixCRSFast< T_Scalar > mat;
    mat = *this;
    mat.saveSpyPlot(path, pointSize, zero, noZero);
  }

  template< class T_Scalar >
  const std::vector< unsigned long long > &MatrixCRS< T_Scalar >::ai() const
  {
    return _ai;
  }

  template< class T_Scalar >
  const std::vector< unsigned long long > &MatrixCRS< T_Scalar >::aj() const
  {
    return _aj;
  }

  template< class T_Scalar >
  const std::vector< T_Scalar > &MatrixCRS< T_Scalar >::a() const
  {
    return _a;
  }

  INSTANTIATE_CLASS(MatrixCRS, 4, TEMPLATE_ARGS(std::complex< double >, std::complex< float >, double, float))


} // namespace gmshfem::algebra
