// GmshFEM - Copyright (C) 2019-2022, A. Royer, E. Béchet, C. Geuzaine, Université de Liège
//
// See the LICENSE.txt file for license information. Please report all
// issues on https://gitlab.onelab.info/gmsh/fem/issues

#include "MatrixCCS.h"

#include "KahanSum.h"
#include "MatrixFactory.h"
#include "OmpInterface.h"
#include "PetscInterface.h"
#include "instantiate.h"
#include "io.h"

namespace gmshfem::algebra
{


  template< class T_Scalar >
  void MatrixCCS< T_Scalar >::_buildPetsc()
  {
    throw common::Exception("PETSc matrix cannot be built from 'MatrixCCS': use a 'MatrixCRS' instead");
  }

  template< class T_Scalar >
  MatrixCCS< T_Scalar >::MatrixCCS() :
    Matrix< T_Scalar >(), _ai(), _aj(), _a()
  {
  }

  template< class T_Scalar >
  MatrixCCS< T_Scalar >::MatrixCCS(const MatrixCCS &other) :
    Matrix< T_Scalar >(other), _ai(other._ai), _aj(other._aj), _a(other._a)
  {
  }

  template< class T_Scalar >
  MatrixCCS< T_Scalar >::MatrixCCS(MatrixCCS &&other) :
    Matrix< T_Scalar >(other), _ai(std::move(other._ai)), _aj(std::move(other._aj)), _a(std::move(other._a))
  {
  }

  template< class T_Scalar >
  MatrixCCS< T_Scalar >::MatrixCCS(std::vector< unsigned long long > &&ai, std::vector< unsigned long long > &&aj, std::vector< T_Scalar > &&a) :
    Matrix< T_Scalar >(*std::max_element(ai.begin(), ai.end()) + 1, aj.size() - 1), _ai(std::move(ai)), _aj(std::move(aj)), _a(std::move(a))
  {
  }

  template< class T_Scalar >
  MatrixCCS< T_Scalar >::~MatrixCCS()
  {
  }

  template< class T_Scalar >
  Eigen::SparseMatrix< T_Scalar, Eigen::ColMajor, long long > MatrixCCS< T_Scalar >::getEigen()
  {
    return Eigen::Map< Eigen::SparseMatrix< T_Scalar, Eigen::ColMajor, long long > >(this->_size[0], this->_size[1], numberOfNonZero(), reinterpret_cast< long long * >(_aj.data()), reinterpret_cast< long long * >(_ai.data()), _a.data());
  }

  template< class T_Scalar >
  MatrixCCS< T_Scalar > &MatrixCCS< T_Scalar >::operator=(const MatrixCCS &other)
  {
    this->copySize(other);
    _ai = other._ai;
    _aj = other._aj;
    _a = other._a;
    return *this;
  }

  template< class T_Scalar >
  MatrixCCS< T_Scalar > &MatrixCCS< T_Scalar >::operator=(MatrixCCS &&other)
  {
    this->copySize(other);
    _ai = std::move(other._ai);
    _aj = std::move(other._aj);
    _a = std::move(other._a);
    return *this;
  }

  template< class T_Scalar >
  MatrixCCS< T_Scalar > &MatrixCCS< T_Scalar >::operator=(const Matrix< T_Scalar > &other)
  {
    clear();
    switch(other.format()) {
    case MatrixFormat::CRSFast: {
      const MatrixCRSFast< T_Scalar > &mat = static_cast< const MatrixCRSFast< T_Scalar > & >(other);
      this->copySize(other);

      _aj.resize(this->_size[1] + 1, 0);
#pragma omp parallel for num_threads(omp::getMaxThreads())
      for(auto i = 0ULL; i < mat.ai()[this->_size[0]]; ++i) {
#pragma omp atomic update
        _aj[mat.aj()[i]]++;
      }

      unsigned long long tmp = _aj[0];
      _aj[0] = 0;
      for(auto i = 1ULL; i < _aj.size(); ++i) {
        const unsigned long long tmp2 = _aj[i];
        _aj[i] = _aj[i - 1] + tmp;
        tmp = tmp2;
      }

      _ai.resize(_aj[this->_size[1]], 0.);
      _a.resize(_aj[this->_size[1]], 0.);

      for(auto i = 0ULL; i < this->_size[0]; ++i) {
        for(auto j = mat.ai()[i]; j < mat.ai()[i + 1]; ++j) {
          _ai[_aj[mat.aj()[j]]] = i;
          _a[_aj[mat.aj()[j]]] = mat.a()[j];
          _aj[mat.aj()[j]]++;
        }
      }

      for(auto i = _aj.size() - 1; i > 0; --i) {
        _aj[i] = _aj[i - 1];
      }
      _aj[0] = 0;
    } break;
    case MatrixFormat::CRS: {
      MatrixCRSFast< T_Scalar > mat;
      mat = other;
      *this = mat;
    } break;
    default:
      throw common::Exception("Matrix conversion not implemented");
      break;
    }
    return *this;
  }

  template< class T_Scalar >
  void MatrixCCS< T_Scalar >::extract(const system::MatrixFactory< T_Scalar > *matrix)
  {
    if(matrix->getModule() == nullptr) {
      return;
    }
    MatrixCRSFast< T_Scalar > mat;
    mat.extract(matrix);
    *this = mat;
  }

  template< class T_Scalar >
  void MatrixCCS< T_Scalar >::clear()
  {
    this->_size[0] = this->_size[1] = 0;
    _ai.clear();
    _aj.clear();
    _a.clear();
  }

  template< class T_Scalar >
  MatrixFormat MatrixCCS< T_Scalar >::format() const
  {
    return MatrixFormat::CCS;
  }

  template< class T_Scalar >
  unsigned long long MatrixCCS< T_Scalar >::numberOfNonZero() const
  {
    return _aj.back();
  }

  template< class T_Scalar >
  bool MatrixCCS< T_Scalar >::isSymmetric(const scalar::Precision< T_Scalar > tolerance) const
  {
    if(!this->isSquare()) return false;
    for(auto j = 0ULL; j < this->_size[1]; ++j) {
      for(auto i = _aj[j]; i < _aj[j + 1]; ++i) {
        const unsigned long long *place = std::lower_bound(&_ai[_aj[_ai[i]]], &_ai[_aj[_ai[i] + 1]], j);
        const std::size_t pos = place - &_ai[0];
        if(std::abs(_a[pos] - _a[i]) > tolerance * std::abs(_a[pos])) {
          return false;
        }
      }
    }

    return true;
  }

  template< class T_Scalar >
  bool MatrixCCS< T_Scalar >::isHermitian(const scalar::Precision< T_Scalar > tolerance) const
  {
    if(!this->isSquare()) return false;
    if(scalar::IsComplex< T_Scalar >::value) {
      for(auto j = 0ULL; j < this->_size[1]; ++j) {
        for(auto i = _aj[j]; i < _aj[j + 1]; ++i) {
          const unsigned long long *place = std::lower_bound(&_ai[_aj[_ai[i]]], &_ai[_aj[_ai[i] + 1]], j);
          const std::size_t pos = place - &_ai[0];
          if(std::abs(_a[pos] - std::conj(_a[i])) > tolerance * std::abs(_a[pos])) {
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
  void MatrixCCS< T_Scalar >::save(const std::string &path) const
  {
    common::CSVio file(path);

    if(!file.isOpen()) {
      throw common::Exception("Cannot open file " + path + ".csv");
    }

    file << csv::precision(scalar::PrecisionDigits< T_Scalar >::value) << csv::scientific;
    for(auto j = 0ULL; j < this->_size[1]; ++j) {
      for(auto i = _aj[j]; i < _aj[j + 1]; ++i) {
        file << _ai[i] << j << _a[i] << csv::endl;
      }
    }

    file.close();
  }

  template< class T_Scalar >
  void MatrixCCS< T_Scalar >::saveSpyPlot(const std::string &path, const unsigned int pointSize, const common::Color &zero, const common::Color &noZero) const
  {
    std::vector< common::Color > line(this->_size[1]); // A line of the picture
    common::PPMio file(path, common::PPMFormat::PPM);
    file.writeHeader(this->_size[1], this->_size[0]);

    for(auto i = 0ULL; i < this->_size[0]; ++i) {
      for(auto j = 0ULL; j < this->_size[1]; ++j) {
        line[j] = zero;
      }

      const unsigned long long iMin = pointSize - 1 > i ? 0 : i - pointSize + 1;
      const unsigned long long iMax = i + pointSize > this->_size[0] ? this->_size[0] : i + pointSize;

      for(auto j = 0ULL; j < this->_size[1]; ++j) {
        const unsigned long long jMin = pointSize - 1 > j ? 0 : j - pointSize + 1;
        const unsigned long long jMax = j + pointSize > this->_size[1] ? this->_size[1] : j + pointSize;

        for(auto ii = _aj[j]; ii < _aj[j + 1]; ++ii) {
          if(_ai[ii] >= iMin && _ai[ii] < iMax) {
            if(_a[ii] != T_Scalar(0.)) {
              const unsigned long long I = _ai[ii];
              for(auto J = jMin; J < jMax; ++J) {
                if(J * J + j * j + I * I + i * i + pointSize <= pointSize * pointSize + 2 * J * j + 2 * I * i + 1) {
                  line[J] = noZero;
                }
              }
            }
          }
          else if(_ai[ii] > iMax) {
            break;
          }
        }
      }
      file.writeLine(line);
    }
    file.close();
  }

  template< class T_Scalar >
  const std::vector< unsigned long long > &MatrixCCS< T_Scalar >::ai() const
  {
    return _ai;
  }

  template< class T_Scalar >
  const std::vector< unsigned long long > &MatrixCCS< T_Scalar >::aj() const
  {
    return _aj;
  }

  template< class T_Scalar >
  const std::vector< T_Scalar > &MatrixCCS< T_Scalar >::a() const
  {
    return _a;
  }

  INSTANTIATE_CLASS(MatrixCCS, 4, TEMPLATE_ARGS(std::complex< double >, std::complex< float >, double, float))


} // namespace gmshfem::algebra
