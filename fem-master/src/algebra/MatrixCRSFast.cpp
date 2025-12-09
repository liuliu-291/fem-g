// GmshFEM - Copyright (C) 2019-2022, A. Royer, E. Béchet, C. Geuzaine, Université de Liège
//
// See the LICENSE.txt file for license information. Please report all
// issues on https://gitlab.onelab.info/gmsh/fem/issues

#include "MatrixCRSFast.h"

#include "KahanSum.h"
#include "MatrixFactory.h"
#include "OmpInterface.h"
#include "PetscInterface.h"
#include "instantiate.h"
#include "io.h"

namespace gmshfem::algebra
{


  template< class T_Scalar >
  void MatrixCRSFast< T_Scalar >::_buildPetsc()
  {
    throw common::Exception("PETSc matrix cannot be built from 'MatrixCRSFast': use a 'MatrixCRS' instead.");
  }

  template< class T_Scalar >
  MatrixCRSFast< T_Scalar >::MatrixCRSFast() :
    Matrix< T_Scalar >(), _ai(nullptr), _aj(nullptr), _a(nullptr)
  {
  }

  template< class T_Scalar >
  MatrixCRSFast< T_Scalar >::MatrixCRSFast(const unsigned long long size0, const unsigned long long size1, const std::vector< unsigned long long > &ai, const std::vector< unsigned long long > &aj, const std::vector< T_Scalar > &a) :
    Matrix< T_Scalar >(size0, size1), _ai(&ai[0]), _aj(&aj[0]), _a(&a[0])
  {
  }

  template< class T_Scalar >
  MatrixCRSFast< T_Scalar >::MatrixCRSFast(const MatrixCRSFast &other) :
    Matrix< T_Scalar >(other), _ai(other._ai), _aj(other._aj), _a(other._a)
  {
  }

  template< class T_Scalar >
  MatrixCRSFast< T_Scalar >::MatrixCRSFast(MatrixCRSFast &&other) :
    Matrix< T_Scalar >(other), _ai(std::move(other._ai)), _aj(std::move(other._aj)), _a(std::move(other._a))
  {
  }

  template< class T_Scalar >
  MatrixCRSFast< T_Scalar >::~MatrixCRSFast()
  {
  }

  template< class T_Scalar >
  const Eigen::SparseMatrix< T_Scalar, Eigen::RowMajor, long long > MatrixCRSFast< T_Scalar >::getEigen() const
  {
    return Eigen::Map< const Eigen::SparseMatrix< T_Scalar, Eigen::RowMajor, long long > >(this->_size[0], this->_size[1], numberOfNonZero(), reinterpret_cast< const long long * >(_ai), reinterpret_cast< const long long * >(_aj), _a);
  }

  template< class T_Scalar >
  MatrixCRSFast< T_Scalar > &MatrixCRSFast< T_Scalar >::operator=(const MatrixCRSFast &other)
  {
    this->copySize(other);
    _ai = other._ai;
    _aj = other._aj;
    _a = other._a;
    return *this;
  }

  template< class T_Scalar >
  MatrixCRSFast< T_Scalar > &MatrixCRSFast< T_Scalar >::operator=(MatrixCRSFast &&other)
  {
    this->copySize(other);
    _ai = other._ai;
    _aj = other._aj;
    _a = other._a;
    other._ai = nullptr;
    other._aj = nullptr;
    other._a = nullptr;
    return *this;
  }

  template< class T_Scalar >
  MatrixCRSFast< T_Scalar > &MatrixCRSFast< T_Scalar >::operator=(const Matrix< T_Scalar > &other)
  {
    clear();
    switch(other.format()) {
    case MatrixFormat::CRS: {
      const MatrixCRS< T_Scalar > &mat = static_cast< const MatrixCRS< T_Scalar > & >(other);
      this->copySize(other);
      _ai = &mat.ai()[0];
      _aj = &mat.aj()[0];
      _a = &mat.a()[0];
    } break;
    default:
      common::Exception("Conversion of matrix not implement.");
      break;
    }
    return *this;
  }

  template< class T_Scalar >
  void MatrixCRSFast< T_Scalar >::extract(const system::MatrixFactory< T_Scalar > *matrix)
  {
    if(matrix->getModule() == nullptr) {
      return;
    }
    this->_size[0] = this->_size[1] = matrix->size();
    _ai = matrix->_ai;
    _aj = matrix->_aj;
    _a = matrix->getModule()->getMatrix();
  }

  template< class T_Scalar >
  void MatrixCRSFast< T_Scalar >::clear()
  {
    this->_size[0] = this->_size[1] = 0;
    _ai = nullptr;
    _aj = nullptr;
    _a = nullptr;
  }

  template< class T_Scalar >
  MatrixFormat MatrixCRSFast< T_Scalar >::format() const
  {
    return MatrixFormat::CRSFast;
  }

  template< class T_Scalar >
  unsigned long long MatrixCRSFast< T_Scalar >::numberOfNonZero() const
  {
    if(_ai == nullptr) {
      return 0;
    }
    return _ai[this->_size[0]];
  }

  template< class T_Scalar >
  bool MatrixCRSFast< T_Scalar >::isSymmetric(const scalar::Precision< T_Scalar > tolerance) const
  {
    if(_ai == nullptr || _aj == nullptr || _a == nullptr) {
      return false;
    }

    if(!this->isSquare()) return false;

    for(auto i = 0ULL; i < this->_size[0]; ++i) {
      for(auto j = _ai[i]; j < _ai[i + 1]; ++j) {
        const unsigned long long *place = std::lower_bound(&_aj[_ai[_aj[j]]], &_aj[_ai[_aj[j] + 1]], i);
        const std::size_t pos = place - &_aj[0];
        if(std::abs(_a[pos] - _a[j]) > tolerance * std::abs(_a[pos])) {
          return false;
        }
      }
    }

    return true;
  }

  template< class T_Scalar >
  bool MatrixCRSFast< T_Scalar >::isHermitian(const scalar::Precision< T_Scalar > tolerance) const
  {
    if(_ai == nullptr || _aj == nullptr || _a == nullptr) {
      return false;
    }

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
  void MatrixCRSFast< T_Scalar >::save(const std::string &path) const
  {
    if(_ai == nullptr || _aj == nullptr || _a == nullptr) {
      return;
    }

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
  void MatrixCRSFast< T_Scalar >::saveSpyPlot(const std::string &path, const unsigned int pointSize, const common::Color &zero, const common::Color &noZero) const
  {
    if(_ai == nullptr || _aj == nullptr || _a == nullptr) {
      return;
    }

    std::vector< common::Color > line(this->_size[1]); // A line of the picture
    common::PPMio file(path, common::PPMFormat::PPM);
    file.writeHeader(this->_size[1], this->_size[0]);

    for(auto i = 0ULL; i < this->_size[0]; ++i) {
      for(auto j = 0ULL; j < this->_size[1]; ++j) {
        line[j] = zero;
      }
      const unsigned long long iMin = pointSize - 1 > i ? 0 : i - pointSize + 1;
      const unsigned long long iMax = i + pointSize > this->_size[0] ? this->_size[0] : i + pointSize;

      for(auto I = iMin; I < iMax; ++I) {
        for(auto jj = _ai[I]; jj < _ai[I + 1]; ++jj) {
          if(_a[jj] != T_Scalar(0.)) {
            const unsigned long long j = _aj[jj];
            const unsigned long long jMin = pointSize - 1 > j ? 0 : j - pointSize + 1;
            const unsigned long long jMax = j + pointSize > this->_size[1] ? this->_size[1] : j + pointSize;
            for(auto J = jMin; J < jMax; ++J) {
              if(J * J + j * j + I * I + i * i + pointSize <= pointSize * pointSize + 2 * J * j + 2 * I * i + 1) {
                line[J] = noZero;
              }
            }
          }
        }
      }
      file.writeLine(line);
    }
    file.close();
  }

  template< class T_Scalar >
  const unsigned long long *MatrixCRSFast< T_Scalar >::ai() const
  {
    return _ai;
  }

  template< class T_Scalar >
  const unsigned long long *MatrixCRSFast< T_Scalar >::aj() const
  {
    return _aj;
  }

  template< class T_Scalar >
  const T_Scalar *MatrixCRSFast< T_Scalar >::a() const
  {
    return _a;
  }

  INSTANTIATE_CLASS(MatrixCRSFast, 4, TEMPLATE_ARGS(std::complex< double >, std::complex< float >, double, float))


} // namespace gmshfem::algebra
