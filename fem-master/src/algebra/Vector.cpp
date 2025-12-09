// GmshFEM - Copyright (C) 2019-2022, A. Royer, E. Béchet, C. Geuzaine, Université de Liège
//
// See the LICENSE.txt file for license information. Please report all
// issues on https://gitlab.onelab.info/gmsh/fem/issues

#include <numeric>

#include "Vector.h"

#include "OmpInterface.h"
#include "PetscInterface.h"
#include "VectorFactory.h"
#include "instantiate.h"
#include "io.h"

namespace gmshfem::algebra
{


  template< class T_Scalar >
  void Vector< T_Scalar >::_buildPetsc()
  {
#ifdef HAVE_PETSC
    system::PetscInterface< T_Scalar, PetscScalar > interface;
    interface(this->_size[0], &_vector[0], &this->_vecPetsc, this->_shouldBeDestroyedWithPetsc);
#endif
  }

  template< class T_Scalar >
  Vector< T_Scalar >::Vector() :
    AlgebraicObject< T_Scalar >(0, 1), _vector(), _havePetsc(false)
  {
  }

  template< class T_Scalar >
  Vector< T_Scalar >::Vector(const unsigned long long size) :
    AlgebraicObject< T_Scalar >(size, 1), _vector(size, 0.), _havePetsc(false)
  {
  }

  template< class T_Scalar >
  Vector< T_Scalar >::Vector(const Vector &other) :
    AlgebraicObject< T_Scalar >(other), _vector(other._vector), _havePetsc(false)
  {
  }

  template< class T_Scalar >
  Vector< T_Scalar >::Vector(Vector &&other) :
    AlgebraicObject< T_Scalar >(other), _vector(std::move(other._vector)), _havePetsc(false)
  {
  }

  template< class T_Scalar >
  Vector< T_Scalar >::Vector(std::vector< T_Scalar > &&vector) :
    AlgebraicObject< T_Scalar >(vector.size(), 1), _vector(std::move(vector)), _havePetsc(false)
  {
  }

  template< class T_Scalar >
  Vector< T_Scalar >::~Vector()
  {
    removePetsc();
  }

  template< class T_Scalar >
  Vector< T_Scalar > &Vector< T_Scalar >::operator=(std::vector< T_Scalar > &&vector)
  {
    this->_size[0] = vector.size();
    _vector = std::move(vector);
    _havePetsc = false;
    return *this;
  }

  template< class T_Scalar >
  Vector< T_Scalar > &Vector< T_Scalar >::operator=(const Vector &other)
  {
    this->_size[0] = other._size[0];
    _vector = other._vector;
    _havePetsc = false;
    return *this;
  }

  template< class T_Scalar >
  Vector< T_Scalar > &Vector< T_Scalar >::operator=(Vector &&other)
  {
    this->_size[0] = other._size[0];
    other._size[0] = 0;
    _vector = std::move(other._vector);
    _havePetsc = false;
    return *this;
  }

  template< class T_Scalar >
  Vector< T_Scalar > &Vector< T_Scalar >::operator+=(const Vector &other)
  {
    eigenop(*this, other, "addition",
            &Eigen::Map< Eigen::VectorX< T_Scalar > >::operator+=);
    return *this;
  }

  template< class T_Scalar >
  Vector< T_Scalar > &Vector< T_Scalar >::operator-=(const Vector &other)
  {
    eigenop(*this, other, "subtraction",
            &Eigen::Map< Eigen::VectorX< T_Scalar > >::operator-=);
    return *this;
  }

  template< class T_Scalar >
  Vector< T_Scalar > &Vector< T_Scalar >::operator*=(const T_Scalar &alpha)
  {
    eigenop(*this, alpha, "multiplication with a scalar",
            &Eigen::Map< Eigen::VectorX< T_Scalar > >::operator*=);
    return *this;
  }

  template< class T_Scalar >
  Vector< T_Scalar > &Vector< T_Scalar >::operator/=(const T_Scalar &alpha)
  {
    eigenop(*this, alpha, "division with a scalar",
            &Eigen::Map< Eigen::VectorX< T_Scalar > >::operator/=);
    return *this;
  }

  template< class T_Scalar >
  Vector< T_Scalar > Vector< T_Scalar >::operator+(const Vector &other) const
  {
    Vector< T_Scalar > copy(*this);
    copy += other;
    return copy;
  }

  template< class T_Scalar >
  Vector< T_Scalar > Vector< T_Scalar >::operator-(const Vector &other) const
  {
    Vector< T_Scalar > copy(*this);
    copy -= other;
    return copy;
  }

  template< class T_Scalar >
  Vector< T_Scalar > Vector< T_Scalar >::operator*(const T_Scalar &alpha) const
  {
    Vector< T_Scalar > copy(*this);
    copy *= alpha;
    return copy;
  }

  template< class T_Scalar >
  Vector< T_Scalar > Vector< T_Scalar >::operator/(const T_Scalar &alpha) const
  {
    Vector< T_Scalar > copy(*this);
    copy /= alpha;
    return copy;
  }

  template< class T_Scalar >
  T_Scalar Vector< T_Scalar >::operator*(const Vector &other) const{
    return eigenop(*this, other, "scalar product",
                   &Eigen::Map<const Eigen::VectorX< T_Scalar > >::dot);
  }

  template< class T_Scalar >
  void Vector< T_Scalar >::extract(const system::VectorFactory< T_Scalar > *vector)
  {
    this->_size[0] = vector->_vec.size();
    _vector = vector->_vec;
  }

  template< class T_Scalar >
  T_Scalar &Vector< T_Scalar >::operator[](const unsigned long long pos)
  {
    if(pos >= this->_size[0]) {
      throw common::Exception("Trying to access element '" + std::to_string(pos) + "' of a vector of size '" + std::to_string(this->_size[0]) + "'");
    }
    return _vector[pos];
  }

  template< class T_Scalar >
  const T_Scalar &Vector< T_Scalar >::operator[](const unsigned long long pos) const
  {
    if(pos >= this->_size[0]) {
      throw common::Exception("Trying to acces element '" + std::to_string(pos) + "' of a vector of size '" + std::to_string(this->_size[0]) + "'");
    }
    return _vector[pos];
  }

  template< class T_Scalar >
  const std::vector< T_Scalar > &Vector< T_Scalar >::getStdVector() const
  {
    return _vector;
  }

  template< class T_Scalar >
  Vec Vector< T_Scalar >::getPetsc()
  {
#ifdef HAVE_PETSC
    if(_havePetsc) {
      return _vecPetsc;
    }
    _buildPetsc();
    _havePetsc = true;
    return _vecPetsc;
#else
    throw common::Exception("GmshFEM is not compiled with PETSc");
#endif
  }

  template< class T_Scalar >
  void Vector< T_Scalar >::removePetsc()
  {
#ifdef HAVE_PETSC
    if(_havePetsc) {
      VecDestroy(&_vecPetsc);
      _havePetsc = false;
    }
#endif
  }

  template< class T_Scalar >
  Eigen::VectorX< T_Scalar > Vector< T_Scalar >::getEigen()
  {
    return Eigen::Map< Eigen::VectorX< T_Scalar > >(_vector.data(), _vector.size());
  }

  template< class T_Scalar >
  AlgebraicObjectType Vector< T_Scalar >::type() const
  {
    return AlgebraicObjectType::Vector;
  }

  template< class T_Scalar >
  unsigned long long Vector< T_Scalar >::size() const
  {
    return this->_size[0];
  }

  template< class T_Scalar >
  void Vector< T_Scalar >::resize(const unsigned long long size)
  {
    this->_size[0] = size;
    _vector.resize(size);
  }

  template< class T_Scalar >
  void Vector< T_Scalar >::clear()
  {
    this->_size[0] = 0;
    _vector.clear();
    _vector.shrink_to_fit();
    removePetsc();
  }

  template< class T_Scalar >
  void Vector< T_Scalar >::save(const std::string &name) const
  {
    common::CSVio file(name);
    if(!file.isOpen()) {
      throw common::Exception("Cannot open file " + name + ".csv");
    }

    file << csv::precision(scalar::PrecisionDigits< T_Scalar >::value) << csv::scientific;
    if(scalar::IsComplex< T_Scalar >::value) {
      for(auto i = 0ULL; i < _vector.size(); ++i) {
        file << std::real(_vector[i]) << std::imag(_vector[i]) << csv::endl;
      }
    }
    else {
      for(auto i = 0ULL; i < _vector.size(); ++i) {
        file << _vector[i] << csv::endl;
      }
    }
    file.close();
  }

  template< class T_Scalar >
  void Vector< T_Scalar >::concatenate(const Vector< T_Scalar > &other)
  {
    unsigned long long firstIndex = _vector.size();
    _vector.resize(_vector.size() + other._vector.size());
#pragma omp parallel for num_threads(omp::getMaxThreads())
    for(auto i = 0ULL; i < other._vector.size(); ++i) {
      _vector[firstIndex + i] = other._vector[i];
    }
    this->_size[0] = _vector.size();
  }

  template< class T_Scalar >
  T_Scalar Vector< T_Scalar >::min() const
  {
    return *min_element(_vector.begin(), _vector.end(), [](const T_Scalar &a, const T_Scalar &b) -> bool { return std::norm(a) < std::norm(b); });
  }

  template< class T_Scalar >
  T_Scalar Vector< T_Scalar >::max() const
  {
    return *max_element(_vector.begin(), _vector.end(), [](const T_Scalar &a, const T_Scalar &b) -> bool { return std::norm(a) < std::norm(b); });
  }

  template< class T_Scalar >
  scalar::Precision< T_Scalar > Vector< T_Scalar >::norm(unsigned int p) const
  {
    return std::pow(std::accumulate(_vector.begin(), _vector.end(), 0.,
                                    [p](const scalar::Precision< T_Scalar > &a,
                                        const T_Scalar &b)
                                    -> scalar::Precision< T_Scalar > {
                                      return a + std::pow(std::abs(b), p);
                                    }),
                    1./p);
  }

  template< class T_Scalar >
  void Vector< T_Scalar >::check(const Vector &lhs, const Vector &rhs, const std::string &name){
    if(lhs.size() != rhs.size()) {
      throw common::Exception("Trying to perform a(n) " + name +
                              "on two vectors of different sizes (" +
                              std::to_string(lhs.size()) + " and " +
                              std::to_string(rhs.size()) + ")");
    }
  }

  template< class T_Scalar >
  void Vector< T_Scalar >::eigenop(Vector &lhs, const Vector &rhs, const std::string &name,
                                   T_EigenOpVVV< Eigen::VectorX< T_Scalar > > op){
    check(lhs, rhs, name);
    Eigen::Map< Eigen::VectorX< T_Scalar > >       mLHS(lhs._vector.data(), lhs.size());
    Eigen::Map< const Eigen::VectorX< T_Scalar > > mRHS(rhs._vector.data(), rhs.size());
    (mLHS.*op)(mRHS);
  }

  template< class T_Scalar >
  T_Scalar Vector< T_Scalar >::eigenop(const Vector &lhs, const Vector &rhs, const std::string &name,
                                       T_EigenOpSVV< Eigen::VectorX< T_Scalar > > op){
    check(lhs, rhs, name);
    Eigen::Map< const Eigen::VectorX< T_Scalar > > mLHS(lhs._vector.data(), lhs.size());
    Eigen::Map< const Eigen::VectorX< T_Scalar > > mRHS(rhs._vector.data(), rhs.size());

    return (mLHS.*op)(mRHS);
  }

  template< class T_Scalar >
  void Vector< T_Scalar >::eigenop(Vector &lhs, const T_Scalar &rhs, const std::string &name,
                                   T_EigenOpVVS< Eigen::VectorX< T_Scalar > > op){
    Eigen::Map< Eigen::VectorX< T_Scalar > > mLHS(lhs._vector.data(), lhs.size());
    (mLHS.*op)(rhs);
  }

  INSTANTIATE_CLASS(Vector, 4, TEMPLATE_ARGS(std::complex< double >, std::complex< float >, double, float))

} // namespace gmshfem::algebra
