// GmshFEM - Copyright (C) 2019-2022, A. Royer, E. Béchet, C. Geuzaine, Université de Liège
//
// See the LICENSE.txt file for license information. Please report all
// issues on https://gitlab.onelab.info/gmsh/fem/issues

#include "Matrix.h"

#include "Exception.h"
#include "PetscInterface.h"
#include "instantiate.h"

namespace gmshfem::algebra
{


  template< class T_Scalar >
  Matrix< T_Scalar >::Matrix() :
    AlgebraicObject< T_Scalar >(), _havePetsc(false)
  {
  }

  template< class T_Scalar >
  Matrix< T_Scalar >::Matrix(const unsigned long long size0) :
    AlgebraicObject< T_Scalar >(size0, size0), _havePetsc(false)
  {
  }

  template< class T_Scalar >
  Matrix< T_Scalar >::Matrix(const unsigned long long size0, const unsigned long long size1) :
    AlgebraicObject< T_Scalar >(size0, size1), _havePetsc(false)
  {
  }

  template< class T_Scalar >
  Matrix< T_Scalar >::Matrix(const Matrix &other) :
    AlgebraicObject< T_Scalar >(other), _havePetsc(false)
  {
  }

  template< class T_Scalar >
  Matrix< T_Scalar >::Matrix(Matrix &&other) :
    AlgebraicObject< T_Scalar >(other), _havePetsc(false)
  {
  }

  template< class T_Scalar >
  Matrix< T_Scalar >::~Matrix()
  {
    removePetsc();
  }

  template< class T_Scalar >
  Mat Matrix< T_Scalar >::getPetsc()
  {
#ifdef HAVE_PETSC
    if(_havePetsc) {
      return _matPetsc;
    }
    _buildPetsc();
    _havePetsc = true;
    return _matPetsc;
#else
    throw common::Exception("GmshFEM is not compiled with PETSc");
#endif
  }

  template< class T_Scalar >
  void Matrix< T_Scalar >::removePetsc()
  {
#ifdef HAVE_PETSC
    if(_havePetsc) {
      MatDestroy(&_matPetsc);
      for(auto i = 0ULL; i < _shouldBeDestroyedWithPetsc.size(); ++i) {
        std::free(_shouldBeDestroyedWithPetsc[i]);
      }
      _shouldBeDestroyedWithPetsc.clear();
    }
#endif
  }

  template< class T_Scalar >
  AlgebraicObjectType Matrix< T_Scalar >::type() const
  {
    return AlgebraicObjectType::Matrix;
  }

  template< class T_Scalar >
  unsigned long long Matrix< T_Scalar >::size(const unsigned int dim) const
  {
    if(dim > 1) {
      throw common::Exception("The 'dim' parameter of Matrix::size should be 0 or 1");
    }
    return this->_size[dim];
  }

  template< class T_Scalar >
  double Matrix< T_Scalar >::sparsity() const
  {
    return 1. - numberOfNonZero() / static_cast< double >(this->_size[0] * this->_size[1]);
  }

  template< class T_Scalar >
  bool Matrix< T_Scalar >::isSquare() const
  {
    return (this->_size[0] == this->_size[1]);
  }

  INSTANTIATE_CLASS(Matrix, 4, TEMPLATE_ARGS(std::complex< double >, std::complex< float >, double, float))


} // namespace gmshfem::algebra
