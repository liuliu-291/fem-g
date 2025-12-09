// GmshFEM - Copyright (C) 2019-2022, A. Royer, E. Béchet, C. Geuzaine, Université de Liège
//
// See the LICENSE.txt file for license information. Please report all
// issues on https://gitlab.onelab.info/gmsh/fem/issues

#include "AlgebraicObject.h"

#include "instantiate.h"

#include <complex>

namespace gmshfem::algebra
{


  template< class T_Scalar >
  AlgebraicObject< T_Scalar >::AlgebraicObject() :
    _size{0, 0}
  {
  }

  template< class T_Scalar >
  AlgebraicObject< T_Scalar >::AlgebraicObject(const unsigned long long size0, const unsigned long long size1) :
    _size{size0, size1}
  {
  }

  template< class T_Scalar >
  AlgebraicObject< T_Scalar >::AlgebraicObject(const AlgebraicObject &other) :
    _size{other._size[0], other._size[1]}
  {
  }

  template< class T_Scalar >
  AlgebraicObject< T_Scalar >::AlgebraicObject(AlgebraicObject &&other) :
    _size{std::move(other._size[0]), std::move(other._size[1])}
  {
  }

  template< class T_Scalar >
  AlgebraicObject< T_Scalar >::~AlgebraicObject()
  {
  }

  template< class T_Scalar >
  void AlgebraicObject< T_Scalar >::copySize(const AlgebraicObject &other)
  {
    _size[0] = other._size[0];
    _size[1] = other._size[1];
  }

  template< class T_Scalar >
  void AlgebraicObject< T_Scalar >::copySize(AlgebraicObject &&other)
  {
    _size[0] = std::move(other._size[0]);
    _size[1] = std::move(other._size[1]);
  }

  INSTANTIATE_CLASS(AlgebraicObject, 4, TEMPLATE_ARGS(std::complex< double >, std::complex< float >, double, float))


} // namespace gmshfem::algebra
