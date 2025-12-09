// GmshFEM - Copyright (C) 2019-2022, A. Royer, E. Béchet, C. Geuzaine, Université de Liège
//
// See the LICENSE.txt file for license information. Please report all
// issues on https://gitlab.onelab.info/gmsh/fem/issues

#ifndef H_GMSHFEM_VECTORFACTORY
#define H_GMSHFEM_VECTORFACTORY

#include "Memory.h"
#include "OmpInterface.h"
#include "Vector.h"
#include "gmshfemDefines.h"

#include <string>
#include <vector>

typedef struct _p_Vec *Vec;

namespace gmshfem::system
{


  template< class T_Scalar >
  class VectorFactory
  {
   private:
    unsigned long long _size;
    std::vector< T_Scalar > _vec;
    std::vector< T_Scalar > _valuesDC;
    std::vector< T_Scalar > _valuesLC;
    std::vector< unsigned long long > _indicesLC;

    mutable std::vector< void * > _shouldBeDestroyedWithPetsc;

   public:
    VectorFactory();
    VectorFactory(const unsigned long long size);
    ~VectorFactory();

    void init(const unsigned long long size);

    void removeVector();
    void removePetscData();
    void setToZero();
    void setToZero(const std::vector< unsigned long long > &indices);

    void setValuesDC(std::vector< T_Scalar > &values);

    void addValue(const unsigned long long n, const T_Scalar &value);
    void addValues(const unsigned long long n, const std::pair< unsigned long long, int > *const index, const T_Scalar *const values);
    void addValuesDC(const unsigned long long m, const unsigned long long n, const std::pair< unsigned long long, int > *const indexI, const std::pair< unsigned long long, int > *const indexJ, T_Scalar *const values);
    void setValuesLC(std::vector< T_Scalar > &values);
    void setIndicesLC(std::vector< unsigned long long > &indices);

#ifdef HAVE_PETSC
    Vec getPetsc() const;
#endif
    common::Memory memory() const;
    unsigned long long size() const;

    friend class algebra::Vector< T_Scalar >;
  };


} // namespace gmshfem::system

#endif // H_GMSHFEM_VECTORFACTORY
