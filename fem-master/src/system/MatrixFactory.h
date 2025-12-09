// GmshFEM - Copyright (C) 2019-2022, A. Royer, E. Béchet, C. Geuzaine, Université de Liège
//
// See the LICENSE.txt file for license information. Please report all
// issues on https://gitlab.onelab.info/gmsh/fem/issues

#ifndef H_GMSHFEM_MATRIXFACTORY
#define H_GMSHFEM_MATRIXFACTORY

#include "DofsManager.h"
#include "Matrix.h"
#include "MatrixFactory.h"
#include "MatrixModule.h"
#include "MatrixOptions.h"
#include "Memory.h"
#include "OmpInterface.h"
#include "gmshfemDefines.h"
#include "scalar.h"

#include <vector>

typedef struct _p_Mat *Mat;

namespace gmshfem::system
{


  template< class T_Scalar >
  class MatrixFactory
  {
   private:
    unsigned long long _size;
    system::MatrixOptions _options;

    unsigned long long *_ai;
    unsigned long long *_aj;
    MatrixModule< T_Scalar > *_module;
    std::vector< T_Scalar > _valuesLC;
    std::vector< unsigned long long > _indicesLC;

    std::vector< unsigned long long > **_pattern;
    std::vector< unsigned long long > *_patternMemory;
    bool _havePattern;
    std::vector< omp::Lock > _locks;

    mutable std::vector< void * > _shouldBeDestroyedWithPetsc;

    void _reorder(const std::vector< unsigned long long > &degree);

   public:
    MatrixFactory(const std::string &options = "");
    MatrixFactory(const unsigned long long size, const std::string &options = "", const unsigned long long nbrBubble = 0);
    ~MatrixFactory();

    void setModule(MatrixModule< T_Scalar > *module);
    const MatrixModule< T_Scalar > *getModule() const;
    void init(const unsigned long long size, const unsigned long long nbrBubble = 0);

    void removeMatrix();
    void removePetscData();
    void setPatternToZero();
    void setToZero();
    void finalizePattern();

    void setValuesLC(std::vector< T_Scalar > &values);
    void setIndicesLC(std::vector< unsigned long long > &indices);

    bool addValues(const unsigned long long m, const unsigned long long n, const std::pair< unsigned long long, int > *const indexI, const std::pair< unsigned long long, int > *const indexJ, T_Scalar *const values);
    void addValue(const unsigned long long i, const unsigned long long j, const T_Scalar &values);
    void addPattern(const unsigned long long m, const unsigned long long n, const std::pair< unsigned long long, int > *const indexI, const std::pair< unsigned long long, int > *const indexJ);
    void addPatternGlobalDof(const unsigned long long i);

    bool havePattern() const;
#ifdef HAVE_PETSC
    Mat getPetsc() const;
#endif

    common::Memory memory() const;
    unsigned long long size() const;
    unsigned long long numberOfNonZeros() const;
    system::MatrixOptions getOptions() const;

    // Useful for the RCM reordering algorithm
    friend void dofs::DofsManager< T_Scalar >::reorderWithRCM(system::MatrixFactory< T_Scalar > *const A);

    friend class algebra::MatrixCCS< T_Scalar >;
    friend class algebra::MatrixCRS< T_Scalar >;
    friend class algebra::MatrixCRSFast< T_Scalar >;
  };


} // namespace gmshfem::system

#endif // H_GMSHFEM_MATRIXFACTORY
