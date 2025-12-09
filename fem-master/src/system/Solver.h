// GmshFEM - Copyright (C) 2019-2022, A. Royer, E. Béchet, C. Geuzaine, Université de Liège
//
// See the LICENSE.txt file for license information. Please report all
// issues on https://gitlab.onelab.info/gmsh/fem/issues

#ifndef H_GMSHFEM_SOLVER
#define H_GMSHFEM_SOLVER

#include "Memory.h"
#include "gmshfemDefines.h"
#include "scalar.h"

#include <vector>

typedef struct _p_KSP *KSP;

namespace gmshfem::system
{
  template< class T_Scalar >
  class MatrixFactory;

  template< class T_Scalar >
  class VectorFactory;
} // namespace gmshfem::system

namespace gmshfem::system
{


  template< class T_Scalar >
  class Solver
  {
   private:
#ifdef HAVE_PETSC
    KSP _ksp;
#endif
    MatrixFactory< T_Scalar > *_A;
    std::vector<VectorFactory< T_Scalar >>& _b;

   public:
    Solver(MatrixFactory< T_Scalar > *const A, std::vector<VectorFactory< T_Scalar >>& b);
    ~Solver();

    common::Memory getEstimatedFactorizationMemoryUsage() const;

    void solve(std::vector< T_Scalar > &values, const bool reusePreconditioner = false);
    void solveAll(std::vector<std::vector< T_Scalar >> &values, const bool reusePreconditioner = false);
    void eigensolve(std::vector< scalar::ComplexPrecision< T_Scalar > > &eigenvalues, std::vector< std::vector< scalar::ComplexPrecision< T_Scalar > > > &eigenvectors, const bool computeEigenvectors, const unsigned int numberOfEigenvalues, const scalar::ComplexPrecision< T_Scalar > target);
  };


} // namespace gmshfem::system

#endif // H_GMSHFEM_SOLVER
