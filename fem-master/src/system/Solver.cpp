// GmshFEM - Copyright (C) 2019-2022, A. Royer, E. Béchet, C. Geuzaine, Université de Liège
//
// See the LICENSE.txt file for license information. Please report all
// issues on https://gitlab.onelab.info/gmsh/fem/issues

#include "Solver.h"

#include "Exception.h"
#include "MatrixFactory.h"
#include "Message.h"
#include "OmpInterface.h"
#include "PetscInterface.h"
#include "VectorFactory.h"
#include "gmshfemDefines.h"
#include "instantiate.h"

#include <complex>

#if(PETSC_VERSION_MAJOR <= 3) && (PETSC_VERSION_MINOR < 9)
#define PCFactorSetMatSolverType PCFactorSetMatSolverPackage
#endif

namespace gmshfem::system
{


  template< class T_Scalar >
  Solver< T_Scalar >::Solver(MatrixFactory< T_Scalar > *const A, std::vector< VectorFactory< T_Scalar > > &b) :
    _A(A), _b(b)
  {
#ifdef HAVE_PETSC
    KSPCreate(PETSC_COMM_SELF, &_ksp);
#endif // HAVE_PETSC
  }

  template< class T_Scalar >
  Solver< T_Scalar >::~Solver()
  {
    _A = nullptr;
#ifdef HAVE_PETSC
    KSPDestroy(&_ksp);
#endif // HAVE_PETSC
  }

  template< class T_Scalar >
  static bool _isCholesky(const bool isSymetric, const bool isHermitian)
  {
    return (isSymetric && scalar::IsReal< T_Scalar >::value) || isHermitian;
  }

  template< class T_Scalar >
  common::Memory Solver< T_Scalar >::getEstimatedFactorizationMemoryUsage() const
  {
#if defined(HAVE_PETSC) && defined(PETSC_HAVE_MUMPS)
    const bool isCholesky = _isCholesky< T_Scalar >(_A->getOptions().symmetric(), _A->getOptions().hermitian());
    Mat mat = _A->getPetsc();
    Mat matFac;
    MatGetFactor(mat, MATSOLVERMUMPS, (_isCholesky< T_Scalar >(_A->getOptions().symmetric(), _A->getOptions().hermitian()) ? MAT_FACTOR_CHOLESKY : MAT_FACTOR_LU), &matFac);
    if(isCholesky) {
      MatCholeskyFactorSymbolic(matFac, mat, NULL, NULL);
    }
    else {
      MatLUFactorSymbolic(matFac, mat, NULL, NULL, NULL);
    }
    PetscInt ival = 0;
    MatMumpsGetInfog(matFac, 16, &ival);
    return common::Memory(ival * 1e6);
#else
#warning "Estimation of factorization memory usage not available without PETSc and MUMPS"
    msg::error << "Unable to get estimated factorization memory usage without PETSc and MUMPS" << msg::endl;
    return common::Memory(0);
#endif
  }

  template< class T_Scalar >
  void Solver< T_Scalar >::solve(std::vector< T_Scalar > &values, const bool reusePreconditioner)
  {
    //*******************
    // Solve :
    //    A x = b for rhs index 0
    //*******************
    try {
#ifdef HAVE_PETSC
      if(_A->getModule()->name() != "A" && _A->getModule()->name() != "AFrequency") {
        throw common::Exception("A separately assembled system can not be solved");
      }

      PC pc;
      KSPSetReusePreconditioner(_ksp, reusePreconditioner ? PETSC_TRUE : PETSC_FALSE);

      Mat A, B;
      KSPGetOperators(_ksp, &A, &B);
      PetscInt m, n;
      MatGetSize(A, &m, &n);
      if((scalar::IsComplex< T_Scalar >::value == true && scalar::IsComplex< PetscScalar >::value == false ? 2 * _A->size() : _A->size()) != unsigned(m)) {
        KSPSetOperators(_ksp, NULL, NULL);
      }

      if((_A->size() != static_cast< unsigned long long >(m) || _A->size() != static_cast< unsigned long long >(n)) || !reusePreconditioner) {
        Mat mat = _A->getPetsc();
        KSPSetOperators(_ksp, mat, mat);
        MatDestroy(&mat);
      }
      KSPGetPC(_ksp, &pc);

      // use direct sparse solver by default
      KSPSetType(_ksp, "preonly");
      PCSetType(pc, (_isCholesky< T_Scalar >(_A->getOptions().symmetric(), _A->getOptions().hermitian()) ? "cholesky" : "lu"));
#if defined(PETSC_HAVE_MUMPS)
      PCFactorSetMatSolverType(pc, "mumps");
#elif defined(PETSC_HAVE_MKL_PARDISO)
      PCFactorSetMatSolverType(pc, "mkl_pardiso");
#elif defined(PETSC_HAVE_UMFPACK) || defined(PETSC_HAVE_SUITESPARSE)
      PCFactorSetMatSolverType(pc, "umfpack");
#endif

      KSPSetFromOptions(_ksp);

      Vec sol;
      Vec b = _b.at(0).getPetsc();
      VecDuplicate(b, &sol);
      VecCopy(b, sol);
      VecScale(sol, -1.0);

      PetscInt size;
      VecGetLocalSize(sol, &size);

      KSPType ksptype;
      KSPGetType(_ksp, &ksptype);
      PCType pctype;
      PCGetType(pc, &pctype);
      MatSolverType stype;
      PCFactorGetMatSolverType(pc, &stype);
      msg::info << "System size: " << size << " - "
                << (ksptype ? ksptype : "none") << " "
                << (pctype ? pctype : "none") << " "
                << (stype ? stype : "none") << msg::endl;

      KSPSolve(_ksp, sol, sol);

      values.resize((scalar::IsComplex< T_Scalar >::value == true && scalar::IsComplex< PetscScalar >::value == false ? size / 2 : size));

      PetscScalar *array;
      VecGetArray(sol, &array);

      PetscInterface< T_Scalar, PetscScalar > interface;
      interface(values, array, size);
      VecRestoreArray(sol, &array);

      VecDestroy(&sol);
      VecDestroy(&b);
#else
      msg::error << "Unable to solve the linear system without PETSc" << msg::endl;
#endif
    }
    catch(const std::exception &exc) {
      values.clear();
      throw;
    }
  }

  template< class T_Scalar >
  void Solver< T_Scalar >::solveAll(std::vector< std::vector< T_Scalar > > &values, const bool reusePreconditioner)
  {
//*******************
// Solve :
//    A X = B for a thin rectangular matrix B (multi RHS)
//*******************
#if !defined(PETSC_HAVE_HPDDM) || !defined(HAVE_PETSC)
    throw gmshfem::common::Exception("Solving of multiple RHS requires PETSc with HPDDM");
#else

    try {
      if(_A->getModule()->name() != "A" && _A->getModule()->name() != "AFrequency") {
        throw common::Exception("A separately assembled system can not be solved");
      }

      PC pc;
      KSPSetReusePreconditioner(_ksp, reusePreconditioner ? PETSC_TRUE : PETSC_FALSE);

      Mat A, B;
      KSPGetOperators(_ksp, &A, &B);
      PetscInt m, n;
      MatGetSize(A, &m, &n);
      if((scalar::IsComplex< T_Scalar >::value == true && scalar::IsComplex< PetscScalar >::value == false ? 2 * _A->size() : _A->size()) != unsigned(m)) {
        KSPSetOperators(_ksp, NULL, NULL);
      }

      if((_A->size() != static_cast< unsigned long long >(m) || _A->size() != static_cast< unsigned long long >(n)) || !reusePreconditioner) {

        Mat mat = _A->getPetsc();
        KSPSetOperators(_ksp, mat, mat);
        MatGetSize(mat, &m, &n);
        MatDestroy(&mat);
      }
      KSPGetPC(_ksp, &pc);

      // use direct sparse solver by default
      KSPSetType(_ksp, KSPHPDDM);
      KSPHPDDMSetType(_ksp, KSP_HPDDM_TYPE_PREONLY);

      PCSetType(pc, (_isCholesky< T_Scalar >(_A->getOptions().symmetric(), _A->getOptions().hermitian()) ? "cholesky" : "lu"));
#if defined(PETSC_HAVE_MUMPS)
      PCFactorSetMatSolverType(pc, "mumps");
#elif defined(PETSC_HAVE_MKL_PARDISO)
      PCFactorSetMatSolverType(pc, "mkl_pardiso");
#elif defined(PETSC_HAVE_UMFPACK) || defined(PETSC_HAVE_SUITESPARSE)
      PCFactorSetMatSolverType(pc, "umfpack");
#endif

      KSPSetFromOptions(_ksp);

      // Define block matrices
      unsigned nRHS = _b.size();
      Mat BlockRHS, BlockY;

      std::vector< PetscInt > rowIndices(m);
      for(PetscInt i = 0; i < m; i++) {
        rowIndices[i] = i;
      }

      MatCreateSeqDense(MPI_COMM_WORLD, m, nRHS, NULL, &BlockRHS);
      MatDuplicate(BlockRHS, MAT_DO_NOT_COPY_VALUES, &BlockY);
      // For each RHS, copy the vector into the column of the block
      for(PetscInt j = 0; j < static_cast< PetscInt >(nRHS); j++) {
        Vec bj = _b.at(j).getPetsc();
        PetscScalar *colValues;
        VecGetArray(bj, &colValues);
        MatSetValues(BlockRHS, m, rowIndices.data(), 1, &j, colValues, INSERT_VALUES);
        VecRestoreArray(bj, &colValues);
      }

      MatAssemblyBegin(BlockRHS, MAT_FINAL_ASSEMBLY);
      MatAssemblyEnd(BlockRHS, MAT_FINAL_ASSEMBLY);

      MatScale(BlockRHS, -1);


      KSPType ksptype;
      KSPGetType(_ksp, &ksptype);
      PCType pctype;
      PCGetType(pc, &pctype);
      MatSolverType stype;
      PCFactorGetMatSolverType(pc, &stype);
      msg::info << "System size: " << m << " - "
                << (ksptype ? ksptype : "none") << " "
                << (pctype ? pctype : "none") << " "
                << (stype ? stype : "none") << msg::endl;
      KSPMatSolve(_ksp, BlockRHS, BlockY);

#warning "Handle FEM complex + real PETSc ?"
      //values.resize((scalar::IsComplex< T_Scalar >::value == true && scalar::IsComplex< PetscScalar >::value == false ? size / 2 : size));


      // Create a single Vec that contains the block solution, column major format.
      PetscScalar *array;
      MatDenseGetArray(BlockY, &array);
      PetscInterface< T_Scalar, PetscScalar > interface;

      for(unsigned iRHS = 0; iRHS < nRHS; ++iRHS) {
        // Get the pointer to the beginning of the iRHS-th column.
        PetscScalar *col_start = array + iRHS * m;

        values[iRHS].resize(m);
        interface(values[iRHS], col_start, m);
      }

      MatDenseRestoreArray(BlockY, &array);

      MatDestroy(&BlockRHS);
      MatDestroy(&BlockY);

      // Don't set fields to the solved system; loadSolution() should be used.
    }

    catch(const std::exception &exc) {
      values.clear();
      throw;
    }
#endif
  }

  template< class T_Scalar >
  void Solver< T_Scalar >::eigensolve(std::vector< scalar::ComplexPrecision< T_Scalar > > &eigenvalues, std::vector< std::vector< scalar::ComplexPrecision< T_Scalar > > > &eigenvectors, const bool computeEigenvectors, const unsigned int numberOfEigenvalues, const scalar::ComplexPrecision< T_Scalar > target)
  {
    //*******************
    // Solve :
    //    A x = \lambda x
    // or
    //    K x = \lambda M x
    // or
    //    K x = \lambda C x
    // or
    //    (\lambda^2 M + \lambda C + K) x = 0
    //*******************
    try {
#ifdef HAVE_SLEPC
      if(scalar::IsComplex< T_Scalar >::value && scalar::IsReal< PetscScalar >::value) {
        throw common::Exception("Eigensolve cannot be used on a complex formulation if PETSc if compiled in real arithmetic");
      }

      if(_A->getModule()->name() == "MCK") {
        msg::debug << "Using PEP solver" << msg::endl;
        PEP pep;
        PEPCreate(PETSC_COMM_SELF, &pep);

        Mat A[3];
        _A->getModule()->activate('K');
        A[0] = _A->getPetsc();
        _A->getModule()->activate('C');
        A[1] = _A->getPetsc();
        _A->getModule()->activate('M');
        A[2] = _A->getPetsc();
        PEPSetOperators(pep, 3, A);
        PEPSetProblemType(pep, (_A->getOptions().symmetric() && scalar::IsReal< T_Scalar >::value) || _A->getOptions().hermitian() ? PEP_HERMITIAN : PEP_GENERAL);
        PEPSetDimensions(pep, (numberOfEigenvalues == 0 ? _A->size() : numberOfEigenvalues), PETSC_DEFAULT, PETSC_DEFAULT);

        // Shift-and-invert spectral transformation toward target
        ST st;
        PEPGetST(pep, &st);
        STSetType(st, "sinvert");
        PEPSetTarget(pep, target);
        PEPSetWhichEigenpairs(pep, PEP_TARGET_MAGNITUDE);

        // Use MUMPS
        KSP ksp;
        PC pc;
        STGetKSP(st, &ksp);
        KSPSetType(ksp, "preonly");
        KSPGetPC(ksp, &pc);
        PCSetType(pc, "lu");
        PCFactorSetMatSolverType(pc, "mumps");

        // Set from options
        PEPSetFromOptions(pep);

        PEPSolve(pep);
        PetscInt nconv = 0;
        PEPGetConverged(pep, &nconv);

        eigenvalues.resize(nconv);
        if(computeEigenvectors) {
          eigenvectors.resize(nconv);
        }

        // Extract eigenpairs
        Vec xr = nullptr, xi = nullptr;
        PetscScalar kr, ki;
        if(computeEigenvectors) {
          MatCreateVecs(_A->getPetsc(), &xr, &xi);
        }
        for(auto i = 0LL; i < nconv; ++i) {
          PEPGetEigenpair(pep, i, &kr, &ki, xr, xi);
#ifdef PETSC_USE_COMPLEX
          eigenvalues[i] = scalar::ComplexPrecision< T_Scalar >(PetscRealPart(kr), PetscImaginaryPart(kr));
          if(computeEigenvectors) {
            eigenvectors[i].resize(_A->size());

            PetscScalar *array;
            VecGetArray(xr, &array);
            for(auto j = 0ULL; j < _A->size(); ++j) {
              eigenvectors[i][j] = scalar::ComplexPrecision< T_Scalar >(PetscRealPart(array[j]), PetscImaginaryPart(array[j]));
            }
            VecRestoreArray(xr, &array);
          }
#else
          eigenvalues[i] = scalar::ComplexPrecision< T_Scalar >(kr, ki);
          if(computeEigenvectors) {
            eigenvectors[i].resize(_A->size());

            PetscScalar *arrayR, *arrayI;
            VecGetArray(xr, &arrayR);
            VecGetArray(xi, &arrayI);
            for(auto j = 0ULL; j < _A->size(); ++j) {
              eigenvectors[i][j] = scalar::ComplexPrecision< T_Scalar >(arrayR[j], arrayI[j]);
            }
            VecRestoreArray(xi, &arrayI);
            VecRestoreArray(xr, &arrayR);
          }
#endif // PETSC_USE_COMPLEX
        }

        VecDestroy(&xr);
        if(xi != nullptr) {
          VecDestroy(&xi);
        }
        MatDestroy(&A[0]);
        MatDestroy(&A[1]);
        MatDestroy(&A[2]);
        _A->removePetscData();

        PEPDestroy(&pep);
      }
      else {
        msg::debug << "Using EPS solver" << msg::endl;
        EPS eps;
        EPSCreate(PETSC_COMM_SELF, &eps);

        Mat A, B;
        EPSProblemType problemType = EPS_HEP;
        if(_A->getModule()->name() == "A") {
          A = _A->getPetsc();
          B = nullptr;
          problemType = ((_A->getOptions().symmetric() && scalar::IsReal< T_Scalar >::value) || _A->getOptions().hermitian() ? EPS_HEP : EPS_NHEP);
        }
        else if(_A->getModule()->name() == "MK" || _A->getModule()->name() == "CK") {
          _A->getModule()->activate(_A->getModule()->name()[1]);
          A = _A->getPetsc();
          _A->getModule()->activate(_A->getModule()->name()[0]);
          B = _A->getPetsc();
          problemType = ((_A->getOptions().symmetric() && scalar::IsReal< T_Scalar >::value) || _A->getOptions().hermitian() ? EPS_GHEP : EPS_GNHEP);
        }
        else {
          throw common::Exception("Unkown matrix module '" + _A->getModule()->name() + "'");
        }
        EPSSetOperators(eps, A, B);
        EPSSetProblemType(eps, problemType);
        EPSSetDimensions(eps, (numberOfEigenvalues == 0 ? _A->size() : numberOfEigenvalues), PETSC_DEFAULT, PETSC_DEFAULT);

        // Shift-and-invert spectral transformation toward target
        ST st;
        EPSGetST(eps, &st);
        STSetType(st, "sinvert");
        EPSSetTarget(eps, target);
        EPSSetWhichEigenpairs(eps, EPS_TARGET_MAGNITUDE);

        // Use MUMPS
        KSP ksp;
        PC pc;
        STGetKSP(st, &ksp);
        KSPSetType(ksp, "preonly");
        KSPGetPC(ksp, &pc);
        PCSetType(pc, "lu");
        PCFactorSetMatSolverType(pc, "mumps");

        // Set from options
        EPSSetFromOptions(eps);

        // Solve
        EPSSolve(eps);
        PetscInt nconv = 0;
        EPSGetConverged(eps, &nconv);

        eigenvalues.resize(nconv);
        if(computeEigenvectors) {
          eigenvectors.resize(nconv);
        }

        // Extract eigenpairs
        Vec xr = nullptr, xi = nullptr;
        PetscScalar kr, ki;
        if(computeEigenvectors) {
          MatCreateVecs(_A->getPetsc(), &xr, &xi);
        }
        for(auto i = 0LL; i < nconv; ++i) {
          EPSGetEigenpair(eps, i, &kr, &ki, xr, xi);
#ifdef PETSC_USE_COMPLEX
          eigenvalues[i] = scalar::ComplexPrecision< T_Scalar >(PetscRealPart(kr), PetscImaginaryPart(kr));
          if(computeEigenvectors) {
            eigenvectors[i].resize(_A->size());

            PetscScalar *array;
            VecGetArray(xr, &array);
            for(auto j = 0ULL; j < _A->size(); ++j) {
              eigenvectors[i][j] = scalar::ComplexPrecision< T_Scalar >(PetscRealPart(array[j]), PetscImaginaryPart(array[j]));
            }
            VecRestoreArray(xr, &array);
          }
#else
          eigenvalues[i] = scalar::ComplexPrecision< T_Scalar >(kr, ki);
          if(computeEigenvectors) {
            eigenvectors[i].resize(_A->size());

            PetscScalar *arrayR, *arrayI;
            VecGetArray(xr, &arrayR);
            VecGetArray(xi, &arrayI);
            for(auto j = 0ULL; j < _A->size(); ++j) {
              eigenvectors[i][j] = scalar::ComplexPrecision< T_Scalar >(arrayR[j], arrayI[j]);
            }
            VecRestoreArray(xi, &arrayI);
            VecRestoreArray(xr, &arrayR);
          }
#endif // PETSC_USE_COMPLEX
        }

        VecDestroy(&xr);
        VecDestroy(&xi);
        if(_A->getModule()->name() == "A") {
          MatDestroy(&A);
        }
        else if(_A->getModule()->name() == "MK" || _A->getModule()->name() == "CK") {
          MatDestroy(&A);
          MatDestroy(&B);
        }
        _A->removePetscData();

        EPSDestroy(&eps);
      }
#else
      msg::error << "Unable to compute eigenvalues without SLEPc" << msg::endl;
#endif
    }
    catch(const std::exception &exc) {
      throw;
    }
  }

  INSTANTIATE_CLASS(Solver, 4, TEMPLATE_ARGS(std::complex< double >, std::complex< float >, double, float))


} // namespace gmshfem::system
