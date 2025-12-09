// GmshFEM - Copyright (C) 2019-2022, A. Royer, E. Béchet, C. Geuzaine, Université de Liège
//
// See the LICENSE.txt file for license information. Please report all
// issues on https://gitlab.onelab.info/gmsh/fem/issues

#include "AlgebraicFunctions.h"

#include "Exception.h"
#include "Message.h"
#include "PetscInterface.h"
#include "instantiate.h"

namespace gmshfem::algebra
{


  template< class T_Scalar >
  void eigenvalues(Vector< scalar::ComplexPrecision< T_Scalar > > &eigenvalue, Matrix< T_Scalar > &matrix, const bool symmetric, const unsigned long long numberOfEigenvalues){
    std::vector< Vector< scalar::ComplexPrecision< T_Scalar > > > dummy;
    eigenvalues(eigenvalue, dummy, matrix, symmetric, numberOfEigenvalues, false);
  }

  INSTANTIATE_FCT(void, , eigenvalues, 0, 4, class, TEMPLATE_ARGS(std::complex< double >, std::complex< float >, double, float), TEMPLATE_PARAMS(Vector< scalar::ComplexPrecision< TEMPLATE_PARAM_1 > > &, Matrix< TEMPLATE_PARAM_1 > &, const bool, const unsigned long long))


  template< class T_Scalar >
  void eigenvalues(Vector< scalar::ComplexPrecision< T_Scalar > > &eigenvalues, std::vector< Vector< scalar::ComplexPrecision< T_Scalar > > > &eigenvectors, Matrix< T_Scalar > &matrix, const bool symmetric, const unsigned long long numberOfEigenvalues, const bool saveeigenvectors)
  {
    if(!matrix.isSquare()) {
      msg::error << "Cannot compute eigenvalues of non-square matrix" << msg::endl;
      return;
    }
    try {
#ifdef HAVE_SLEPC
      Mat A = matrix.getPetsc();

      EPS eps;
      EPSCreate(PETSC_COMM_SELF, &eps);

      EPSSetOperators(eps, A, nullptr);
      EPSSetProblemType(eps, (symmetric ? EPS_HEP : EPS_NHEP));
      EPSSetDimensions(eps, (numberOfEigenvalues == 0 ? matrix.size(0) : numberOfEigenvalues), PETSC_DEFAULT, PETSC_DEFAULT);
      EPSSetFromOptions(eps);

      EPSSolve(eps);
      PetscInt nconv = 0;
      EPSGetConverged(eps, &nconv);
      msg::info << "Found " << nconv << " eigenvalues" << msg::endl;
      std::vector< scalar::ComplexPrecision< T_Scalar > > eig(nconv);
      Vec xr = nullptr, xi = nullptr;
      PetscScalar kr, ki;
      for(auto i = 0LL; i < nconv; ++i) {
        EPSGetEigenpair(eps, i, &kr, &ki, xr, xi);
#ifdef PETSC_USE_COMPLEX
        eig[i] = scalar::ComplexPrecision< T_Scalar >(PetscRealPart(kr), PetscImaginaryPart(kr));
#else
        eig[i] = scalar::ComplexPrecision< T_Scalar >(kr, ki);
#endif // PETSC_USE_COMPLEX
      }

      if(saveeigenvectors){
        eigenvectors.resize(nconv);
        MatCreateVecs(A, nullptr, &xr);
        MatCreateVecs(A, nullptr, &xi);
        for(auto i = 0LL; i < nconv; ++i) {
          const PetscScalar *xrp, *xip;
          EPSGetEigenpair(eps, i, &kr, &ki, xr, xi);
          VecGetArrayRead(xr, &xrp);
          VecGetArrayRead(xi, &xip);

          PetscInt size;
          VecGetLocalSize(xr, &size);
          std::vector< scalar::ComplexPrecision< T_Scalar > > buf(size);

#ifdef PETSC_USE_COMPLEX
          for(size_t v = 0; v < buf.size(); v++) {
            buf[v] = xrp[v];
          }
#else
          for(size_t v = 0; v < buf.size(); v++) {
            buf[v] = scalar::ComplexPrecision< T_Scalar >(xrp[v], xip[v]);
          }
#endif // PETSC_USE_COMPLEX

          VecRestoreArrayRead(xr, &xrp);
          VecRestoreArrayRead(xi, &xip);
          eigenvectors[i] = std::move(buf);
        }
      }

      EPSDestroy(&eps);

      eigenvalues = std::move(eig);
#else
      msg::error << "Cannot compute eigenvalues without SLEPc" << msg::endl;
#endif
    }
    catch(const std::exception &exc) {
      throw;
    }
  }

  INSTANTIATE_FCT(void, , eigenvalues, 1, 4, class, TEMPLATE_ARGS(std::complex< double >, std::complex< float >, double, float), TEMPLATE_PARAMS(Vector< scalar::ComplexPrecision< TEMPLATE_PARAM_1 > > &, std::vector< Vector< scalar::ComplexPrecision< TEMPLATE_PARAM_1 > > > &, Matrix< TEMPLATE_PARAM_1 > &, const bool, const unsigned long long, const bool))


  template< class T_Scalar >
  void eigenvalues(Vector< scalar::ComplexPrecision< T_Scalar > > &eigenvalues, Matrix< T_Scalar > &left, Matrix< T_Scalar > &right, const bool symmetric, const unsigned long long numberOfEigenvalues)
  {
    if(!left.isSquare() || !right.isSquare()) {
      msg::error << "Cannot compute eigenvalues of non-square matrix" << msg::endl;
      return;
    }
    if(left.size(0) != right.size(0) || left.size(1) != right.size(1)) {
      msg::error << "Left and right matrices are not of the same size" << msg::endl;
      return;
    }
    try {
#ifdef HAVE_SLEPC
      EPS eps;
      EPSCreate(PETSC_COMM_SELF, &eps);

      EPSSetOperators(eps, left.getPetsc(), right.getPetsc());
      EPSSetProblemType(eps, (symmetric ? EPS_HEP : EPS_NHEP));
      EPSSetDimensions(eps, (numberOfEigenvalues == 0 ? left.size(0) : numberOfEigenvalues), PETSC_DEFAULT, PETSC_DEFAULT);
      EPSSetFromOptions(eps);

      EPSSolve(eps);
      PetscInt nconv = 0;
      EPSGetConverged(eps, &nconv);
      msg::info << "Found " << nconv << " eigenvalues" << msg::endl;
      std::vector< scalar::ComplexPrecision< T_Scalar > > eig(nconv);
      Vec xr = nullptr, xi = nullptr;
      PetscScalar kr, ki;
      for(auto i = 0LL; i < nconv; ++i) {
        EPSGetEigenpair(eps, i, &kr, &ki, xr, xi);
#ifdef PETSC_USE_COMPLEX
        eig[i] = scalar::ComplexPrecision< T_Scalar >(PetscRealPart(kr), PetscImaginaryPart(kr));
#else
        eig[i] = scalar::ComplexPrecision< T_Scalar >(kr, ki);
#endif // PETSC_USE_COMPLEX
      }

      EPSDestroy(&eps);

      eigenvalues = std::move(eig);
#else
      msg::error << "Cannot compute eigenvalues without SLEPc" << msg::endl;
#endif
    }
    catch(const std::exception &exc) {
      throw;
    }
  }

  INSTANTIATE_FCT(void, , eigenvalues, 2, 4, class, TEMPLATE_ARGS(std::complex< double >, std::complex< float >, double, float), TEMPLATE_PARAMS(Vector< scalar::ComplexPrecision< TEMPLATE_PARAM_1 > > &, Matrix< TEMPLATE_PARAM_1 > &, Matrix< TEMPLATE_PARAM_1 > &, const bool, const unsigned long long))


  template< class T_Scalar >
  void singularValues(Vector< scalar::Precision< T_Scalar > > &singularValues, Matrix< T_Scalar > &matrix, const bool smallest, const unsigned long long numberOfSingularValues)
  {
    try {
#ifdef HAVE_SLEPC
      SVD svd;
      SVDCreate(PETSC_COMM_SELF, &svd);

      //        SVDSetOperator(svd, matrix.getPetsc());
      SVDSetFromOptions(svd);
      SVDSetDimensions(svd, (numberOfSingularValues == 0 ? std::min(matrix.size(0), matrix.size(1)) : numberOfSingularValues), PETSC_DEFAULT, PETSC_DEFAULT);
      SVDSetWhichSingularTriplets(svd, (smallest ? SVD_SMALLEST : SVD_LARGEST));
      SVDSetFromOptions(svd);

      SVDSolve(svd);
      PetscInt nconv = 0;
      SVDGetConverged(svd, &nconv);
      msg::info << "Found " << nconv << " singular values" << msg::endl;

      std::vector< scalar::Precision< T_Scalar > > sigma(nconv);
      PetscReal singularValue = 0.;
      Vec u = nullptr, v = nullptr;
      for(auto i = 0LL; i < nconv; ++i) {
        SVDGetSingularTriplet(svd, i, &singularValue, u, v);
        sigma[i] = singularValue;
      }

      SVDDestroy(&svd);

      singularValues = std::move(sigma);
#else
      msg::error << "Cannot compute singular values without SLEPc" << msg::endl;
#endif
    }
    catch(const std::exception &exc) {
      throw;
    }
  }

  INSTANTIATE_FCT(void, , singularValues, 0, 4, class, TEMPLATE_ARGS(std::complex< double >, std::complex< float >, double, float), TEMPLATE_PARAMS(Vector< scalar::Precision< TEMPLATE_PARAM_1 > > &, Matrix< TEMPLATE_PARAM_1 > &, const bool, const unsigned long long))


  template< class T_Scalar >
  T_Scalar bilinear(const Vector< T_Scalar > &x, const Matrix< T_Scalar > &A, const Vector< T_Scalar > &y)
  {
    if(x.size() != A.size(0) || A.size(1) != y.size() || !A.isSquare()) {
      throw common::Exception("Size of objects are incompatible in 'bilinearProduct' x = (" + std::to_string(x.size()) + "), A = (" + std::to_string(A.size(0)) + ", " + std::to_string(A.size(1)) + "), y = (" + std::to_string(y.size()) + ")");
    }

    T_Scalar value = 0.;
    unsigned long long size = A.size(0);

    if(A.format() == MatrixFormat::CRS) {
      const MatrixCRS< T_Scalar > &tmp = static_cast< const MatrixCRS< T_Scalar > & >(A);
#pragma omp parallel num_threads(omp::getMaxThreads())
      {
        T_Scalar product = 0.;
#pragma omp for
        for(auto i = 0ULL; i < size; ++i) {
          T_Scalar line = 0.;
          for(auto jj = tmp.ai()[i]; jj < tmp.ai()[i + 1]; ++jj) {
            const unsigned long long j = tmp.aj()[jj];
            line += tmp.a()[jj] * y[j];
          }
          product += x[i] * line;
        }

#pragma omp critical
        value += product;
      }
    }
    else if(A.format() == MatrixFormat::CRSFast) {
      const MatrixCRSFast< T_Scalar > &tmp = static_cast< const MatrixCRSFast< T_Scalar > & >(A);
#pragma omp parallel num_threads(omp::getMaxThreads())
      {
        T_Scalar product = 0.;
#pragma omp for
        for(auto i = 0ULL; i < size; ++i) {
          T_Scalar line = 0.;
          for(auto jj = tmp.ai()[i]; jj < tmp.ai()[i + 1]; ++jj) {
            const unsigned long long j = tmp.aj()[jj];
            line += tmp.a()[jj] * y[j];
          }
          product += x[i] * line;
        }

#pragma omp critical
        value += product;
      }
    }
    else {
      msg::error << "'bilinearProduct' is not implemented for this format of matrix" << msg::endl;
    }

    return value;
  }

  INSTANTIATE_FCT(TEMPLATE_PARAM_1, , bilinear, 0, 4, class, TEMPLATE_ARGS(std::complex< double >, std::complex< float >, double, float), TEMPLATE_PARAMS(const Vector< TEMPLATE_PARAM_1 > &, const Matrix< TEMPLATE_PARAM_1 > &, const Vector< TEMPLATE_PARAM_1 > &))


  template< class T_Scalar >
  void linear(Vector< T_Scalar > &b, const Matrix< T_Scalar > &A, const Vector< T_Scalar > &x)
  {
    if(x.size() != A.size(1)) {
      throw common::Exception("Size of objects are incompatible in 'linear' A = (" + std::to_string(A.size(0)) + ", " + std::to_string(A.size(1)) + "), x = (" + std::to_string(x.size()) + ")");
    }

    b.resize(A.size(0));
    unsigned long long size = A.size(0);

    if(A.format() == MatrixFormat::CRS) {
      const MatrixCRS< T_Scalar > &tmp = static_cast< const MatrixCRS< T_Scalar > & >(A);
#pragma omp parallel num_threads(omp::getMaxThreads())
      {
#pragma omp for
        for(auto i = 0ULL; i < size; ++i) {
          T_Scalar line = 0.;
          for(auto jj = tmp.ai()[i]; jj < tmp.ai()[i + 1]; ++jj) {
            const unsigned long long j = tmp.aj()[jj];
            line += tmp.a()[jj] * x[j];
          }
          b[i] = line;
        }
      }
    }
    else if(A.format() == MatrixFormat::CRSFast) {
      const MatrixCRSFast< T_Scalar > &tmp = static_cast< const MatrixCRSFast< T_Scalar > & >(A);
#pragma omp parallel num_threads(omp::getMaxThreads())
      {
#pragma omp for
        for(auto i = 0ULL; i < size; ++i) {
          T_Scalar line = 0.;
          for(auto jj = tmp.ai()[i]; jj < tmp.ai()[i + 1]; ++jj) {
            const unsigned long long j = tmp.aj()[jj];
            line += tmp.a()[jj] * x[j];
          }
          b[i] = line;
        }
      }
    }
    else {
      msg::error << "'linear' is not implemented for this format of matrix" << msg::endl;
    }
  }

  INSTANTIATE_FCT(void, , linear, 0, 4, class, TEMPLATE_ARGS(std::complex< double >, std::complex< float >, double, float), TEMPLATE_PARAMS(Vector< TEMPLATE_PARAM_1 > &, const Matrix< TEMPLATE_PARAM_1 > &, const Vector< TEMPLATE_PARAM_1 > &))


  template< class T_Scalar >
  scalar::Precision< T_Scalar > residual(const Vector< T_Scalar > &b, const Matrix< T_Scalar > &A, const Vector< T_Scalar > &x)
  {
    if(x.size() != A.size(1) || A.size(0) != b.size()) {
      throw common::Exception("Size of objects are incompatible in 'residual' b = (" + std::to_string(b.size()) + "), A = (" + std::to_string(A.size(0)) + ", " + std::to_string(A.size(1)) + "), x = (" + std::to_string(x.size()) + ")");
    }

    scalar::Precision< T_Scalar > value = 0.;
    unsigned long long size = A.size(0);

    if(A.format() == MatrixFormat::CRS) {
      const MatrixCRS< T_Scalar > &tmp = static_cast< const MatrixCRS< T_Scalar > & >(A);
#pragma omp parallel num_threads(omp::getMaxThreads())
      {
        scalar::Precision< T_Scalar > product = 0.;
#pragma omp for
        for(auto i = 0ULL; i < size; ++i) {
          T_Scalar line = 0.;
          for(auto jj = tmp.ai()[i]; jj < tmp.ai()[i + 1]; ++jj) {
            const unsigned long long j = tmp.aj()[jj];
            line += tmp.a()[jj] * x[j];
          }
          product += std::abs((b[i] + line) * (b[i] + line));
        }

#pragma omp critical
        value += product;
      }
    }
    else if(A.format() == MatrixFormat::CRSFast) {
      const MatrixCRSFast< T_Scalar > &tmp = static_cast< const MatrixCRSFast< T_Scalar > & >(A);
#pragma omp parallel num_threads(omp::getMaxThreads())
      {
        scalar::Precision< T_Scalar > product = 0.;
#pragma omp for
        for(auto i = 0ULL; i < size; ++i) {
          T_Scalar line = 0.;
          for(auto jj = tmp.ai()[i]; jj < tmp.ai()[i + 1]; ++jj) {
            const unsigned long long j = tmp.aj()[jj];
            line += tmp.a()[jj] * x[j];
          }
          product += std::abs((b[i] + line) * (b[i] + line));
        }

#pragma omp critical
        value += product;
      }
    }
    else {
      msg::error << "'residual' is not implemented for this format of matrix" << msg::endl;
    }

    value = std::sqrt(value);

    return value;
  }

  INSTANTIATE_FCT(scalar::Precision< TEMPLATE_PARAM_1 >, , residual, 0, 4, class, TEMPLATE_ARGS(std::complex< double >, std::complex< float >, double, float), TEMPLATE_PARAMS(const Vector< TEMPLATE_PARAM_1 > &, const Matrix< TEMPLATE_PARAM_1 > &, const Vector< TEMPLATE_PARAM_1 > &))


} // namespace gmshfem::algebra
