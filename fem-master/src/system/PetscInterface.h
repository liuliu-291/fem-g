// GmshFEM - Copyright (C) 2019-2022, A. Royer, E. Béchet, C. Geuzaine, Université de Liège
//
// See the LICENSE.txt file for license information. Please report all
// issues on https://gitlab.onelab.info/gmsh/fem/issues

#ifndef H_GMSHFEM_PETSCINTERFACE
#define H_GMSHFEM_PETSCINTERFACE

#include "gmshfemDefines.h"

#ifdef HAVE_PETSC

#include "Message.h"
#include "OmpInterface.h"
#include "scalar.h"

#include <petsc.h>

#ifdef HAVE_SLEPC
#include <slepc.h>
#endif // HAVE_SLEPC

#include <vector>

namespace gmshfem::system
{


  // T_Scalar, PetscScalar
  template< class T_Scalar1, class T_Scalar2, bool T_isComplex1 = scalar::IsComplex< T_Scalar1 >::value, bool T_isComplex2 = scalar::IsComplex< T_Scalar2 >::value >
  struct PetscInterface {
    // For matrix
    void operator()(const unsigned long long size0, const unsigned long long size1, const unsigned long long *const row, const unsigned long long *const indices, T_Scalar1 *const values, Mat *matPetsc, std::vector< void * > &shouldBeDestroyedWithPetsc, const bool isSymmetric)
    {
    }

    // For vector
    void operator()(const unsigned long long size, const T_Scalar1 *const values, Vec *vecPetsc, std::vector< void * > &shouldBeDestroyedWithPetsc)
    {
    }

    // For solution
    void operator()(std::vector< T_Scalar1 > &values, T_Scalar2 *array, unsigned long long size)
    {
    }
  };

  template< class T_Scalar1 >
  struct PetscInterface< T_Scalar1, T_Scalar1, true, true > {
    // For matrix
    void operator()(const unsigned long long size0, const unsigned long long size1, const unsigned long long *const row, const unsigned long long *const indices, T_Scalar1 *const values, Mat *matPetsc, std::vector< void * > &shouldBeDestroyedWithPetsc, const bool isSymmetric)
    {
      PetscInt *petscRow;
      PetscInt *petscIndices;
      if(sizeof(PetscInt) == sizeof(unsigned long long)) {
        petscRow = const_cast< PetscInt * >(reinterpret_cast< const PetscInt * >(row));
        petscIndices = const_cast< PetscInt * >(reinterpret_cast< const PetscInt * >(indices));
      }
      else {
        petscRow = (PetscInt *)std::malloc((size0 + 1) * sizeof(PetscInt));
        petscIndices = (PetscInt *)std::malloc(row[size0] * sizeof(PetscInt));
#pragma omp parallel for num_threads(omp::getMaxThreads())
        for(auto i = 0ULL; i < size0 + 1; ++i) {
          petscRow[i] = row[i];
        }
#pragma omp parallel for num_threads(omp::getMaxThreads())
        for(auto i = 0ULL; i < row[size0]; ++i) {
          petscIndices[i] = indices[i];
        }
        shouldBeDestroyedWithPetsc.push_back(petscRow);
        shouldBeDestroyedWithPetsc.push_back(petscIndices);
      }
      if(isSymmetric) {
        MatCreateSeqSBAIJWithArrays(PETSC_COMM_SELF, 1, size0, size1, petscRow, petscIndices, values, matPetsc);
      }
      else {
        MatCreateSeqAIJWithArrays(PETSC_COMM_SELF, size0, size1, petscRow, petscIndices, values, matPetsc);
      }
    }

    // For vector
    void operator()(const unsigned long long size, const T_Scalar1 *const values, Vec *vecPetsc, std::vector< void * > &shouldBeDestroyedWithPetsc)
    {
      VecCreateSeqWithArray(PETSC_COMM_SELF, 1, size, values, vecPetsc);
    }

    // For solution
    void operator()(std::vector< T_Scalar1 > &values, T_Scalar1 *array, unsigned long long size)
    {
#pragma omp parallel for num_threads(omp::getMaxThreads())
      for(auto i = 0ULL; i < size; ++i) {
        values[i] = array[i];
      }
    }
  };

  template< class T_Scalar1 >
  struct PetscInterface< T_Scalar1, T_Scalar1, false, false > {
    // For matrix
    void operator()(const unsigned long long size0, const unsigned long long size1, const unsigned long long *const row, const unsigned long long *const indices, T_Scalar1 *const values, Mat *matPetsc, std::vector< void * > &shouldBeDestroyedWithPetsc, const bool isSymmetric)
    {
      PetscInt *petscRow;
      PetscInt *petscIndices;
      if(sizeof(PetscInt) == sizeof(unsigned long long)) {
        petscRow = const_cast< PetscInt * >(reinterpret_cast< const PetscInt * >(row));
        petscIndices = const_cast< PetscInt * >(reinterpret_cast< const PetscInt * >(indices));
      }
      else {
        petscRow = (PetscInt *)std::malloc((size0 + 1) * sizeof(PetscInt));
        petscIndices = (PetscInt *)std::malloc(row[size0] * sizeof(PetscInt));
#pragma omp parallel for num_threads(omp::getMaxThreads())
        for(auto i = 0ULL; i < size0 + 1; ++i) {
          petscRow[i] = row[i];
        }
#pragma omp parallel for num_threads(omp::getMaxThreads())
        for(auto i = 0ULL; i < row[size0]; ++i) {
          petscIndices[i] = indices[i];
        }
        shouldBeDestroyedWithPetsc.push_back(petscRow);
        shouldBeDestroyedWithPetsc.push_back(petscIndices);
      }

      if(isSymmetric) {
        MatCreateSeqSBAIJWithArrays(PETSC_COMM_SELF, 1, size0, size1, petscRow, petscIndices, values, matPetsc);
      }
      else {
        MatCreateSeqAIJWithArrays(PETSC_COMM_SELF, size0, size1, petscRow, petscIndices, values, matPetsc);
      }
    }

    // For vector
    void operator()(const unsigned long long size, const T_Scalar1 *const values, Vec *vecPetsc, std::vector< void * > &shouldBeDestroyedWithPetsc)
    {
      VecCreateSeqWithArray(PETSC_COMM_SELF, 1, size, values, vecPetsc);
    }

    // For solution
    void operator()(std::vector< T_Scalar1 > &values, T_Scalar1 *array, unsigned long long size)
    {
#pragma omp parallel for num_threads(omp::getMaxThreads())
      for(auto i = 0ULL; i < size; ++i) {
        values[i] = array[i];
      }
    }
  };

  template< class T_Scalar1, class T_Scalar2 >
  struct PetscInterface< T_Scalar1, T_Scalar2, true, true > {
    // For matrix
    void operator()(const unsigned long long size0, const unsigned long long size1, const unsigned long long *const row, const unsigned long long *const indices, T_Scalar1 *const values, Mat *matPetsc, std::vector< void * > &shouldBeDestroyedWithPetsc, const bool isSymmetric)
    {
      msg::warning << "Loss of efficiency can be noticed if you use a '" << scalar::Name< T_Scalar1 >::name << "' formulation with a '" << scalar::Name< T_Scalar2 >::name << "' PETSc" << msg::endl;
      T_Scalar2 *petscValues = (T_Scalar2 *)std::malloc(row[size0] * sizeof(T_Scalar2));
#pragma omp parallel for num_threads(omp::getMaxThreads())
      for(auto i = 0ULL; i < row[size0]; ++i) {
        petscValues[i] = T_Scalar2(std::real(values[i]), std::imag(values[i]));
      }

      PetscInt *petscRow;
      PetscInt *petscIndices;
      if(sizeof(PetscInt) == sizeof(unsigned long long)) {
        petscRow = const_cast< PetscInt * >(reinterpret_cast< const PetscInt * >(row));
        petscIndices = const_cast< PetscInt * >(reinterpret_cast< const PetscInt * >(indices));
      }
      else {
        petscRow = (PetscInt *)std::malloc((size0 + 1) * sizeof(PetscInt));
        petscIndices = (PetscInt *)std::malloc(row[size0] * sizeof(PetscInt));
#pragma omp parallel for num_threads(omp::getMaxThreads())
        for(auto i = 0ULL; i < size0 + 1; ++i) {
          petscRow[i] = row[i];
        }
#pragma omp parallel for num_threads(omp::getMaxThreads())
        for(auto i = 0ULL; i < row[size0]; ++i) {
          petscIndices[i] = indices[i];
        }
        shouldBeDestroyedWithPetsc.push_back(petscRow);
        shouldBeDestroyedWithPetsc.push_back(petscIndices);
      }

      if(isSymmetric) {
        MatCreateSeqSBAIJWithArrays(PETSC_COMM_SELF, 1, size0, size1, petscRow, petscIndices, petscValues, matPetsc);
      }
      else {
        MatCreateSeqAIJWithArrays(PETSC_COMM_SELF, size0, size1, petscRow, petscIndices, petscValues, matPetsc);
      }
      shouldBeDestroyedWithPetsc.push_back(petscValues);
    }

    // For vector
    void operator()(const unsigned long long size, const T_Scalar1 *const values, Vec *vecPetsc, std::vector< void * > &shouldBeDestroyedWithPetsc)
    {
      msg::warning << "Loss of efficiency can be noticed if you use a '" << scalar::Name< T_Scalar1 >::name << "' formulation with a '" << scalar::Name< T_Scalar2 >::name << "' PETSc" << msg::endl;
      T_Scalar2 *petscValues = (T_Scalar2 *)std::malloc(size * sizeof(T_Scalar2));
#pragma omp parallel for num_threads(omp::getMaxThreads())
      for(auto i = 0ULL; i < size; ++i) {
        petscValues[i] = T_Scalar2(std::real(values[i]), std::imag(values[i]));
      }
      VecCreateSeqWithArray(PETSC_COMM_SELF, 1, size, petscValues, vecPetsc);
      shouldBeDestroyedWithPetsc.push_back(petscValues);
    }

    // For solution
    void operator()(std::vector< T_Scalar1 > &values, T_Scalar2 *array, unsigned long long size)
    {
      msg::warning << "The solution is cast because you use a '" << scalar::Name< T_Scalar1 >::name << "' formulation with a '" << scalar::Name< T_Scalar2 >::name << "' PETSc" << msg::endl;
#pragma omp parallel for num_threads(omp::getMaxThreads())
      for(auto i = 0ULL; i < size; ++i) {
        values[i] = T_Scalar1(array[i]);
      }
    }
  };

  template< class T_Scalar1, class T_Scalar2 >
  struct PetscInterface< T_Scalar1, T_Scalar2, false, false > {
    // For matrix
    void operator()(const unsigned long long size0, const unsigned long long size1, const unsigned long long *const row, const unsigned long long *const indices, T_Scalar1 *const values, Mat *matPetsc, std::vector< void * > &shouldBeDestroyedWithPetsc, const bool isSymmetric)
    {
      msg::warning << "Loss of efficiency can be noticed if you use a '" << scalar::Name< T_Scalar1 >::name << "' formulation with a '" << scalar::Name< T_Scalar2 >::name << "' PETSc" << msg::endl;
      T_Scalar2 *petscValues = (T_Scalar2 *)std::malloc(row[size0] * sizeof(T_Scalar2));
#pragma omp parallel for num_threads(omp::getMaxThreads())
      for(auto i = 0ULL; i < row[size0]; ++i) {
        petscValues[i] = T_Scalar2(values[i]);
      }

      PetscInt *petscRow;
      PetscInt *petscIndices;
      if(sizeof(PetscInt) == sizeof(unsigned long long)) {
        petscRow = const_cast< PetscInt * >(reinterpret_cast< const PetscInt * >(row));
        petscIndices = const_cast< PetscInt * >(reinterpret_cast< const PetscInt * >(indices));
      }
      else {
        petscRow = (PetscInt *)std::malloc((size0 + 1) * sizeof(PetscInt));
        petscIndices = (PetscInt *)std::malloc(row[size0] * sizeof(PetscInt));
#pragma omp parallel for num_threads(omp::getMaxThreads())
        for(auto i = 0ULL; i < size0 + 1; ++i) {
          petscRow[i] = row[i];
        }
#pragma omp parallel for num_threads(omp::getMaxThreads())
        for(auto i = 0ULL; i < row[size0]; ++i) {
          petscIndices[i] = indices[i];
        }
        shouldBeDestroyedWithPetsc.push_back(petscRow);
        shouldBeDestroyedWithPetsc.push_back(petscIndices);
      }

      if(isSymmetric) {
        MatCreateSeqSBAIJWithArrays(PETSC_COMM_SELF, 1, size0, size1, petscRow, petscIndices, petscValues, matPetsc);
      }
      else {
        MatCreateSeqAIJWithArrays(PETSC_COMM_SELF, size0, size1, petscRow, petscIndices, petscValues, matPetsc);
      }
      shouldBeDestroyedWithPetsc.push_back(petscValues);
    }

    // For vector
    void operator()(const unsigned long long size, const T_Scalar1 *const values, Vec *vecPetsc, std::vector< void * > &shouldBeDestroyedWithPetsc)
    {
      msg::warning << "Loss of efficiency can be noticed if you use a '" << scalar::Name< T_Scalar1 >::name << "' formulation with a '" << scalar::Name< T_Scalar2 >::name << "' PETSc" << msg::endl;
      T_Scalar2 *petscValues = (T_Scalar2 *)std::malloc(size * sizeof(T_Scalar2));
#pragma omp parallel for num_threads(omp::getMaxThreads())
      for(auto i = 0ULL; i < size; ++i) {
        petscValues[i] = T_Scalar2(values[i]);
      }
      VecCreateSeqWithArray(PETSC_COMM_SELF, 1, size, petscValues, vecPetsc);
      shouldBeDestroyedWithPetsc.push_back(petscValues);
    }

    // For solution
    void operator()(std::vector< T_Scalar1 > &values, T_Scalar2 *array, unsigned long long size)
    {
      msg::warning << "The solution is cast because you use a '" << scalar::Name< T_Scalar1 >::name << "' formulation with a '" << scalar::Name< T_Scalar2 >::name << "' PETSc" << msg::endl;
#pragma omp parallel for num_threads(omp::getMaxThreads())
      for(auto i = 0ULL; i < size; ++i) {
        values[i] = T_Scalar1(array[i]);
      }
    }
  };

  template< class T_Scalar1, class T_Scalar2 >
  struct PetscInterface< T_Scalar1, T_Scalar2, false, true > {
    // For matrix
    void operator()(const unsigned long long size0, const unsigned long long size1, const unsigned long long *const row, const unsigned long long *const indices, T_Scalar1 *const values, Mat *matPetsc, std::vector< void * > &shouldBeDestroyedWithPetsc, const bool isSymmetric)
    {
      msg::warning << "Loss of efficiency can be noticed if you use a '" << scalar::Name< T_Scalar1 >::name << "' formulation with a '" << scalar::Name< T_Scalar2 >::name << "' PETSc" << msg::endl;
      T_Scalar2 *petscValues = (T_Scalar2 *)std::malloc(row[size0] * sizeof(T_Scalar2));
#pragma omp parallel for num_threads(omp::getMaxThreads())
      for(auto i = 0ULL; i < row[size0]; ++i) {
        petscValues[i] = T_Scalar2(values[i]);
      }

      PetscInt *petscRow;
      PetscInt *petscIndices;
      if(sizeof(PetscInt) == sizeof(unsigned long long)) {
        petscRow = const_cast< PetscInt * >(reinterpret_cast< const PetscInt * >(row));
        petscIndices = const_cast< PetscInt * >(reinterpret_cast< const PetscInt * >(indices));
      }
      else {
        petscRow = (PetscInt *)std::malloc((size0 + 1) * sizeof(PetscInt));
        petscIndices = (PetscInt *)std::malloc(row[size0] * sizeof(PetscInt));
#pragma omp parallel for num_threads(omp::getMaxThreads())
        for(auto i = 0ULL; i < size0 + 1; ++i) {
          petscRow[i] = row[i];
        }
#pragma omp parallel for num_threads(omp::getMaxThreads())
        for(auto i = 0ULL; i < row[size0]; ++i) {
          petscIndices[i] = indices[i];
        }
        shouldBeDestroyedWithPetsc.push_back(petscRow);
        shouldBeDestroyedWithPetsc.push_back(petscIndices);
      }

      if(isSymmetric) {
        MatCreateSeqSBAIJWithArrays(PETSC_COMM_SELF, 1, size0, size1, petscRow, petscIndices, petscValues, matPetsc);
      }
      else {
        MatCreateSeqAIJWithArrays(PETSC_COMM_SELF, size0, size1, petscRow, petscIndices, petscValues, matPetsc);
      }
      shouldBeDestroyedWithPetsc.push_back(petscValues);
    }

    // For vector
    void operator()(const unsigned long long size, const T_Scalar1 *const values, Vec *vecPetsc, std::vector< void * > &shouldBeDestroyedWithPetsc)
    {
      msg::warning << "Loss of efficiency can be noticed if you use a '" << scalar::Name< T_Scalar1 >::name << "' formulation with a '" << scalar::Name< T_Scalar2 >::name << "' PETSc" << msg::endl;
      T_Scalar2 *petscValues = (T_Scalar2 *)std::malloc(size * sizeof(T_Scalar2));
#pragma omp parallel for num_threads(omp::getMaxThreads())
      for(auto i = 0ULL; i < size; ++i) {
        petscValues[i] = T_Scalar2(values[i]);
      }
      VecCreateSeqWithArray(PETSC_COMM_SELF, 1, size, petscValues, vecPetsc);
      shouldBeDestroyedWithPetsc.push_back(petscValues);
    }

    // For solution
    void operator()(std::vector< T_Scalar1 > &values, T_Scalar2 *array, unsigned long long size)
    {
      msg::warning << "The solution is cast because you use a '" << scalar::Name< T_Scalar1 >::name << "' formulation with a '" << scalar::Name< T_Scalar2 >::name << "' PETSc" << msg::endl;
#pragma omp parallel for num_threads(omp::getMaxThreads())
      for(auto i = 0ULL; i < size; ++i) {
        values[i] = T_Scalar1(std::real(array[i]));
      }
    }
  };

  template< class T_Scalar1, class T_Scalar2 >
  struct PetscInterface< T_Scalar1, T_Scalar2, true, false > {
    // For matrix
    void operator()(const unsigned long long size0, const unsigned long long size1, const unsigned long long *const row, const unsigned long long *const indices, T_Scalar1 *const values, Mat *matPetsc, std::vector< void * > &shouldBeDestroyedWithPetsc, const bool isSymmetric)
    {
      msg::warning << "Loss of efficiency can be noticed if you use a '" << scalar::Name< T_Scalar1 >::name << "' formulation with a '" << scalar::Name< T_Scalar2 >::name << "' PETSc" << msg::endl;

      PetscInt *petscRow = (PetscInt *)std::malloc((size0 * 2 + 1) * sizeof(PetscInt));
      PetscInt *petscIndices = (PetscInt *)std::malloc((row[size0] * 4) * sizeof(PetscInt));
      T_Scalar2 *petscValues = (T_Scalar2 *)std::malloc((row[size0] * 4) * sizeof(T_Scalar2));

#pragma omp parallel num_threads(omp::getMaxThreads())
      {
        // row
#pragma omp for
        for(auto i = 0ULL; i <= size0; ++i) {
          petscRow[i] = 2 * row[i];
        }
#pragma omp for
        for(auto i = 1ULL; i <= size0; ++i) {
          petscRow[i + size0] = petscRow[size0] + 2 * row[i];
        }

        // indices
#pragma omp for
        for(auto i = 0ULL; i < size0; ++i) {
          for(auto j = row[i]; j < row[i + 1]; ++j) {
            petscIndices[petscRow[i] + j - row[i]] = indices[j];
          }
          for(auto j = row[i]; j < row[i + 1]; ++j) {
            petscIndices[petscRow[i] + j + row[i + 1] - 2 * row[i]] = indices[j] + size1;
          }
        }
#pragma omp for
        for(auto i = 0ULL; i < size0; ++i) {
          for(auto j = row[i]; j < row[i + 1]; ++j) {
            petscIndices[petscRow[i + size0] + j - row[i]] = indices[j];
          }
          for(auto j = row[i]; j < row[i + 1]; ++j) {
            petscIndices[petscRow[i + size0] + j + row[i + 1] - 2 * row[i]] = indices[j] + size1;
          }
        }

        // values
#pragma omp for
        for(auto i = 0ULL; i < size0; ++i) {
          for(auto j = row[i]; j < row[i + 1]; ++j) {
            petscValues[petscRow[i] + j - row[i]] = T_Scalar2(std::real(values[j]));
          }
          for(auto j = row[i]; j < row[i + 1]; ++j) {
            petscValues[petscRow[i] + j + row[i + 1] - 2 * row[i]] = (T_Scalar2(std::imag(values[j])) == 0. ? 0. : -T_Scalar2(std::imag(values[j])));
          }
        }
#pragma omp for
        for(auto i = 0ULL; i < size0; ++i) {
          for(auto j = row[i]; j < row[i + 1]; ++j) {
            petscValues[petscRow[i + size0] + j - row[i]] = T_Scalar2(std::imag(values[j]));
          }
          for(auto j = row[i]; j < row[i + 1]; ++j) {
            petscValues[petscRow[i + size0] + j + row[i + 1] - 2 * row[i]] = T_Scalar2(std::real(values[j]));
          }
        }
      }

      if(isSymmetric) {
        MatCreateSeqSBAIJWithArrays(PETSC_COMM_SELF, 1, 2 * size0, 2 * size1, petscRow, petscIndices, petscValues, matPetsc);
      }
      else {
        MatCreateSeqAIJWithArrays(PETSC_COMM_SELF, 2 * size0, 2 * size1, petscRow, petscIndices, petscValues, matPetsc);
      }
      shouldBeDestroyedWithPetsc.push_back(petscRow);
      shouldBeDestroyedWithPetsc.push_back(petscIndices);
      shouldBeDestroyedWithPetsc.push_back(petscValues);
    }

    // For vector
    void operator()(const unsigned long long size, const T_Scalar1 *const values, Vec *vecPetsc, std::vector< void * > &shouldBeDestroyedWithPetsc)
    {
      msg::warning << "Loss of efficiency can be noticed if you use a '" << scalar::Name< T_Scalar1 >::name << "' formulation with a '" << scalar::Name< T_Scalar2 >::name << "' PETSc" << msg::endl;

      T_Scalar2 *petscValues = (T_Scalar2 *)std::malloc((size * 2) * sizeof(T_Scalar2));
#pragma omp parallel for num_threads(omp::getMaxThreads())
      for(auto i = 0ULL; i < size; ++i) {
        petscValues[i] = std::real(values[i]);
        petscValues[size + i] = std::imag(values[i]);
      }

      VecCreateSeqWithArray(PETSC_COMM_SELF, 1, 2 * size, petscValues, vecPetsc);
      shouldBeDestroyedWithPetsc.push_back(petscValues);
    }

    // For solution
    void operator()(std::vector< T_Scalar1 > &values, T_Scalar2 *array, unsigned long long size)
    {
      msg::warning << "The solution is cast because you use a '" << scalar::Name< T_Scalar1 >::name << "' formulation with a '" << scalar::Name< T_Scalar2 >::name << "' PETSc" << msg::endl;
#pragma omp parallel for num_threads(omp::getMaxThreads())
      for(auto i = 0ULL; i < size / 2; ++i) {
        values[i] = T_Scalar1(array[i], array[size / 2 + i]);
      }
    }
  };


} // namespace gmshfem::system

#endif

#endif // H_GMSHFEM_PETSCINTERFACE
