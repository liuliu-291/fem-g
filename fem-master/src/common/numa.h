// GmshFEM - Copyright (C) 2019-2022, A. Royer, E. Béchet, C. Geuzaine, Université de Liège
//
// See the LICENSE.txt file for license information. Please report all
// issues on https://gitlab.onelab.info/gmsh/fem/issues

#ifndef H_GMSHFEM_NUMA
#define H_GMSHFEM_NUMA

#include "OmpInterface.h"
#include "gmshfemDefines.h"
#include "scalar.h"

#include <memory>
#include <vector>

namespace gmshfem::numa
{


#if GMSHFEM_NUMA_ALLOCATOR == 1
  template< class T >
  class allocator : private std::allocator< T >
  {
   public:
    typedef T value_type;
    typedef std::size_t size_type;
    typedef std::ptrdiff_t difference_type;
    typedef std::true_type propagate_on_container_move_assignment;
    typedef std::true_type is_always_equal;

   public:
    allocator() noexcept {}
    allocator(const allocator &other) noexcept {}
    ~allocator() {}

    T *allocate(size_type n)
    {
      T *p = std::allocator< T >::allocate(n);

      // first touch
#pragma omp parallel for num_threads(omp::getMaxThreads())
      for(unsigned long long i = 0; i < n; ++i) {
        *reinterpret_cast< char * >(p + i) = 0.;
      }

      return p;
    }

    void deallocate(T *p, size_type n)
    {
      std::allocator< T >::deallocate(p, n);
    }
  };

  template< class T_Scalar1, class T_Scalar2 >
  void copy(std::vector< T_Scalar1, allocator< T_Scalar1 > > &v1, std::vector< T_Scalar2 > &v2)
  {
    v1.resize(v2.size());
#pragma omp parallel for num_threads(omp::getMaxThreads())
    for(auto i = 0ULL; i < v2.size(); ++i) {
      v1[i] = static_cast< T_Scalar1 >(v2[i]);
    }
  }
#else
  template< class T >
  using allocator = std::allocator< T >;

  template< class T_Scalar1, class T_Scalar2 >
  void copy(std::vector< T_Scalar1, allocator< T_Scalar1 > > &v1, std::vector< T_Scalar2 > &v2)
  {
    scalar::move(v1, v2);
  }
#endif // GMSHFEM_NUMA_ALLOCATOR

} // namespace gmshfem::numa

#endif // H_GMSHFEM_NUMA
