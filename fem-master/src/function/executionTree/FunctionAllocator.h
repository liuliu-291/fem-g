// GmshFEM - Copyright (C) 2019-2022, A. Royer, E. Béchet, C. Geuzaine, Université de Liège
//
// See the LICENSE.txt file for license information. Please report all
// issues on https://gitlab.onelab.info/gmsh/fem/issues

#ifndef H_GMSHFEM_FUNCTIONALLOCATOR
#define H_GMSHFEM_FUNCTIONALLOCATOR

#include "MathObject.h"
#include "gmshfemDefines.h"

#include <vector>

namespace gmshfem::function
{


  class Pool
  {
   private:
    Pool *_next;
    const std::size_t _size;
    void *const _buffer;
    bool _isFree;

   public:
    Pool(const std::size_t size);
    ~Pool();

    Pool *getNext() const;
    void setNext(Pool *pool);

    std::size_t getSize() const;

    void *getBuffer() const;

    bool getState() const;
    void setState(const bool state);

    void printInfo() const;
  };

  // Singleton
  class MemoryPoolAllocator
  {
   private:
    static MemoryPoolAllocator *_instance;
    Pool *_pool;
    std::size_t _size;
    bool _smallBufferState[GMSHFEM_FUNCTION_MEMORY_SMALLBUFFER_SIZE];
    unsigned char *_smallBuffer[GMSHFEM_FUNCTION_MEMORY_SMALLBUFFER_SIZE];
    unsigned char *_smallBufferPtr;

   protected:
    MemoryPoolAllocator();
    ~MemoryPoolAllocator();
    MemoryPoolAllocator(const MemoryPoolAllocator &other) = delete;
    MemoryPoolAllocator(MemoryPoolAllocator &&other) = delete;

    MemoryPoolAllocator &operator=(const MemoryPoolAllocator &other) = delete;
    MemoryPoolAllocator &operator=(MemoryPoolAllocator &&other) = delete;

   public:
    static MemoryPoolAllocator *instance();
    static void destroy();

    void *allocate(std::size_t n);
    void *allocateSmall(std::size_t n);
    void deallocate(void *p, std::size_t n);
    void deallocateSmall(void *p, std::size_t n);

    void shrinkTo(const std::size_t size);
  };


  template< class T >
  class FunctionAllocator
  {
   public:
    typedef T value_type;
    typedef std::size_t size_type;
    typedef std::ptrdiff_t difference_type;
    typedef std::true_type propagate_on_container_move_assignment;
    typedef std::true_type is_always_equal;

   private:
    T *const _wrapper;

   public:
    FunctionAllocator() noexcept :
      _wrapper(nullptr) {}
    FunctionAllocator(T *data) noexcept :
      _wrapper(data) {}
    FunctionAllocator(const FunctionAllocator &other) noexcept :
      _wrapper(other._wrapper) {}
    ~FunctionAllocator() {}

    T *allocate(size_type n)
    {
      if(_wrapper) {
        return _wrapper;
      }

      if(std::numeric_limits< size_type >::max() / sizeof(value_type) < n) {
        throw std::bad_array_new_length();
      }
      T *p = nullptr;
      if(n * sizeof(value_type) <= GMSHFEM_FUNCTION_MEMORY_ALIGNMENT) {
        p = static_cast< T * >(MemoryPoolAllocator::instance()->allocateSmall(n * sizeof(value_type)));
      }
      else {
        p = static_cast< T * >(MemoryPoolAllocator::instance()->allocate(n * sizeof(value_type)));
      }

      return p;
    }

    void deallocate(T *p, size_type n)
    {
      if(_wrapper) {
        return;
      }

      if(n * sizeof(value_type) <= GMSHFEM_FUNCTION_MEMORY_ALIGNMENT) {
        MemoryPoolAllocator::instance()->deallocateSmall(p, n * sizeof(value_type));
      }
      else {
        MemoryPoolAllocator::instance()->deallocate(p, n * sizeof(value_type));
      }
    }

    static void initalizeNUMAMemory(size_type n)
    {
      if(std::numeric_limits< size_type >::max() / sizeof(value_type) < n) {
        return;
      }

      if(n * sizeof(value_type) <= GMSHFEM_FUNCTION_MEMORY_ALIGNMENT) {
        return;
      }
      else {
        static T *s_p = nullptr;
#pragma omp master
        s_p = static_cast< T * >(MemoryPoolAllocator::instance()->allocate(n * sizeof(value_type)));
#pragma omp barrier
#pragma omp for
        for(size_type i = 0; i < n; ++i) {
          *reinterpret_cast< char * >(s_p + i) = 0.;
        }
#pragma omp barrier
#pragma omp master
        MemoryPoolAllocator::instance()->deallocate(s_p, n * sizeof(value_type));
#pragma omp barrier
      }
    }
  };

  template< class T1, class T2 >
  bool operator==(const FunctionAllocator< T1 > &lhs, const FunctionAllocator< T2 > &rhs) noexcept
  {
    return true;
  }

  template< class T1, class T2 >
  bool operator!=(const FunctionAllocator< T1 > &lhs, const FunctionAllocator< T2 > &rhs) noexcept
  {
    return false;
  }


  template< class T_Scalar, Degree T_Degree >
  using OutputVector = std::vector< typename MathObject< T_Scalar, T_Degree >::Object, FunctionAllocator< typename MathObject< T_Scalar, T_Degree >::Object > >;


} // namespace gmshfem::function


#endif // H_GMSHFEM_FUNCTIONALLOCATOR
