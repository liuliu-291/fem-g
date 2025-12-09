// GmshFEM - Copyright (C) 2019-2022, A. Royer, E. Béchet, C. Geuzaine, Université de Liège
//
// See the LICENSE.txt file for license information. Please report all
// issues on https://gitlab.onelab.info/gmsh/fem/issues

#include "FunctionAllocator.h"

#include "Message.h"
#include "gmshfemDefines.h"
#include "instantiate.h"

#include <cstdlib>

#ifdef HAVE_WIN32
#define ALIGNED_ALLOC(ALIGNMENT, SIZE) (_aligned_malloc(SIZE, ALIGNMENT))
#define ALIGNED_FREE(PTR) (_aligned_free(PTR))
#else
#define ALIGNED_ALLOC(ALIGNMENT, SIZE) (std::aligned_alloc(ALIGNMENT, SIZE))
#define ALIGNED_FREE(PTR) (std::free(PTR))
#endif

namespace gmshfem::function
{


  //
  // class Pool
  //

  Pool::Pool(const std::size_t size) :
    _next(nullptr), _size(size), _buffer(ALIGNED_ALLOC(GMSHFEM_FUNCTION_MEMORY_ALIGNMENT, size)), _isFree(false)
  {
    if(_buffer == nullptr) {
      *(const_cast< void ** >(&_buffer)) = std::malloc(size);
    }

    if(_buffer == nullptr) {
      throw std::bad_alloc();
    }
  }

  Pool::~Pool()
  {
    ALIGNED_FREE(_buffer);
  }

  Pool *Pool::getNext() const
  {
    return _next;
  }

  void Pool::setNext(Pool *pool)
  {
    _next = pool;
  }

  std::size_t Pool::getSize() const
  {
    return _size;
  }

  void *Pool::getBuffer() const
  {
    return _buffer;
  }

  bool Pool::getState() const
  {
    return _isFree;
  }

  void Pool::setState(const bool state)
  {
    _isFree = state;
  }

  void Pool::printInfo() const
  {
    msg::info << "Pool: size = " << _size << ", state = " << (_isFree ? "free" : "not free") << msg::endl;
  }


  //
  // class MemoryPoolAllocator
  //

  MemoryPoolAllocator *MemoryPoolAllocator::_instance = nullptr;

  MemoryPoolAllocator::MemoryPoolAllocator() :
    _pool(nullptr), _size(0)
  {
    _smallBufferPtr = static_cast< unsigned char * >(ALIGNED_ALLOC(GMSHFEM_FUNCTION_MEMORY_ALIGNMENT, GMSHFEM_FUNCTION_MEMORY_SMALLBUFFER_SIZE * GMSHFEM_FUNCTION_MEMORY_ALIGNMENT));
    for(auto i = 0; i < GMSHFEM_FUNCTION_MEMORY_SMALLBUFFER_SIZE; ++i) {
      _smallBufferState[i] = true;
      _smallBuffer[i] = _smallBufferPtr + i * GMSHFEM_FUNCTION_MEMORY_ALIGNMENT;
    }
  }

  MemoryPoolAllocator::~MemoryPoolAllocator()
  {
    std::vector< Pool * > allPools;
    Pool *currentPool = _pool;

    while(currentPool != nullptr) {
      allPools.push_back(currentPool);
      currentPool = currentPool->getNext();
    }

    for(auto it = allPools.begin(); it != allPools.end(); ++it) {
      delete *it;
    }

    ALIGNED_FREE(_smallBufferPtr);
  }

  MemoryPoolAllocator *MemoryPoolAllocator::instance()
  {
    if(_instance == nullptr) _instance = new MemoryPoolAllocator();

    return _instance;
  }

  void MemoryPoolAllocator::destroy()
  {
    if(_instance != nullptr) {
      delete _instance;
      _instance = nullptr;
    }
  }

  void *MemoryPoolAllocator::allocate(std::size_t n)
  {
    Pool *currentPool = _pool;
    Pool *lastOfLessSize = nullptr;

    while(currentPool != nullptr) {
      if(currentPool->getSize() > GMSHFEM_FUNCTION_MEMORY_OVERSIZE_TOLERENCE * n) {
        break;
      }
      else {
        if(currentPool->getSize() >= n) {
          if(currentPool->getState()) {
            currentPool->setState(false);
            return currentPool->getBuffer();
          }
        }
        else {
          lastOfLessSize = currentPool;
        }
      }

      currentPool = currentPool->getNext();
    }

    Pool *newPool = new Pool(n);
    _size += n;
    if(lastOfLessSize == nullptr) {
      newPool->setNext(_pool);
      _pool = newPool;
    }
    else {
      newPool->setNext(lastOfLessSize->getNext());
      lastOfLessSize->setNext(newPool);
    }

    return newPool->getBuffer();
  }

  void *MemoryPoolAllocator::allocateSmall(std::size_t n)
  {
    for(auto i = 0; i < GMSHFEM_FUNCTION_MEMORY_SMALLBUFFER_SIZE; ++i) {
      if(_smallBufferState[i]) {
        _smallBufferState[i] = false;
        return _smallBuffer[i];
      }
    }
    return std::malloc(n);
  }

  void MemoryPoolAllocator::deallocate(void *p, std::size_t n)
  {
    Pool *currentPool = _pool;
    while(currentPool != nullptr) {
      if(currentPool->getBuffer() == p) {
        currentPool->setState(true);
        break;
      }
      currentPool = currentPool->getNext();
    }
  }

  void MemoryPoolAllocator::deallocateSmall(void *p, std::size_t n)
  {
    for(auto i = 0; i < GMSHFEM_FUNCTION_MEMORY_SMALLBUFFER_SIZE; ++i) {
      if(_smallBuffer[i] == p) {
        _smallBufferState[i] = true;
        return;
      }
    }
    std::free(p);
  }

  void MemoryPoolAllocator::shrinkTo(const std::size_t size)
  {
    if(_size > size) {
      std::vector< Pool * > allPools;
      Pool *currentPool = _pool;
      while(currentPool != nullptr) {
        allPools.push_back(currentPool);
        currentPool = currentPool->getNext();
      }

      // step 1
      unsigned int count = 0;
      std::size_t currentSize = 0;
      for(auto it = allPools.begin(); it != allPools.end(); ++it) {
        if(currentSize != (*it)->getSize()) {
          if(count > 2) {
            auto it2 = it - 1;
            while(count > 2 && it2 != allPools.begin()) {
              if((*it2)->getState()) {
                _size -= (*it2)->getSize();
                delete *it2;
                *it2 = nullptr;
              }
              --count;
              --it2;
            }
            count = 0;
            if(_size <= size) {
              break;
            }
          }
          currentSize = (*it)->getSize();
        }
        count++;
      }

      // step 2
      if(_size > size) {
        for(auto it = allPools.rbegin(); it != allPools.rend(); ++it) {
          if((*it) != nullptr) {
            if((*it)->getState()) {
              _size -= (*it)->getSize();
              delete *it;
              *it = nullptr;
              if(_size <= size) {
                break;
              }
            }
          }
        }
      }

      // rebuild list
      _pool = nullptr;
      for(auto it = allPools.begin(); it != allPools.end(); ++it) {
        if(_pool == nullptr && (*it) != nullptr) {
          _pool = *it;
        }
        if((*it) != nullptr) {
          (*it)->setNext(nullptr);
        }
      }
      for(auto it = allPools.begin(); it != allPools.end(); ++it) {
        for(auto it2 = it; it2 != allPools.end(); ++it2) {
          if((*it2) != nullptr && it2 != it) {
            (*it)->setNext((*it2));
            it = --it2;
            break;
          }
        }
      }
    }
  }


} // namespace gmshfem::function
