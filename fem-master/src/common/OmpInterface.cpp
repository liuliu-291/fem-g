// GmshFEM - Copyright (C) 2019-2022, A. Royer, E. Béchet, C. Geuzaine, Université de Liège
//
// See the LICENSE.txt file for license information. Please report all
// issues on https://gitlab.onelab.info/gmsh/fem/issues

#include "OmpInterface.h"

#include "Options.h"

#ifndef HAVE_OPENMP
#include <chrono>
#endif


namespace gmshfem::omp
{


#ifdef HAVE_OPENMP
  Lock::Lock() :
    _lock()
  {
    omp_init_lock(&_lock);
#else
  Lock::Lock(){
#endif
  }

  Lock::~Lock()
  {
#ifdef HAVE_OPENMP
    omp_destroy_lock(&_lock);
#endif
  }

  void Lock::lock()
  {
#ifdef HAVE_OPENMP
    omp_set_lock(&_lock);
#endif
  }

  void Lock::unlock()
  {
#ifdef HAVE_OPENMP
    omp_unset_lock(&_lock);
#endif
  }

  bool Lock::trySetLock()
  {
#ifdef HAVE_OPENMP
    return omp_test_lock(&_lock);
#else
    return true;
#endif
  }

  unsigned int getMaxThreads()
  {
#ifdef HAVE_OPENMP
    if(common::Options::instance()->maxThreads == 0) {
      return omp_get_max_threads();
    }
    return common::Options::instance()->maxThreads;
#endif
    return 1;
  }

  unsigned int getNumThreads()
  {
#ifdef HAVE_OPENMP
    return omp_get_num_threads();
#endif
    return 1;
  }

  unsigned int getThreadNum()
  {
#ifdef HAVE_OPENMP
    return omp_get_thread_num();
#endif
    return 0;
  }

  double getTime()
  {
#ifdef HAVE_OPENMP
    return omp_get_wtime();
#else
    auto currentTime = std::chrono::system_clock::now();
    auto durationInSeconds = std::chrono::duration< double >(currentTime.time_since_epoch());
    return durationInSeconds.count();
#endif
  }

  void setMaxThreads(const unsigned int maxThreads)
  {
#ifdef HAVE_OPENMP
    omp_set_num_threads(maxThreads == 0 ? 1 : maxThreads);
#endif
  }

  bool isInParallel()
  {
#ifdef HAVE_OPENMP
    return omp_in_parallel();
#endif
    return false;
  }


} // namespace gmshfem::omp
