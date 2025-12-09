// GmshFEM - Copyright (C) 2019-2022, A. Royer, E. Béchet, C. Geuzaine, Université de Liège
//
// See the LICENSE.txt file for license information. Please report all
// issues on https://gitlab.onelab.info/gmsh/fem/issues

#ifndef H_GMSHFEM_OMPINTERFACE
#define H_GMSHFEM_OMPINTERFACE

#include "gmshfemDefines.h"

#ifdef HAVE_OPENMP
#include <omp.h>
#endif

namespace gmshfem::omp
{


  class Lock
  {
   private:
#ifdef HAVE_OPENMP
    omp_lock_t _lock;
#endif

   public:
    Lock();
    ~Lock();

    void lock();
    void unlock();

    bool trySetLock();
  };

  unsigned int getMaxThreads();
  unsigned int getNumThreads();
  unsigned int getThreadNum();

  double getTime();

  void setMaxThreads(const unsigned int maxThreads);

  bool isInParallel();


} // namespace gmshfem::omp

#endif // H_GMSHFEM_OMPINTERFACE
