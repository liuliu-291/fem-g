// GmshFEM - Copyright (C) 2019-2022, A. Royer, E. Béchet, C. Geuzaine, Université de Liège
//
// See the LICENSE.txt file for license information. Please report all
// issues on https://gitlab.onelab.info/gmsh/fem/issues

#ifndef H_GMSHFEM_OPTIONS
#define H_GMSHFEM_OPTIONS

#include "gmshfemDefines.h"
#include "optionsEnums.h"

#include <string>

namespace gmshfem::common
{


  // Singleton
  class Options
  {
   private:
    static Options *_instance;

   protected:
    Options();
    ~Options();
    Options(const Options &other) = delete;
    Options(Options &&other) = delete;

    Options &operator=(const Options &other) = delete;
    Options &operator=(Options &&other) = delete;

    void _setDefault();

   public:
    bool debug;
    problem::DofsSort::Algorithm dofsSortAlgorithm;
    problem::ElementsSort::Algorithm elementsSortAlgorithm;
    bool interface;
    unsigned int maxThreads;
    bool memory;
    bool forceHermitian;
    bool forceSymmetric;
    int verbose;

   public:
    static Options *instance();
    static void destroy();
  };


} // namespace gmshfem::common

#endif // H_GMSHFEM_OPTIONS
