// GmshFEM - Copyright (C) 2019-2022, A. Royer, E. Béchet, C. Geuzaine, Université de Liège
//
// See the LICENSE.txt file for license information. Please report all
// issues on https://gitlab.onelab.info/gmsh/fem/issues

#include "Options.h"

namespace gmshfem::common
{


  // protected
  Options::Options()
  {
    _setDefault();
  }

  // protected
  Options::~Options()
  {
  }

  // protected
  void Options::_setDefault()
  {
    debug = false;
    dofsSortAlgorithm = problem::DofsSort::Algorithm::Default;
    elementsSortAlgorithm = problem::ElementsSort::Algorithm::None;
    interface = false;
    maxThreads = 0;
    memory = false;
    forceHermitian = false;
    forceSymmetric = false;
    verbose = 3;
  }

  Options *Options::_instance = nullptr;

  Options *Options::instance()
  {
    if(_instance == nullptr) _instance = new Options();

    return _instance;
  }

  void Options::destroy()
  {
    if(_instance != nullptr) {
      delete _instance;
      _instance = nullptr;
    }
  }


} // namespace gmshfem::common
