// GmshFEM - Copyright (C) 2019-2022, A. Royer, E. Béchet, C. Geuzaine, Université de Liège
//
// See the LICENSE.txt file for license information. Please report all
// issues on https://gitlab.onelab.info/gmsh/fem/issues

#include "GmshFem.h"

#include "Exception.h"
#include "FunctionAllocator.h"
#include "Message.h"
#include "OmpInterface.h"
#include "Options.h"
#include "Timer.h"
#include "gmshfemDefines.h"

#include <gmsh.h>

#ifdef HAVE_PETSC
#include <petsc.h>
#ifdef HAVE_SLEPC
#include <slepc.h>
#endif // HAVE_SLEPC
#endif // HAVE_PETSC

namespace gmshfem::common
{


  static GmshFem *s_gmshfem = nullptr;

  GmshFem::GmshFem(const std::string &userLibName) :
    _argsManager(this), _userLibName(userLibName)
  {
    if(s_gmshfem == nullptr) {
      s_gmshfem = this;
    }
    else {
      throw common::Exception("Only one object of class GmshFem can be instantiated");
    }
  }

  GmshFem::GmshFem(int argc, char **argv) :
    _argsManager(this), _userLibName("")
  {
    if(s_gmshfem == nullptr) {
      s_gmshfem = this;
      finalizeInit(argc, argv);
    }
    else {
      throw common::Exception("Only one object of class GmshFem can be instantiated");
    }
  }

  GmshFem::~GmshFem()
  {
    _argsManager.post();

    if(common::Options::instance()->interface) {
      gmsh::fltk::run();
    }

    common::Options::destroy();
    function::MemoryPoolAllocator::destroy();
    gmsh::finalize();

#ifdef HAVE_SLEPC
    //Finalize SLEPc and PETSc
    SlepcFinalize();
#else
#ifdef HAVE_PETSC
    //Finalize PETSc
    PetscFinalize();
#endif // HAVE_PETSC
#endif // HAVE_SLEPC
  }

  void GmshFem::finalizeInit(int argc, char **argv)
  {
    _argsManager.setArgcArgv(argc, argv);
    _argsManager.processArgcArgv();
    if(_argsManager.haveKillParam()) _argsManager.pre();

    //Initialize gmsh
    gmsh::initialize();

#ifdef HAVE_SLEPC
    //Initialize SLEPc and PETSc
    SlepcInitialize(&argc, &argv, nullptr, nullptr);
#else
#ifdef HAVE_PETSC
    //Initialize PETSc
    PetscInitialize(&argc, &argv, nullptr, nullptr);
#endif // HAVE_PETSC
#endif // HAVE_SLEPC

    Eigen::initParallel();

    _argsManager.pre();

    if(common::Options::instance()->interface) {
      gmsh::fltk::initialize();
      gmsh::option::setNumber("General.Terminal", 0);
    }
    else {
      gmsh::option::setNumber("General.Terminal", 1);
    }
    function::MemoryPoolAllocator::instance();

    std::string name;
    if(_userLibName == "") {
      name = "GmshFEM(";
    }
    else {
      name = _userLibName + " (GmshFEM) (";
    }
    msg::print << name << omp::getMaxThreads() << (omp::getMaxThreads() == 1 ? " thread):" : " threads):") << msg::endl;
    msg::print << "# " << today() << ", " << hour() << msg::endl;
    msg::print << "# ";
    for(auto i = 0; i < argc; ++i) {
      msg::print << argv[i] << " ";
    }
    msg::print << msg::endl;
  }

  void GmshFem::setVerbosity(const int verbose)
  {
    common::Options::instance()->verbose = verbose;
    switch(verbose) {
    case 0: gmsh::option::setNumber("General.Verbosity", 0.); break;
    case 1: gmsh::option::setNumber("General.Verbosity", 1.); break;
    case 2: gmsh::option::setNumber("General.Verbosity", 2.); break;
    case 3: gmsh::option::setNumber("General.Verbosity", 4.); break;
    case 4: gmsh::option::setNumber("General.Verbosity", 99.); break;
    default: gmsh::option::setNumber("General.Verbosity", 4.); break;
    }
  }

  void GmshFem::setMaxThreads(const int maxThreads)
  {
    common::Options::instance()->maxThreads = maxThreads;
  }

  bool GmshFem::userDefinedParameter(int &value, const std::string &name) const
  {
    return _argsManager.userDefinedParameter(value, name);
  }

  bool GmshFem::userDefinedParameter(unsigned int &value, const std::string &name) const
  {
    int valueInt = static_cast< int >(value);
    bool parameterFound = _argsManager.userDefinedParameter(valueInt, name);
    value = static_cast< unsigned int >(valueInt);
    return parameterFound;
  }

  bool GmshFem::userDefinedParameter(double &value, const std::string &name) const
  {
    return _argsManager.userDefinedParameter(value, name);
  }

  bool GmshFem::userDefinedParameter(float &value, const std::string &name) const
  {
    double valueDouble = static_cast< double >(value);
    bool parameterFound = _argsManager.userDefinedParameter(valueDouble, name);
    value = static_cast< float >(valueDouble);
    return parameterFound;
  }

  bool GmshFem::userDefinedParameter(std::string &value, const std::string &name) const
  {
    return _argsManager.userDefinedParameter(value, name);
  }

  bool GmshFem::userDefinedParameter(char *value, const std::string &name) const
  {
    std::string tmp(value);
    bool parameterFound = _argsManager.userDefinedParameter(tmp, name);
    std::strcpy(value, tmp.c_str());
    return parameterFound;
  }

  bool GmshFem::userDefinedParameter(bool &value, const std::string &name) const
  {
    return _argsManager.userDefinedParameter(value, name);
  }

  std::string GmshFem::version() const
  {
    std::string version = std::to_string(GMSHFEM_MAJOR_VERSION) + "." + std::to_string(GMSHFEM_MINOR_VERSION) + "." + std::to_string(GMSHFEM_PATCH_VERSION);
    return version;
  }

  std::string GmshFem::date() const
  {
    return GMSHFEM_DATE;
  }

  std::string GmshFem::host() const
  {
    return GMSHFEM_HOST;
  }

  std::string GmshFem::packager() const
  {
    return GMSHFEM_PACKAGER;
  }

  std::string GmshFem::os() const
  {
    return GMSHFEM_OS;
  }

  std::string GmshFem::configOptions() const
  {
    return GMSHFEM_CONFIG_OPTIONS;
  }

  std::string GmshFem::cxxFalgs() const
  {
    return GMSHFEM_CXX_FLAGS;
  }

  ArgsManager *GmshFem::getArgsManager()
  {
    return &_argsManager;
  }

  // static method
  GmshFem *GmshFem::currentInstance()
  {
    return s_gmshfem;
  }


} // namespace gmshfem::common
