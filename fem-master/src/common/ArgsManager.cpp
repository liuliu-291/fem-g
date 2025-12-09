// GmshFEM - Copyright (C) 2019-2022, A. Royer, E. Béchet, C. Geuzaine, Université de Liège
//
// See the LICENSE.txt file for license information. Please report all
// issues on https://gitlab.onelab.info/gmsh/fem/issues

#include "ArgsManager.h"

#include "GmshFem.h"
#include "Memory.h"
#include "Message.h"
#include "OmpInterface.h"
#include "Options.h"
#include "gmshfemDefines.h"
#include "scalar.h"

#ifdef HAVE_PETSC
#include <petsc.h>
#endif // HAVE_PETSC
#include <fstream>
#include <gmsh.h>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

namespace gmshfem::common
{


  //
  // function
  //

  void help(const std::set< CommandLineInterface * > &commands)
  {
    msg::print << "./programName [ParameterFile ...] [-Command [value(s) for command ...] ...]" << msg::endl;
    msg::print << "Commands:" << msg::endl;
    for(auto it = commands.begin(); it != commands.end(); ++it) {
      //Help command
      if((*it)->type() == CommandLine::Type::Help) {
        msg::print << " -" << (*it)->name() << msg::fill(50) << ": " << (*it)->help() << msg::endl;
      }
      //Unary commands
      else if((*it)->type() == CommandLine::Type::Unary) {
        msg::print << " -" << (*it)->name() << msg::fill(50) << ": " << (*it)->help() << msg::endl;
      }
      //Binary commands
      else if((*it)->type() == CommandLine::Type::Binary) {
        msg::print << " -" << (*it)->name() << " param1(" << (*it)->getParameterTypeName() << ")" << msg::fill(50) << ": " << (*it)->help() << msg::endl;
      }
      //Ternary commands
      else if((*it)->type() == CommandLine::Type::Ternary) {
        msg::print << " -" << (*it)->name() << " param1(" << (*it)->getParameterTypeName(0) << ") param2(" << (*it)->getParameterTypeName(1) << ")" << msg::fill(50) << ": " << (*it)->help() << msg::endl;
      }
    }
  }

  //Unary commands
  void info(ArgsManager *argsManager)
  {
    msg::print << msg::endl;
    msg::print << "###################" << msg::endl;
    msg::print << "#  G m s h F E M  #" << msg::endl;
    msg::print << "###################" << msg::endl;
    msg::print << "# Version: " << GMSHFEM_MAJOR_VERSION << "." << GMSHFEM_MINOR_VERSION << "." << GMSHFEM_PATCH_VERSION << msg::endl;
    msg::print << "# Build date: " << GMSHFEM_DATE << msg::endl;
    msg::print << "# Build host: " << GMSHFEM_HOST << msg::endl;
    msg::print << "# Packager: " << GMSHFEM_PACKAGER << msg::endl;
    msg::print << "# Build OS: " << GMSHFEM_OS << msg::endl;
    msg::print << "# Build type: " << GMSHFEM_CMAKE_BUILD_TYPE << msg::endl;
    msg::print << "# Build options: " << GMSHFEM_CONFIG_OPTIONS << msg::endl;
    msg::print << "# CXX flags: " << GMSHFEM_CXX_FLAGS << msg::endl;
#ifdef HAVE_GMSH
    msg::print << "# Gmsh API version: " << GMSH_API_VERSION << msg::endl;
#endif
#ifdef HAVE_PETSC
    msg::print << "# PETSc arithmetic: " << scalar::Name< PetscScalar >::name << msg::endl;
    msg::print << "# PETSc integer: " << sizeof(PetscInt) << " bytes" << msg::endl;
#endif
    msg::print << "# Constants: " << msg::endl;
    msg::print << "# \t GMSHFEM_DOF_FIELD_OFFSET = " << GMSHFEM_DOF_FIELD_OFFSET << msg::endl;
    msg::print << "# \t GMSHFEM_DOF_MEMORY_POOL_SIZE = " << GMSHFEM_DOF_MEMORY_POOL_SIZE << msg::endl;
    msg::print << "# \t GMSHFEM_FUNCTION_MEMORY_ALIGNMENT = " << GMSHFEM_FUNCTION_MEMORY_ALIGNMENT << msg::endl;
    msg::print << "# \t GMSHFEM_FUNCTION_MEMORY_OVERSIZE_TOLERENCE = " << GMSHFEM_FUNCTION_MEMORY_OVERSIZE_TOLERENCE << msg::endl;
    msg::print << "# \t GMSHFEM_FUNCTION_MEMORY_RESIDUAL = " << Memory(GMSHFEM_FUNCTION_MEMORY_RESIDUAL) << " ( = " << GMSHFEM_FUNCTION_MEMORY_RESIDUAL << " [B])" << msg::endl;
    msg::print << "# \t GMSHFEM_FUNCTION_MEMORY_SMALLBUFFER_SIZE = " << GMSHFEM_FUNCTION_MEMORY_SMALLBUFFER_SIZE << msg::endl;
    msg::print << "# \t GMSHFEM_NUMA_ALLOCATOR = " << GMSHFEM_NUMA_ALLOCATOR << msg::endl;
    msg::print << msg::endl;
  }
  
  void license(ArgsManager *argsManager)
  {
    const std::string copyright = "Copyright (C) 2019-2022, A. Royer, E. Béchet, C. Geuzaine, Université de Liège";
    msg::print << msg::endl;
    msg::print << "GmshFEM is an efficient finite element library based on Gmsh" << msg::endl;
    msg::print << copyright << msg::endl;
    msg::print << msg::endl;
    msg::print << "This program is free software: you can redistribute it and/or modify" << msg::endl;
    msg::print << "it under the terms of the GNU Affero General Public License as" << msg::endl;
    msg::print << "published by the Free Software Foundation, either version 3 of the" << msg::endl;
    msg::print << "License, or (at your option) any later version." << msg::endl;
    msg::print << msg::endl;
    msg::print << "This program is distributed in the hope that it will be useful," << msg::endl;
    msg::print << "but WITHOUT ANY WARRANTY; without even the implied warranty of" << msg::endl;
    msg::print << "MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the" << msg::endl;
    msg::print << "GNU Affero General Public License for more details." << msg::endl;
    msg::print << msg::endl;
    msg::print << "You should have received a copy of the GNU Affero General Public License" << msg::endl;
    msg::print << "along with this program.  If not, see <https://www.gnu.org/licenses/>." << msg::endl;
  }

  void version(ArgsManager *argsManager)
  {
    msg::print << "GsmhFEM " << GMSHFEM_MAJOR_VERSION << "." << GMSHFEM_MINOR_VERSION << "." << GMSHFEM_PATCH_VERSION << msg::endl;
  }

  //Binary commands
  void dofsSort(ArgsManager *argsManager, const CommandLine::Parameter &parameters)
  {
    switch(parameters.i) {
    case 1:
      common::Options::instance()->dofsSortAlgorithm = problem::DofsSort::Algorithm::None;
      break;
    case 2:
      common::Options::instance()->dofsSortAlgorithm = problem::DofsSort::Algorithm::Hilbert;
      break;
    case 3:
      common::Options::instance()->dofsSortAlgorithm = problem::DofsSort::Algorithm::RCM;
      break;
    default: break;
    }
  }

  void elmsSort(ArgsManager *argsManager, const CommandLine::Parameter &parameters)
  {
    switch(parameters.i) {
    case 1:
      common::Options::instance()->elementsSortAlgorithm = problem::ElementsSort::Algorithm::None;
      break;
    case 2:
      common::Options::instance()->elementsSortAlgorithm = problem::ElementsSort::Algorithm::Hilbert;
      break;
    default: break;
    }
  }

#ifdef HAVE_PROFILE
  void p_lockEvents(ArgsManager *argsManager)
  {
    common::Options::instance()->ProfilerOptions.lockEvents = true;
  }
#endif

  //
  // ArgsManager
  //

  static void loadCommand(std::set< CommandLineInterface * > &commands)
  {
    UnaryHelpCL *unaryHelpCL = new UnaryHelpCL();
    commands.insert(unaryHelpCL);

    // CommandLine -> name  /  help  /  kill  /  pre  /  function

    //Unary commands
    UnaryCL *unaryCL = nullptr;
    unaryCL = new UnaryCL("debug", "Start in debug mode", false, true, [](ArgsManager *argsManager) { common::Options::instance()->debug = true; });
    commands.insert(unaryCL);
    unaryCL = new UnaryCL("info", "Get compile information about GmshFEM", true, true, info);
    commands.insert(unaryCL);
    unaryCL = new UnaryCL("interface", "Switch on Gmsh graphical interface", false, true, [](ArgsManager *argsManager) { common::Options::instance()->interface = true; });
    commands.insert(unaryCL);
    unaryCL = new UnaryCL("license", "Print the license notive", true, true, license);
    commands.insert(unaryCL);
    unaryCL = new UnaryCL("memory", "Print memory information", false, true, [](ArgsManager *argsManager) { common::Options::instance()->memory = true; });
    commands.insert(unaryCL);
    unaryCL = new UnaryCL("forceSymmetric", "Force to use symmetric matrix", false, true, [](ArgsManager *argsManager) { common::Options::instance()->forceSymmetric = true; });
    commands.insert(unaryCL);
    unaryCL = new UnaryCL("forceHermitian", "Force to use hermitian matrix", false, true, [](ArgsManager *argsManager) { common::Options::instance()->forceHermitian = true; });
    commands.insert(unaryCL);
    unaryCL = new UnaryCL("version", "Get GmshFEM version", true, true, version);
    commands.insert(unaryCL);

    //Binary commands
    BinaryCL *binaryCL = nullptr;
    binaryCL = new BinaryCL("dofsSort", "Choose the dof reordering algorithm (1:none, 2:Hilbert, 3:RCM)", false, true, CommandLine::Parameter::Type::Integer, dofsSort);
    commands.insert(binaryCL);
    binaryCL = new BinaryCL("elmsSort", "Choose the element reordering algorithm (1:none, 2:Hilbert)", false, true, CommandLine::Parameter::Type::Integer, elmsSort);
    commands.insert(binaryCL);
    binaryCL = new BinaryCL("maxThreads", "Set the maximum number of threads that can be used (given by 'param1')", false, true, CommandLine::Parameter::Type::Integer, [](ArgsManager *argsManager, const CommandLine::Parameter &parameters) { common::Options::instance()->maxThreads = parameters.i; });
    commands.insert(binaryCL);
    binaryCL = new BinaryCL("maxThreadsGmsh", "Set the maximum number of threads used by Gmsh (given by 'param1')", false, true, CommandLine::Parameter::Type::Integer, [](ArgsManager *argsManager, const CommandLine::Parameter &parameters) { gmsh::option::setNumber("General.NumThreads", parameters.i); });
    commands.insert(binaryCL);
    binaryCL = new BinaryCL("verbose", "Set the verbose level (0:nothing, 1:errors, 2:add warnings, 3:add infos, 4:add debugs)", false, true, CommandLine::Parameter::Type::Integer, [](ArgsManager *argsManager, const CommandLine::Parameter &parameters) { argsManager->gmshFem()->setVerbosity(parameters.i); });
    commands.insert(binaryCL);
  }

  static void loadCommandProfile(std::set< CommandLineInterface * > &commands)
  {
#ifdef HAVE_PROFILE
    UnaryCL *unaryCL = nullptr;
    unaryCL = new UnaryCL("p_lockEvents", "Switch on lock events Profiler", false, true, p_lockEvents);
    commands.insert(unaryCL);
#endif
  }

  void ArgsManager::_parameterFile(const std::string &path)
  {
    std::ifstream file(path);
    if(file.is_open()) {
      std::string line;
      while(std::getline(file, line)) {
        if(line[0] == '#') {
          continue;
        }
        else {
          std::string parameterName;
          std::string parameterValue;
          bool equalFound = false;
          bool firstNonWhiteFound = false;
          for(auto i = 0ULL; i < line.size(); ++i) {
            if(line[i] == '=') {
              equalFound = true;
              continue;
            }

            if(equalFound) {
              if(line[i] != ' ' && line[i] != '\t' && firstNonWhiteFound == false) {
                firstNonWhiteFound = true;
              }
              if(firstNonWhiteFound) {
                parameterValue += line[i];
              }
            }
            else {
              if(line[i] != ' ' && line[i] != '\t') {
                parameterName += line[i];
              }
            }
          }
          msg::debug << "Read parameter " << parameterName << " = " << parameterValue << msg::endl;
          if(parameterValue.size() == 0) {
            msg::warning << "Missing value for parameter '" << parameterName << "'" << msg::endl;
          }
          else {
            addUserDefinedParameter(parameterName, parameterValue);
          }
        }
      }
      file.close();
    }
    else {
      msg::warning << "Unable to open parameter file '" << path << "'" << msg::endl;
    }
  }

  ArgsManager::ArgsManager(GmshFem *gmshFem) :
    _gmshFem(gmshFem), _commands(), _command(), _userDefinedParameters(), _haveKillParam(false)
  {
    loadCommand(_commands);
    loadCommandProfile(_commands);
  }

  ArgsManager::~ArgsManager()
  {
    for(auto it = _commands.begin(); it != _commands.end(); ++it) {
      delete *it;
    }

    for(auto it = _command.begin(); it != _command.end(); ++it) {
      delete *it;
    }
  }

  void ArgsManager::setArgcArgv(int argc, char **argv)
  {
    _argc = argc;
    _argv = argv;
  }

  void ArgsManager::processArgcArgv()
  {
    bool firstParameterFound = false;

    for(auto i = 1; i < _argc; ++i) {
      if(_argv[i][0] == '-') {
        firstParameterFound = true;
        SearchCL searchCL(&_argv[i][1]);
        auto it = _commands.find(&searchCL);
        if(it != _commands.end()) {
          CommandLineInterface *currentCL = (*it)->copy();
          bool CLignore = false;
          for(auto j = 0U; j < currentCL->nbrParameters(); ++j) {
            CommandLine::Parameter parameter;
            if(currentCL->getParameterType(j) == CommandLine::Parameter::Type::Real) {
              if(i + j + 1 < static_cast< unsigned int >(_argc)) {
                if(_argv[i + j + 1][0] == '-' && (_argv[i + j + 1][1] < 47 || _argv[i + j + 1][1] > 57)) {
                  msg::warning << "Missing parameter for command '" << _argv[i] << "'" << msg::endl;
                  msg::warning << "Command ignored" << msg::endl;
                  CLignore = true;
                }
                else {
                  parameter.d = std::stod(std::string(_argv[i + j + 1]));
                }
              }
              else {
                msg::warning << "Missing parameter for command '" << _argv[i] << "'" << msg::endl;
                msg::warning << "Command ignored" << msg::endl;
                CLignore = true;
              }
            }
            else if(currentCL->getParameterType(j) == CommandLine::Parameter::Type::Integer) {
              if(i + j + 1 < static_cast< unsigned int >(_argc)) {
                if(_argv[i + j + 1][0] == '-' && (_argv[i + j + 1][1] < 47 || _argv[i + j + 1][1] > 57)) {
                  msg::warning << "Missing parameter for command '" << _argv[i] << "'" << msg::endl;
                  msg::warning << "Command ignored" << msg::endl;
                  CLignore = true;
                }
                else {
                  parameter.i = std::stoi(std::string(_argv[i + j + 1]));
                }
              }
              else {
                msg::warning << "Missing parameter for command '" << _argv[i] << "'" << msg::endl;
                msg::warning << "Command ignored" << msg::endl;
                CLignore = true;
              }
            }
            else if(currentCL->getParameterType(j) == CommandLine::Parameter::Type::String) {
              if(i + j + 1 < static_cast< unsigned int >(_argc)) {
                if(_argv[i + j + 1][0] == '-') {
                  msg::warning << "Missing parameter for command '" << _argv[i] << "'" << msg::endl;
                  msg::warning << "Command ignored" << msg::endl;
                  CLignore = true;
                }
                else {
                  parameter.s = std::string(_argv[i + j + 1]);
                }
              }
              else {
                msg::warning << "Missing parameter for command '" << _argv[i] << "'" << msg::endl;
                msg::warning << "Command ignored" << msg::endl;
                CLignore = true;
              }
            }
            currentCL->setParameter(parameter, j);
          }
          if(!CLignore) {
            _command.push_back(currentCL);
            _haveKillParam |= currentCL->kill();
            i += currentCL->nbrParameters();
          }
          else {
            i++;
          }
        }
        else {
          addUserDefinedParameter(std::string(&_argv[i][1]), (i + 1 < _argc ? std::string(&_argv[i + 1][0]) : ""));
        }
      }
      else {
        if(!firstParameterFound) {
          _parameterFile(&_argv[i][0]);
        }
      }
    }
  }

  void ArgsManager::insert(CommandLineInterface *newCL)
  {
    SearchCL searchCL(newCL->name());
    auto it = _commands.find(&searchCL);
    if(it == _commands.end()) {
      _commands.insert(newCL);
    }
    else {
      throw common::Exception("Command '" + newCL->name() + "' already exists: please choose another name");
    }
  }

  void ArgsManager::pre()
  {
    for(auto i = 0ULL; i < _command.size(); ++i) {
      if(_command[i]->type() == CommandLine::Type::Help) {
        help(_commands);
        exit(0);
      }
    }

    for(auto i = 0ULL; i < _command.size(); ++i) {
      if((_command[i]->kill() && _haveKillParam) || !_haveKillParam) {
        if(_command[i]->pre()) {
          _command[i]->run(this);
          if(_haveKillParam) exit(0);
        }
      }
    }
  }

  void ArgsManager::post()
  {
    for(auto i = 0ULL; i < _command.size(); ++i) {
      if((_command[i]->kill() && _haveKillParam) || !_haveKillParam) {
        if(_command[i]->post()) {
          _command[i]->run(this);
          if(_haveKillParam) exit(0);
        }
      }
    }
  }

  bool ArgsManager::haveKillParam() const
  {
    return _haveKillParam;
  }

  GmshFem *ArgsManager::gmshFem() const
  {
    return _gmshFem;
  }

  static bool checkParamameterName(const std::string &name)
  {
    if(name[0] == '#') {
      msg::warning << "A user-defined parameter name cannot starts with a '#'" << msg::endl;
      return false;
    }
    for(auto i = 0ULL; i < name.size(); ++i) {
      if(name[i] == '=') {
        msg::warning << "A user-defined parameter name cannot contains an equal character ('=')" << msg::endl;
        return false;
      }
      else if(name[i] == ' ') {
        msg::warning << "A user-defined parameter name cannot contains a whitespace character" << msg::endl;
        return false;
      }
      else if(name[i] == '\t') {
        msg::warning << "A user-defined parameter name cannot contains a tab character" << msg::endl;
        return false;
      }
    }

    return true;
  }

  void ArgsManager::addUserDefinedParameter(const std::string &name, const std::string &value)
  {
    auto it = _userDefinedParameters.find(name);
    if(it != _userDefinedParameters.end()) {
      msg::warning << "User defined parameter '" << name << "' is already defined: ignoring" << msg::endl;
    }
    else {
      _userDefinedParameters.insert(std::make_pair(name, value));
    }
  }

  bool ArgsManager::userDefinedParameter(int &value, const std::string &name) const
  {
    if(!checkParamameterName(name)) {
      return false;
    }

    auto it = _userDefinedParameters.find(name);
    bool parameterFound;
    if(it != _userDefinedParameters.end()) {
      parameterFound = true;
      value = std::stoi(it->second);
    }
    else {
      parameterFound = false;
      msg::debug << "User-defined parameter: " << name << " cannot be found. The default value: " << value << " is used" << msg::endl;
    }
    return parameterFound;
  }

  bool ArgsManager::userDefinedParameter(double &value, const std::string &name) const
  {
    if(!checkParamameterName(name)) {
      return false;
    }

    auto it = _userDefinedParameters.find(name);
    bool parameterFound;
    if(it != _userDefinedParameters.end()) {
      parameterFound = true;
      value = std::stod(it->second);
    }
    else {
      parameterFound = false;
      msg::debug << "User-defined parameter: " << name << " cannot be found. The default value: " << value << " is used" << msg::endl;
    }
    return parameterFound;
  }

  bool ArgsManager::userDefinedParameter(std::string &value, const std::string &name) const
  {
    if(!checkParamameterName(name)) {
      return false;
    }

    auto it = _userDefinedParameters.find(name);
    bool parameterFound;
    if(it != _userDefinedParameters.end()) {
      parameterFound = true;
      value = it->second;
    }
    else {
      parameterFound = false;
      msg::debug << "User-defined parameter: " << name << " cannot be found. The default value: " << value << " is used" << msg::endl;
    }
    return parameterFound;
  }

  bool ArgsManager::userDefinedParameter(bool &value, const std::string &name) const
  {
    if(!checkParamameterName(name)) {
      return false;
    }

    auto it = _userDefinedParameters.find(name);
    bool parameterFound;
    if(it != _userDefinedParameters.end()) {
      parameterFound = true;
      value = !value;
    }
    else {
      parameterFound = false;
      msg::debug << "User-defined parameter: " << name << " cannot be found. The default value: " << value << " is used" << msg::endl;
    }
    return parameterFound;
  }

  std::vector< CommandLineInterface * >::const_iterator ArgsManager::begin() const
  {
    return _command.begin();
  }

  std::vector< CommandLineInterface * >::const_iterator ArgsManager::end() const
  {
    return _command.end();
  }


} // namespace gmshfem::common
