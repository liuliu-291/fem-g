// GmshFEM - Copyright (C) 2019-2022, A. Royer, E. Béchet, C. Geuzaine, Université de Liège
//
// See the LICENSE.txt file for license information. Please report all
// issues on https://gitlab.onelab.info/gmsh/fem/issues

#include "CommandLine.h"

#include "Exception.h"

namespace gmshfem::common
{


  //
  // class CommandLineInterface
  //

  CommandLineInterface::CommandLineInterface(const std::string &name, const std::string &help, const bool kill, const bool pre, const std::vector< CommandLine::Parameter::Type > &parametersType) :
    _name(name), _help(help), _kill(kill), _pre(pre), _parametersType(parametersType), _parameters()
  {
  }

  CommandLineInterface::~CommandLineInterface() {}

  const std::string &CommandLineInterface::name() const
  {
    return _name;
  }

  const std::string &CommandLineInterface::help() const
  {
    return _help;
  }

  bool CommandLineInterface::kill() const
  {
    return _kill;
  }

  bool CommandLineInterface::pre() const
  {
    return _pre;
  }

  bool CommandLineInterface::post() const
  {
    return !_pre;
  }

  CommandLine::Parameter::Type CommandLineInterface::getParameterType(const unsigned int i)
  {
    if(_parametersType.size() <= i) {
      throw common::Exception("Trying to access parameter type number " + std::to_string(i) + " while there are only " + std::to_string(_parametersType.size()));
    }

    return _parametersType[i];
  }

  std::string CommandLineInterface::getParameterTypeName(const unsigned int i)
  {
    if(_parametersType.size() <= i) {
      throw common::Exception("Trying to access parameter type number " + std::to_string(i) + " while there are only " + std::to_string(_parametersType.size()));
    }

    switch(_parametersType[i]) {
    case CommandLine::Parameter::Type::Real: return "real"; break;
    case CommandLine::Parameter::Type::Integer: return "integer"; break;
    case CommandLine::Parameter::Type::String: return "string"; break;
    default: break;
    }

    return "";
  }

  void CommandLineInterface::setParameter(const CommandLine::Parameter &parameter,
                                          const unsigned int i)
  {
    if(_parametersType.size() <= i) {
      throw common::Exception("Trying to define parameter number " + std::to_string(i) + " while there are only " + std::to_string(_parametersType.size()));
    }
    _parameters[i] = parameter;
  }

  unsigned int CommandLineInterface::nbrParameters() const
  {
    return _parametersType.size();
  }

  //
  // class SearchCL : public CommandLineInterface
  //

  SearchCL::SearchCL(const std::string &name) :
    CommandLineInterface(name, "", false, false, {})
  {
  }

  SearchCL::~SearchCL() {}

  CommandLine::Type SearchCL::type() const
  {
    return CommandLine::Type::Search;
  }

  void SearchCL::run(ArgsManager *argsManager) const
  {
  }

  CommandLineInterface *SearchCL::copy() const
  {
    SearchCL *cl = new SearchCL(_name);
    return cl;
  }

  //
  // class UnaryCL : public CommandLineInterface
  //

  UnaryCL::UnaryCL(const std::string &name, const std::string &help,
                   const bool kill, const bool pre,
                   void (*action)(ArgsManager *argsManager)) :
    CommandLineInterface(name, help, kill, pre, {}),
    _action(action)
  {
  }

  UnaryCL::~UnaryCL() {}

  CommandLine::Type UnaryCL::type() const
  {
    return CommandLine::Type::Unary;
  }

  void UnaryCL::run(ArgsManager *argsManager) const
  {
    _action(argsManager);
  }

  CommandLineInterface *UnaryCL::copy() const
  {
    UnaryCL *cl = new UnaryCL(_name, _help, _kill, _pre, _action);
    return cl;
  }

  void UnaryCL::action(void (*action)(ArgsManager *argsManager))
  {
    _action = action;
  }

  //
  // class UnaryHelpCL : public UnaryCL
  //

  UnaryHelpCL::UnaryHelpCL() :
    UnaryCL("help", "Print this help message", true, true, nullptr)
  {
  }

  UnaryHelpCL::~UnaryHelpCL() {}

  CommandLine::Type UnaryHelpCL::type() const
  {
    return CommandLine::Type::Help;
  }

  void UnaryHelpCL::run(ArgsManager *argsManager) const {}

  CommandLineInterface *UnaryHelpCL::copy() const
  {
    UnaryHelpCL *cl = new UnaryHelpCL();
    return cl;
  }

  //
  // class BinaryCL : public CommandLineInterface
  //

  BinaryCL::BinaryCL(const std::string &name, const std::string &help,
                     const bool kill, const bool pre,
                     const CommandLine::Parameter::Type &ParameterType,
                     void (*action)(ArgsManager *argsManager, const CommandLine::Parameter &parameter)) :
    CommandLineInterface(name, help, kill, pre, {ParameterType}),
    _action(action)
  {
    _parameters.resize(1);
  }

  BinaryCL::~BinaryCL() {}

  CommandLine::Type BinaryCL::type() const
  {
    return CommandLine::Type::Binary;
  }

  void BinaryCL::run(ArgsManager *argsManager) const
  {
    _action(argsManager, _parameters[0]);
  }

  CommandLineInterface *BinaryCL::copy() const
  {
    BinaryCL *cl = new BinaryCL(_name, _help, _kill, _pre, _parametersType[0], _action);
    return cl;
  }

  void BinaryCL::action(void (*action)(ArgsManager *argsManager, const CommandLine::Parameter &parameter))
  {
    _action = action;
  }


} // namespace gmshfem::common
