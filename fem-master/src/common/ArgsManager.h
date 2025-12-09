// GmshFEM - Copyright (C) 2019-2022, A. Royer, E. Béchet, C. Geuzaine, Université de Liège
//
// See the LICENSE.txt file for license information. Please report all
// issues on https://gitlab.onelab.info/gmsh/fem/issues

#ifndef H_GMSHFEM_ARGSMANAGER
#define H_GMSHFEM_ARGSMANAGER

#include "CommandLine.h"

#include <map>
#include <set>
#include <string>
#include <vector>

namespace gmshfem::common
{
  class GmshFem;
}

namespace gmshfem::common
{


  class ArgsManager
  {
   protected:
    int _argc;
    char **_argv;
    GmshFem *_gmshFem;
    std::set< CommandLineInterface * > _commands;
    std::vector< CommandLineInterface * > _command;
    std::map< std::string, std::string > _userDefinedParameters;
    bool _haveKillParam;

    void _parameterFile(const std::string &path);

   public:
    ArgsManager(GmshFem *gmshFem);
    virtual ~ArgsManager();

    void setArgcArgv(int argc, char **argv);
    void processArgcArgv();

    void insert(CommandLineInterface *newCL);

    void pre();
    void post();

    bool haveKillParam() const;

    GmshFem *gmshFem() const;

    void addUserDefinedParameter(const std::string &name, const std::string &value);

    bool userDefinedParameter(int &value, const std::string &name) const;
    bool userDefinedParameter(double &value, const std::string &name) const;
    bool userDefinedParameter(std::string &value, const std::string &name) const;
    bool userDefinedParameter(bool &value, const std::string &name) const;

    std::vector< CommandLineInterface * >::const_iterator begin() const;
    std::vector< CommandLineInterface * >::const_iterator end() const;
  };


} // namespace gmshfem::common

#endif // H_GMSHFEM_ARGSMANAGER
