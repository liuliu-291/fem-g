// GmshFEM - Copyright (C) 2019-2022, A. Royer, E. Béchet, C. Geuzaine, Université de Liège
//
// See the LICENSE.txt file for license information. Please report all
// issues on https://gitlab.onelab.info/gmsh/fem/issues

#ifndef H_GMSHFEM_GMSHFEM
#define H_GMSHFEM_GMSHFEM

#include "ArgsManager.h"

#include <map>
#include <string>
#include <vector>

namespace gmshfem::common
{


  class GmshFem
  {
   private:
    ArgsManager _argsManager;
    const std::string _userLibName;

   public:
    GmshFem(const std::string &userLibName);
    GmshFem(int argc, char **argv);
    ~GmshFem();

    void finalizeInit(int argc, char **argv);

    void setVerbosity(const int verbose);
    void setMaxThreads(const int maxThreads);

    bool userDefinedParameter(int &value, const std::string &name) const;
    bool userDefinedParameter(unsigned int &value, const std::string &name) const;
    bool userDefinedParameter(double &value, const std::string &name) const;
    bool userDefinedParameter(float &value, const std::string &name) const;
    bool userDefinedParameter(std::string &value, const std::string &name) const;
    bool userDefinedParameter(char *value, const std::string &name) const;
    bool userDefinedParameter(bool &value, const std::string &name) const;

    std::string version() const;
    std::string date() const;
    std::string host() const;
    std::string packager() const;
    std::string os() const;
    std::string configOptions() const;
    std::string cxxFalgs() const;
    ArgsManager *getArgsManager();

    static GmshFem *currentInstance();
  };


} // namespace gmshfem::common

#endif // H_GMSHFEM_GMSHFEM
