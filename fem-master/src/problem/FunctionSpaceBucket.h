// GmshFEM - Copyright (C) 2019-2022, A. Royer, E. Béchet, C. Geuzaine, Université de Liège
//
// See the LICENSE.txt file for license information. Please report all
// issues on https://gitlab.onelab.info/gmsh/fem/issues

#ifndef H_GMSHFEM_FUNCTIONSPACEBUCKET
#define H_GMSHFEM_FUNCTIONSPACEBUCKET

#include "Memory.h"

#include <string>
#include <unordered_map>
#include <vector>

namespace gmshfem::problem
{


  template< class T_PScalar >
  class FunctionSpaceBucket
  {
   private:
    std::unordered_map< std::string, std::vector< T_PScalar > > _functionSpace;
    std::unordered_map< std::string, unsigned int > _nbrOfDofsByElement;
    std::unordered_map< std::string, std::vector< int > > _functionSpaceOrientation;

   public:
    FunctionSpaceBucket();
    ~FunctionSpaceBucket();

    void clear();

    FunctionSpaceBucket(const FunctionSpaceBucket &other) = delete;
    FunctionSpaceBucket(FunctionSpaceBucket &&other);

    FunctionSpaceBucket &operator=(const FunctionSpaceBucket &other) = delete;
    FunctionSpaceBucket &operator=(FunctionSpaceBucket &&other);

    void functionSpace(const std::string &name, const std::string &gauss, std::vector< T_PScalar > &functionSpace);
    void nbrOfDofsByElement(const std::string &name, const std::string &gauss, const unsigned int nbrOfDofsByElement);
    void functionSpaceOrientation(const std::string &name, std::vector< int > &functionSpaceOrientation);

    const std::vector< T_PScalar > *functionSpace(const std::string &name, const std::string &gauss) const;
    unsigned int nbrOfDofsByElement(const std::string &name, const std::string &gauss) const;
    const std::vector< int > *functionSpaceOrientation(const std::string &name) const;

    bool have(const std::string &name, const std::string &gauss) const;
    bool haveOrientation(const std::string &name) const;

    common::Memory memory() const;
  };


} // namespace gmshfem::problem

#endif // H_GMSHFEM_ELEMENTBUCKET
