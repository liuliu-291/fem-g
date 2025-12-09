// GmshFEM - Copyright (C) 2019-2022, A. Royer, E. Béchet, C. Geuzaine, Université de Liège
//
// See the LICENSE.txt file for license information. Please report all
// issues on https://gitlab.onelab.info/gmsh/fem/issues

#include "FunctionSpaceBucket.h"

#include "instantiate.h"

namespace gmshfem::problem
{


  template< class T_PScalar >
  FunctionSpaceBucket< T_PScalar >::FunctionSpaceBucket() :
    _functionSpace(), _nbrOfDofsByElement(), _functionSpaceOrientation()
  {
  }

  template< class T_PScalar >
  FunctionSpaceBucket< T_PScalar >::~FunctionSpaceBucket()
  {
  }

  template< class T_PScalar >
  void FunctionSpaceBucket< T_PScalar >::clear()
  {
    _functionSpace.clear();
    _nbrOfDofsByElement.clear();
    _functionSpaceOrientation.clear();
  }

  template< class T_PScalar >
  FunctionSpaceBucket< T_PScalar >::FunctionSpaceBucket(FunctionSpaceBucket< T_PScalar > &&other) :
    _functionSpace(std::move(other._functionSpace)), _nbrOfDofsByElement(std::move(other._nbrOfDofsByElement)), _functionSpaceOrientation(std::move(other._functionSpaceOrientation))
  {
  }

  template< class T_PScalar >
  FunctionSpaceBucket< T_PScalar > &FunctionSpaceBucket< T_PScalar >::operator=(FunctionSpaceBucket< T_PScalar > &&other)
  {
    _functionSpace = std::move(other._functionSpace);
    _nbrOfDofsByElement = std::move(other._nbrOfDofsByElement);
    _functionSpaceOrientation = std::move(other._functionSpaceOrientation);

    return *this;
  }

  template< class T_PScalar >
  void FunctionSpaceBucket< T_PScalar >::functionSpace(const std::string &name, const std::string &gauss, std::vector< T_PScalar > &functionSpace)
  {
    _functionSpace[name + gauss] = std::move(functionSpace);
  }

  template< class T_PScalar >
  void FunctionSpaceBucket< T_PScalar >::nbrOfDofsByElement(const std::string &name, const std::string &gauss, const unsigned int nbrOfDofsByElement)
  {
    _nbrOfDofsByElement[name + gauss] = nbrOfDofsByElement;
  }

  template< class T_PScalar >
  void FunctionSpaceBucket< T_PScalar >::functionSpaceOrientation(const std::string &name, std::vector< int > &functionSpaceOrientation)
  {
    _functionSpaceOrientation[name] = std::move(functionSpaceOrientation);
  }

  template< class T_PScalar >
  const std::vector< T_PScalar > *FunctionSpaceBucket< T_PScalar >::functionSpace(const std::string &name, const std::string &gauss) const
  {
    auto it = _functionSpace.find(name + gauss);
    if(it != _functionSpace.end()) {
      return &(it->second);
    }
    return nullptr;
  }

  template< class T_PScalar >
  unsigned int FunctionSpaceBucket< T_PScalar >::nbrOfDofsByElement(const std::string &name, const std::string &gauss) const
  {
    auto it = _nbrOfDofsByElement.find(name + gauss);
    if(it != _nbrOfDofsByElement.end()) {
      return it->second;
    }
    return 0;
  }

  template< class T_PScalar >
  const std::vector< int > *FunctionSpaceBucket< T_PScalar >::functionSpaceOrientation(const std::string &name) const
  {
    auto it = _functionSpaceOrientation.find(name);
    if(it != _functionSpaceOrientation.end()) {
      return &(it->second);
    }
    return nullptr;
  }

  template< class T_PScalar >
  bool FunctionSpaceBucket< T_PScalar >::have(const std::string &name, const std::string &gauss) const
  {
    auto it = _functionSpace.find(name + gauss);
    if(it != _functionSpace.end()) {
      return true;
    }
    return false;
  }

  template< class T_PScalar >
  bool FunctionSpaceBucket< T_PScalar >::haveOrientation(const std::string &name) const
  {
    auto it = _functionSpaceOrientation.find(name);
    if(it != _functionSpaceOrientation.end()) {
      return true;
    }
    return false;
  }

  template< class T_PScalar >
  common::Memory FunctionSpaceBucket< T_PScalar >::memory() const
  {
    common::Memory memory;
    for(auto it = _functionSpace.begin(); it != _functionSpace.end(); ++it) {
      memory += it->second.size() * sizeof(T_PScalar);
    }

    memory += _nbrOfDofsByElement.size() * sizeof(unsigned int);

    for(auto it = _functionSpaceOrientation.begin(); it != _functionSpaceOrientation.end(); ++it) {
      memory += it->second.size() * sizeof(int);
    }

    return memory;
  }

  INSTANTIATE_CLASS(FunctionSpaceBucket, 2, TEMPLATE_ARGS(double, float))


} // namespace gmshfem::problem
