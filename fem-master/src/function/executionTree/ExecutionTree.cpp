// GmshFEM - Copyright (C) 2019-2022, A. Royer, E. Béchet, C. Geuzaine, Université de Liège
//
// See the LICENSE.txt file for license information. Please report all
// issues on https://gitlab.onelab.info/gmsh/fem/issues

#include "ExecutionTree.h"

namespace gmshfem::function
{


  static unsigned int s_nbrExecutionTree = 0;

  ExecutionTreeInterface::ExecutionTreeInterface(const ExecutionTreeInterface &other) :
    _tag(other._tag)
  {
  }

  ExecutionTreeInterface::ExecutionTreeInterface() :
    _tag(++s_nbrExecutionTree)
  {
  }

  ExecutionTreeInterface::ExecutionTreeInterface(const unsigned int tag) :
    _tag(tag)
  {
  }

  ExecutionTreeInterface::~ExecutionTreeInterface()
  {
  }

  unsigned int ExecutionTreeInterface::tag() const
  {
    return _tag;
  }

  ExecutionTreeIterator ExecutionTreeInterface::begin() const
  {
    return ExecutionTreeIterator(this);
  }

  ExecutionTreeIterator ExecutionTreeInterface::end() const
  {
    return ExecutionTreeIterator(nullptr);
  }


} // namespace gmshfem::function
