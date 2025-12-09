// GmshFEM - Copyright (C) 2019-2022, A. Royer, E. Béchet, C. Geuzaine, Université de Liège
//
// See the LICENSE.txt file for license information. Please report all
// issues on https://gitlab.onelab.info/gmsh/fem/issues

#include "ExecutionTreeIterator.h"

#include "ExecutionTree.h"

namespace gmshfem::function
{


  ExecutionTreeIterator::ExecutionTreeIterator(const ExecutionTreeInterface *tree) :
    _tree(tree), _path()
  {
  }

  ExecutionTreeIterator::~ExecutionTreeIterator()
  {
  }

  bool ExecutionTreeIterator::operator==(const ExecutionTreeIterator &other) const
  {
    return (_tree == other._tree && _path == other._path);
  }

  bool ExecutionTreeIterator::operator!=(const ExecutionTreeIterator &other) const
  {
    return (_tree != other._tree || _path != other._path);
  }

  ExecutionTreeIterator &ExecutionTreeIterator::operator++()
  {
    const ExecutionTreeInterface *pathTree = _tree;
    for(auto i = 0ULL; i < _path.size(); ++i) {
      pathTree = pathTree->getLeaf(_path[i]);
    }

    if(pathTree->numberOfLeaves() != 0) {
      _path.push_back(0);
    }
    else {
      while(_path.size() != 0) {
        unsigned int lastLeave = _path.back();
        _path.pop_back();
        pathTree = _tree;
        for(auto i = 0ULL; i < _path.size(); ++i) {
          pathTree = pathTree->getLeaf(_path[i]);
        }
        if(pathTree->numberOfLeaves() != lastLeave + 1) {
          _path.push_back(lastLeave + 1);
          return *this;
        }
      }
      _tree = nullptr;
    }

    return *this;
  }

  const ExecutionTreeInterface *ExecutionTreeIterator::operator*()
  {
    const ExecutionTreeInterface *pathTree = _tree;
    for(auto i = 0ULL; i < _path.size(); ++i) {
      pathTree = pathTree->getLeaf(_path[i]);
    }
    return pathTree;
  }


} // namespace gmshfem::function
