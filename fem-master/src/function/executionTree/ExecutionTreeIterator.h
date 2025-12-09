// GmshFEM - Copyright (C) 2019-2022, A. Royer, E. Béchet, C. Geuzaine, Université de Liège
//
// See the LICENSE.txt file for license information. Please report all
// issues on https://gitlab.onelab.info/gmsh/fem/issues

#ifndef H_GMSHFEM_EXECUTIONTREEITERATOR
#define H_GMSHFEM_EXECUTIONTREEITERATOR

#include "MathObject.h"

namespace gmshfem::function
{
  class ExecutionTreeInterface;
}

namespace gmshfem::function
{


  class ExecutionTreeIterator
  {
   private:
    const ExecutionTreeInterface *_tree;
    std::vector< unsigned int > _path;

   public:
    ExecutionTreeIterator(const ExecutionTreeInterface *tree);
    ~ExecutionTreeIterator();

    bool operator==(const ExecutionTreeIterator &other) const;
    bool operator!=(const ExecutionTreeIterator &other) const;

    ExecutionTreeIterator &operator++();
    const ExecutionTreeInterface *operator*();
  };


} // namespace gmshfem::function


#endif // H_GMSHFEM_EXECUTIONTREEITERATOR
