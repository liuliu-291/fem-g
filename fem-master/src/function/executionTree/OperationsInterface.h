// GmshFEM - Copyright (C) 2019-2022, A. Royer, E. Béchet, C. Geuzaine, Université de Liège
//
// See the LICENSE.txt file for license information. Please report all
// issues on https://gitlab.onelab.info/gmsh/fem/issues

#ifndef H_GMSHFEM_OPERATIONSINTERFACE
#define H_GMSHFEM_OPERATIONSINTERFACE

#include "InputVector.h"
#include "MathObject.h"
#include "scalar.h"

#include <utility>
#include <vector>

namespace gmshfem::function
{


  template< class T_Scalar, Degree T_Degree >
  class OperationsInterface
  {
   protected:
    OperationsInterface();
    OperationsInterface(const OperationsInterface &other);
    virtual ~OperationsInterface();

   public:
    virtual bool canUseSameVectorsForOutputAndInputs() const;
    virtual std::string name() const;
  };


} // namespace gmshfem::function


#endif // H_GMSHFEM_OPERATIONSINTERFACE
