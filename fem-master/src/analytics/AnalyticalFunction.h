// GmshFEM - Copyright (C) 2019-2022, A. Royer, E. Béchet, C. Geuzaine, Université de Liège
//
// See the LICENSE.txt file for license information. Please report all
// issues on https://gitlab.onelab.info/gmsh/fem/issues

#ifndef H_GMSHFEM_ANALYTICALFUNCTION
#define H_GMSHFEM_ANALYTICALFUNCTION

#include "Function.h"
#include "analyticalOperations.h"

namespace gmshfem::analytics
{


  template< class T_Op >
  class AnalyticalFunction final : public function::Function< typename T_Op::D_Scalar, T_Op::D_Degree >
  {
   public:
    AnalyticalFunction() :
      function::Function< typename T_Op::D_Scalar, T_Op::D_Degree >(new function::AnalyticalNode< T_Op >())
    {
    }

    template< class... T_Params >
    AnalyticalFunction(T_Params... params) :
      function::Function< typename T_Op::D_Scalar, T_Op::D_Degree >(new function::AnalyticalNode< T_Op >(params...))
    {
    }

    virtual ~AnalyticalFunction()
    {
    }
  };


} // namespace gmshfem::analytics

#endif // H_GMSHFEM_ANALYTICALFUNCTION
