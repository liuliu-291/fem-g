// GmshFEM - Copyright (C) 2019-2022, A. Royer, E. Béchet, C. Geuzaine, Université de Liège
//
// See the LICENSE.txt file for license information. Please report all
// issues on https://gitlab.onelab.info/gmsh/fem/issues

#ifndef H_GMSHFEM_PIECEWISEFUNCTION
#define H_GMSHFEM_PIECEWISEFUNCTION

#include "Domain.h"
#include "Function.h"

namespace gmshfem::function
{


  // *****************************
  // PiecewiseFunction< T_Scalar, T_Degree >
  // *****************************

  template< class T_Scalar, Degree T_Degree >
  class PiecewiseFunction : public Function< T_Scalar, T_Degree >
  {
   public:
    PiecewiseFunction();
    virtual ~PiecewiseFunction();

    void addFunction(const Function< T_Scalar, T_Degree > &f, const domain::Domain &domain) const;
  };

  template< class T_Scalar >
  using ScalarPiecewiseFunction = PiecewiseFunction< T_Scalar, Degree::Degree0 >;
  template< class T_Scalar >
  using VectorPiecewiseFunction = PiecewiseFunction< T_Scalar, Degree::Degree1 >;
  template< class T_Scalar, unsigned int T_Rank = 2 >
  using TensorPiecewiseFunction = PiecewiseFunction< T_Scalar, Degree(TensorRank< T_Rank >::value) >;


} // namespace gmshfem::function

#endif // H_GMSHFEM_PIECEWISEFUNCTION
