// GmshFEM - Copyright (C) 2019-2022, A. Royer, E. Béchet, C. Geuzaine, Université de Liège
//
// See the LICENSE.txt file for license information. Please report all
// issues on https://gitlab.onelab.info/gmsh/fem/issues

#include "OperationsInterface.h"

#include "instantiate.h"

namespace gmshfem::function
{


  template< class T_Scalar, Degree T_Degree >
  OperationsInterface< T_Scalar, T_Degree >::OperationsInterface()
  {
  }

  template< class T_Scalar, Degree T_Degree >
  OperationsInterface< T_Scalar, T_Degree >::OperationsInterface(const OperationsInterface &other)
  {
  }

  template< class T_Scalar, Degree T_Degree >
  OperationsInterface< T_Scalar, T_Degree >::~OperationsInterface()
  {
  }

  template< class T_Scalar, Degree T_Degree >
  bool OperationsInterface< T_Scalar, T_Degree >::canUseSameVectorsForOutputAndInputs() const
  {
    return false;
  }

  template< class T_Scalar, Degree T_Degree >
  std::string OperationsInterface< T_Scalar, T_Degree >::name() const
  {
    return "#UnknownName";
  }

  INSTANTIATE_CLASS_2(OperationsInterface, 4, 4, TEMPLATE_ARGS(std::complex< double >, std::complex< float >, double, float), TEMPLATE_ARGS(Degree::Degree0, Degree::Degree1, Degree::Degree2, Degree::Degree4))


} // namespace gmshfem::function
