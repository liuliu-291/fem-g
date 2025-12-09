// GmshFEM - Copyright (C) 2019-2022, A. Royer, E. Béchet, C. Geuzaine, Université de Liège
//
// See the LICENSE.txt file for license information. Please report all
// issues on https://gitlab.onelab.info/gmsh/fem/issues

#ifndef H_GMSHFEM_OPTIONSENUMS
#define H_GMSHFEM_OPTIONSENUMS

namespace gmshfem::problem::DofsSort
{


  enum class Algorithm {
    Default = 0,
    None = 1,
    Hilbert = 2,
    RCM = 3
  };


} // namespace gmshfem::problem::DofsSort

namespace gmshfem::problem::ElementsSort
{


  enum class Algorithm {
    None = 1,
    Hilbert = 2
  };


} // namespace gmshfem::problem::ElementsSort

#endif // H_GMSHFEM_OPTIONSENUMS
