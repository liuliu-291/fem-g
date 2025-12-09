// GmshFEM - Copyright (C) 2019-2022, A. Royer, E. Béchet, C. Geuzaine, Université de Liège
//
// See the LICENSE.txt file for license information. Please report all
// issues on https://gitlab.onelab.info/gmsh/fem/issues

#ifndef H_GMSHFEM_RCM
#define H_GMSHFEM_RCM

#include <cstdint>
#include <vector>

namespace gmshfem::reorder
{


  class RCM
  {
   private:
   public:
    RCM();
    ~RCM();

    void apply(std::vector< unsigned long long > &sorted, const unsigned long long *const row, const unsigned long long *const indices);
  };


} // namespace gmshfem::reorder

#endif // H_GMSHFEM_RCM
