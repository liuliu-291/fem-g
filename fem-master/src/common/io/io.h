// GmshFEM - Copyright (C) 2019-2022, A. Royer, E. Béchet, C. Geuzaine, Université de Liège
//
// See the LICENSE.txt file for license information. Please report all
// issues on https://gitlab.onelab.info/gmsh/fem/issues

#ifndef H_GMSHFEM_IO
#define H_GMSHFEM_IO

namespace gmshfem::common
{


  enum class OpeningMode {
    NewFile,
    Append,
    Reading
  };


} // namespace gmshfem::common

#include "CSVio.h"
#include "PGFPlotsio.h"
#include "PPMio.h"

#endif // H_GMSHFEM_IO
