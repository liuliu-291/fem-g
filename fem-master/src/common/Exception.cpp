// GmshFEM - Copyright (C) 2019-2022, A. Royer, E. Béchet, C. Geuzaine, Université de Liège
//
// See the LICENSE.txt file for license information. Please report all
// issues on https://gitlab.onelab.info/gmsh/fem/issues

#include "Exception.h"

namespace gmshfem::common
{


  Exception::Exception(const std::string &sentence) noexcept :
    _sentence(sentence)
  {
  }

  Exception::~Exception() noexcept
  {
  }

  const char *Exception::what() const noexcept
  {
    return _sentence.c_str();
  }


} // namespace gmshfem::common
