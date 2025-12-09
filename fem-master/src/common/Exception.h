// GmshFEM - Copyright (C) 2019-2022, A. Royer, E. Béchet, C. Geuzaine, Université de Liège
//
// See the LICENSE.txt file for license information. Please report all
// issues on https://gitlab.onelab.info/gmsh/fem/issues

#ifndef H_GMSHFEM_EXCEPTION
#define H_GMSHFEM_EXCEPTION

#include <exception>
#include <string>

namespace gmshfem::common
{


  class Exception : public std::exception
  {
   private:
    const std::string _sentence;

   public:
    Exception(const std::string &sentence = "") noexcept;
    ~Exception() noexcept;

    const char *what() const noexcept override;
  };


} // namespace gmshfem::common

#endif // H_GMSHFEM_EXCEPTION
