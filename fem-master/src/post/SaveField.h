// GmshFEM - Copyright (C) 2019-2022, A. Royer, E. Béchet, C. Geuzaine, Université de Liège
//
// See the LICENSE.txt file for license information. Please report all
// issues on https://gitlab.onelab.info/gmsh/fem/issues

#ifndef H_GMSHFEM_SAVEFIELD
#define H_GMSHFEM_SAVEFIELD

#include "Post.h"

#include <string>

namespace gmshfem::field
{
  template< class T_Scalar, field::Form T_Form >
  class Field;
}

namespace gmshfem::post
{


  template< class T_Scalar, field::Form T_Form >
  class SaveField final : public PostInterface
  {
   private:
    const field::Field< T_Scalar, T_Form > *const _field;
    const std::string _format;
    const std::string _path;
    const bool _memorize;
    const int _step;
    const double _time;
    const int _partition;

   public:
    SaveField(const field::Field< T_Scalar, T_Form > &field, const std::string &format, const std::string &path, const double memorize, const int step, const double time, const int partition);
    virtual ~SaveField();

    virtual int run() override;
  };


} // namespace gmshfem::post


#endif // H_GMSHFEM_SAVEFIELD
