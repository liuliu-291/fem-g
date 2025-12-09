// GmshFEM - Copyright (C) 2019-2022, A. Royer, E. Béchet, C. Geuzaine, Université de Liège
//
// See the LICENSE.txt file for license information. Please report all
// issues on https://gitlab.onelab.info/gmsh/fem/issues

#ifndef H_GMSHFEM_SAVEFUNCTION
#define H_GMSHFEM_SAVEFUNCTION

#include "Domain.h"
#include "Function.h"
#include "MathObject.h"
#include "Post.h"
#include "scalar.h"

#include <string>

namespace gmshfem::post
{


  template< class T_Scalar, Degree T_Degree >
  class SaveFunction final : public PostInterface
  {
   private:
    const function::Function< T_Scalar, T_Degree > _function;
    domain::GeometricObject _domain;
    const std::string _name;
    const std::string _format;
    const std::string _path;
    const bool _memorize;
    const int _step;
    const double _time;
    const int _partition;

   public:
    SaveFunction(const function::Function< T_Scalar, T_Degree > &function, const domain::GeometricObject &domain, const std::string &name, const std::string &format, const std::string &path, const double memorize, const int step, const double time, const int partition);
    virtual ~SaveFunction();

    virtual int run() override;
  };


} // namespace gmshfem::post


#endif // H_GMSHFEM_SAVEFUNCTION
