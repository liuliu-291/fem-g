// GmshFEM - Copyright (C) 2019-2022, A. Royer, E. Béchet, C. Geuzaine, Université de Liège
//
// See the LICENSE.txt file for license information. Please report all
// issues on https://gitlab.onelab.info/gmsh/fem/issues

#ifndef H_GMSHFEM_POSTPROMAP
#define H_GMSHFEM_POSTPROMAP

#include "Domain.h"
#include "Function.h"
#include "MathObject.h"
#include "scalar.h"

#include <string>

namespace gmshfem::post
{


  template< class T_Scalar, Degree T_Degree >
  class PostproMap
  {
   private:
    const std::string _name;
    int _tag;

   public:
    PostproMap(const std::string &name, const int tag = 0);
    ~PostproMap();

    void append(const function::Function< T_Scalar, T_Degree > &function, const domain::GeometricObject &domain, const int step, const double time, const int partition = 0);
    template< field::Form T_Form, class = typename std::enable_if< field::DegreeOfForm< T_Form >::value == T_Degree >::type >
    void append(const field::Field< T_Scalar, T_Form > &field, const domain::GeometricObject &domain, const int step, const double time, const int partition = 0);
    template< field::Form T_Form, unsigned int T_NumFields, class = typename std::enable_if< field::DegreeOfCompoundForm< T_Form >::value == T_Degree >::type >
    void append(const field::CompoundField< T_Scalar, T_Form, T_NumFields > &field, const domain::GeometricObject &domain, const int step, const double time, const int partition = 0);
    void write(const std::string &format = "msh", const std::string &path = "") const;
  };


} // namespace gmshfem::post


#endif // H_GMSHFEM_POSTPROMAP
