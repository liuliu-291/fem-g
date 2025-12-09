// GmshFEM - Copyright (C) 2019-2022, A. Royer, E. Béchet, C. Geuzaine, Université de Liège
//
// See the LICENSE.txt file for license information. Please report all
// issues on https://gitlab.onelab.info/gmsh/fem/issues

#ifndef H_GMSHFEM_INPUTVECTOR
#define H_GMSHFEM_INPUTVECTOR

#include "FunctionAllocator.h"
#include "Message.h"

namespace gmshfem::function
{


  template< class T_Scalar, Degree T_Degree, class T_IsConstant = std::false_type >
  class InputVector
  {
   private:
    const OutputVector< T_Scalar, T_Degree > &_vector;

   public:
    InputVector(const OutputVector< T_Scalar, T_Degree > &vector) :
      _vector(vector)
    {
    }

    ~InputVector()
    {
    }

    inline const typename MathObject< T_Scalar, T_Degree >::Object &operator[](const unsigned long long i) const
    {
      return _vector[i];
    }
  };

  template< class T_Scalar, Degree T_Degree >
  class InputVector< T_Scalar, T_Degree, std::true_type >
  {
   private:
    const OutputVector< T_Scalar, T_Degree > &_vector;

   public:
    InputVector(const OutputVector< T_Scalar, T_Degree > &vector) :
      _vector(vector)
    {
    }

    ~InputVector()
    {
    }

    inline const typename MathObject< T_Scalar, T_Degree >::Object &operator[](const unsigned long long i) const
    {
      return _vector[0];
    }
  };


} // namespace gmshfem::function


#endif // H_GMSHFEM_INPUTVECTOR
