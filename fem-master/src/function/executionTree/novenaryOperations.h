// GmshFEM - Copyright (C) 2019-2022, A. Royer, E. Béchet, C. Geuzaine, Université de Liège
//
// See the LICENSE.txt file for license information. Please report all
// issues on https://gitlab.onelab.info/gmsh/fem/issues

#ifndef H_GMSHFEM_NOVENARYOPERATIONS
#define H_GMSHFEM_NOVENARYOPERATIONS

#include "NovenaryNode.h"
#include "OmpInterface.h"

namespace gmshfem::function
{


  // Available novenary operations:
  //  Tensor
  //  Tensor4

  //
  // Tensor
  //

  template< class T_Scalar >
  class Tensor final : public NovenaryOperation< T_Scalar, Degree::Degree2, Degree::Degree0, Degree::Degree0, Degree::Degree0, Degree::Degree0, Degree::Degree0, Degree::Degree0, Degree::Degree0, Degree::Degree0, Degree::Degree0 >
  {
   public:
    Tensor()
    {
    }

    Tensor(const Tensor &other) :
      NovenaryOperation< T_Scalar, Degree::Degree2, Degree::Degree0, Degree::Degree0, Degree::Degree0, Degree::Degree0, Degree::Degree0, Degree::Degree0, Degree::Degree0, Degree::Degree0, Degree::Degree0 >(other)
    {
    }

    template< class T_A, class T_B, class T_C, class T_D, class T_E, class T_F, class T_G, class T_H, class T_I >
    void operator()(OutputVector< T_Scalar, Degree::Degree2 > &values,
                    const InputVector< T_Scalar, Degree::Degree0, T_A > &a,
                    const InputVector< T_Scalar, Degree::Degree0, T_B > &b,
                    const InputVector< T_Scalar, Degree::Degree0, T_C > &c,
                    const InputVector< T_Scalar, Degree::Degree0, T_D > &d,
                    const InputVector< T_Scalar, Degree::Degree0, T_E > &e,
                    const InputVector< T_Scalar, Degree::Degree0, T_F > &f,
                    const InputVector< T_Scalar, Degree::Degree0, T_G > &g,
                    const InputVector< T_Scalar, Degree::Degree0, T_H > &h,
                    const InputVector< T_Scalar, Degree::Degree0, T_I > &i,
                    const std::vector< scalar::Precision< T_Scalar >, numa::allocator< scalar::Precision< T_Scalar > > > &points, const std::vector< scalar::Precision< T_Scalar > > &gaussPoints, const int elementType, const std::pair< int, int > &entity) const
    {
#pragma omp for
      for(auto j = 0ULL; j < values.size(); ++j) {
        values[j] << a[j], b[j], c[j],
          d[j], e[j], f[j],
          g[j], h[j], i[j];
      }
    }

    std::string name() const override
    {
      return "create tensor [a b c, d e f, g h i]";
    }
  };


  template< class T_Scalar >
  class Tensor4 final : public NovenaryOperation< T_Scalar, Degree::Degree4, Degree::Degree2, Degree::Degree2, Degree::Degree2, Degree::Degree2, Degree::Degree2, Degree::Degree2, Degree::Degree2, Degree::Degree2, Degree::Degree2 >
  {
   public:
    Tensor4()
    {
    }

    Tensor4(const Tensor4 &other) :
      NovenaryOperation< T_Scalar, Degree::Degree4, Degree::Degree2, Degree::Degree2, Degree::Degree2, Degree::Degree2, Degree::Degree2, Degree::Degree2, Degree::Degree2, Degree::Degree2, Degree::Degree2 >(other)
    {
    }

    template< class T_A, class T_B, class T_C, class T_D, class T_E, class T_F, class T_G, class T_H, class T_I >
    void operator()(OutputVector< T_Scalar, Degree::Degree4 > &values,
                    const InputVector< T_Scalar, Degree::Degree2, T_A > &a,
                    const InputVector< T_Scalar, Degree::Degree2, T_B > &b,
                    const InputVector< T_Scalar, Degree::Degree2, T_C > &c,
                    const InputVector< T_Scalar, Degree::Degree2, T_D > &d,
                    const InputVector< T_Scalar, Degree::Degree2, T_E > &e,
                    const InputVector< T_Scalar, Degree::Degree2, T_F > &f,
                    const InputVector< T_Scalar, Degree::Degree2, T_G > &g,
                    const InputVector< T_Scalar, Degree::Degree2, T_H > &h,
                    const InputVector< T_Scalar, Degree::Degree2, T_I > &i,
                    const std::vector< scalar::Precision< T_Scalar >, numa::allocator< scalar::Precision< T_Scalar > > > &points, const std::vector< scalar::Precision< T_Scalar > > &gaussPoints, const int elementType, const std::pair< int, int > &entity) const
    {
#pragma omp for
      for(auto j = 0ULL; j < values.size(); ++j) {
        values[j](0, 0) = a[j];
        values[j](0, 1) = b[j];
        values[j](0, 2) = c[j];
        values[j](1, 0) = d[j];
        values[j](1, 1) = e[j];
        values[j](1, 2) = f[j];
        values[j](2, 0) = g[j];
        values[j](2, 1) = h[j];
        values[j](2, 2) = i[j];
      }
    }

    std::string name() const override
    {
      return "create tensor<4> [a b c, d e f, g h i]";
    }
  };


} // namespace gmshfem::function


#endif // H_GMSHFEM_NOVENARYOPERATIONS
