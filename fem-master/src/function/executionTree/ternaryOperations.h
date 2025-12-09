// GmshFEM - Copyright (C) 2019-2022, A. Royer, E. Béchet, C. Geuzaine, Université de Liège
//
// See the LICENSE.txt file for license information. Please report all
// issues on https://gitlab.onelab.info/gmsh/fem/issues

#ifndef H_GMSHFEM_TERNARYOPERATIONS
#define H_GMSHFEM_TERNARYOPERATIONS

#include "OmpInterface.h"
#include "TernaryNode.h"

namespace gmshfem::function
{


  // Available ternary operations:
  //  TensorDiag
  //  Vector

  //
  // TensorDiag
  //

  template< class T_Scalar >
  class TensorDiag final : public TernaryOperation< T_Scalar, Degree::Degree2, Degree::Degree0, Degree::Degree0, Degree::Degree0 >
  {
   public:
    TensorDiag()
    {
    }

    TensorDiag(const TensorDiag &other) :
      TernaryOperation< T_Scalar, Degree::Degree2, Degree::Degree0, Degree::Degree0, Degree::Degree0 >(other)
    {
    }

    template< class T_A, class T_B, class T_C >
    void operator()(OutputVector< T_Scalar, Degree::Degree2 > &values,
                    const InputVector< T_Scalar, Degree::Degree0, T_A > &a,
                    const InputVector< T_Scalar, Degree::Degree0, T_B > &b,
                    const InputVector< T_Scalar, Degree::Degree0, T_C > &c,
                    const std::vector< scalar::Precision< T_Scalar >, numa::allocator< scalar::Precision< T_Scalar > > > &points, const std::vector< scalar::Precision< T_Scalar > > &gaussPoints, const int elementType, const std::pair< int, int > &entity) const
    {
#pragma omp for
      for(auto j = 0ULL; j < values.size(); ++j) {
        values[j] = typename MathObject< T_Scalar, Degree::Degree1 >::Object(a[j], b[j], c[j]).asDiagonal();
      }
    }

    std::string name() const override
    {
      return "create diagonal tensor [a 0 0, 0 b 0, 0 0 c]";
    }
  };

  //
  // Vector
  //

  template< class T_Scalar >
  class Vector final : public TernaryOperation< T_Scalar, Degree::Degree1, Degree::Degree0, Degree::Degree0, Degree::Degree0 >
  {
   public:
    Vector()
    {
    }

    Vector(const Vector &other) :
      TernaryOperation< T_Scalar, Degree::Degree1, Degree::Degree0, Degree::Degree0, Degree::Degree0 >(other)
    {
    }

    template< class T_A, class T_B, class T_C >
    void operator()(OutputVector< T_Scalar, Degree::Degree1 > &values,
                    const InputVector< T_Scalar, Degree::Degree0, T_A > &a,
                    const InputVector< T_Scalar, Degree::Degree0, T_B > &b,
                    const InputVector< T_Scalar, Degree::Degree0, T_C > &c,
                    const std::vector< scalar::Precision< T_Scalar >, numa::allocator< scalar::Precision< T_Scalar > > > &points, const std::vector< scalar::Precision< T_Scalar > > &gaussPoints, const int elementType, const std::pair< int, int > &entity) const
    {
#pragma omp for
      for(auto j = 0ULL; j < values.size(); ++j) {
        values[j] << a[j], b[j], c[j];
      }
    }

    std::string name() const override
    {
      return "create vector [a b c]";
    }
  };


} // namespace gmshfem::function


#endif // H_GMSHFEM_TERNARYOPERATIONS
