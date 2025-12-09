// GmshFEM - Copyright (C) 2019-2022, A. Royer, E. Béchet, C. Geuzaine, Université de Liège
//
// See the LICENSE.txt file for license information. Please report all
// issues on https://gitlab.onelab.info/gmsh/fem/issues

#ifndef H_GMSHFEM_PRODUCT
#define H_GMSHFEM_PRODUCT


namespace gmshfem::equation
{


  enum class Product {
    Empty,
    ScalarProduct,
    VectorProduct
  };


  template< Degree T_DegreeLhs, Product T_Product, Degree T_DegreeRhs >
  struct degreeOfProduct {
  };

  //
  // degreeOfProduct
  //

  template< Degree T_DegreeLhs, Degree T_DegreeRhs >
  struct degreeOfProduct< T_DegreeLhs, Product::ScalarProduct, T_DegreeRhs > {
    constexpr static Degree value = MathProduct< T_DegreeLhs, T_DegreeRhs >::ProductDegree;
  };

  template<>
  struct degreeOfProduct< Degree::Degree1, Product::VectorProduct, Degree::Degree1 > {
    constexpr static Degree value = Degree::Degree1;
  };


  template< Degree T_DegreeLhs, Product T_Product, Degree T_DegreeRhs >
  using DegreeOfProduct = degreeOfProduct< T_DegreeLhs, T_Product, T_DegreeRhs >;


  //
  // ComputeProduct
  //

  template< class T_Scalar, Product T_Product, Degree T_DegreeLhs, Degree T_DegreeRhs >
  struct ComputeProduct {
    static function::Function< T_Scalar, DegreeOfProduct< T_DegreeLhs, T_Product, T_DegreeRhs >::value > compute(const function::Function< T_Scalar, T_DegreeRhs > &rhs, const function::Function< T_Scalar, T_DegreeLhs > &lhs)
    {
    }
  };

  template< class T_Scalar, Degree T_DegreeLhs, Degree T_DegreeRhs >
  struct ComputeProduct< T_Scalar, Product::ScalarProduct, T_DegreeLhs, T_DegreeRhs > {
    static function::Function< T_Scalar, DegreeOfProduct< T_DegreeLhs, Product::ScalarProduct, T_DegreeRhs >::value > compute(const function::Function< T_Scalar, T_DegreeLhs > &lhs, const function::Function< T_Scalar, T_DegreeRhs > &rhs)
    {
      return lhs * rhs;
    }
  };

  template< class T_Scalar, Degree T_DegreeLhs, Degree T_DegreeRhs >
  struct ComputeProduct< T_Scalar, Product::VectorProduct, T_DegreeLhs, T_DegreeRhs > {
    static function::Function< T_Scalar, DegreeOfProduct< T_DegreeLhs, Product::VectorProduct, T_DegreeRhs >::value > compute(const function::Function< T_Scalar, T_DegreeLhs > &lhs, const function::Function< T_Scalar, T_DegreeRhs > &rhs)
    {
      return lhs % rhs;
    }
  };


} // namespace gmshfem::equation


#endif // H_GMSHFEM_PRODUCT
