// GmshFEM - Copyright (C) 2019-2022, A. Royer, E. Béchet, C. Geuzaine, Université de Liège
//
// See the LICENSE.txt file for license information. Please report all
// issues on https://gitlab.onelab.info/gmsh/fem/issues

#ifndef H_GMSHFEM_SCALAROPERATIONS
#define H_GMSHFEM_SCALAROPERATIONS

#include "Message.h"
#include "OmpInterface.h"
#include "ScalarTypeNode.h"

namespace gmshfem::function
{


  // Available scalar type operations:
  //  ComplexImaginaryPart
  //  ComplexRealPart
  //  RealToComplex


  //
  // ComplexImaginaryPart
  //

  template< class T_Scalar, Degree T_Degree >
  class ComplexImaginaryPart final : public ScalarTypeOperation< scalar::Precision< T_Scalar >, T_Scalar, T_Degree >
  {
   public:
    ComplexImaginaryPart()
    {
    }

    ComplexImaginaryPart(const ComplexImaginaryPart &other) :
      ScalarTypeOperation< scalar::Precision< T_Scalar >, T_Scalar, T_Degree >(other)
    {
    }

    void operator()(OutputVector< scalar::Precision< T_Scalar >, T_Degree > &values, OutputVector< T_Scalar, T_Degree > &a, const std::vector< scalar::Precision< T_Scalar >, numa::allocator< scalar::Precision< T_Scalar > > > &points, const std::vector< scalar::Precision< T_Scalar > > &gaussPoints, const int elementType, const std::pair< int, int > &entity) const
    {
#pragma omp for
      for(auto i = 0ULL; i < values.size(); ++i) {
        values[i] = scalar::Precision< T_Scalar >(std::imag(a[i]));
      }
    }

    std::string name() const override
    {
      return "imaginary part of a";
    }
  };

  template< class T_Scalar >
  class ComplexImaginaryPart< T_Scalar, Degree::Degree1 > final : public ScalarTypeOperation< scalar::Precision< T_Scalar >, T_Scalar, Degree::Degree1 >
  {
   public:
    ComplexImaginaryPart()
    {
    }

    ComplexImaginaryPart(const ComplexImaginaryPart &other) :
      ScalarTypeOperation< scalar::Precision< T_Scalar >, T_Scalar, Degree::Degree1 >(other)
    {
    }

    void operator()(OutputVector< scalar::Precision< T_Scalar >, Degree::Degree1 > &values, OutputVector< T_Scalar, Degree::Degree1 > &a, const std::vector< scalar::Precision< T_Scalar >, numa::allocator< scalar::Precision< T_Scalar > > > &points, const std::vector< scalar::Precision< T_Scalar > > &gaussPoints, const int elementType, const std::pair< int, int > &entity) const
    {
#pragma omp for
      for(auto i = 0ULL; i < values.size(); ++i) {
        values[i] << scalar::Precision< T_Scalar >(std::imag(a[i](0))), scalar::Precision< T_Scalar >(std::imag(a[i](1))), scalar::Precision< T_Scalar >(std::imag(a[i](2)));
      }
    }

    std::string name() const override
    {
      return "imaginary part of a";
    }
  };

  template< class T_Scalar >
  class ComplexImaginaryPart< T_Scalar, Degree::Degree2 > final : public ScalarTypeOperation< scalar::Precision< T_Scalar >, T_Scalar, Degree::Degree2 >
  {
   public:
    ComplexImaginaryPart()
    {
    }

    ComplexImaginaryPart(const ComplexImaginaryPart &other) :
      ScalarTypeOperation< scalar::Precision< T_Scalar >, T_Scalar, Degree::Degree2 >(other)
    {
    }

    void operator()(OutputVector< scalar::Precision< T_Scalar >, Degree::Degree2 > &values, OutputVector< T_Scalar, Degree::Degree2 > &a, const std::vector< scalar::Precision< T_Scalar >, numa::allocator< scalar::Precision< T_Scalar > > > &points, const std::vector< scalar::Precision< T_Scalar > > &gaussPoints, const int elementType, const std::pair< int, int > &entity) const
    {
#pragma omp for
      for(auto i = 0ULL; i < values.size(); ++i) {
        values[i] << scalar::Precision< T_Scalar >(std::imag(a[i](0, 0))), scalar::Precision< T_Scalar >(std::imag(a[i](0, 1))), scalar::Precision< T_Scalar >(std::imag(a[i](0, 2))), scalar::Precision< T_Scalar >(std::imag(a[i](1, 0))), scalar::Precision< T_Scalar >(std::imag(a[i](1, 1))), scalar::Precision< T_Scalar >(std::imag(a[i](1, 2))), scalar::Precision< T_Scalar >(std::imag(a[i](2, 0))), scalar::Precision< T_Scalar >(std::imag(a[i](2, 1))), scalar::Precision< T_Scalar >(std::imag(a[i](2, 2)));
      }
    }

    std::string name() const override
    {
      return "imaginary part of a";
    }
  };

  template< class T_Scalar >
  class ComplexImaginaryPart< T_Scalar, Degree::Degree4 > final : public ScalarTypeOperation< scalar::Precision< T_Scalar >, T_Scalar, Degree::Degree4 >
  {
   public:
    ComplexImaginaryPart()
    {
    }

    ComplexImaginaryPart(const ComplexImaginaryPart &other) :
      ScalarTypeOperation< scalar::Precision< T_Scalar >, T_Scalar, Degree::Degree4 >(other)
    {
    }

    void operator()(OutputVector< scalar::Precision< T_Scalar >, Degree::Degree4 > &values, OutputVector< T_Scalar, Degree::Degree4 > &a, const std::vector< scalar::Precision< T_Scalar >, numa::allocator< scalar::Precision< T_Scalar > > > &points, const std::vector< scalar::Precision< T_Scalar > > &gaussPoints, const int elementType, const std::pair< int, int > &entity) const
    {
#pragma omp for
      for(auto i = 0ULL; i < values.size(); ++i) {
        for(unsigned int I = 0; I < 3; ++I) {
          for(unsigned int J = 0; J < 3; ++J) {
            values[i](I, J) << scalar::Precision< T_Scalar >(std::imag(a[i](I, J)(0, 0))), scalar::Precision< T_Scalar >(std::imag(a[i](I, J)(0, 1))), scalar::Precision< T_Scalar >(std::imag(a[i](I, J)(0, 2))), scalar::Precision< T_Scalar >(std::imag(a[i](I, J)(1, 0))), scalar::Precision< T_Scalar >(std::imag(a[i](I, J)(1, 1))), scalar::Precision< T_Scalar >(std::imag(a[i](I, J)(1, 2))), scalar::Precision< T_Scalar >(std::imag(a[i](I, J)(2, 0))), scalar::Precision< T_Scalar >(std::imag(a[i](I, J)(2, 1))), scalar::Precision< T_Scalar >(std::imag(a[i](I, J)(2, 2)));
          }
        }
      }
    }

    std::string name() const override
    {
      return "imaginary part of a";
    }
  };

  //
  // ComplexRealPart
  //

  template< class T_Scalar, Degree T_Degree >
  class ComplexRealPart final : public ScalarTypeOperation< scalar::Precision< T_Scalar >, T_Scalar, T_Degree >
  {
   public:
    ComplexRealPart()
    {
    }

    ComplexRealPart(const ComplexRealPart &other) :
      ScalarTypeOperation< scalar::Precision< T_Scalar >, T_Scalar, T_Degree >(other)
    {
    }

    void operator()(OutputVector< scalar::Precision< T_Scalar >, T_Degree > &values, OutputVector< T_Scalar, T_Degree > &a, const std::vector< scalar::Precision< T_Scalar >, numa::allocator< scalar::Precision< T_Scalar > > > &points, const std::vector< scalar::Precision< T_Scalar > > &gaussPoints, const int elementType, const std::pair< int, int > &entity) const
    {
#pragma omp for
      for(auto i = 0ULL; i < values.size(); ++i) {
        values[i] = scalar::Precision< T_Scalar >(std::real(a[i]));
      }
    }

    std::string name() const override
    {
      return "real part of a";
    }
  };

  template< class T_Scalar >
  class ComplexRealPart< T_Scalar, Degree::Degree1 > final : public ScalarTypeOperation< scalar::Precision< T_Scalar >, T_Scalar, Degree::Degree1 >
  {
   public:
    ComplexRealPart()
    {
    }

    ComplexRealPart(const ComplexRealPart &other) :
      ScalarTypeOperation< scalar::Precision< T_Scalar >, T_Scalar, Degree::Degree1 >(other)
    {
    }

    void operator()(OutputVector< scalar::Precision< T_Scalar >, Degree::Degree1 > &values, OutputVector< T_Scalar, Degree::Degree1 > &a, const std::vector< scalar::Precision< T_Scalar >, numa::allocator< scalar::Precision< T_Scalar > > > &points, const std::vector< scalar::Precision< T_Scalar > > &gaussPoints, const int elementType, const std::pair< int, int > &entity) const
    {
#pragma omp for
      for(auto i = 0ULL; i < values.size(); ++i) {
        values[i] << scalar::Precision< T_Scalar >(std::real(a[i](0))), scalar::Precision< T_Scalar >(std::real(a[i](1))), scalar::Precision< T_Scalar >(std::real(a[i](2)));
      }
    }

    std::string name() const override
    {
      return "real part of a";
    }
  };

  template< class T_Scalar >
  class ComplexRealPart< T_Scalar, Degree::Degree2 > final : public ScalarTypeOperation< scalar::Precision< T_Scalar >, T_Scalar, Degree::Degree2 >
  {
   public:
    ComplexRealPart()
    {
    }

    ComplexRealPart(const ComplexRealPart &other) :
      ScalarTypeOperation< scalar::Precision< T_Scalar >, T_Scalar, Degree::Degree2 >(other)
    {
    }

    void operator()(OutputVector< scalar::Precision< T_Scalar >, Degree::Degree2 > &values, OutputVector< T_Scalar, Degree::Degree2 > &a, const std::vector< scalar::Precision< T_Scalar >, numa::allocator< scalar::Precision< T_Scalar > > > &points, const std::vector< scalar::Precision< T_Scalar > > &gaussPoints, const int elementType, const std::pair< int, int > &entity) const
    {
#pragma omp for
      for(auto i = 0ULL; i < values.size(); ++i) {
        values[i] << scalar::Precision< T_Scalar >(std::real(a[i](0, 0))), scalar::Precision< T_Scalar >(std::real(a[i](0, 1))), scalar::Precision< T_Scalar >(std::real(a[i](0, 2))), scalar::Precision< T_Scalar >(std::real(a[i](1, 0))), scalar::Precision< T_Scalar >(std::real(a[i](1, 1))), scalar::Precision< T_Scalar >(std::real(a[i](1, 2))), scalar::Precision< T_Scalar >(std::real(a[i](2, 0))), scalar::Precision< T_Scalar >(std::real(a[i](2, 1))), scalar::Precision< T_Scalar >(std::real(a[i](2, 2)));
      }
    }

    std::string name() const override
    {
      return "real part of a";
    }
  };

  template< class T_Scalar >
  class ComplexRealPart< T_Scalar, Degree::Degree4 > final : public ScalarTypeOperation< scalar::Precision< T_Scalar >, T_Scalar, Degree::Degree4 >
  {
   public:
    ComplexRealPart()
    {
    }

    ComplexRealPart(const ComplexRealPart &other) :
      ScalarTypeOperation< scalar::Precision< T_Scalar >, T_Scalar, Degree::Degree4 >(other)
    {
    }

    void operator()(OutputVector< scalar::Precision< T_Scalar >, Degree::Degree4 > &values, OutputVector< T_Scalar, Degree::Degree4 > &a, const std::vector< scalar::Precision< T_Scalar >, numa::allocator< scalar::Precision< T_Scalar > > > &points, const std::vector< scalar::Precision< T_Scalar > > &gaussPoints, const int elementType, const std::pair< int, int > &entity) const
    {
#pragma omp for
      for(auto i = 0ULL; i < values.size(); ++i) {
        for(unsigned int I = 0; I < 3; ++I) {
          for(unsigned int J = 0; J < 3; ++J) {
            values[i](I, J) << scalar::Precision< T_Scalar >(std::real(a[i](I, J)(0, 0))), scalar::Precision< T_Scalar >(std::real(a[i](I, J)(0, 1))), scalar::Precision< T_Scalar >(std::real(a[i](I, J)(0, 2))), scalar::Precision< T_Scalar >(std::real(a[i](I, J)(1, 0))), scalar::Precision< T_Scalar >(std::real(a[i](I, J)(1, 1))), scalar::Precision< T_Scalar >(std::real(a[i](I, J)(1, 2))), scalar::Precision< T_Scalar >(std::real(a[i](I, J)(2, 0))), scalar::Precision< T_Scalar >(std::real(a[i](I, J)(2, 1))), scalar::Precision< T_Scalar >(std::real(a[i](I, J)(2, 2)));
          }
        }
      }
    }

    std::string name() const override
    {
      return "real part of a";
    }
  };

  //
  // RealToComplex
  //

  template< class T_Scalar, Degree T_Degree >
  class RealToComplex final : public ScalarTypeOperation< scalar::ComplexPrecision< T_Scalar >, T_Scalar, T_Degree >
  {
   public:
    RealToComplex()
    {
    }

    RealToComplex(const RealToComplex &other) :
      ScalarTypeOperation< scalar::ComplexPrecision< T_Scalar >, T_Scalar, T_Degree >(other)
    {
    }

    void operator()(OutputVector< scalar::ComplexPrecision< T_Scalar >, T_Degree > &values, OutputVector< T_Scalar, T_Degree > &a, const std::vector< scalar::Precision< T_Scalar >, numa::allocator< scalar::Precision< T_Scalar > > > &points, const std::vector< scalar::Precision< T_Scalar > > &gaussPoints, const int elementType, const std::pair< int, int > &entity) const
    {
#pragma omp for
      for(auto i = 0ULL; i < values.size(); ++i) {
        values[i] = scalar::ComplexPrecision< T_Scalar >(a[i]);
      }
    }

    std::string name() const override
    {
      return "build complex number a";
    }
  };

  template< class T_Scalar >
  class RealToComplex< T_Scalar, Degree::Degree1 > final : public ScalarTypeOperation< scalar::ComplexPrecision< T_Scalar >, T_Scalar, Degree::Degree1 >
  {
   public:
    RealToComplex()
    {
    }

    RealToComplex(const RealToComplex &other) :
      ScalarTypeOperation< scalar::ComplexPrecision< T_Scalar >, T_Scalar, Degree::Degree1 >(other)
    {
    }

    void operator()(OutputVector< scalar::ComplexPrecision< T_Scalar >, Degree::Degree1 > &values, OutputVector< T_Scalar, Degree::Degree1 > &a, const std::vector< scalar::Precision< T_Scalar >, numa::allocator< scalar::Precision< T_Scalar > > > &points, const std::vector< scalar::Precision< T_Scalar > > &gaussPoints, const int elementType, const std::pair< int, int > &entity) const
    {
#pragma omp for
      for(auto i = 0ULL; i < values.size(); ++i) {
        values[i] << scalar::ComplexPrecision< T_Scalar >(a[i][0]), scalar::ComplexPrecision< T_Scalar >(a[i][1]), scalar::ComplexPrecision< T_Scalar >(a[i][2]);
      }
    }

    std::string name() const override
    {
      return "build complex number a";
    }
  };

  template< class T_Scalar >
  class RealToComplex< T_Scalar, Degree::Degree2 > final : public ScalarTypeOperation< scalar::ComplexPrecision< T_Scalar >, T_Scalar, Degree::Degree2 >
  {
   public:
    RealToComplex()
    {
    }

    RealToComplex(const RealToComplex &other) :
      ScalarTypeOperation< scalar::ComplexPrecision< T_Scalar >, T_Scalar, Degree::Degree2 >(other)
    {
    }

    void operator()(OutputVector< scalar::ComplexPrecision< T_Scalar >, Degree::Degree2 > &values, OutputVector< T_Scalar, Degree::Degree2 > &a, const std::vector< scalar::Precision< T_Scalar >, numa::allocator< scalar::Precision< T_Scalar > > > &points, const std::vector< scalar::Precision< T_Scalar > > &gaussPoints, const int elementType, const std::pair< int, int > &entity) const
    {
#pragma omp for
      for(auto i = 0ULL; i < values.size(); ++i) {
        values[i] << scalar::ComplexPrecision< T_Scalar >(a[i](0, 0)), scalar::ComplexPrecision< T_Scalar >(a[i](0, 1)), scalar::ComplexPrecision< T_Scalar >(a[i](0, 2)), scalar::ComplexPrecision< T_Scalar >(a[i](1, 0)), scalar::ComplexPrecision< T_Scalar >(a[i](1, 1)), scalar::ComplexPrecision< T_Scalar >(a[i](1, 2)), scalar::ComplexPrecision< T_Scalar >(a[i](2, 0)), scalar::ComplexPrecision< T_Scalar >(a[i](2, 1)), scalar::ComplexPrecision< T_Scalar >(a[i](2, 2));
      }
    }

    std::string name() const override
    {
      return "build complex number a";
    }
  };

  template< class T_Scalar >
  class RealToComplex< T_Scalar, Degree::Degree4 > final : public ScalarTypeOperation< scalar::ComplexPrecision< T_Scalar >, T_Scalar, Degree::Degree4 >
  {
   public:
    RealToComplex()
    {
    }

    RealToComplex(const RealToComplex &other) :
      ScalarTypeOperation< scalar::ComplexPrecision< T_Scalar >, T_Scalar, Degree::Degree4 >(other)
    {
    }

    void operator()(OutputVector< scalar::ComplexPrecision< T_Scalar >, Degree::Degree4 > &values, OutputVector< T_Scalar, Degree::Degree4 > &a, const std::vector< scalar::Precision< T_Scalar >, numa::allocator< scalar::Precision< T_Scalar > > > &points, const std::vector< scalar::Precision< T_Scalar > > &gaussPoints, const int elementType, const std::pair< int, int > &entity) const
    {
#pragma omp for
      for(auto i = 0ULL; i < values.size(); ++i) {
        for(unsigned int I = 0; I < 3; ++I) {
          for(unsigned int J = 0; J < 3; ++J) {
            values[i](I, J) << scalar::ComplexPrecision< T_Scalar >(a[i](I, J)(0, 0)), scalar::ComplexPrecision< T_Scalar >(a[i](I, J)(0, 1)), scalar::ComplexPrecision< T_Scalar >(a[i](I, J)(0, 2)), scalar::ComplexPrecision< T_Scalar >(a[i](I, J)(1, 0)), scalar::ComplexPrecision< T_Scalar >(a[i](I, J)(1, 1)), scalar::ComplexPrecision< T_Scalar >(a[i](I, J)(1, 2)), scalar::ComplexPrecision< T_Scalar >(a[i](I, J)(2, 0)), scalar::ComplexPrecision< T_Scalar >(a[i](I, J)(2, 1)), scalar::ComplexPrecision< T_Scalar >(a[i](I, J)(2, 2));
          }
        }
      }
    }

    std::string name() const override
    {
      return "build complex number a";
    }
  };


} // namespace gmshfem::function


#endif // H_GMSHFEM_SCALAROPERATIONS
