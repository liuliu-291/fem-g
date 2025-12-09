// GmshFEM - Copyright (C) 2019-2022, A. Royer, E. Béchet, C. Geuzaine, Université de Liège
//
// See the LICENSE.txt file for license information. Please report all
// issues on https://gitlab.onelab.info/gmsh/fem/issues

#ifndef H_GMSHFEM_BINARYOPERATIONS
#define H_GMSHFEM_BINARYOPERATIONS

#include "BinaryNode.h"
#include "OmpInterface.h"
#include "numa.h"

namespace gmshfem::function
{


  // Available binary opertions:
  //  Add
  //  Atan2
  //  Cross
  //  Div
  //  Dyadic
  //  Hadamard
  //  Mul
  //  Sub
  //

  //
  // Add
  //

  template< class T_Scalar, Degree T_Degree >
  class Add final : public BinaryOperation< T_Scalar, T_Degree, T_Degree, T_Degree >
  {
   public:
    Add()
    {
    }

    Add(const Add &other) :
      BinaryOperation< T_Scalar, T_Degree, T_Degree, T_Degree >(other)
    {
    }

    template< class T_A, class T_B >
    void operator()(OutputVector< T_Scalar, T_Degree > &values, const InputVector< T_Scalar, T_Degree, T_A > &a, const InputVector< T_Scalar, T_Degree, T_B > &b, const std::vector< scalar::Precision< T_Scalar >, numa::allocator< scalar::Precision< T_Scalar > > > &points, const std::vector< scalar::Precision< T_Scalar > > &gaussPoints, const int elementType, const std::pair< int, int > &entity) const
    {
#pragma omp for
      for(auto i = 0ULL; i < values.size(); ++i) {
        values[i] = a[i] + b[i];
      }
    }

    bool canUseSameVectorsForOutputAndInputs() const override
    {
      return true;
    }

    std::string name() const override
    {
      return "a + b";
    }
  };

  //
  // Atan2
  //

  template< class T_Scalar, Degree T_Degree >
  class Atan2 final : public BinaryOperation< T_Scalar, T_Degree, T_Degree, T_Degree >
  {
   public:
    Atan2()
    {
    }

    Atan2(const Atan2 &other) :
      BinaryOperation< T_Scalar, T_Degree, T_Degree, T_Degree >(other)
    {
    }

    template< class T_A, class T_B >
    void operator()(OutputVector< T_Scalar, T_Degree > &values, const InputVector< T_Scalar, T_Degree, T_A > &a, const InputVector< T_Scalar, T_Degree, T_B > &b, const std::vector< scalar::Precision< T_Scalar >, numa::allocator< scalar::Precision< T_Scalar > > > &points, const std::vector< scalar::Precision< T_Scalar > > &gaussPoints, const int elementType, const std::pair< int, int > &entity) const
    {
#pragma omp for
      for(auto i = 0ULL; i < values.size(); ++i) {
        values[i] = std::atan2(a[i], b[i]);
      }
    }

    bool canUseSameVectorsForOutputAndInputs() const override
    {
      return true;
    }

    std::string name() const override
    {
      return "atan2(a, b)";
    }
  };

  //
  // Cross
  //

  template< class T_Scalar >
  class Cross final : public BinaryOperation< T_Scalar, Degree::Degree1, Degree::Degree1, Degree::Degree1 >
  {
   public:
    Cross()
    {
    }

    Cross(const Cross &other) :
      BinaryOperation< T_Scalar, Degree::Degree1, Degree::Degree1, Degree::Degree1 >(other)
    {
    }

    template< class T_A, class T_B >
    void operator()(OutputVector< T_Scalar, Degree::Degree1 > &values, const InputVector< T_Scalar, Degree::Degree1, T_A > &a, const InputVector< T_Scalar, Degree::Degree1, T_B > &b, const std::vector< scalar::Precision< T_Scalar >, numa::allocator< scalar::Precision< T_Scalar > > > &points, const std::vector< scalar::Precision< T_Scalar > > &gaussPoints, const int elementType, const std::pair< int, int > &entity) const
    {
#pragma omp for
      for(auto i = 0ULL; i < values.size(); ++i) {
        values[i] = (a[i].cross(b[i])).eval();
      }
    }

    bool canUseSameVectorsForOutputAndInputs() const override
    {
      return true;
    }

    std::string name() const override
    {
      return "cross product(a, b)";
    }
  };

  template< class T_ComplexScalar >
  class Cross< std::complex< T_ComplexScalar > > final : public BinaryOperation< std::complex< T_ComplexScalar >, Degree::Degree1, Degree::Degree1, Degree::Degree1 >
  {
   public:
    Cross()
    {
    }

    Cross(const Cross &other) :
      BinaryOperation< std::complex< T_ComplexScalar >, Degree::Degree1, Degree::Degree1, Degree::Degree1 >(other)
    {
    }

    template< class T_A, class T_B >
    void operator()(OutputVector< std::complex< T_ComplexScalar >, Degree::Degree1 > &values, const InputVector< std::complex< T_ComplexScalar >, Degree::Degree1, T_A > &a, const InputVector< std::complex< T_ComplexScalar >, Degree::Degree1, T_B > &b, const std::vector< scalar::Precision< std::complex< T_ComplexScalar > >, numa::allocator< scalar::Precision< std::complex< T_ComplexScalar > > > > &points, const std::vector< scalar::Precision< std::complex< T_ComplexScalar > > > &gaussPoints, const int elementType, const std::pair< int, int > &entity) const
    {
#pragma omp for
      for(auto i = 0ULL; i < values.size(); ++i) {
        values[i] = (a[i].cross(b[i]).conjugate()).eval();
      }
    }

    bool canUseSameVectorsForOutputAndInputs() const override
    {
      return true;
    }

    std::string name() const override
    {
      return "cross product(a, b)";
    }
  };

  //
  // Div
  //

  template< class T_Scalar, Degree T_Degree >
  class Div final : public BinaryOperation< T_Scalar, T_Degree, T_Degree, Degree::Degree0 >
  {
   public:
    Div()
    {
    }

    Div(const Div &other) :
      BinaryOperation< T_Scalar, T_Degree, T_Degree, Degree::Degree0 >(other)
    {
    }

    template< class T_A, class T_B >
    void operator()(OutputVector< T_Scalar, T_Degree > &values, const InputVector< T_Scalar, T_Degree, T_A > &a, const InputVector< T_Scalar, Degree::Degree0, T_B > &b, const std::vector< scalar::Precision< T_Scalar >, numa::allocator< scalar::Precision< T_Scalar > > > &points, const std::vector< scalar::Precision< T_Scalar > > &gaussPoints, const int elementType, const std::pair< int, int > &entity) const
    {
#pragma omp for
      for(auto i = 0ULL; i < values.size(); ++i) {
        values[i] = a[i] / b[i];
      }
    }

    bool canUseSameVectorsForOutputAndInputs() const override
    {
      return true;
    }

    std::string name() const override
    {
      return "a / b";
    }
  };

  template< class T_Scalar >
  class Div< T_Scalar, Degree::Degree4 > final : public BinaryOperation< T_Scalar, Degree::Degree4, Degree::Degree4, Degree::Degree0 >
  {
   public:
    Div()
    {
    }

    Div(const Div &other) :
      BinaryOperation< T_Scalar, Degree::Degree4, Degree::Degree4, Degree::Degree0 >(other)
    {
    }

    template< class T_A, class T_B >
    void operator()(OutputVector< T_Scalar, Degree::Degree4 > &values, const InputVector< T_Scalar, Degree::Degree4, T_A > &a, const InputVector< T_Scalar, Degree::Degree0, T_B > &b, const std::vector< scalar::Precision< T_Scalar >, numa::allocator< scalar::Precision< T_Scalar > > > &points, const std::vector< scalar::Precision< T_Scalar > > &gaussPoints, const int elementType, const std::pair< int, int > &entity) const
    {
#pragma omp for
      for(auto i = 0ULL; i < values.size(); ++i) {
        for(unsigned int I = 0; I < 3; ++I) {
          for(unsigned int J = 0; J < 3; ++J) {
            values[i](I, J) = a[i](I, J) / b[i];
          }
        }
      }
    }

    bool canUseSameVectorsForOutputAndInputs() const override
    {
      return true;
    }

    std::string name() const override
    {
      return "a / b";
    }
  };

  //
  // Dyadic
  //

  template< class T_Scalar, Degree T_Degree >
  class Dyadic;
  
  template< class T_Scalar >
  class Dyadic< T_Scalar, Degree::Degree1 > final : public BinaryOperation< T_Scalar, Degree::Degree2, Degree::Degree1, Degree::Degree1 >
  {
   public:
    Dyadic()
    {
    }

    Dyadic(const Dyadic &other) :
      BinaryOperation< T_Scalar, Degree::Degree2, Degree::Degree1, Degree::Degree1 >(other)
    {
    }

    template< class T_A, class T_B >
    void operator()(OutputVector< T_Scalar, Degree::Degree2 > &values, const InputVector< T_Scalar, Degree::Degree1, T_A > &a, const InputVector< T_Scalar, Degree::Degree1, T_B > &b, const std::vector< scalar::Precision< T_Scalar >, numa::allocator< scalar::Precision< T_Scalar > > > &points, const std::vector< scalar::Precision< T_Scalar > > &gaussPoints, const int elementType, const std::pair< int, int > &entity) const
    {
#pragma omp for
      for(auto i = 0ULL; i < values.size(); ++i) {
        values[i] = a[i] * b[i].transpose();
      }
    }

    bool canUseSameVectorsForOutputAndInputs() const override
    {
      return false;
    }

    std::string name() const override
    {
      return "dyadic(a, b)";
    }
  };
  
  template< class T_Scalar >
  class Dyadic< T_Scalar, Degree::Degree2 > final : public BinaryOperation< T_Scalar, Degree::Degree4, Degree::Degree2, Degree::Degree2 >
  {
   public:
    Dyadic()
    {
    }

    Dyadic(const Dyadic &other) :
      BinaryOperation< T_Scalar, Degree::Degree4, Degree::Degree2, Degree::Degree2 >(other)
    {
    }

    template< class T_A, class T_B >
    void operator()(OutputVector< T_Scalar, Degree::Degree4 > &values, const InputVector< T_Scalar, Degree::Degree2, T_A > &a, const InputVector< T_Scalar, Degree::Degree2, T_B > &b, const std::vector< scalar::Precision< T_Scalar >, numa::allocator< scalar::Precision< T_Scalar > > > &points, const std::vector< scalar::Precision< T_Scalar > > &gaussPoints, const int elementType, const std::pair< int, int > &entity) const
    {
#pragma omp for
      for(auto i = 0ULL; i < values.size(); ++i) {
        for(unsigned int I = 0; I < 3; ++I) {
          for(unsigned int J = 0; J < 3; ++J) {
            for(unsigned int K = 0; K < 3; ++K) {
              for(unsigned int L = 0; L < 3; ++L) {
                values[i](I, J)(K, L) = a[i](I, J) * b[i](K, L);
              }
            }
          }
        }
      }
    }

    bool canUseSameVectorsForOutputAndInputs() const override
    {
      return false;
    }

    std::string name() const override
    {
      return "dyadic(a, b)";
    }
  };

  //
  // Hadamard
  //

  template< class T_Scalar, Degree T_Degree >
  class Hadamard final : public BinaryOperation< T_Scalar, T_Degree, T_Degree, T_Degree >
  {
   public:
    Hadamard()
    {
    }

    Hadamard(const Hadamard &other) :
      BinaryOperation< T_Scalar, T_Degree, T_Degree, T_Degree >(other)
    {
    }

    template< class T_A, class T_B >
    void operator()(OutputVector< T_Scalar, T_Degree > &values, const InputVector< T_Scalar, T_Degree, T_A > &a, const InputVector< T_Scalar, T_Degree, T_B > &b, const std::vector< scalar::Precision< T_Scalar >, numa::allocator< scalar::Precision< T_Scalar > > > &points, const std::vector< scalar::Precision< T_Scalar > > &gaussPoints, const int elementType, const std::pair< int, int > &entity) const
    {
#pragma omp for
      for(auto i = 0ULL; i < values.size(); ++i) {
        values[i] = a[i].cwiseProduct(b[i]);
      }
    }

    bool canUseSameVectorsForOutputAndInputs() const override
    {
      return false;
    }

    std::string name() const override
    {
      return "hadamard(a, b)";
    }
  };

  template< class T_Scalar >
  class Hadamard< T_Scalar, Degree::Degree0 > final : public BinaryOperation< T_Scalar, Degree::Degree0, Degree::Degree0, Degree::Degree0 >
  {
   public:
    Hadamard()
    {
    }

    Hadamard(const Hadamard &other) :
      BinaryOperation< T_Scalar, Degree::Degree0, Degree::Degree0, Degree::Degree0 >(other)
    {
    }

    template< class T_A, class T_B >
    void operator()(OutputVector< T_Scalar, Degree::Degree0 > &values, const InputVector< T_Scalar, Degree::Degree0, T_A > &a, const InputVector< T_Scalar, Degree::Degree0, T_B > &b, const std::vector< scalar::Precision< T_Scalar >, numa::allocator< scalar::Precision< T_Scalar > > > &points, const std::vector< scalar::Precision< T_Scalar > > &gaussPoints, const int elementType, const std::pair< int, int > &entity) const
    {
#pragma omp for
      for(auto i = 0ULL; i < values.size(); ++i) {
        values[i] = a[i] * b[i];
      }
    }

    bool canUseSameVectorsForOutputAndInputs() const override
    {
      return true;
    }

    std::string name() const override
    {
      return "hadamard(a, b)";
    }
  };

  template< class T_Scalar >
  class Hadamard< T_Scalar, Degree::Degree4 > final : public BinaryOperation< T_Scalar, Degree::Degree4, Degree::Degree4, Degree::Degree4 >
  {
   public:
    Hadamard()
    {
    }

    Hadamard(const Hadamard &other) :
      BinaryOperation< T_Scalar, Degree::Degree4, Degree::Degree4, Degree::Degree4 >(other)
    {
    }

    template< class T_A, class T_B >
    void operator()(OutputVector< T_Scalar, Degree::Degree4 > &values, const InputVector< T_Scalar, Degree::Degree4, T_A > &a, const InputVector< T_Scalar, Degree::Degree4, T_B > &b, const std::vector< scalar::Precision< T_Scalar >, numa::allocator< scalar::Precision< T_Scalar > > > &points, const std::vector< scalar::Precision< T_Scalar > > &gaussPoints, const int elementType, const std::pair< int, int > &entity) const
    {
#pragma omp for
      for(auto i = 0ULL; i < values.size(); ++i) {
        for(unsigned int J = 0; J < 3; ++J) {
          for(unsigned int I = 0; I < 3; ++I) {
            values[i](I, J) = a[i](I, J).cwiseProduct(b[i](I, J));
          }
        }
      }
    }

    bool canUseSameVectorsForOutputAndInputs() const override
    {
      return false;
    }

    std::string name() const override
    {
      return "hadamard(a, b)";
    }
  };

  //
  // Mul
  //

  template< class T_Scalar, Degree T_DegreeA, Degree T_DegreeB >
  class Mul final : public BinaryOperation< T_Scalar, MathProduct< T_DegreeA, T_DegreeB >::ProductDegree, T_DegreeA, T_DegreeB >
  {
   public:
    Mul()
    {
    }

    Mul(const Mul &other) :
      BinaryOperation< T_Scalar, MathProduct< T_DegreeA, T_DegreeB >::ProductDegree, T_DegreeA, T_DegreeB >(other)
    {
    }

    template< class T_A, class T_B >
    void operator()(OutputVector< T_Scalar, MathProduct< T_DegreeA, T_DegreeB >::ProductDegree > &values, const InputVector< T_Scalar, T_DegreeA, T_A > &a, const InputVector< T_Scalar, T_DegreeB, T_B > &b, const std::vector< scalar::Precision< T_Scalar >, numa::allocator< scalar::Precision< T_Scalar > > > &points, const std::vector< scalar::Precision< T_Scalar > > &gaussPoints, const int elementType, const std::pair< int, int > &entity) const
    {
#pragma omp for
      for(auto i = 0ULL; i < values.size(); ++i) {
        values[i] = a[i] * b[i];
      }
    }

    bool canUseSameVectorsForOutputAndInputs() const override
    {
      return true;
    }

    std::string name() const override
    {
      return "a * b";
    }
  };

  template< class T_Scalar >
  class Mul< T_Scalar, Degree::Degree0, Degree::Degree4 > final : public BinaryOperation< T_Scalar, MathProduct< Degree::Degree0, Degree::Degree4 >::ProductDegree, Degree::Degree0, Degree::Degree4 >
  {
   public:
    Mul()
    {
    }

    Mul(const Mul &other) :
      BinaryOperation< T_Scalar, MathProduct< Degree::Degree0, Degree::Degree4 >::ProductDegree, Degree::Degree0, Degree::Degree4 >(other)
    {
    }

    template< class T_A, class T_B >
    void operator()(OutputVector< T_Scalar, MathProduct< Degree::Degree0, Degree::Degree4 >::ProductDegree > &values, const InputVector< T_Scalar, Degree::Degree0, T_A > &a, const InputVector< T_Scalar, Degree::Degree4, T_B > &b, const std::vector< scalar::Precision< T_Scalar >, numa::allocator< scalar::Precision< T_Scalar > > > &points, const std::vector< scalar::Precision< T_Scalar > > &gaussPoints, const int elementType, const std::pair< int, int > &entity) const
    {
#pragma omp for
      for(auto i = 0ULL; i < values.size(); ++i) {
        for(unsigned int J = 0; J < 3; ++J) {
          for(unsigned int I = 0; I < 3; ++I) {
            values[i](I, J) = a[i] * b[i](I, J);
          }
        }
      }
    }

    bool canUseSameVectorsForOutputAndInputs() const override
    {
      return true;
    }

    std::string name() const override
    {
      return "a * b";
    }
  };

  template< class T_Scalar >
  class Mul< T_Scalar, Degree::Degree1, Degree::Degree1 > final : public BinaryOperation< T_Scalar, MathProduct< Degree::Degree1, Degree::Degree1 >::ProductDegree, Degree::Degree1, Degree::Degree1 >
  {
   public:
    Mul()
    {
    }

    Mul(const Mul &other) :
      BinaryOperation< T_Scalar, MathProduct< Degree::Degree1, Degree::Degree1 >::ProductDegree, Degree::Degree1, Degree::Degree1 >(other)
    {
    }

    template< class T_A, class T_B >
    void operator()(OutputVector< T_Scalar, MathProduct< Degree::Degree1, Degree::Degree1 >::ProductDegree > &values, const InputVector< T_Scalar, Degree::Degree1, T_A > &a, const InputVector< T_Scalar, Degree::Degree1, T_B > &b, const std::vector< scalar::Precision< T_Scalar >, numa::allocator< scalar::Precision< T_Scalar > > > &points, const std::vector< scalar::Precision< T_Scalar > > &gaussPoints, const int elementType, const std::pair< int, int > &entity) const
    {
#pragma omp for
      for(auto i = 0ULL; i < values.size(); ++i) {
        values[i] = (a[i].transpose() * b[i])(0);
      }
    }

    std::string name() const override
    {
      return "a * b";
    }
  };

  template< class T_Scalar >
  class Mul< T_Scalar, Degree::Degree1, Degree::Degree2 > final : public BinaryOperation< T_Scalar, MathProduct< Degree::Degree1, Degree::Degree2 >::ProductDegree, Degree::Degree1, Degree::Degree2 >
  {
   public:
    Mul()
    {
    }

    Mul(const Mul &other) :
      BinaryOperation< T_Scalar, MathProduct< Degree::Degree1, Degree::Degree2 >::ProductDegree, Degree::Degree1, Degree::Degree2 >(other)
    {
    }

    template< class T_A, class T_B >
    void operator()(OutputVector< T_Scalar, MathProduct< Degree::Degree1, Degree::Degree2 >::ProductDegree > &values, const InputVector< T_Scalar, Degree::Degree1, T_A > &a, const InputVector< T_Scalar, Degree::Degree2, T_B > &b, const std::vector< scalar::Precision< T_Scalar >, numa::allocator< scalar::Precision< T_Scalar > > > &points, const std::vector< scalar::Precision< T_Scalar > > &gaussPoints, const int elementType, const std::pair< int, int > &entity) const
    {
#pragma omp for
      for(auto i = 0ULL; i < values.size(); ++i) {
        values[i] = (a[i].transpose() * b[i]).eval();
      }
    }

    bool canUseSameVectorsForOutputAndInputs() const override
    {
      return true;
    }

    std::string name() const override
    {
      return "a * b";
    }
  };

  template< class T_Scalar >
  class Mul< T_Scalar, Degree::Degree2, Degree::Degree2 > final : public BinaryOperation< T_Scalar, MathProduct< Degree::Degree2, Degree::Degree2 >::ProductDegree, Degree::Degree2, Degree::Degree2 >
  {
   public:
    Mul()
    {
    }

    Mul(const Mul &other) :
      BinaryOperation< T_Scalar, MathProduct< Degree::Degree2, Degree::Degree2 >::ProductDegree, Degree::Degree2, Degree::Degree2 >(other)
    {
    }

    template< class T_A, class T_B >
    void operator()(OutputVector< T_Scalar, MathProduct< Degree::Degree2, Degree::Degree2 >::ProductDegree > &values, const InputVector< T_Scalar, Degree::Degree2, T_A > &a, const InputVector< T_Scalar, Degree::Degree2, T_B > &b, const std::vector< scalar::Precision< T_Scalar >, numa::allocator< scalar::Precision< T_Scalar > > > &points, const std::vector< scalar::Precision< T_Scalar > > &gaussPoints, const int elementType, const std::pair< int, int > &entity) const
    {
#pragma omp for
      for(auto i = 0ULL; i < values.size(); ++i) {
        values[i] = (a[i] * b[i]).eval();
      }
    }

    bool canUseSameVectorsForOutputAndInputs() const override
    {
      return true;
    }

    std::string name() const override
    {
      return "a * b";
    }
  };

  template< class T_Scalar >
  class Mul< T_Scalar, Degree::Degree4, Degree::Degree0 > final : public BinaryOperation< T_Scalar, MathProduct< Degree::Degree4, Degree::Degree0 >::ProductDegree, Degree::Degree4, Degree::Degree0 >
  {
   public:
    Mul()
    {
    }

    Mul(const Mul &other) :
      BinaryOperation< T_Scalar, MathProduct< Degree::Degree4, Degree::Degree0 >::ProductDegree, Degree::Degree4, Degree::Degree0 >(other)
    {
    }

    template< class T_A, class T_B >
    void operator()(OutputVector< T_Scalar, MathProduct< Degree::Degree4, Degree::Degree0 >::ProductDegree > &values, const InputVector< T_Scalar, Degree::Degree4, T_A > &a, const InputVector< T_Scalar, Degree::Degree0, T_B > &b, const std::vector< scalar::Precision< T_Scalar >, numa::allocator< scalar::Precision< T_Scalar > > > &points, const std::vector< scalar::Precision< T_Scalar > > &gaussPoints, const int elementType, const std::pair< int, int > &entity) const
    {
#pragma omp for
      for(auto i = 0ULL; i < values.size(); ++i) {
        typename MathObject< T_Scalar, Degree::Degree4 >::Object tmp;
        for(unsigned int J = 0; J < 3; ++J) {
          for(unsigned int I = 0; I < 3; ++I) {
            tmp(I, J) = a[i](I, J) * b[i];
          }
        }
        values[i] = tmp;
      }
    }

    bool canUseSameVectorsForOutputAndInputs() const override
    {
      return true;
    }

    std::string name() const override
    {
      return "a * b";
    }
  };

  template< class T_Scalar >
  class Mul< T_Scalar, Degree::Degree4, Degree::Degree2 > final : public BinaryOperation< T_Scalar, MathProduct< Degree::Degree4, Degree::Degree2 >::ProductDegree, Degree::Degree4, Degree::Degree2 >
  {
   public:
    Mul()
    {
    }

    Mul(const Mul &other) :
      BinaryOperation< T_Scalar, MathProduct< Degree::Degree4, Degree::Degree2 >::ProductDegree, Degree::Degree4, Degree::Degree2 >(other)
    {
    }

    template< class T_A, class T_B >
    void operator()(OutputVector< T_Scalar, MathProduct< Degree::Degree4, Degree::Degree2 >::ProductDegree > &values, const InputVector< T_Scalar, Degree::Degree4, T_A > &a, const InputVector< T_Scalar, Degree::Degree2, T_B > &b, const std::vector< scalar::Precision< T_Scalar >, numa::allocator< scalar::Precision< T_Scalar > > > &points, const std::vector< scalar::Precision< T_Scalar > > &gaussPoints, const int elementType, const std::pair< int, int > &entity) const
    {
#pragma omp for
      for(auto i = 0ULL; i < values.size(); ++i) {
        typename MathObject< T_Scalar, Degree::Degree2 >::Object tmp;
        for(unsigned int J = 0; J < 3; ++J) {
          for(unsigned int I = 0; I < 3; ++I) {
            tmp(I, J) = (a[i](I, J) * b[i].transpose()).trace();
          }
        }
        values[i] = tmp;
      }
    }

    bool canUseSameVectorsForOutputAndInputs() const override
    {
      return true;
    }

    std::string name() const override
    {
      return "a * b";
    }
  };
  
  template< class T_Scalar >
  class Mul< T_Scalar, Degree::Degree4, Degree::Degree4 > final : public BinaryOperation< T_Scalar, MathProduct< Degree::Degree4, Degree::Degree4 >::ProductDegree, Degree::Degree4, Degree::Degree4 >
  {
   public:
    Mul()
    {
    }

    Mul(const Mul &other) :
      BinaryOperation< T_Scalar, MathProduct< Degree::Degree4, Degree::Degree4 >::ProductDegree, Degree::Degree4, Degree::Degree4 >(other)
    {
    }

    template< class T_A, class T_B >
    void operator()(OutputVector< T_Scalar, MathProduct< Degree::Degree4, Degree::Degree4 >::ProductDegree > &values, const InputVector< T_Scalar, Degree::Degree4, T_A > &a, const InputVector< T_Scalar, Degree::Degree4, T_B > &b, const std::vector< scalar::Precision< T_Scalar >, numa::allocator< scalar::Precision< T_Scalar > > > &points, const std::vector< scalar::Precision< T_Scalar > > &gaussPoints, const int elementType, const std::pair< int, int > &entity) const
    {
#pragma omp for
      for(auto i = 0ULL; i < values.size(); ++i) {
        typename MathObject< T_Scalar, Degree::Degree4 >::Object tmp;
        
        for(unsigned int I = 0; I < 3; ++I) {
          for(unsigned int J = 0; J < 3; ++J) {
            for(unsigned int M = 0; M < 3; ++M) {
              for(unsigned int N = 0; N < 3; ++N) {
                tmp(I, J)(M, N) = 0.;
                for(unsigned int K = 0; K < 3; ++K) {
                  for(unsigned int L = 0; L < 3; ++L) {
                    tmp(I, J)(M, N) += a[i](I, J)(K, L) * b[i](K, L)(M, N);
                  }
                }
              }
            }
          }
        }

        values[i] = tmp;
      }
    }

    bool canUseSameVectorsForOutputAndInputs() const override
    {
      return true;
    }

    std::string name() const override
    {
      return "a * b";
    }
  };

  //
  // Sub
  //

  template< class T_Scalar, Degree T_Degree >
  class Sub final : public BinaryOperation< T_Scalar, T_Degree, T_Degree, T_Degree >
  {
   public:
    Sub()
    {
    }

    Sub(const Sub &other) :
      BinaryOperation< T_Scalar, T_Degree, T_Degree, T_Degree >(other)
    {
    }

    template< class T_A, class T_B >
    void operator()(OutputVector< T_Scalar, T_Degree > &values, const InputVector< T_Scalar, T_Degree, T_A > &a, const InputVector< T_Scalar, T_Degree, T_B > &b, const std::vector< scalar::Precision< T_Scalar >, numa::allocator< scalar::Precision< T_Scalar > > > &points, const std::vector< scalar::Precision< T_Scalar > > &gaussPoints, const int elementType, const std::pair< int, int > &entity) const
    {
#pragma omp for
      for(auto i = 0ULL; i < values.size(); ++i) {
        values[i] = a[i] - b[i];
      }
    }

    bool canUseSameVectorsForOutputAndInputs() const override
    {
      return true;
    }

    std::string name() const override
    {
      return "a - b";
    }
  };


} // namespace gmshfem::function


#endif // H_GMSHFEM_BINARYOPERATIONS
