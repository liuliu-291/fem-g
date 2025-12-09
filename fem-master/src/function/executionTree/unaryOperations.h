// GmshFEM - Copyright (C) 2019-2022, A. Royer, E. Béchet, C. Geuzaine, Université de Liège
//
// See the LICENSE.txt file for license information. Please report all
// issues on https://gitlab.onelab.info/gmsh/fem/issues

#ifndef H_GMSHFEM_UNARYOPERATIONS
#define H_GMSHFEM_UNARYOPERATIONS

#include "Message.h"
#include "OmpInterface.h"
#include "UnaryNode.h"

#ifdef HAVE_BESSEL
#include "Bessel.h"
#endif // HAVE_BESSEL

namespace gmshfem::function
{


  // Available unary operations:
  //  Abs
  //  Acos
  //  Acosh
  //  AngularComp
  //  Asin
  //  Asinh
  //  Atan
  //  Atanh
  //  Cbrt
  //  Commutator
  //  Conj
  //  Cos
  //  Cosh
  //  CylBesselJ (need ENABLE_BESSEL)
  //  CylNeumann (need ENABLE_BESSEL)
  //  Exp
  //  Heaviside
  //  Imag
  //  Inv
  //  Log
  //  Ln
  //  Minus
  //  Norm
  //  Plus
  //  Pow
  //  R2Dcomp
  //  Real
  //  Sin
  //  Sinh
  //  Sqrt
  //  Tan
  //  Tanh
  //  Trace
  //  Xcomp
  //  XXcomp
  //  XYcomp
  //  XZcomp
  //  Ycomp
  //  YXcomp
  //  YYcomp
  //  YZcomp
  //  Zcomp
  //  ZXcomp
  //  ZYcomp
  //  ZZcomp
  //

  //
  // Abs
  //

  template< class T_Scalar, Degree T_Degree >
  class Abs final : public UnaryOperation< T_Scalar, T_Degree, T_Degree >
  {
   public:
    Abs()
    {
    }

    Abs(const Abs &other) :
      UnaryOperation< T_Scalar, T_Degree, T_Degree >(other)
    {
    }

    template< class T_A >
    void operator()(OutputVector< T_Scalar, T_Degree > &values, const InputVector< T_Scalar, T_Degree, T_A > &a, const std::vector< scalar::Precision< T_Scalar >, numa::allocator< scalar::Precision< T_Scalar > > > &points, const std::vector< scalar::Precision< T_Scalar > > &gaussPoints, const int elementType, const std::pair< int, int > &entity) const
    {
#pragma omp for
      for(auto i = 0ULL; i < values.size(); ++i) {
        values[i] = std::abs(a[i]);
      }
    }

    bool canUseSameVectorsForOutputAndInputs() const override
    {
      return true;
    }

    std::string name() const override
    {
      return "abs(a)";
    }
  };

  //
  // Acos
  //

  template< class T_Scalar, Degree T_Degree >
  class Acos final : public UnaryOperation< T_Scalar, T_Degree, T_Degree >
  {
   public:
    Acos()
    {
    }

    Acos(const Acos &other) :
      UnaryOperation< T_Scalar, T_Degree, T_Degree >(other)
    {
    }

    template< class T_A >
    void operator()(OutputVector< T_Scalar, T_Degree > &values, const InputVector< T_Scalar, T_Degree, T_A > &a, const std::vector< scalar::Precision< T_Scalar >, numa::allocator< scalar::Precision< T_Scalar > > > &points, const std::vector< scalar::Precision< T_Scalar > > &gaussPoints, const int elementType, const std::pair< int, int > &entity) const
    {
#pragma omp for
      for(auto i = 0ULL; i < values.size(); ++i) {
        values[i] = std::acos(a[i]);
      }
    }

    bool canUseSameVectorsForOutputAndInputs() const override
    {
      return true;
    }

    std::string name() const override
    {
      return "acos(a)";
    }
  };

  //
  // Acosh
  //

  template< class T_Scalar, Degree T_Degree >
  class Acosh final : public UnaryOperation< T_Scalar, T_Degree, T_Degree >
  {
   public:
    Acosh()
    {
    }

    Acosh(const Acosh &other) :
      UnaryOperation< T_Scalar, T_Degree, T_Degree >(other)
    {
    }

    template< class T_A >
    void operator()(OutputVector< T_Scalar, T_Degree > &values, const InputVector< T_Scalar, T_Degree, T_A > &a, const std::vector< scalar::Precision< T_Scalar >, numa::allocator< scalar::Precision< T_Scalar > > > &points, const std::vector< scalar::Precision< T_Scalar > > &gaussPoints, const int elementType, const std::pair< int, int > &entity) const
    {
#pragma omp for
      for(auto i = 0ULL; i < values.size(); ++i) {
        values[i] = std::acosh(a[i]);
      }
    }

    bool canUseSameVectorsForOutputAndInputs() const override
    {
      return true;
    }

    std::string name() const override
    {
      return "acosh(a)";
    }
  };

  //
  // AngularComp
  //

  template< class T_Scalar >
  class AngularComp final : public UnaryOperation< T_Scalar, Degree::Degree0, Degree::Degree1 >
  {
   public:
    AngularComp()
    {
    }

    AngularComp(const AngularComp &other) :
      UnaryOperation< T_Scalar, Degree::Degree0, Degree::Degree1 >(other)
    {
    }

    template< class T_A >
    void operator()(OutputVector< T_Scalar, Degree::Degree0 > &values, const InputVector< T_Scalar, Degree::Degree1, T_A > &a, const std::vector< scalar::Precision< T_Scalar >, numa::allocator< scalar::Precision< T_Scalar > > > &points, const std::vector< scalar::Precision< T_Scalar > > &gaussPoints, const int elementType, const std::pair< int, int > &entity) const
    {
#pragma omp for
      for(auto i = 0ULL; i < values.size(); ++i) {
        values[i] = std::atan2(a[i](1), a[i](0));
      }
    }

    std::string name() const override
    {
      return "angular component(a)";
    }
  };

  template< class T_ComplexScalar >
  class AngularComp< std::complex< T_ComplexScalar > > final : public UnaryOperation< std::complex< T_ComplexScalar >, Degree::Degree0, Degree::Degree1 >
  {
   public:
    AngularComp()
    {
    }

    AngularComp(const AngularComp &other) :
      UnaryOperation< std::complex< T_ComplexScalar >, Degree::Degree0, Degree::Degree1 >(other)
    {
    }

    template< class T_A >
    void operator()(OutputVector< std::complex< T_ComplexScalar >, Degree::Degree0 > &values, const InputVector< std::complex< T_ComplexScalar >, Degree::Degree1, T_A > &a, const std::vector< scalar::Precision< std::complex< T_ComplexScalar > >, numa::allocator< scalar::Precision< std::complex< T_ComplexScalar > > > > &points, const std::vector< scalar::Precision< std::complex< T_ComplexScalar > > > &gaussPoints, const int elementType, const std::pair< int, int > &entity) const
    {
#pragma omp for
      for(auto i = 0ULL; i < values.size(); ++i) {
        values[i] = std::complex< T_ComplexScalar >(std::atan2(std::real(a[i](1)), std::real(a[i](0))), std::atan2(std::imag(a[i](1)), std::imag(a[i](0))));
      }
    }

    std::string name() const override
    {
      return "angular component(a)";
    }
  };

  //
  // Asin
  //

  template< class T_Scalar, Degree T_Degree >
  class Asin final : public UnaryOperation< T_Scalar, T_Degree, T_Degree >
  {
   public:
    Asin()
    {
    }

    Asin(const Asin &other) :
      UnaryOperation< T_Scalar, T_Degree, T_Degree >(other)
    {
    }

    template< class T_A >
    void operator()(OutputVector< T_Scalar, T_Degree > &values, const InputVector< T_Scalar, T_Degree, T_A > &a, const std::vector< scalar::Precision< T_Scalar >, numa::allocator< scalar::Precision< T_Scalar > > > &points, const std::vector< scalar::Precision< T_Scalar > > &gaussPoints, const int elementType, const std::pair< int, int > &entity) const
    {
#pragma omp for
      for(auto i = 0ULL; i < values.size(); ++i) {
        values[i] = std::asin(a[i]);
      }
    }

    bool canUseSameVectorsForOutputAndInputs() const override
    {
      return true;
    }

    std::string name() const override
    {
      return "asin(a)";
    }
  };

  //
  // Asinh
  //

  template< class T_Scalar, Degree T_Degree >
  class Asinh final : public UnaryOperation< T_Scalar, T_Degree, T_Degree >
  {
   public:
    Asinh()
    {
    }

    Asinh(const Asinh &other) :
      UnaryOperation< T_Scalar, T_Degree, T_Degree >(other)
    {
    }

    template< class T_A >
    void operator()(OutputVector< T_Scalar, T_Degree > &values, const InputVector< T_Scalar, T_Degree, T_A > &a, const std::vector< scalar::Precision< T_Scalar >, numa::allocator< scalar::Precision< T_Scalar > > > &points, const std::vector< scalar::Precision< T_Scalar > > &gaussPoints, const int elementType, const std::pair< int, int > &entity) const
    {
#pragma omp for
      for(auto i = 0ULL; i < values.size(); ++i) {
        values[i] = std::asinh(a[i]);
      }
    }

    bool canUseSameVectorsForOutputAndInputs() const override
    {
      return true;
    }

    std::string name() const override
    {
      return "asinh(a)";
    }
  };

  //
  // Atan
  //

  template< class T_Scalar, Degree T_Degree >
  class Atan final : public UnaryOperation< T_Scalar, T_Degree, T_Degree >
  {
   public:
    Atan()
    {
    }

    Atan(const Atan &other) :
      UnaryOperation< T_Scalar, T_Degree, T_Degree >(other)
    {
    }

    template< class T_A >
    void operator()(OutputVector< T_Scalar, T_Degree > &values, const InputVector< T_Scalar, T_Degree, T_A > &a, const std::vector< scalar::Precision< T_Scalar >, numa::allocator< scalar::Precision< T_Scalar > > > &points, const std::vector< scalar::Precision< T_Scalar > > &gaussPoints, const int elementType, const std::pair< int, int > &entity) const
    {
#pragma omp for
      for(auto i = 0ULL; i < values.size(); ++i) {
        values[i] = std::atan(a[i]);
      }
    }

    bool canUseSameVectorsForOutputAndInputs() const override
    {
      return true;
    }

    std::string name() const override
    {
      return "atan(a)";
    }
  };

  //
  // Atanh
  //

  template< class T_Scalar, Degree T_Degree >
  class Atanh final : public UnaryOperation< T_Scalar, T_Degree, T_Degree >
  {
   public:
    Atanh()
    {
    }

    Atanh(const Atanh &other) :
      UnaryOperation< T_Scalar, T_Degree, T_Degree >(other)
    {
    }

    template< class T_A >
    void operator()(OutputVector< T_Scalar, T_Degree > &values, const InputVector< T_Scalar, T_Degree, T_A > &a, const std::vector< scalar::Precision< T_Scalar >, numa::allocator< scalar::Precision< T_Scalar > > > &points, const std::vector< scalar::Precision< T_Scalar > > &gaussPoints, const int elementType, const std::pair< int, int > &entity) const
    {
#pragma omp for
      for(auto i = 0ULL; i < values.size(); ++i) {
        values[i] = std::atanh(a[i]);
      }
    }

    bool canUseSameVectorsForOutputAndInputs() const override
    {
      return true;
    }

    std::string name() const override
    {
      return "atanh(a)";
    }
  };

  //
  // Cbrt
  //

  template< class T_Scalar, Degree T_Degree >
  class Cbrt final : public UnaryOperation< T_Scalar, T_Degree, T_Degree >
  {
   public:
    Cbrt()
    {
    }

    Cbrt(const Cbrt &other) :
      UnaryOperation< T_Scalar, T_Degree, T_Degree >(other)
    {
    }

    template< class T_A >
    void operator()(OutputVector< T_Scalar, T_Degree > &values, const InputVector< T_Scalar, T_Degree, T_A > &a, const std::vector< scalar::Precision< T_Scalar >, numa::allocator< scalar::Precision< T_Scalar > > > &points, const std::vector< scalar::Precision< T_Scalar > > &gaussPoints, const int elementType, const std::pair< int, int > &entity) const
    {
#pragma omp for
      for(auto i = 0ULL; i < values.size(); ++i) {
        values[i] = std::cbrt(a[i]);
      }
    }

    bool canUseSameVectorsForOutputAndInputs() const override
    {
      return true;
    }

    std::string name() const override
    {
      return "cbrt(a)";
    }
  };

  template< class T_ComplexScalar, Degree T_Degree >
  class Cbrt< std::complex< T_ComplexScalar >, T_Degree > final : public UnaryOperation< std::complex< T_ComplexScalar >, T_Degree, T_Degree >
  {
   public:
    Cbrt()
    {
    }

    Cbrt(const Cbrt &other) :
      UnaryOperation< std::complex< T_ComplexScalar >, T_Degree, T_Degree >(other)
    {
    }

    template< class T_A >
    void operator()(OutputVector< std::complex< T_ComplexScalar >, T_Degree > &values, const InputVector< std::complex< T_ComplexScalar >, T_Degree, T_A > &a, const std::vector< scalar::Precision< std::complex< T_ComplexScalar > >, numa::allocator< scalar::Precision< std::complex< T_ComplexScalar > > > > &points, const std::vector< scalar::Precision< std::complex< T_ComplexScalar > > > &gaussPoints, const int elementType, const std::pair< int, int > &entity) const
    {
#pragma omp for
      for(auto i = 0ULL; i < values.size(); ++i) {
        values[i] = std::pow(a[i], 1. / 3.);
      }
    }

    bool canUseSameVectorsForOutputAndInputs() const override
    {
      return true;
    }

    std::string name() const override
    {
      return "cbrt(a)";
    }
  };

  //
  // Commutator
  //

  template< class T_Scalar, Degree T_Degree >
  class Commutator final : public UnaryOperation< T_Scalar, T_Degree, T_Degree >
  {
   private:
    bool *_state;
    unsigned int _instance;

   public:
    Commutator(bool *state) :
      _state(state)
    {
    }

    Commutator(const Commutator &other) :
      UnaryOperation< T_Scalar, T_Degree, T_Degree >(other), _state(other._state)
    {
    }

    template< class T_A >
    void operator()(OutputVector< T_Scalar, T_Degree > &values, const InputVector< T_Scalar, T_Degree, T_A > &a, const std::vector< scalar::Precision< T_Scalar >, numa::allocator< scalar::Precision< T_Scalar > > > &points, const std::vector< scalar::Precision< T_Scalar > > &gaussPoints, const int elementType, const std::pair< int, int > &entity) const
    {
      if(*_state) {
#pragma omp for
        for(auto i = 0ULL; i < values.size(); ++i) {
          values[i] = a[i];
        }
      }
      else {
#pragma omp for
        for(auto i = 0ULL; i < values.size(); ++i) {
          MathObject< T_Scalar, T_Degree >::zero(values[i]);
        }
      }
    }

    bool canUseSameVectorsForOutputAndInputs() const override
    {
      return true;
    }

    std::string name() const override
    {
      return "commutator(a)";
    }

    bool operator==(const UnaryOperation< T_Scalar, T_Degree, T_Degree > &other) const override
    {
      if(_state == static_cast< const Commutator & >(other)._state) {
        return true;
      }
      return false;
    }
  };

  //
  // Conj
  //

  template< class T_Scalar, Degree T_Degree >
  class Conj final : public UnaryOperation< T_Scalar, T_Degree, T_Degree >
  {
   public:
    Conj()
    {
    }

    Conj(const Conj &other) :
      UnaryOperation< T_Scalar, T_Degree, T_Degree >(other)
    {
    }

    template< class T_A >
    void operator()(OutputVector< T_Scalar, T_Degree > &values, const InputVector< T_Scalar, T_Degree, T_A > &a, const std::vector< scalar::Precision< T_Scalar >, numa::allocator< scalar::Precision< T_Scalar > > > &points, const std::vector< scalar::Precision< T_Scalar > > &gaussPoints, const int elementType, const std::pair< int, int > &entity) const
    {
#pragma omp for
      for(auto i = 0ULL; i < values.size(); ++i) {
        values[i] = a[i];
      }
    }

    bool canUseSameVectorsForOutputAndInputs() const override
    {
      return true;
    }

    std::string name() const override
    {
      return "conj(a)";
    }
  };

  template< class T_ComplexScalar >
  class Conj< std::complex< T_ComplexScalar >, Degree::Degree0 > final : public UnaryOperation< std::complex< T_ComplexScalar >, Degree::Degree0, Degree::Degree0 >
  {
   public:
    Conj()
    {
    }

    Conj(const Conj &other) :
      UnaryOperation< std::complex< T_ComplexScalar >, Degree::Degree0, Degree::Degree0 >(other)
    {
    }

    template< class T_A >
    void operator()(OutputVector< std::complex< T_ComplexScalar >, Degree::Degree0 > &values, const InputVector< std::complex< T_ComplexScalar >, Degree::Degree0, T_A > &a, const std::vector< scalar::Precision< std::complex< T_ComplexScalar > >, numa::allocator< scalar::Precision< std::complex< T_ComplexScalar > > > > &points, const std::vector< scalar::Precision< std::complex< T_ComplexScalar > > > &gaussPoints, const int elementType, const std::pair< int, int > &entity) const
    {
#pragma omp for
      for(auto i = 0ULL; i < values.size(); ++i) {
        values[i] = std::conj(a[i]);
      }
    }

    bool canUseSameVectorsForOutputAndInputs() const override
    {
      return true;
    }

    std::string name() const override
    {
      return "conj(a)";
    }
  };

  template< class T_ComplexScalar >
  class Conj< std::complex< T_ComplexScalar >, Degree::Degree1 > final : public UnaryOperation< std::complex< T_ComplexScalar >, Degree::Degree1, Degree::Degree1 >
  {
   public:
    Conj()
    {
    }

    Conj(const Conj &other) :
      UnaryOperation< std::complex< T_ComplexScalar >, Degree::Degree1, Degree::Degree1 >(other)
    {
    }

    template< class T_A >
    void operator()(OutputVector< std::complex< T_ComplexScalar >, Degree::Degree1 > &values, const InputVector< std::complex< T_ComplexScalar >, Degree::Degree1, T_A > &a, const std::vector< scalar::Precision< std::complex< T_ComplexScalar > >, numa::allocator< scalar::Precision< std::complex< T_ComplexScalar > > > > &points, const std::vector< scalar::Precision< std::complex< T_ComplexScalar > > > &gaussPoints, const int elementType, const std::pair< int, int > &entity) const
    {
#pragma omp for
      for(auto i = 0ULL; i < values.size(); ++i) {
        values[i] = (a[i].conjugate()).eval();
      }
    }

    bool canUseSameVectorsForOutputAndInputs() const override
    {
      return true;
    }

    std::string name() const override
    {
      return "conj(a)";
    }
  };

  template< class T_ComplexScalar >
  class Conj< std::complex< T_ComplexScalar >, Degree::Degree2 > final : public UnaryOperation< std::complex< T_ComplexScalar >, Degree::Degree2, Degree::Degree2 >
  {
   public:
    Conj()
    {
    }

    Conj(const Conj &other) :
      UnaryOperation< std::complex< T_ComplexScalar >, Degree::Degree2, Degree::Degree2 >(other)
    {
    }

    template< class T_A >
    void operator()(OutputVector< std::complex< T_ComplexScalar >, Degree::Degree2 > &values, const InputVector< std::complex< T_ComplexScalar >, Degree::Degree2, T_A > &a, const std::vector< scalar::Precision< std::complex< T_ComplexScalar > >, numa::allocator< scalar::Precision< std::complex< T_ComplexScalar > > > > &points, const std::vector< scalar::Precision< std::complex< T_ComplexScalar > > > &gaussPoints, const int elementType, const std::pair< int, int > &entity) const
    {
#pragma omp for
      for(auto i = 0ULL; i < values.size(); ++i) {
        values[i] = (a[i].conjugate()).eval();
      }
    }

    bool canUseSameVectorsForOutputAndInputs() const override
    {
      return true;
    }

    std::string name() const override
    {
      return "conj(a)";
    }
  };

  //
  // Cos
  //

  template< class T_Scalar, Degree T_Degree >
  class Cos final : public UnaryOperation< T_Scalar, T_Degree, T_Degree >
  {
   public:
    Cos()
    {
    }

    Cos(const Cos &other) :
      UnaryOperation< T_Scalar, T_Degree, T_Degree >(other)
    {
    }

    template< class T_A >
    void operator()(OutputVector< T_Scalar, T_Degree > &values, const InputVector< T_Scalar, T_Degree, T_A > &a, const std::vector< scalar::Precision< T_Scalar >, numa::allocator< scalar::Precision< T_Scalar > > > &points, const std::vector< scalar::Precision< T_Scalar > > &gaussPoints, const int elementType, const std::pair< int, int > &entity) const
    {
#pragma omp for
      for(auto i = 0ULL; i < values.size(); ++i) {
        values[i] = std::cos(a[i]);
      }
    }

    bool canUseSameVectorsForOutputAndInputs() const override
    {
      return true;
    }

    std::string name() const override
    {
      return "cos(a)";
    }
  };

  //
  // Cosh
  //

  template< class T_Scalar, Degree T_Degree >
  class Cosh final : public UnaryOperation< T_Scalar, T_Degree, T_Degree >
  {
   public:
    Cosh()
    {
    }

    Cosh(const Cosh &other) :
      UnaryOperation< T_Scalar, T_Degree, T_Degree >(other)
    {
    }

    template< class T_A >
    void operator()(OutputVector< T_Scalar, T_Degree > &values, const InputVector< T_Scalar, T_Degree, T_A > &a, const std::vector< scalar::Precision< T_Scalar >, numa::allocator< scalar::Precision< T_Scalar > > > &points, const std::vector< scalar::Precision< T_Scalar > > &gaussPoints, const int elementType, const std::pair< int, int > &entity) const
    {
#pragma omp for
      for(auto i = 0ULL; i < values.size(); ++i) {
        values[i] = std::cosh(a[i]);
      }
    }

    bool canUseSameVectorsForOutputAndInputs() const override
    {
      return true;
    }

    std::string name() const override
    {
      return "cosh(a)";
    }
  };

  //
  // CylBesselJ
  //

  template< class T_Scalar >
  class CylBesselJ final : public UnaryOperation< T_Scalar, Degree::Degree0, Degree::Degree0 >
  {
   private:
    const double _v;

   public:
    CylBesselJ(const double v) :
      _v(v)
    {
    }

    CylBesselJ(const CylBesselJ &other) :
      UnaryOperation< T_Scalar, Degree::Degree0, Degree::Degree0 >(other), _v(other._v)
    {
    }

    template< class T_A >
    void operator()(OutputVector< T_Scalar, Degree::Degree0 > &values, const InputVector< T_Scalar, Degree::Degree0, T_A > &a, const std::vector< scalar::Precision< T_Scalar >, numa::allocator< scalar::Precision< T_Scalar > > > &points, const std::vector< scalar::Precision< T_Scalar > > &gaussPoints, const int elementType, const std::pair< int, int > &entity) const
    {
#ifdef HAVE_BESSEL
      double valr, vali;
#pragma omp for
      for(auto i = 0ULL; i < values.size(); ++i) {
        int err = BesselJnComplex(_v, 1, double(std::real(a[i])), double(std::imag(a[i])), &valr, &vali);
        if(err != 0) {
#pragma omp critical
          msg::error << "Issue with complex bessel function, error output: " << err << msg::endl;
        }
        values[i] = T_Scalar(valr);
      }
#else
      throw common::Exception("Cannot use cylindrical Bessel functions j (of the first kind) without Bessel functions support (ENABLE_BESSEL)");
#endif // HAVE_BESSEL
    }

    bool canUseSameVectorsForOutputAndInputs() const override
    {
      return true;
    }

    std::string name() const override
    {
      return "cylBesselJ(" + std::to_string(_v) + ", a)";
    }

    bool operator==(const UnaryOperation< T_Scalar, Degree::Degree0, Degree::Degree0 > &other) const override
    {
      if(_v == static_cast< const CylBesselJ & >(other)._v) {
        return true;
      }
      return false;
    }
  };

  template< class T_ComplexScalar >
  class CylBesselJ< std::complex< T_ComplexScalar > > final : public UnaryOperation< std::complex< T_ComplexScalar >, Degree::Degree0, Degree::Degree0 >
  {
   private:
    const double _v;

   public:
    CylBesselJ(const double v) :
      _v(v)
    {
    }

    CylBesselJ(const CylBesselJ &other) :
      UnaryOperation< std::complex< T_ComplexScalar >, Degree::Degree0, Degree::Degree0 >(other), _v(other._v)
    {
    }

    template< class T_A >
    void operator()(OutputVector< std::complex< T_ComplexScalar >, Degree::Degree0 > &values, const InputVector< std::complex< T_ComplexScalar >, Degree::Degree0, T_A > &a, const std::vector< scalar::Precision< std::complex< T_ComplexScalar > >, numa::allocator< scalar::Precision< std::complex< T_ComplexScalar > > > > &points, const std::vector< scalar::Precision< std::complex< T_ComplexScalar > > > &gaussPoints, const int elementType, const std::pair< int, int > &entity) const
    {
#ifdef HAVE_BESSEL
      double valr, vali;
#pragma omp for
      for(auto i = 0ULL; i < values.size(); ++i) {
        int err = BesselJnComplex(_v, 1, double(std::real(a[i])), double(std::imag(a[i])), &valr, &vali);
        if(err != 0) {
#pragma omp critical
          msg::error << "Issue with complex bessel function, error output: " << err << msg::endl;
        }
        values[i] = std::complex< T_ComplexScalar >(T_ComplexScalar(valr), T_ComplexScalar(vali));
      }
#else
      throw common::Exception("Cannot use cylindrical Bessel functions j (of the first kind) without Bessel functions support (ENABLE_BESSEL)");
#endif // HAVE_BESSEL
    }

    bool canUseSameVectorsForOutputAndInputs() const override
    {
      return true;
    }

    std::string name() const override
    {
      return "cylBesselJ(" + std::to_string(_v) + ",a)";
    }

    bool operator==(const UnaryOperation< std::complex< T_ComplexScalar >, Degree::Degree0, Degree::Degree0 > &other) const override
    {
      if(_v == static_cast< const CylBesselJ & >(other)._v) {
        return true;
      }
      return false;
    }
  };

  //
  // CylNeumann
  //

  template< class T_Scalar >
  class CylNeumann final : public UnaryOperation< T_Scalar, Degree::Degree0, Degree::Degree0 >
  {
   private:
    const double _v;

   public:
    CylNeumann(const double v) :
      _v(v)
    {
    }

    CylNeumann(const CylNeumann &other) :
      UnaryOperation< T_Scalar, Degree::Degree0, Degree::Degree0 >(other), _v(other._v)
    {
    }

    template< class T_A >
    void operator()(OutputVector< T_Scalar, Degree::Degree0 > &values, const InputVector< T_Scalar, Degree::Degree0, T_A > &a, const std::vector< scalar::Precision< T_Scalar >, numa::allocator< scalar::Precision< T_Scalar > > > &points, const std::vector< scalar::Precision< T_Scalar > > &gaussPoints, const int elementType, const std::pair< int, int > &entity) const
    {
#ifdef HAVE_BESSEL
      double valr, vali;
#pragma omp for
      for(auto i = 0ULL; i < values.size(); ++i) {
        int err = BesselYnComplex(_v, 1, double(std::real(a[i])), double(std::imag(a[i])), &valr, &vali);
        if(err != 0) {
#pragma omp critical
          msg::error << "Issue with complex bessel function, error output: " << err << msg::endl;
        }
        values[i] = T_Scalar(valr);
      }
#else
      throw common::Exception("Cannot use cylindrical Neumann functions without Bessel functions support (ENABLE_BESSEL)");
#endif // HAVE_BESSEL
    }

    bool canUseSameVectorsForOutputAndInputs() const override
    {
      return true;
    }

    std::string name() const override
    {
      return "cylNeumann(" + std::to_string(_v) + ",a)";
    }

    bool operator==(const UnaryOperation< T_Scalar, Degree::Degree0, Degree::Degree0 > &other) const override
    {
      if(_v == static_cast< const CylNeumann & >(other)._v) {
        return true;
      }
      return false;
    }
  };

  template< class T_ComplexScalar >
  class CylNeumann< std::complex< T_ComplexScalar > > final : public UnaryOperation< std::complex< T_ComplexScalar >, Degree::Degree0, Degree::Degree0 >
  {
   private:
    mutable double _v;

   public:
    CylNeumann(const double v) :
      _v(v)
    {
    }

    CylNeumann(const CylNeumann &other) :
      UnaryOperation< std::complex< T_ComplexScalar >, Degree::Degree0, Degree::Degree0 >(other), _v(other._v)
    {
    }

    template< class T_A >
    void operator()(OutputVector< std::complex< T_ComplexScalar >, Degree::Degree0 > &values, const InputVector< std::complex< T_ComplexScalar >, Degree::Degree0, T_A > &a, const std::vector< scalar::Precision< std::complex< T_ComplexScalar > >, numa::allocator< scalar::Precision< std::complex< T_ComplexScalar > > > > &points, const std::vector< scalar::Precision< std::complex< T_ComplexScalar > > > &gaussPoints, const int elementType, const std::pair< int, int > &entity) const
    {
#ifdef HAVE_BESSEL
      double valr, vali;
#pragma omp for
      for(auto i = 0ULL; i < values.size(); ++i) {
        int err = BesselYnComplex(_v, 1, double(std::real(a[i])), double(std::imag(a[i])), &valr, &vali);
        if(err != 0) {
#pragma omp critical
          msg::error << "Issue with complex bessel function, error output:" << err << msg::endl;
        }
        values[i] = std::complex< T_ComplexScalar >(valr, vali);
      }
#else
      throw common::Exception("Cannot use cylindrical Neumann functions without Bessel functions support (ENABLE_BESSEL)");
#endif // HAVE_BESSEL
    }

    bool canUseSameVectorsForOutputAndInputs() const override
    {
      return true;
    }

    std::string name() const override
    {
      return "cylNeumann(" + std::to_string(_v) + ",a)";
    }

    bool operator==(const UnaryOperation< std::complex< T_ComplexScalar >, Degree::Degree0, Degree::Degree0 > &other) const override
    {
      if(_v == static_cast< const CylNeumann & >(other)._v) {
        return true;
      }
      return false;
    }
  };

  //
  // Exp
  //

  template< class T_Scalar, Degree T_Degree >
  class Exp final : public UnaryOperation< T_Scalar, T_Degree, T_Degree >
  {
   public:
    Exp()
    {
    }

    Exp(const Exp &other) :
      UnaryOperation< T_Scalar, T_Degree, T_Degree >(other)
    {
    }

    template< class T_A >
    void operator()(OutputVector< T_Scalar, T_Degree > &values, const InputVector< T_Scalar, T_Degree, T_A > &a, const std::vector< scalar::Precision< T_Scalar >, numa::allocator< scalar::Precision< T_Scalar > > > &points, const std::vector< scalar::Precision< T_Scalar > > &gaussPoints, const int elementType, const std::pair< int, int > &entity) const
    {
#pragma omp for
      for(auto i = 0ULL; i < values.size(); ++i) {
        values[i] = std::exp(a[i]);
      }
    }

    bool canUseSameVectorsForOutputAndInputs() const override
    {
      return true;
    }

    std::string name() const override
    {
      return "exp(a)";
    }
  };

  //
  // Heaviside
  //

  template< class T_Scalar, Degree T_Degree >
  class Heaviside final : public UnaryOperation< T_Scalar, T_Degree, T_Degree >
  {
   public:
    Heaviside()
    {
    }

    Heaviside(const Heaviside &other) :
      UnaryOperation< T_Scalar, T_Degree, T_Degree >(other)
    {
    }

    template< class T_A >
    void operator()(OutputVector< T_Scalar, T_Degree > &values, const InputVector< T_Scalar, T_Degree, T_A > &a, const std::vector< scalar::Precision< T_Scalar >, numa::allocator< scalar::Precision< T_Scalar > > > &points, const std::vector< scalar::Precision< T_Scalar > > &gaussPoints, const int elementType, const std::pair< int, int > &entity) const
    {
#pragma omp for
      for(auto i = 0ULL; i < values.size(); ++i) {
        if(std::real(a[i]) > 0.) {
          values[i] = 1.;
        }
        else if(std::real(a[i]) < 0.) {
          values[i] = 0.;
        }
        else {
          values[i] = 0.5;
        }
      }
    }

    bool canUseSameVectorsForOutputAndInputs() const override
    {
      return true;
    }

    std::string name() const override
    {
      return "heaviside(a)";
    }
  };

  //
  // Imag
  //

  template< class T_Scalar, Degree T_Degree >
  class Imag final : public UnaryOperation< T_Scalar, T_Degree, T_Degree >
  {
   public:
    Imag()
    {
    }

    Imag(const Imag &other) :
      UnaryOperation< T_Scalar, T_Degree, T_Degree >(other)
    {
    }

    template< class T_A >
    void operator()(OutputVector< T_Scalar, T_Degree > &values, const InputVector< T_Scalar, T_Degree, T_A > &a, const std::vector< scalar::Precision< T_Scalar >, numa::allocator< scalar::Precision< T_Scalar > > > &points, const std::vector< scalar::Precision< T_Scalar > > &gaussPoints, const int elementType, const std::pair< int, int > &entity) const
    {
#pragma omp for
      for(auto i = 0ULL; i < values.size(); ++i) {
        values[i] = 0.;
      }
    }

    bool canUseSameVectorsForOutputAndInputs() const override
    {
      return true;
    }

    std::string name() const override
    {
      return "imag(a)";
    }
  };

  template< class T_ComplexScalar, Degree T_Degree >
  class Imag< std::complex< T_ComplexScalar >, T_Degree > final : public UnaryOperation< std::complex< T_ComplexScalar >, T_Degree, T_Degree >
  {
   public:
    Imag()
    {
    }

    Imag(const Imag &other) :
      UnaryOperation< std::complex< T_ComplexScalar >, T_Degree, T_Degree >(other)
    {
    }

    template< class T_A >
    void operator()(OutputVector< std::complex< T_ComplexScalar >, T_Degree > &values, const InputVector< std::complex< T_ComplexScalar >, T_Degree, T_A > &a, const std::vector< scalar::Precision< std::complex< T_ComplexScalar > >, numa::allocator< scalar::Precision< std::complex< T_ComplexScalar > > > > &points, const std::vector< scalar::Precision< std::complex< T_ComplexScalar > > > &gaussPoints, const int elementType, const std::pair< int, int > &entity) const
    {
#pragma omp for
      for(auto i = 0ULL; i < values.size(); ++i) {
        values[i] = std::imag(a[i]);
      }
    }

    bool canUseSameVectorsForOutputAndInputs() const override
    {
      return true;
    }

    std::string name() const override
    {
      return "imag(a)";
    }
  };

  //
  // Inv
  //

  template< class T_Scalar, Degree T_Degree >
  class Inv final : public UnaryOperation< T_Scalar, T_Degree, T_Degree >
  {
   public:
    Inv()
    {
    }

    Inv(const Inv &other) :
      UnaryOperation< T_Scalar, T_Degree, T_Degree >(other)
    {
    }

    template< class T_A >
    void operator()(OutputVector< T_Scalar, T_Degree > &values, const InputVector< T_Scalar, T_Degree, T_A > &a, const std::vector< scalar::Precision< T_Scalar >, numa::allocator< scalar::Precision< T_Scalar > > > &points, const std::vector< scalar::Precision< T_Scalar > > &gaussPoints, const int elementType, const std::pair< int, int > &entity) const
    {
#pragma omp for
      for(auto i = 0ULL; i < values.size(); ++i) {
        values[i] = T_Scalar(1.) / a[i];
      }
    }

    bool canUseSameVectorsForOutputAndInputs() const override
    {
      return true;
    }

    std::string name() const override
    {
      return "inv(a)";
    }
  };

  template< class T_Scalar >
  class Inv< T_Scalar, Degree::Degree2 > final : public UnaryOperation< T_Scalar, Degree::Degree2, Degree::Degree2 >
  {
   public:
    Inv()
    {
    }

    Inv(const Inv &other) :
      UnaryOperation< T_Scalar, Degree::Degree2, Degree::Degree2 >(other)
    {
    }

    template< class T_A >
    void operator()(OutputVector< T_Scalar, Degree::Degree2 > &values, const InputVector< T_Scalar, Degree::Degree2, T_A > &a, const std::vector< scalar::Precision< T_Scalar >, numa::allocator< scalar::Precision< T_Scalar > > > &points, const std::vector< scalar::Precision< T_Scalar > > &gaussPoints, const int elementType, const std::pair< int, int > &entity) const
    {
#pragma omp for
      for(auto i = 0ULL; i < values.size(); ++i) {
        values[i] = (a[i].inverse()).eval();
      }
    }

    bool canUseSameVectorsForOutputAndInputs() const override
    {
      return true;
    }

    std::string name() const override
    {
      return "inv(a)";
    }
  };

  //
  // Log
  //

  template< class T_Scalar, Degree T_Degree >
  class Log final : public UnaryOperation< T_Scalar, T_Degree, T_Degree >
  {
   public:
    Log()
    {
    }

    Log(const Log &other) :
      UnaryOperation< T_Scalar, T_Degree, T_Degree >(other)
    {
    }

    template< class T_A >
    void operator()(OutputVector< T_Scalar, T_Degree > &values, const InputVector< T_Scalar, T_Degree, T_A > &a, const std::vector< scalar::Precision< T_Scalar >, numa::allocator< scalar::Precision< T_Scalar > > > &points, const std::vector< scalar::Precision< T_Scalar > > &gaussPoints, const int elementType, const std::pair< int, int > &entity) const
    {
#pragma omp for
      for(auto i = 0ULL; i < values.size(); ++i) {
        values[i] = std::log10(a[i]);
      }
    }

    bool canUseSameVectorsForOutputAndInputs() const override
    {
      return true;
    }

    std::string name() const override
    {
      return "log(a)";
    }
  };

  //
  // Ln
  //

  template< class T_Scalar, Degree T_Degree >
  class Ln final : public UnaryOperation< T_Scalar, T_Degree, T_Degree >
  {
   public:
    Ln()
    {
    }

    Ln(const Ln &other) :
      UnaryOperation< T_Scalar, T_Degree, T_Degree >(other)
    {
    }

    template< class T_A >
    void operator()(OutputVector< T_Scalar, T_Degree > &values, const InputVector< T_Scalar, T_Degree, T_A > &a, const std::vector< scalar::Precision< T_Scalar >, numa::allocator< scalar::Precision< T_Scalar > > > &points, const std::vector< scalar::Precision< T_Scalar > > &gaussPoints, const int elementType, const std::pair< int, int > &entity) const
    {
#pragma omp for
      for(auto i = 0ULL; i < values.size(); ++i) {
        values[i] = std::log(a[i]);
      }
    }

    bool canUseSameVectorsForOutputAndInputs() const override
    {
      return true;
    }

    std::string name() const override
    {
      return "ln(a)";
    }
  };

  //
  // Minus
  //

  template< class T_Scalar, Degree T_Degree >
  class Minus final : public UnaryOperation< T_Scalar, T_Degree, T_Degree >
  {
   public:
    Minus()
    {
    }

    Minus(const Minus &other) :
      UnaryOperation< T_Scalar, T_Degree, T_Degree >(other)
    {
    }

    template< class T_A >
    void operator()(OutputVector< T_Scalar, T_Degree > &values, const InputVector< T_Scalar, T_Degree, T_A > &a, const std::vector< scalar::Precision< T_Scalar >, numa::allocator< scalar::Precision< T_Scalar > > > &points, const std::vector< scalar::Precision< T_Scalar > > &gaussPoints, const int elementType, const std::pair< int, int > &entity) const
    {
#pragma omp for
      for(auto i = 0ULL; i < values.size(); ++i) {
        values[i] = -a[i];
      }
    }

    bool canUseSameVectorsForOutputAndInputs() const override
    {
      return true;
    }

    std::string name() const override
    {
      return "- a";
    }
  };

  //
  // Norm
  //

  template< class T_Scalar, Degree T_DegreeA >
  class Norm final : public UnaryOperation< T_Scalar, Degree::Degree0, T_DegreeA >
  {
   private:
    mutable int _order;

   public:
    Norm(const int order) :
      _order(order)
    {
    }

    Norm(const Norm &other) :
      UnaryOperation< T_Scalar, Degree::Degree0, T_DegreeA >(other), _order(other._order)
    {
    }

    template< class T_A >
    void operator()(OutputVector< T_Scalar, Degree::Degree0 > &values, const InputVector< T_Scalar, T_DegreeA, T_A > &a, const std::vector< scalar::Precision< T_Scalar >, numa::allocator< scalar::Precision< T_Scalar > > > &points, const std::vector< scalar::Precision< T_Scalar > > &gaussPoints, const int elementType, const std::pair< int, int > &entity) const
    {
#pragma omp for
      for(auto i = 0ULL; i < values.size(); ++i) {
        values[i] = std::abs(a[i]);
      }
    }

    bool canUseSameVectorsForOutputAndInputs() const override
    {
      return true;
    }

    std::string name() const override
    {
      return "norm(a, " + std::to_string(_order) + ")";
    }

    bool operator==(const UnaryOperation< T_Scalar, Degree::Degree0, T_DegreeA > &other) const override
    {
      if(_order == static_cast< const Norm & >(other)._order) {
        return true;
      }
      return false;
    }
  };

  template< class T_Scalar >
  class Norm< T_Scalar, Degree::Degree1 > final : public UnaryOperation< T_Scalar, Degree::Degree0, Degree::Degree1 >
  {
   private:
    mutable int _order;

   public:
    Norm(const int order) :
      _order(order)
    {
    }

    Norm(const Norm &other) :
      UnaryOperation< T_Scalar, Degree::Degree0, Degree::Degree1 >(other), _order(other._order)
    {
    }

    template< class T_A >
    void operator()(OutputVector< T_Scalar, Degree::Degree0 > &values, const InputVector< T_Scalar, Degree::Degree1, T_A > &a, const std::vector< scalar::Precision< T_Scalar >, numa::allocator< scalar::Precision< T_Scalar > > > &points, const std::vector< scalar::Precision< T_Scalar > > &gaussPoints, const int elementType, const std::pair< int, int > &entity) const
    {
      if(_order == -1) { // Infinit norm
#pragma omp for
        for(auto i = 0ULL; i < values.size(); ++i) {
          values[i] = a[i].template lpNorm< Eigen::Infinity >();
        }
      }
      else if(_order == 1) {
#pragma omp for
        for(auto i = 0ULL; i < values.size(); ++i) {
          values[i] = a[i].template lpNorm< 1 >();
        }
      }
      else if(_order == 2) {
#pragma omp for
        for(auto i = 0ULL; i < values.size(); ++i) {
          values[i] = a[i].norm();
        }
      }
      else {
        msg::warning << "Cannot apply the vectorial " + std::to_string(_order) + "-norm. The 2-norm is used" << msg::endl;
#pragma omp for
        for(auto i = 0ULL; i < values.size(); ++i) {
          values[i] = a[i].norm();
        }
      }
    }

    std::string name() const override
    {
      return "norm(a, " + std::to_string(_order) + ")";
    }

    bool operator==(const UnaryOperation< T_Scalar, Degree::Degree0, Degree::Degree1 > &other) const override
    {
      if(_order == static_cast< const Norm & >(other)._order) {
        return true;
      }
      return false;
    }
  };

  template< class T_Scalar >
  class Norm< T_Scalar, Degree::Degree2 > final : public UnaryOperation< T_Scalar, Degree::Degree0, Degree::Degree2 >
  {
   private:
    mutable int _order;

   public:
    Norm(const int order) :
      _order(order)
    {
    }

    Norm(const Norm &other) :
      UnaryOperation< T_Scalar, Degree::Degree0, Degree::Degree2 >(other), _order(other._order)
    {
    }

    template< class T_A >
    void operator()(OutputVector< T_Scalar, Degree::Degree0 > &values, const InputVector< T_Scalar, Degree::Degree2, T_A > &a, const std::vector< scalar::Precision< T_Scalar >, numa::allocator< scalar::Precision< T_Scalar > > > &points, const std::vector< scalar::Precision< T_Scalar > > &gaussPoints, const int elementType, const std::pair< int, int > &entity) const
    {
      if(_order == -1) { // Infinit norm
#pragma omp for
        for(auto i = 0ULL; i < values.size(); ++i) {
          values[i] = a[i].cwiseAbs().rowwise().sum().maxCoeff();
        }
      }
      else if(_order == 1) {
#pragma omp for
        for(auto i = 0ULL; i < values.size(); ++i) {
          values[i] = a[i].cwiseAbs().colwise().sum().maxCoeff();
        }
      }
      else if(_order == 2) { // Frobenius norm
#pragma omp for
        for(auto i = 0ULL; i < values.size(); ++i) {
          values[i] = a[i].norm();
        }
      }
      else {
        msg::warning << "Cannot apply the tensorial " + std::to_string(_order) + "-norm. The 2-norm is used" << msg::endl;
#pragma omp for
        for(auto i = 0ULL; i < values.size(); ++i) {
          values[i] = a[i].norm();
        }
      }
    }

    std::string name() const override
    {
      return "norm(a, " + std::to_string(_order) + ")";
    }

    bool operator==(const UnaryOperation< T_Scalar, Degree::Degree0, Degree::Degree2 > &other) const override
    {
      if(_order == static_cast< const Norm & >(other)._order) {
        return true;
      }
      return false;
    }
  };

  //
  // Plus
  //

  template< class T_Scalar, Degree T_Degree >
  class Plus final : public UnaryOperation< T_Scalar, T_Degree, T_Degree >
  {
   public:
    Plus()
    {
    }

    Plus(const Plus &other) :
      UnaryOperation< T_Scalar, T_Degree, T_Degree >(other)
    {
    }

    template< class T_A >
    void operator()(OutputVector< T_Scalar, T_Degree > &values, const InputVector< T_Scalar, T_Degree, T_A > &a, const std::vector< scalar::Precision< T_Scalar >, numa::allocator< scalar::Precision< T_Scalar > > > &points, const std::vector< scalar::Precision< T_Scalar > > &gaussPoints, const int elementType, const std::pair< int, int > &entity) const
    {
#pragma omp for
      for(auto i = 0ULL; i < values.size(); ++i) {
        values[i] = a[i];
      }
    }

    bool canUseSameVectorsForOutputAndInputs() const override
    {
      return true;
    }

    std::string name() const override
    {
      return "+ a";
    }
  };

  //
  // Pow
  //

  template< class T_Scalar, Degree T_Degree >
  class Pow final : public UnaryOperation< T_Scalar, T_Degree, T_Degree >
  {
   private:
    mutable unsigned int _pow;

   public:
    Pow(const unsigned int pow) :
      _pow(pow)
    {
    }

    Pow(const Pow &other) :
      UnaryOperation< T_Scalar, T_Degree, T_Degree >(other), _pow(other._pow)
    {
    }

    template< class T_A >
    void operator()(OutputVector< T_Scalar, T_Degree > &values, const InputVector< T_Scalar, T_Degree, T_A > &a, const std::vector< scalar::Precision< T_Scalar >, numa::allocator< scalar::Precision< T_Scalar > > > &points, const std::vector< scalar::Precision< T_Scalar > > &gaussPoints, const int elementType, const std::pair< int, int > &entity) const
    {
      if(_pow == 1) {
#pragma omp for
        for(auto i = 0ULL; i < values.size(); ++i) {
          values[i] = a[i];
        }
      }
      else if(_pow == 2) {
#pragma omp for
        for(auto i = 0ULL; i < values.size(); ++i) {
          values[i] = a[i] * a[i];
        }
      }
      else if(_pow == 3) {
#pragma omp for
        for(auto i = 0ULL; i < values.size(); ++i) {
          values[i] = a[i] * a[i] * a[i];
        }
      }
      else if(_pow == 4) {
#pragma omp for
        for(auto i = 0ULL; i < values.size(); ++i) {
          values[i] = a[i] * a[i];
          values[i] *= values[i];
        }
      }
      else {
#pragma omp for
        for(auto i = 0ULL; i < values.size(); ++i) {
          values[i] = std::pow(a[i], _pow);
        }
      }
    }

    bool canUseSameVectorsForOutputAndInputs() const override
    {
      return true;
    }

    std::string name() const override
    {
      return "pow(a, " + std::to_string(_pow) + ")";
    }

    bool operator==(const UnaryOperation< T_Scalar, T_Degree, T_Degree > &other) const override
    {
      if(_pow == static_cast< const Pow & >(other)._pow) {
        return true;
      }
      return false;
    }
  };

  //
  // R2Dcomp
  //

  template< class T_Scalar >
  class R2Dcomp final : public UnaryOperation< T_Scalar, Degree::Degree0, Degree::Degree1 >
  {
   public:
    R2Dcomp()
    {
    }

    R2Dcomp(const R2Dcomp &other) :
      UnaryOperation< T_Scalar, Degree::Degree0, Degree::Degree1 >(other)
    {
    }

    template< class T_A >
    void operator()(OutputVector< T_Scalar, Degree::Degree0 > &values, const InputVector< T_Scalar, Degree::Degree1, T_A > &a, const std::vector< scalar::Precision< T_Scalar >, numa::allocator< scalar::Precision< T_Scalar > > > &points, const std::vector< scalar::Precision< T_Scalar > > &gaussPoints, const int elementType, const std::pair< int, int > &entity) const
    {
#pragma omp for
      for(auto i = 0ULL; i < values.size(); ++i) {
        values[i] = std::sqrt(a[i](0) * a[i](0) + a[i](1) * a[i](1));
      }
    }

    std::string name() const override
    {
      return "r2D component(a)";
    }
  };

  template< class T_ComplexScalar >
  class R2Dcomp< std::complex< T_ComplexScalar > > final : public UnaryOperation< std::complex< T_ComplexScalar >, Degree::Degree0, Degree::Degree1 >
  {
   public:
    R2Dcomp()
    {
    }

    R2Dcomp(const R2Dcomp &other) :
      UnaryOperation< std::complex< T_ComplexScalar >, Degree::Degree0, Degree::Degree1 >(other)
    {
    }

    template< class T_A >
    void operator()(OutputVector< std::complex< T_ComplexScalar >, Degree::Degree0 > &values, const InputVector< std::complex< T_ComplexScalar >, Degree::Degree1, T_A > &a, const std::vector< scalar::Precision< std::complex< T_ComplexScalar > >, numa::allocator< scalar::Precision< std::complex< T_ComplexScalar > > > > &points, const std::vector< scalar::Precision< std::complex< T_ComplexScalar > > > &gaussPoints, const int elementType, const std::pair< int, int > &entity) const
    {
#pragma omp for
      for(auto i = 0ULL; i < values.size(); ++i) {
        values[i] = std::complex< T_ComplexScalar >(std::sqrt(std::real(a[i](0)) * std::real(a[i](0)) + std::real(a[i](1)) * std::real(a[i](1))), std::sqrt(std::imag(a[i](0)) * std::imag(a[i](0)) + std::imag(a[i](1)) * std::imag(a[i](1))));
      }
    }

    std::string name() const override
    {
      return "r2D component(a)";
    }
  };

  //
  // Real
  //

  template< class T_Scalar, Degree T_Degree >
  class Real final : public UnaryOperation< T_Scalar, T_Degree, T_Degree >
  {
   public:
    Real()
    {
    }

    Real(const Real &other) :
      UnaryOperation< T_Scalar, T_Degree, T_Degree >(other)
    {
    }

    template< class T_A >
    void operator()(OutputVector< T_Scalar, T_Degree > &values, const InputVector< T_Scalar, T_Degree, T_A > &a, const std::vector< scalar::Precision< T_Scalar >, numa::allocator< scalar::Precision< T_Scalar > > > &points, const std::vector< scalar::Precision< T_Scalar > > &gaussPoints, const int elementType, const std::pair< int, int > &entity) const
    {
#pragma omp for
      for(auto i = 0ULL; i < values.size(); ++i) {
        values[i] = a[i];
      }
    }

    bool canUseSameVectorsForOutputAndInputs() const override
    {
      return true;
    }

    std::string name() const override
    {
      return "real(a)";
    }
  };

  template< class T_ComplexScalar, Degree T_Degree >
  class Real< std::complex< T_ComplexScalar >, T_Degree > final : public UnaryOperation< std::complex< T_ComplexScalar >, T_Degree, T_Degree >
  {
   public:
    Real()
    {
    }

    Real(const Real &other) :
      UnaryOperation< std::complex< T_ComplexScalar >, T_Degree, T_Degree >(other)
    {
    }

    template< class T_A >
    void operator()(OutputVector< std::complex< T_ComplexScalar >, T_Degree > &values, const InputVector< std::complex< T_ComplexScalar >, T_Degree, T_A > &a, const std::vector< scalar::Precision< std::complex< T_ComplexScalar > >, numa::allocator< scalar::Precision< std::complex< T_ComplexScalar > > > > &points, const std::vector< scalar::Precision< std::complex< T_ComplexScalar > > > &gaussPoints, const int elementType, const std::pair< int, int > &entity) const
    {
#pragma omp for
      for(auto i = 0ULL; i < values.size(); ++i) {
        values[i] = std::real(a[i]);
      }
    }

    bool canUseSameVectorsForOutputAndInputs() const override
    {
      return true;
    }

    std::string name() const override
    {
      return "real(a)";
    }
  };

  //
  // Sin
  //

  template< class T_Scalar, Degree T_Degree >
  class Sin final : public UnaryOperation< T_Scalar, T_Degree, T_Degree >
  {
   public:
    Sin()
    {
    }

    Sin(const Sin &other) :
      UnaryOperation< T_Scalar, T_Degree, T_Degree >(other)
    {
    }

    template< class T_A >
    void operator()(OutputVector< T_Scalar, T_Degree > &values, const InputVector< T_Scalar, T_Degree, T_A > &a, const std::vector< scalar::Precision< T_Scalar >, numa::allocator< scalar::Precision< T_Scalar > > > &points, const std::vector< scalar::Precision< T_Scalar > > &gaussPoints, const int elementType, const std::pair< int, int > &entity) const
    {
#pragma omp for
      for(auto i = 0ULL; i < values.size(); ++i) {
        values[i] = std::sin(a[i]);
      }
    }

    bool canUseSameVectorsForOutputAndInputs() const override
    {
      return true;
    }

    std::string name() const override
    {
      return "sin(a)";
    }
  };

  //
  // Sinh
  //

  template< class T_Scalar, Degree T_Degree >
  class Sinh final : public UnaryOperation< T_Scalar, T_Degree, T_Degree >
  {
   public:
    Sinh()
    {
    }

    Sinh(const Sinh &other) :
      UnaryOperation< T_Scalar, T_Degree, T_Degree >(other)
    {
    }

    template< class T_A >
    void operator()(OutputVector< T_Scalar, T_Degree > &values, const InputVector< T_Scalar, T_Degree, T_A > &a, const std::vector< scalar::Precision< T_Scalar >, numa::allocator< scalar::Precision< T_Scalar > > > &points, const std::vector< scalar::Precision< T_Scalar > > &gaussPoints, const int elementType, const std::pair< int, int > &entity) const
    {
#pragma omp for
      for(auto i = 0ULL; i < values.size(); ++i) {
        values[i] = std::sinh(a[i]);
      }
    }

    bool canUseSameVectorsForOutputAndInputs() const override
    {
      return true;
    }

    std::string name() const override
    {
      return "sinh(a)";
    }
  };

  //
  // Sqrt
  //

  template< class T_Scalar, Degree T_Degree >
  class Sqrt final : public UnaryOperation< T_Scalar, T_Degree, T_Degree >
  {
   public:
    Sqrt()
    {
    }

    Sqrt(const Sqrt &other) :
      UnaryOperation< T_Scalar, T_Degree, T_Degree >(other)
    {
    }

    template< class T_A >
    void operator()(OutputVector< T_Scalar, T_Degree > &values, const InputVector< T_Scalar, T_Degree, T_A > &a, const std::vector< scalar::Precision< T_Scalar >, numa::allocator< scalar::Precision< T_Scalar > > > &points, const std::vector< scalar::Precision< T_Scalar > > &gaussPoints, const int elementType, const std::pair< int, int > &entity) const
    {
#pragma omp for
      for(auto i = 0ULL; i < values.size(); ++i) {
        values[i] = std::sqrt(a[i]);
      }
    }

    bool canUseSameVectorsForOutputAndInputs() const override
    {
      return true;
    }

    std::string name() const override
    {
      return "sqrt(a)";
    }
  };

  //
  // Tan
  //

  template< class T_Scalar, Degree T_Degree >
  class Tan final : public UnaryOperation< T_Scalar, T_Degree, T_Degree >
  {
   public:
    Tan()
    {
    }

    Tan(const Tan &other) :
      UnaryOperation< T_Scalar, T_Degree, T_Degree >(other)
    {
    }

    template< class T_A >
    void operator()(OutputVector< T_Scalar, T_Degree > &values, const InputVector< T_Scalar, T_Degree, T_A > &a, const std::vector< scalar::Precision< T_Scalar >, numa::allocator< scalar::Precision< T_Scalar > > > &points, const std::vector< scalar::Precision< T_Scalar > > &gaussPoints, const int elementType, const std::pair< int, int > &entity) const
    {
#pragma omp for
      for(auto i = 0ULL; i < values.size(); ++i) {
        values[i] = std::tan(a[i]);
      }
    }

    bool canUseSameVectorsForOutputAndInputs() const override
    {
      return true;
    }

    std::string name() const override
    {
      return "tan(a)";
    }
  };

  //
  // Tanh
  //

  template< class T_Scalar, Degree T_Degree >
  class Tanh final : public UnaryOperation< T_Scalar, T_Degree, T_Degree >
  {
   public:
    Tanh()
    {
    }

    Tanh(const Tanh &other) :
      UnaryOperation< T_Scalar, T_Degree, T_Degree >(other)
    {
    }

    template< class T_A >
    void operator()(OutputVector< T_Scalar, T_Degree > &values, const InputVector< T_Scalar, T_Degree, T_A > &a, const std::vector< scalar::Precision< T_Scalar >, numa::allocator< scalar::Precision< T_Scalar > > > &points, const std::vector< scalar::Precision< T_Scalar > > &gaussPoints, const int elementType, const std::pair< int, int > &entity) const
    {
#pragma omp for
      for(auto i = 0ULL; i < values.size(); ++i) {
        values[i] = std::tanh(a[i]);
      }
    }

    bool canUseSameVectorsForOutputAndInputs() const override
    {
      return true;
    }

    std::string name() const override
    {
      return "tanh(a)";
    }
  };

  //
  // Trace
  //

  template< class T_Scalar >
  class Trace final : public UnaryOperation< T_Scalar, Degree::Degree0, Degree::Degree2 >
  {
   public:
    Trace()
    {
    }

    Trace(const Trace &other) :
      UnaryOperation< T_Scalar, Degree::Degree0, Degree::Degree2 >(other)
    {
    }

    template< class T_A >
    void operator()(OutputVector< T_Scalar, Degree::Degree0 > &values, const InputVector< T_Scalar, Degree::Degree2, T_A > &a, const std::vector< scalar::Precision< T_Scalar >, numa::allocator< scalar::Precision< T_Scalar > > > &points, const std::vector< scalar::Precision< T_Scalar > > &gaussPoints, const int elementType, const std::pair< int, int > &entity) const
    {
#pragma omp for
      for(auto i = 0ULL; i < values.size(); ++i) {
        values[i] = a[i].trace();
      }
    }

    bool canUseSameVectorsForOutputAndInputs() const override
    {
      return false;
    }

    std::string name() const override
    {
      return "tr(a)";
    }
  };

  //
  // Xcomp
  //

  template< class T_Scalar >
  class Xcomp final : public UnaryOperation< T_Scalar, Degree::Degree0, Degree::Degree1 >
  {
   public:
    Xcomp()
    {
    }

    Xcomp(const Xcomp &other) :
      UnaryOperation< T_Scalar, Degree::Degree0, Degree::Degree1 >(other)
    {
    }

    template< class T_A >
    void operator()(OutputVector< T_Scalar, Degree::Degree0 > &values, const InputVector< T_Scalar, Degree::Degree1, T_A > &a, const std::vector< scalar::Precision< T_Scalar >, numa::allocator< scalar::Precision< T_Scalar > > > &points, const std::vector< scalar::Precision< T_Scalar > > &gaussPoints, const int elementType, const std::pair< int, int > &entity) const
    {
#pragma omp for
      for(auto i = 0ULL; i < values.size(); ++i) {
        values[i] = a[i](0);
      }
    }

    std::string name() const override
    {
      return "x component(a)";
    }
  };

  //
  // XXcomp
  //

  template< class T_Scalar >
  class XXcomp final : public UnaryOperation< T_Scalar, Degree::Degree0, Degree::Degree2 >
  {
   public:
    XXcomp()
    {
    }

    XXcomp(const XXcomp &other) :
      UnaryOperation< T_Scalar, Degree::Degree0, Degree::Degree2 >(other)
    {
    }

    template< class T_A >
    void operator()(OutputVector< T_Scalar, Degree::Degree0 > &values, const InputVector< T_Scalar, Degree::Degree2, T_A > &a, const std::vector< scalar::Precision< T_Scalar >, numa::allocator< scalar::Precision< T_Scalar > > > &points, const std::vector< scalar::Precision< T_Scalar > > &gaussPoints, const int elementType, const std::pair< int, int > &entity) const
    {
#pragma omp for
      for(auto i = 0ULL; i < values.size(); ++i) {
        values[i] = a[i](0, 0);
      }
    }

    std::string name() const override
    {
      return "xx component(a)";
    }
  };

  //
  // XYcomp
  //

  template< class T_Scalar >
  class XYcomp final : public UnaryOperation< T_Scalar, Degree::Degree0, Degree::Degree2 >
  {
   public:
    XYcomp()
    {
    }

    XYcomp(const XYcomp &other) :
      UnaryOperation< T_Scalar, Degree::Degree0, Degree::Degree2 >(other)
    {
    }

    template< class T_A >
    void operator()(OutputVector< T_Scalar, Degree::Degree0 > &values, const InputVector< T_Scalar, Degree::Degree2, T_A > &a, const std::vector< scalar::Precision< T_Scalar >, numa::allocator< scalar::Precision< T_Scalar > > > &points, const std::vector< scalar::Precision< T_Scalar > > &gaussPoints, const int elementType, const std::pair< int, int > &entity) const
    {
#pragma omp for
      for(auto i = 0ULL; i < values.size(); ++i) {
        values[i] = a[i](0, 1);
      }
    }

    std::string name() const override
    {
      return "xy component(a)";
    }
  };

  //
  // XZcomp
  //

  template< class T_Scalar >
  class XZcomp final : public UnaryOperation< T_Scalar, Degree::Degree0, Degree::Degree2 >
  {
   public:
    XZcomp()
    {
    }

    XZcomp(const XZcomp &other) :
      UnaryOperation< T_Scalar, Degree::Degree0, Degree::Degree2 >(other)
    {
    }

    template< class T_A >
    void operator()(OutputVector< T_Scalar, Degree::Degree0 > &values, const InputVector< T_Scalar, Degree::Degree2, T_A > &a, const std::vector< scalar::Precision< T_Scalar >, numa::allocator< scalar::Precision< T_Scalar > > > &points, const std::vector< scalar::Precision< T_Scalar > > &gaussPoints, const int elementType, const std::pair< int, int > &entity) const
    {
#pragma omp for
      for(auto i = 0ULL; i < values.size(); ++i) {
        values[i] = a[i](0, 2);
      }
    }

    std::string name() const override
    {
      return "xz component(a)";
    }
  };

  //
  // Ycomp
  //

  template< class T_Scalar >
  class Ycomp final : public UnaryOperation< T_Scalar, Degree::Degree0, Degree::Degree1 >
  {
   public:
    Ycomp()
    {
    }

    Ycomp(const Ycomp &other) :
      UnaryOperation< T_Scalar, Degree::Degree0, Degree::Degree1 >(other)
    {
    }

    template< class T_A >
    void operator()(OutputVector< T_Scalar, Degree::Degree0 > &values, const InputVector< T_Scalar, Degree::Degree1, T_A > &a, const std::vector< scalar::Precision< T_Scalar >, numa::allocator< scalar::Precision< T_Scalar > > > &points, const std::vector< scalar::Precision< T_Scalar > > &gaussPoints, const int elementType, const std::pair< int, int > &entity) const
    {
#pragma omp for
      for(auto i = 0ULL; i < values.size(); ++i) {
        values[i] = a[i](1);
      }
    }

    std::string name() const override
    {
      return "y component(a)";
    }
  };

  //
  // YXcomp
  //

  template< class T_Scalar >
  class YXcomp final : public UnaryOperation< T_Scalar, Degree::Degree0, Degree::Degree2 >
  {
   public:
    YXcomp()
    {
    }

    YXcomp(const YXcomp &other) :
      UnaryOperation< T_Scalar, Degree::Degree0, Degree::Degree2 >(other)
    {
    }

    template< class T_A >
    void operator()(OutputVector< T_Scalar, Degree::Degree0 > &values, const InputVector< T_Scalar, Degree::Degree2, T_A > &a, const std::vector< scalar::Precision< T_Scalar >, numa::allocator< scalar::Precision< T_Scalar > > > &points, const std::vector< scalar::Precision< T_Scalar > > &gaussPoints, const int elementType, const std::pair< int, int > &entity) const
    {
#pragma omp for
      for(auto i = 0ULL; i < values.size(); ++i) {
        values[i] = a[i](1, 0);
      }
    }

    std::string name() const override
    {
      return "yx component(a)";
    }
  };

  //
  // YYcomp
  //

  template< class T_Scalar >
  class YYcomp final : public UnaryOperation< T_Scalar, Degree::Degree0, Degree::Degree2 >
  {
   public:
    YYcomp()
    {
    }

    YYcomp(const YYcomp &other) :
      UnaryOperation< T_Scalar, Degree::Degree0, Degree::Degree2 >(other)
    {
    }

    template< class T_A >
    void operator()(OutputVector< T_Scalar, Degree::Degree0 > &values, const InputVector< T_Scalar, Degree::Degree2, T_A > &a, const std::vector< scalar::Precision< T_Scalar >, numa::allocator< scalar::Precision< T_Scalar > > > &points, const std::vector< scalar::Precision< T_Scalar > > &gaussPoints, const int elementType, const std::pair< int, int > &entity) const
    {
#pragma omp for
      for(auto i = 0ULL; i < values.size(); ++i) {
        values[i] = a[i](1, 1);
      }
    }

    std::string name() const override
    {
      return "yy component(a)";
    }
  };

  //
  // YZcomp
  //

  template< class T_Scalar >
  class YZcomp final : public UnaryOperation< T_Scalar, Degree::Degree0, Degree::Degree2 >
  {
   public:
    YZcomp()
    {
    }

    YZcomp(const YZcomp &other) :
      UnaryOperation< T_Scalar, Degree::Degree0, Degree::Degree2 >(other)
    {
    }

    template< class T_A >
    void operator()(OutputVector< T_Scalar, Degree::Degree0 > &values, const InputVector< T_Scalar, Degree::Degree2, T_A > &a, const std::vector< scalar::Precision< T_Scalar >, numa::allocator< scalar::Precision< T_Scalar > > > &points, const std::vector< scalar::Precision< T_Scalar > > &gaussPoints, const int elementType, const std::pair< int, int > &entity) const
    {
#pragma omp for
      for(auto i = 0ULL; i < values.size(); ++i) {
        values[i] = a[i](1, 2);
      }
    }

    std::string name() const override
    {
      return "yz component(a)";
    }
  };

  //
  // Zcomp
  //

  template< class T_Scalar >
  class Zcomp final : public UnaryOperation< T_Scalar, Degree::Degree0, Degree::Degree1 >
  {
   public:
    Zcomp()
    {
    }

    Zcomp(const Zcomp &other) :
      UnaryOperation< T_Scalar, Degree::Degree0, Degree::Degree1 >(other)
    {
    }

    template< class T_A >
    void operator()(OutputVector< T_Scalar, Degree::Degree0 > &values, const InputVector< T_Scalar, Degree::Degree1, T_A > &a, const std::vector< scalar::Precision< T_Scalar >, numa::allocator< scalar::Precision< T_Scalar > > > &points, const std::vector< scalar::Precision< T_Scalar > > &gaussPoints, const int elementType, const std::pair< int, int > &entity) const
    {
#pragma omp for
      for(auto i = 0ULL; i < values.size(); ++i) {
        values[i] = a[i](2);
      }
    }

    std::string name() const override
    {
      return "z component(a)";
    }
  };

  //
  // ZXcomp
  //

  template< class T_Scalar >
  class ZXcomp final : public UnaryOperation< T_Scalar, Degree::Degree0, Degree::Degree2 >
  {
   public:
    ZXcomp()
    {
    }

    ZXcomp(const ZXcomp &other) :
      UnaryOperation< T_Scalar, Degree::Degree0, Degree::Degree2 >(other)
    {
    }

    template< class T_A >
    void operator()(OutputVector< T_Scalar, Degree::Degree0 > &values, const InputVector< T_Scalar, Degree::Degree2, T_A > &a, const std::vector< scalar::Precision< T_Scalar >, numa::allocator< scalar::Precision< T_Scalar > > > &points, const std::vector< scalar::Precision< T_Scalar > > &gaussPoints, const int elementType, const std::pair< int, int > &entity) const
    {
#pragma omp for
      for(auto i = 0ULL; i < values.size(); ++i) {
        values[i] = a[i](2, 0);
      }
    }

    std::string name() const override
    {
      return "zx component(a)";
    }
  };

  //
  // ZYcomp
  //

  template< class T_Scalar >
  class ZYcomp final : public UnaryOperation< T_Scalar, Degree::Degree0, Degree::Degree2 >
  {
   public:
    ZYcomp()
    {
    }

    ZYcomp(const ZYcomp &other) :
      UnaryOperation< T_Scalar, Degree::Degree0, Degree::Degree2 >(other)
    {
    }

    template< class T_A >
    void operator()(OutputVector< T_Scalar, Degree::Degree0 > &values, const InputVector< T_Scalar, Degree::Degree2, T_A > &a, const std::vector< scalar::Precision< T_Scalar >, numa::allocator< scalar::Precision< T_Scalar > > > &points, const std::vector< scalar::Precision< T_Scalar > > &gaussPoints, const int elementType, const std::pair< int, int > &entity) const
    {
#pragma omp for
      for(auto i = 0ULL; i < values.size(); ++i) {
        values[i] = a[i](2, 1);
      }
    }

    std::string name() const override
    {
      return "zy component(a)";
    }
  };

  //
  // ZZcomp
  //

  template< class T_Scalar >
  class ZZcomp final : public UnaryOperation< T_Scalar, Degree::Degree0, Degree::Degree2 >
  {
   public:
    ZZcomp()
    {
    }

    ZZcomp(const ZZcomp &other) :
      UnaryOperation< T_Scalar, Degree::Degree0, Degree::Degree2 >(other)
    {
    }

    template< class T_A >
    void operator()(OutputVector< T_Scalar, Degree::Degree0 > &values, const InputVector< T_Scalar, Degree::Degree2, T_A > &a, const std::vector< scalar::Precision< T_Scalar >, numa::allocator< scalar::Precision< T_Scalar > > > &points, const std::vector< scalar::Precision< T_Scalar > > &gaussPoints, const int elementType, const std::pair< int, int > &entity) const
    {
#pragma omp for
      for(auto i = 0ULL; i < values.size(); ++i) {
        values[i] = a[i](2, 2);
      }
    }

    std::string name() const override
    {
      return "zz component(a)";
    }
  };


} // namespace gmshfem::function


#endif // H_GMSHFEM_UNARYOPERATIONS
