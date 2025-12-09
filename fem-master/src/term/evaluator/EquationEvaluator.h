// GmshFEM - Copyright (C) 2019-2022, A. Royer, E. Béchet, C. Geuzaine, Université de Liège
//
// See the LICENSE.txt file for license information. Please report all
// issues on https://gitlab.onelab.info/gmsh/fem/issues

#ifndef H_GMSHFEM_EQUATIONEVALUATOR
#define H_GMSHFEM_EQUATIONEVALUATOR

#include "FieldObject.h"
#include "Function.h"
#include "Memory.h"
#include "OmpInterface.h"
#include "Product.h"
#include "numa.h"

#include <Eigen/Dense>
#include <iostream>

namespace gmshfem::term
{
  enum class ProductType {
    Hermitian,
    Scalar
  };
}

namespace gmshfem::term::evaluator
{


  enum class EvaluationOrder {
    Left,
    Right
  };

  template< class T_Scalar >
  static inline Eigen::MatrixX< T_Scalar > cross(const typename MathObject< T_Scalar, Degree::Degree1 >::Object &a, const Eigen::MatrixX< scalar::Precision< T_Scalar > > &b)
  {
    Eigen::MatrixX< T_Scalar > result = b;
    const unsigned int size = b.cols();
    for(auto i = 0U; i < size; ++i) {
      typename MathObject< T_Scalar, Degree::Degree1 >::Object vec = b.col(i);
      result.col(i) = a.cross(vec).conjugate();
    }
    return result;
  }

  template< class T_Scalar >
  static inline Eigen::MatrixX< T_Scalar > crossT(const typename MathObject< T_Scalar, Degree::Degree1 >::Object &a, const Eigen::MatrixX< T_Scalar > &b)
  {
    Eigen::MatrixX< T_Scalar > result = b;
    const unsigned int size = b.cols();
    for(auto i = 0U; i < size; ++i) {
      typename MathObject< T_Scalar, Degree::Degree1 >::Object vec = b.col(i);
      result.col(i) = a.cross(vec).conjugate();
    }
    return result;
  }

  template< class T_Scalar >
  static inline Eigen::MatrixX< T_Scalar > cross(const Eigen::MatrixX< scalar::Precision< T_Scalar > > &a, const typename MathObject< T_Scalar, Degree::Degree1 >::Object &b)
  {
    Eigen::MatrixX< T_Scalar > result = a;
    const unsigned int size = a.cols();
    for(auto i = 0U; i < size; ++i) {
      typename MathObject< T_Scalar, Degree::Degree1 >::Object vec = a.col(i);
      result.col(i) = vec.cross(b).conjugate();
    }
    return result;
  }

  template< class T_Scalar >
  static Eigen::MatrixX< T_Scalar > crossT(const Eigen::MatrixX< T_Scalar > &a, const typename MathObject< T_Scalar, Degree::Degree1 >::Object &b)
  {
    Eigen::MatrixX< T_Scalar > result = a;
    const unsigned int size = a.cols();
    for(auto i = 0U; i < size; ++i) {
      typename MathObject< T_Scalar, Degree::Degree1 >::Object vec = a.col(i);
      result.col(i) = vec.cross(b).conjugate();
    }
    return result;
  }

  template< class T_Scalar, field::Form T_Form >
  class EquationEvaluatorInterface
  {
   public:
    EquationEvaluatorInterface() {}
    virtual ~EquationEvaluatorInterface() {}

    virtual bool isConstant(const std::pair< int, int > &entity) const = 0;

    virtual void initialize(const std::vector< scalar::Precision< T_Scalar >, numa::allocator< scalar::Precision< T_Scalar > > > &points, const std::vector< scalar::Precision< T_Scalar > > &gaussPoints, const int elementType, const std::pair< int, int > &entity) const = 0;
    virtual void operator()(Eigen::Ref< Eigen::MatrixX< T_Scalar > > expression, const Eigen::Ref< const Eigen::MatrixX< scalar::Precision< T_Scalar > > > field, const unsigned int offset) const = 0;

    virtual unsigned int size1() const = 0;
    virtual unsigned int size2() const = 0;
    virtual void clear() const = 0;
    virtual common::Memory memory() const = 0;
  };

  // ############################
  // Field
  // ############################

  template< class T_Scalar, field::Form T_Form, Degree T_DegreeLhs, equation::Product T_ProductLhs, Degree T_DegreeRhs, equation::Product T_ProductRhs, EvaluationOrder T_EvaluationOrder > // T_EvaluationOrder = Left -> (T_DegreeLhs * dof) * T_DegreeRhs, = Right - > T_DegreeLhs * (dof * T_DegreeRhs)
  class EquationEvaluator : public EquationEvaluatorInterface< T_Scalar, T_Form >
  {
   private:
    function::Function< T_Scalar, T_DegreeLhs > _left;
    function::Function< T_Scalar, T_DegreeLhs > _right;

   public:
    EquationEvaluator(const function::Function< T_Scalar, T_DegreeLhs > &left, const function::Function< T_Scalar, T_DegreeLhs > &right) :
      _left(left), _right(right) {}
    ~EquationEvaluator() {}

    virtual bool isConstant(const std::pair< int, int > &entity) const = 0;

    virtual void initialize(const std::vector< scalar::Precision< T_Scalar >, numa::allocator< scalar::Precision< T_Scalar > > > &points, const std::vector< scalar::Precision< T_Scalar > > &gaussPoints, const int elementType, const std::pair< int, int > &entity) const = 0;
    virtual void operator()(Eigen::Ref< Eigen::MatrixX< T_Scalar > > expression, const Eigen::Ref< const Eigen::MatrixX< scalar::Precision< T_Scalar > > > field, const unsigned int offset) const {}

    virtual unsigned int size1() const = 0;
    virtual unsigned int size2() const = 0;
    virtual void clear() const = 0;
    virtual common::Memory memory() const = 0;
  };

  //
  // Left
  //

  template< class T_Scalar, field::Form T_Form, Degree T_DegreeLhs >
  class EquationEvaluator< T_Scalar, T_Form, T_DegreeLhs, equation::Product::ScalarProduct, Degree::Empty, equation::Product::Empty, EvaluationOrder::Left > final : public EquationEvaluatorInterface< T_Scalar, T_Form >
  {
   private:
    function::Function< T_Scalar, T_DegreeLhs > _left;
    mutable std::vector< typename MathObject< T_Scalar, T_DegreeLhs >::Object, numa::allocator< typename MathObject< T_Scalar, T_DegreeLhs >::Object > > _values;
    mutable bool _isConstant;

   public:
    EquationEvaluator(const function::Function< T_Scalar, T_DegreeLhs > &left) :
      _left(left), _values(), _isConstant(false)
    {
    }
    ~EquationEvaluator()
    {
    }

    bool isConstant(const std::pair< int, int > &entity) const override
    {
      return _left.isConstant(entity);
    }

    void initialize(const std::vector< scalar::Precision< T_Scalar >, numa::allocator< scalar::Precision< T_Scalar > > > &points, const std::vector< scalar::Precision< T_Scalar > > &gaussPoints, const int elementType, const std::pair< int, int > &entity) const override
    {
      _isConstant = isConstant(entity);
      if(_isConstant) {
        _values.resize(1);
        _left.evaluate(_values[0], 0., 0., 0., entity);
      }
      else {
#pragma omp parallel num_threads(omp::getMaxThreads())
        {
          _left.evaluate(_values, points, gaussPoints, elementType, entity);
        }
      }
    }

    void operator()(Eigen::Ref< Eigen::MatrixX< T_Scalar > > expression, const Eigen::Ref< const Eigen::MatrixX< scalar::Precision< T_Scalar > > > field, const unsigned int offset) const override
    {
      expression = _values[_isConstant ? 0 : offset] * field;
    }

    unsigned int size1() const override
    {
      return 4; // means field size
    }

    unsigned int size2() const override
    {
      return 0; // means number of dofs by element
    }

    void clear() const override
    {
      _values.clear();
      _isConstant = false;
    }

    common::Memory memory() const override
    {
      return common::Memory(_values.size() * sizeof(typename MathObject< T_Scalar, T_DegreeLhs >::Object));
    }
  };

  template< class T_Scalar, field::Form T_Form >
  class EquationEvaluator< T_Scalar, T_Form, Degree::Degree1, equation::Product::ScalarProduct, Degree::Empty, equation::Product::Empty, EvaluationOrder::Left > final : public EquationEvaluatorInterface< T_Scalar, T_Form >
  {
   private:
    function::Function< T_Scalar, Degree::Degree1 > _left;
    mutable std::vector< typename MathObject< T_Scalar, Degree::Degree1 >::Object, numa::allocator< typename MathObject< T_Scalar, Degree::Degree1 >::Object > > _values;
    mutable bool _isConstant;

   public:
    EquationEvaluator(const function::Function< T_Scalar, Degree::Degree1 > &left) :
      _left(left), _values(), _isConstant(false)
    {
    }
    ~EquationEvaluator()
    {
    }

    bool isConstant(const std::pair< int, int > &entity) const override
    {
      return _left.isConstant(entity);
    }

    void initialize(const std::vector< scalar::Precision< T_Scalar >, numa::allocator< scalar::Precision< T_Scalar > > > &points, const std::vector< scalar::Precision< T_Scalar > > &gaussPoints, const int elementType, const std::pair< int, int > &entity) const override
    {
      _isConstant = isConstant(entity);
      if(_isConstant) {
        _values.resize(1);
        _left.evaluate(_values[0], 0., 0., 0., entity);
      }
      else {
#pragma omp parallel num_threads(omp::getMaxThreads())
        {
          _left.evaluate(_values, points, gaussPoints, elementType, entity);
        }
      }
    }

    void operator()(Eigen::Ref< Eigen::MatrixX< T_Scalar > > expression, const Eigen::Ref< const Eigen::MatrixX< scalar::Precision< T_Scalar > > > field, const unsigned int offset) const override
    {
      expression = _values[_isConstant ? 0 : offset].transpose() * field;
    }

    unsigned int size1() const override
    {
      return 1;
    }

    unsigned int size2() const override
    {
      return 0; // means number of dofs by element
    }

    void clear() const override
    {
      _values.clear();
      _isConstant = false;
    }

    common::Memory memory() const override
    {
      return common::Memory(_values.size() * sizeof(typename MathObject< T_Scalar, Degree::Degree1 >::Object));
    }
  };

  template< class T_Scalar >
  class EquationEvaluator< T_Scalar, field::Form::Form0, Degree::Degree1, equation::Product::ScalarProduct, Degree::Empty, equation::Product::Empty, EvaluationOrder::Left > final : public EquationEvaluatorInterface< T_Scalar, field::Form::Form0 >
  {
   private:
    function::Function< T_Scalar, Degree::Degree1 > _left;
    mutable std::vector< typename MathObject< T_Scalar, Degree::Degree1 >::Object, numa::allocator< typename MathObject< T_Scalar, Degree::Degree1 >::Object > > _values;
    mutable bool _isConstant;

   public:
    EquationEvaluator(const function::Function< T_Scalar, Degree::Degree1 > &left) :
      _left(left), _values(), _isConstant(false)
    {
    }
    ~EquationEvaluator()
    {
    }

    bool isConstant(const std::pair< int, int > &entity) const override
    {
      return _left.isConstant(entity);
    }

    void initialize(const std::vector< scalar::Precision< T_Scalar >, numa::allocator< scalar::Precision< T_Scalar > > > &points, const std::vector< scalar::Precision< T_Scalar > > &gaussPoints, const int elementType, const std::pair< int, int > &entity) const override
    {
      _isConstant = isConstant(entity);
      if(_isConstant) {
        _values.resize(1);
        _left.evaluate(_values[0], 0., 0., 0., entity);
      }
      else {
#pragma omp parallel num_threads(omp::getMaxThreads())
        {
          _left.evaluate(_values, points, gaussPoints, elementType, entity);
        }
      }
    }

    void operator()(Eigen::Ref< Eigen::MatrixX< T_Scalar > > expression, const Eigen::Ref< const Eigen::MatrixX< scalar::Precision< T_Scalar > > > field, const unsigned int offset) const override
    {
      expression = _values[_isConstant ? 0 : offset] * field;
    }

    unsigned int size1() const override
    {
      return 3;
    }

    unsigned int size2() const override
    {
      return 0; // means number of dofs by element
    }

    void clear() const override
    {
      _values.clear();
      _isConstant = false;
    }

    common::Memory memory() const override
    {
      return common::Memory(_values.size() * sizeof(typename MathObject< T_Scalar, Degree::Degree1 >::Object));
    }
  };

  template< class T_Scalar >
  class EquationEvaluator< T_Scalar, field::Form::Form3, Degree::Degree1, equation::Product::ScalarProduct, Degree::Empty, equation::Product::Empty, EvaluationOrder::Left > final : public EquationEvaluatorInterface< T_Scalar, field::Form::Form3 >
  {
   private:
    function::Function< T_Scalar, Degree::Degree1 > _left;
    mutable std::vector< typename MathObject< T_Scalar, Degree::Degree1 >::Object, numa::allocator< typename MathObject< T_Scalar, Degree::Degree1 >::Object > > _values;
    mutable bool _isConstant;

   public:
    EquationEvaluator(const function::Function< T_Scalar, Degree::Degree1 > &left) :
      _left(left), _values(), _isConstant(false)
    {
    }
    ~EquationEvaluator()
    {
    }

    bool isConstant(const std::pair< int, int > &entity) const override
    {
      return _left.isConstant(entity);
    }

    void initialize(const std::vector< scalar::Precision< T_Scalar >, numa::allocator< scalar::Precision< T_Scalar > > > &points, const std::vector< scalar::Precision< T_Scalar > > &gaussPoints, const int elementType, const std::pair< int, int > &entity) const override
    {
      _isConstant = isConstant(entity);
      if(_isConstant) {
        _values.resize(1);
        _left.evaluate(_values[0], 0., 0., 0., entity);
      }
      else {
#pragma omp parallel num_threads(omp::getMaxThreads())
        {
          _left.evaluate(_values, points, gaussPoints, elementType, entity);
        }
      }
    }

    void operator()(Eigen::Ref< Eigen::MatrixX< T_Scalar > > expression, const Eigen::Ref< const Eigen::MatrixX< scalar::Precision< T_Scalar > > > field, const unsigned int offset) const override
    {
      expression = _values[_isConstant ? 0 : offset] * field;
    }

    unsigned int size1() const override
    {
      return 3;
    }

    unsigned int size2() const override
    {
      return 0; // means number of dofs by element
    }

    void clear() const override
    {
      _values.clear();
      _isConstant = false;
    }

    common::Memory memory() const override
    {
      return common::Memory(_values.size() * sizeof(typename MathObject< T_Scalar, Degree::Degree1 >::Object));
    }
  };

  template< class T_Scalar, field::Form T_Form >
  class EquationEvaluator< T_Scalar, T_Form, Degree::Degree1, equation::Product::VectorProduct, Degree::Empty, equation::Product::Empty, EvaluationOrder::Left > final : public EquationEvaluatorInterface< T_Scalar, T_Form >
  {
   private:
    function::Function< T_Scalar, Degree::Degree1 > _left;
    mutable std::vector< typename MathObject< T_Scalar, Degree::Degree1 >::Object, numa::allocator< typename MathObject< T_Scalar, Degree::Degree1 >::Object > > _values;
    mutable bool _isConstant;

   public:
    EquationEvaluator(const function::Function< T_Scalar, Degree::Degree1 > &left) :
      _left(left), _values(), _isConstant(false)
    {
    }
    ~EquationEvaluator()
    {
    }

    bool isConstant(const std::pair< int, int > &entity) const override
    {
      return _left.isConstant(entity);
    }

    void initialize(const std::vector< scalar::Precision< T_Scalar >, numa::allocator< scalar::Precision< T_Scalar > > > &points, const std::vector< scalar::Precision< T_Scalar > > &gaussPoints, const int elementType, const std::pair< int, int > &entity) const override
    {
      _isConstant = isConstant(entity);
      if(_isConstant) {
        _values.resize(1);
        _left.evaluate(_values[0], 0., 0., 0., entity);
      }
      else {
#pragma omp parallel num_threads(omp::getMaxThreads())
        {
          _left.evaluate(_values, points, gaussPoints, elementType, entity);
        }
      }
    }

    void operator()(Eigen::Ref< Eigen::MatrixX< T_Scalar > > expression, const Eigen::Ref< const Eigen::MatrixX< scalar::Precision< T_Scalar > > > field, const unsigned int offset) const override
    {
      expression = cross< T_Scalar >(_values[_isConstant ? 0 : offset], field);
    }

    unsigned int size1() const override
    {
      return 3;
    }

    unsigned int size2() const override
    {
      return 0; // means number of dofs by element
    }

    void clear() const override
    {
      _values.clear();
      _isConstant = false;
    }

    common::Memory memory() const override
    {
      return common::Memory(_values.size() * sizeof(typename MathObject< T_Scalar, Degree::Degree1 >::Object));
    }
  };

  //
  // Right
  //

  template< class T_Scalar, field::Form T_Form, Degree T_DegreeRhs >
  class EquationEvaluator< T_Scalar, T_Form, Degree::Empty, equation::Product::Empty, T_DegreeRhs, equation::Product::ScalarProduct, EvaluationOrder::Right > final : public EquationEvaluatorInterface< T_Scalar, T_Form >
  {
   private:
    function::Function< T_Scalar, T_DegreeRhs > _right;
    mutable std::vector< typename MathObject< T_Scalar, T_DegreeRhs >::Object, numa::allocator< typename MathObject< T_Scalar, T_DegreeRhs >::Object > > _values;
    mutable bool _isConstant;

   public:
    EquationEvaluator(const function::Function< T_Scalar, T_DegreeRhs > &right) :
      _right(right), _values(), _isConstant(false)
    {
    }
    ~EquationEvaluator()
    {
    }

    bool isConstant(const std::pair< int, int > &entity) const override
    {
      return _right.isConstant(entity);
    }

    void initialize(const std::vector< scalar::Precision< T_Scalar >, numa::allocator< scalar::Precision< T_Scalar > > > &points, const std::vector< scalar::Precision< T_Scalar > > &gaussPoints, const int elementType, const std::pair< int, int > &entity) const override
    {
      _isConstant = isConstant(entity);
      if(_isConstant) {
        _values.resize(1);
        _right.evaluate(_values[0], 0., 0., 0., entity);
      }
      else {
#pragma omp parallel num_threads(omp::getMaxThreads())
        {
          _right.evaluate(_values, points, gaussPoints, elementType, entity);
        }
      }
    }

    void operator()(Eigen::Ref< Eigen::MatrixX< T_Scalar > > expression, const Eigen::Ref< const Eigen::MatrixX< scalar::Precision< T_Scalar > > > field, const unsigned int offset) const override
    {
      expression = field.transpose() * _values[_isConstant ? 0 : offset];
    }

    unsigned int size1() const override
    {
      return 0; // means number of dofs by element
    }

    unsigned int size2() const override
    {
      return 4; // means field size
    }

    void clear() const override
    {
      _values.clear();
      _isConstant = false;
    }

    common::Memory memory() const override
    {
      return common::Memory(_values.size() * sizeof(typename MathObject< T_Scalar, T_DegreeRhs >::Object));
    }
  };

  template< class T_Scalar >
  class EquationEvaluator< T_Scalar, field::Form::Form0, Degree::Empty, equation::Product::Empty, Degree::Degree1, equation::Product::ScalarProduct, EvaluationOrder::Right > final : public EquationEvaluatorInterface< T_Scalar, field::Form::Form0 >
  {
   private:
    function::Function< T_Scalar, Degree::Degree1 > _right;
    mutable std::vector< typename MathObject< T_Scalar, Degree::Degree1 >::Object, numa::allocator< typename MathObject< T_Scalar, Degree::Degree1 >::Object > > _values;
    mutable bool _isConstant;

   public:
    EquationEvaluator(const function::Function< T_Scalar, Degree::Degree1 > &right) :
      _right(right), _values(), _isConstant(false)
    {
    }
    ~EquationEvaluator()
    {
    }

    bool isConstant(const std::pair< int, int > &entity) const override
    {
      return _right.isConstant(entity);
    }

    void initialize(const std::vector< scalar::Precision< T_Scalar >, numa::allocator< scalar::Precision< T_Scalar > > > &points, const std::vector< scalar::Precision< T_Scalar > > &gaussPoints, const int elementType, const std::pair< int, int > &entity) const override
    {
      _isConstant = isConstant(entity);
      if(_isConstant) {
        _values.resize(1);
        _right.evaluate(_values[0], 0., 0., 0., entity);
      }
      else {
#pragma omp parallel num_threads(omp::getMaxThreads())
        {
          _right.evaluate(_values, points, gaussPoints, elementType, entity);
        }
      }
    }

    void operator()(Eigen::Ref< Eigen::MatrixX< T_Scalar > > expression, const Eigen::Ref< const Eigen::MatrixX< scalar::Precision< T_Scalar > > > field, const unsigned int offset) const override
    {
      expression = field.transpose() * _values[_isConstant ? 0 : offset].transpose();
    }

    unsigned int size1() const override
    {
      return 0; // means number of dofs by element
    }

    unsigned int size2() const override
    {
      return 3;
    }

    void clear() const override
    {
      _values.clear();
      _isConstant = false;
    }

    common::Memory memory() const override
    {
      return common::Memory(_values.size() * sizeof(typename MathObject< T_Scalar, Degree::Degree1 >::Object));
    }
  };

  template< class T_Scalar >
  class EquationEvaluator< T_Scalar, field::Form::Form3, Degree::Empty, equation::Product::Empty, Degree::Degree1, equation::Product::ScalarProduct, EvaluationOrder::Right > final : public EquationEvaluatorInterface< T_Scalar, field::Form::Form3 >
  {
   private:
    function::Function< T_Scalar, Degree::Degree1 > _right;
    mutable std::vector< typename MathObject< T_Scalar, Degree::Degree1 >::Object, numa::allocator< typename MathObject< T_Scalar, Degree::Degree1 >::Object > > _values;
    mutable bool _isConstant;

   public:
    EquationEvaluator(const function::Function< T_Scalar, Degree::Degree1 > &right) :
      _right(right), _values(), _isConstant(false)
    {
    }
    ~EquationEvaluator()
    {
    }

    bool isConstant(const std::pair< int, int > &entity) const override
    {
      return _right.isConstant(entity);
    }

    void initialize(const std::vector< scalar::Precision< T_Scalar >, numa::allocator< scalar::Precision< T_Scalar > > > &points, const std::vector< scalar::Precision< T_Scalar > > &gaussPoints, const int elementType, const std::pair< int, int > &entity) const override
    {
      _isConstant = isConstant(entity);
      if(_isConstant) {
        _values.resize(1);
        _right.evaluate(_values[0], 0., 0., 0., entity);
      }
      else {
#pragma omp parallel num_threads(omp::getMaxThreads())
        {
          _right.evaluate(_values, points, gaussPoints, elementType, entity);
        }
      }
    }

    void operator()(Eigen::Ref< Eigen::MatrixX< T_Scalar > > expression, const Eigen::Ref< const Eigen::MatrixX< scalar::Precision< T_Scalar > > > field, const unsigned int offset) const override
    {
      expression = field.transpose() * _values[_isConstant ? 0 : offset].transpose();
    }

    unsigned int size1() const override
    {
      return 0; // means number of dofs by element
    }

    unsigned int size2() const override
    {
      return 3;
    }

    void clear() const override
    {
      _values.clear();
      _isConstant = false;
    }

    common::Memory memory() const override
    {
      return common::Memory(_values.size() * sizeof(typename MathObject< T_Scalar, Degree::Degree1 >::Object));
    }
  };

  template< class T_Scalar, field::Form T_Form, Degree T_DegreeRhs >
  class EquationEvaluator< T_Scalar, T_Form, Degree::Empty, equation::Product::Empty, T_DegreeRhs, equation::Product::VectorProduct, EvaluationOrder::Right > final : public EquationEvaluatorInterface< T_Scalar, T_Form >
  {
   private:
    function::Function< T_Scalar, T_DegreeRhs > _right;
    mutable std::vector< typename MathObject< T_Scalar, T_DegreeRhs >::Object, numa::allocator< typename MathObject< T_Scalar, T_DegreeRhs >::Object > > _values;
    mutable bool _isConstant;

   public:
    EquationEvaluator(const function::Function< T_Scalar, T_DegreeRhs > &right) :
      _right(right), _values(), _isConstant(false)
    {
    }
    ~EquationEvaluator()
    {
    }

    bool isConstant(const std::pair< int, int > &entity) const override
    {
      return _right.isConstant(entity);
    }

    void initialize(const std::vector< scalar::Precision< T_Scalar >, numa::allocator< scalar::Precision< T_Scalar > > > &points, const std::vector< scalar::Precision< T_Scalar > > &gaussPoints, const int elementType, const std::pair< int, int > &entity) const override
    {
      _isConstant = isConstant(entity);
      if(_isConstant) {
        _values.resize(1);
        _right.evaluate(_values[0], 0., 0., 0., entity);
      }
      else {
#pragma omp parallel num_threads(omp::getMaxThreads())
        {
          _right.evaluate(_values, points, gaussPoints, elementType, entity);
        }
      }
    }

    void operator()(Eigen::Ref< Eigen::MatrixX< T_Scalar > > expression, const Eigen::Ref< const Eigen::MatrixX< scalar::Precision< T_Scalar > > > field, const unsigned int offset) const override
    {
      expression = cross< T_Scalar >(field, _values[_isConstant ? 0 : offset]);
    }

    void clear() const override
    {
      _values.clear();
      _isConstant = false;
    }

    unsigned int size1() const override
    {
      return 3;
    }

    unsigned int size2() const override
    {
      return 0; // means number of dofs by element
    }

    common::Memory memory() const override
    {
      return common::Memory(_values.size() * sizeof(typename MathObject< T_Scalar, T_DegreeRhs >::Object));
    }
  };


  //
  // Left and Right
  //

  template< class T_Scalar, field::Form T_Form >
  class EquationEvaluator< T_Scalar, T_Form, Degree::Degree1, equation::Product::VectorProduct, Degree::Degree1, equation::Product::VectorProduct, EvaluationOrder::Right > final : public EquationEvaluatorInterface< T_Scalar, T_Form >
  {
   private:
    function::Function< T_Scalar, Degree::Degree1 > _left;
    function::Function< T_Scalar, Degree::Degree1 > _right;
    mutable std::vector< typename MathObject< T_Scalar, Degree::Degree1 >::Object, numa::allocator< typename MathObject< T_Scalar, Degree::Degree1 >::Object > > _valuesLeft;
    mutable std::vector< typename MathObject< T_Scalar, Degree::Degree1 >::Object, numa::allocator< typename MathObject< T_Scalar, Degree::Degree1 >::Object > > _valuesRight;
    mutable bool _isConstantLeft;
    mutable bool _isConstantRight;

   public:
    EquationEvaluator(const function::Function< T_Scalar, Degree::Degree1 > &left, const function::Function< T_Scalar, Degree::Degree1 > &right) :
      _left(left), _right(right), _valuesLeft(), _valuesRight(), _isConstantLeft(false), _isConstantRight(false)
    {
    }
    ~EquationEvaluator()
    {
    }

    bool isConstant(const std::pair< int, int > &entity) const override
    {
      return _left.isConstant(entity) && _right.isConstant(entity);
    }

    void initialize(const std::vector< scalar::Precision< T_Scalar >, numa::allocator< scalar::Precision< T_Scalar > > > &points, const std::vector< scalar::Precision< T_Scalar > > &gaussPoints, const int elementType, const std::pair< int, int > &entity) const override
    {
      _isConstantLeft = _left.isConstant(entity);
      _isConstantRight = _right.isConstant(entity);

      if(_isConstantLeft) {
        _valuesLeft.resize(1);
        _left.evaluate(_valuesLeft[0], 0., 0., 0., entity);
      }
      else {
#pragma omp parallel num_threads(omp::getMaxThreads())
        {
          _left.evaluate(_valuesLeft, points, gaussPoints, elementType, entity);
        }
      }

      if(_isConstantRight) {
        _valuesRight.resize(1);
        _right.evaluate(_valuesRight[0], 0., 0., 0., entity);
      }
      else {
#pragma omp parallel num_threads(omp::getMaxThreads())
        {
          _right.evaluate(_valuesRight, points, gaussPoints, elementType, entity);
        }
      }
    }

    void operator()(Eigen::Ref< Eigen::MatrixX< T_Scalar > > expression, const Eigen::Ref< const Eigen::MatrixX< scalar::Precision< T_Scalar > > > field, const unsigned int offset) const override
    {
      expression = crossT< T_Scalar >(_valuesLeft[_isConstantLeft ? 0 : offset], cross< T_Scalar >(field, _valuesRight[_isConstantRight ? 0 : offset]));
    }

    unsigned int size1() const override
    {
      return 3;
    }

    unsigned int size2() const override
    {
      return 0; // means number of dofs by element
    }

    void clear() const override
    {
      _valuesLeft.clear();
      _valuesRight.clear();
      _isConstantLeft = false;
      _isConstantRight = false;
    }

    common::Memory memory() const override
    {
      return common::Memory(_valuesLeft.size() * sizeof(typename MathObject< T_Scalar, Degree::Degree1 >::Object) + _valuesRight.size() * sizeof(typename MathObject< T_Scalar, Degree::Degree1 >::Object));
    }
  };

  // ############################
  // Compound Field
  // ############################

  template< class T_Scalar, field::Form T_Form, Degree T_DegreeLhs, equation::Product T_ProductLhs, Degree T_DegreeRhs, equation::Product T_ProductRhs, EvaluationOrder T_EvaluationOrder > // T_EvaluationOrder = Left -> (T_DegreeLhs * dof) * T_DegreeRhs, = Right - > T_DegreeLhs * (dof * T_DegreeRhs)
  class EquationEvaluatorCompound : public EquationEvaluatorInterface< T_Scalar, T_Form >
  {
   private:
    function::Function< T_Scalar, T_DegreeLhs > _left;
    function::Function< T_Scalar, T_DegreeLhs > _right;

   public:
    EquationEvaluatorCompound(const function::Function< T_Scalar, T_DegreeLhs > &left, const function::Function< T_Scalar, T_DegreeLhs > &right) :
      _left(left), _right(right) {}
    ~EquationEvaluatorCompound() {}

    virtual bool isConstant(const std::pair< int, int > &entity) const = 0;

    virtual void initialize(const std::vector< scalar::Precision< T_Scalar >, numa::allocator< scalar::Precision< T_Scalar > > > &points, const std::vector< scalar::Precision< T_Scalar > > &gaussPoints, const int elementType, const std::pair< int, int > &entity) const = 0;
    virtual void operator()(Eigen::Ref< Eigen::MatrixX< T_Scalar > > expression, const Eigen::Ref< const Eigen::MatrixX< scalar::Precision< T_Scalar > > > field, const unsigned int offset) const {}

    virtual unsigned int size1() const = 0;
    virtual unsigned int size2() const = 0;
    virtual void clear() const = 0;
    virtual common::Memory memory() const = 0;
  };

  //
  // Left
  //

  template< class T_Scalar, field::Form T_Form, Degree T_DegreeLhs >
  class EquationEvaluatorCompound< T_Scalar, T_Form, T_DegreeLhs, equation::Product::ScalarProduct, Degree::Empty, equation::Product::Empty, EvaluationOrder::Left > final : public EquationEvaluatorInterface< T_Scalar, T_Form >
  {
   private:
    function::Function< T_Scalar, T_DegreeLhs > _left;
    mutable std::vector< typename MathObject< T_Scalar, T_DegreeLhs >::Object, numa::allocator< typename MathObject< T_Scalar, T_DegreeLhs >::Object > > _values;
    mutable bool _isConstant;

   public:
    EquationEvaluatorCompound(const function::Function< T_Scalar, T_DegreeLhs > &left) :
      _left(left), _values(), _isConstant(false)
    {
    }
    ~EquationEvaluatorCompound()
    {
    }

    bool isConstant(const std::pair< int, int > &entity) const override
    {
      return _left.isConstant(entity);
    }

    void initialize(const std::vector< scalar::Precision< T_Scalar >, numa::allocator< scalar::Precision< T_Scalar > > > &points, const std::vector< scalar::Precision< T_Scalar > > &gaussPoints, const int elementType, const std::pair< int, int > &entity) const override
    {
      _isConstant = isConstant(entity);
      if(_isConstant) {
        _values.resize(1);
        _left.evaluate(_values[0], 0., 0., 0., entity);
      }
      else {
#pragma omp parallel num_threads(omp::getMaxThreads())
        {
          _left.evaluate(_values, points, gaussPoints, elementType, entity);
        }
      }
    }

    void operator()(Eigen::Ref< Eigen::MatrixX< T_Scalar > > expression, const Eigen::Ref< const Eigen::MatrixX< scalar::Precision< T_Scalar > > > field, const unsigned int offset) const override
    {
      expression = _values[_isConstant ? 0 : offset] * field;
    }

    unsigned int size1() const override
    {
      return 4; // means field size
    }

    unsigned int size2() const override
    {
      return 0; // means number of dofs by element
    }

    void clear() const override
    {
      _values.clear();
      _isConstant = false;
    }

    common::Memory memory() const override
    {
      return common::Memory(_values.size() * sizeof(typename MathObject< T_Scalar, T_DegreeLhs >::Object));
    }
  };

  template< class T_Scalar >
  class EquationEvaluatorCompound< T_Scalar, field::Form::Form0, Degree::Degree1, equation::Product::ScalarProduct, Degree::Empty, equation::Product::Empty, EvaluationOrder::Left > final : public EquationEvaluatorInterface< T_Scalar, field::Form::Form0 >
  {
   private:
    function::Function< T_Scalar, Degree::Degree1 > _left;
    mutable std::vector< typename MathObject< T_Scalar, Degree::Degree1 >::Object, numa::allocator< typename MathObject< T_Scalar, Degree::Degree1 >::Object > > _values;
    mutable bool _isConstant;

   public:
    EquationEvaluatorCompound(const function::Function< T_Scalar, Degree::Degree1 > &left) :
      _left(left), _values(), _isConstant(false)
    {
    }
    ~EquationEvaluatorCompound()
    {
    }

    bool isConstant(const std::pair< int, int > &entity) const override
    {
      return _left.isConstant(entity);
    }

    void initialize(const std::vector< scalar::Precision< T_Scalar >, numa::allocator< scalar::Precision< T_Scalar > > > &points, const std::vector< scalar::Precision< T_Scalar > > &gaussPoints, const int elementType, const std::pair< int, int > &entity) const override
    {
      _isConstant = isConstant(entity);
      if(_isConstant) {
        _values.resize(1);
        _left.evaluate(_values[0], 0., 0., 0., entity);
      }
      else {
#pragma omp parallel num_threads(omp::getMaxThreads())
        {
          _left.evaluate(_values, points, gaussPoints, elementType, entity);
        }
      }
    }

    void operator()(Eigen::Ref< Eigen::MatrixX< T_Scalar > > expression, const Eigen::Ref< const Eigen::MatrixX< scalar::Precision< T_Scalar > > > field, const unsigned int offset) const override
    {
      expression = _values[_isConstant ? 0 : offset].transpose() * field;
    }

    unsigned int size1() const override
    {
      return 1;
    }

    unsigned int size2() const override
    {
      return 0; // means number of dofs by element
    }

    void clear() const override
    {
      _values.clear();
      _isConstant = false;
    }

    common::Memory memory() const override
    {
      return common::Memory(_values.size() * sizeof(typename MathObject< T_Scalar, Degree::Degree1 >::Object));
    }
  };

  template< class T_Scalar >
  class EquationEvaluatorCompound< T_Scalar, field::Form::Form3, Degree::Degree1, equation::Product::ScalarProduct, Degree::Empty, equation::Product::Empty, EvaluationOrder::Left > final : public EquationEvaluatorInterface< T_Scalar, field::Form::Form3 >
  {
   private:
    function::Function< T_Scalar, Degree::Degree1 > _left;
    mutable std::vector< typename MathObject< T_Scalar, Degree::Degree1 >::Object, numa::allocator< typename MathObject< T_Scalar, Degree::Degree1 >::Object > > _values;
    mutable bool _isConstant;

   public:
    EquationEvaluatorCompound(const function::Function< T_Scalar, Degree::Degree1 > &left) :
      _left(left), _values(), _isConstant(false)
    {
    }
    ~EquationEvaluatorCompound()
    {
    }

    bool isConstant(const std::pair< int, int > &entity) const override
    {
      return _left.isConstant(entity);
    }

    void initialize(const std::vector< scalar::Precision< T_Scalar >, numa::allocator< scalar::Precision< T_Scalar > > > &points, const std::vector< scalar::Precision< T_Scalar > > &gaussPoints, const int elementType, const std::pair< int, int > &entity) const override
    {
      _isConstant = isConstant(entity);
      if(_isConstant) {
        _values.resize(1);
        _left.evaluate(_values[0], 0., 0., 0., entity);
      }
      else {
#pragma omp parallel num_threads(omp::getMaxThreads())
        {
          _left.evaluate(_values, points, gaussPoints, elementType, entity);
        }
      }
    }

    void operator()(Eigen::Ref< Eigen::MatrixX< T_Scalar > > expression, const Eigen::Ref< const Eigen::MatrixX< scalar::Precision< T_Scalar > > > field, const unsigned int offset) const override
    {
      expression = _values[_isConstant ? 0 : offset].transpose() * field;
    }

    unsigned int size1() const override
    {
      return 1;
    }

    unsigned int size2() const override
    {
      return 0; // means number of dofs by element
    }

    void clear() const override
    {
      _values.clear();
      _isConstant = false;
    }

    common::Memory memory() const override
    {
      return common::Memory(_values.size() * sizeof(typename MathObject< T_Scalar, Degree::Degree1 >::Object));
    }
  };

  template< class T_Scalar, field::Form T_Form >
  class EquationEvaluatorCompound< T_Scalar, T_Form, Degree::Degree2, equation::Product::ScalarProduct, Degree::Empty, equation::Product::Empty, EvaluationOrder::Left > final : public EquationEvaluatorInterface< T_Scalar, T_Form >
  {
   private:
    function::Function< T_Scalar, Degree::Degree2 > _left;
    mutable std::vector< typename MathObject< T_Scalar, Degree::Degree2 >::Object, numa::allocator< typename MathObject< T_Scalar, Degree::Degree2 >::Object > > _values;
    mutable bool _isConstant;

   public:
    EquationEvaluatorCompound(const function::Function< T_Scalar, Degree::Degree2 > &left) :
      _left(left), _values(), _isConstant(false)
    {
    }
    ~EquationEvaluatorCompound()
    {
    }

    bool isConstant(const std::pair< int, int > &entity) const override
    {
      return _left.isConstant(entity);
    }

    void initialize(const std::vector< scalar::Precision< T_Scalar >, numa::allocator< scalar::Precision< T_Scalar > > > &points, const std::vector< scalar::Precision< T_Scalar > > &gaussPoints, const int elementType, const std::pair< int, int > &entity) const override
    {
      _isConstant = isConstant(entity);
      if(_isConstant) {
        _values.resize(1);
        _left.evaluate(_values[0], 0., 0., 0., entity);
      }
      else {
#pragma omp parallel num_threads(omp::getMaxThreads())
        {
          _left.evaluate(_values, points, gaussPoints, elementType, entity);
        }
      }
    }

    void operator()(Eigen::Ref< Eigen::MatrixX< T_Scalar > > expression, const Eigen::Ref< const Eigen::MatrixX< scalar::Precision< T_Scalar > > > field, const unsigned int offset) const override
    {
      for(unsigned int d = 0; d < field.cols(); ++d) {
        expression.col(d) = Eigen::Map< const Eigen::VectorX< T_Scalar > >((Eigen::Map< const Eigen::Matrix3< scalar::Precision< T_Scalar > > >(field.col(d).data(), 3, 3) * _values[_isConstant ? 0 : offset]).eval().data(), 9, 1);
      }
    }

    unsigned int size1() const override
    {
      return 4; // means field size
    }

    unsigned int size2() const override
    {
      return 0; // means number of dofs by element
    }

    void clear() const override
    {
      _values.clear();
      _isConstant = false;
    }

    common::Memory memory() const override
    {
      return common::Memory(_values.size() * sizeof(typename MathObject< T_Scalar, Degree::Degree2 >::Object));
    }
  };

  template< class T_Scalar >
  class EquationEvaluatorCompound< T_Scalar, field::Form::Form0, Degree::Degree2, equation::Product::ScalarProduct, Degree::Empty, equation::Product::Empty, EvaluationOrder::Left > final : public EquationEvaluatorInterface< T_Scalar, field::Form::Form0 >
  {
   private:
    function::Function< T_Scalar, Degree::Degree2 > _left;
    mutable std::vector< typename MathObject< T_Scalar, Degree::Degree2 >::Object, numa::allocator< typename MathObject< T_Scalar, Degree::Degree2 >::Object > > _values;
    mutable bool _isConstant;

   public:
    EquationEvaluatorCompound(const function::Function< T_Scalar, Degree::Degree2 > &left) :
      _left(left), _values(), _isConstant(false)
    {
    }
    ~EquationEvaluatorCompound()
    {
    }

    bool isConstant(const std::pair< int, int > &entity) const override
    {
      return _left.isConstant(entity);
    }

    void initialize(const std::vector< scalar::Precision< T_Scalar >, numa::allocator< scalar::Precision< T_Scalar > > > &points, const std::vector< scalar::Precision< T_Scalar > > &gaussPoints, const int elementType, const std::pair< int, int > &entity) const override
    {
      _isConstant = isConstant(entity);
      if(_isConstant) {
        _values.resize(1);
        _left.evaluate(_values[0], 0., 0., 0., entity);
      }
      else {
#pragma omp parallel num_threads(omp::getMaxThreads())
        {
          _left.evaluate(_values, points, gaussPoints, elementType, entity);
        }
      }
    }

    void operator()(Eigen::Ref< Eigen::MatrixX< T_Scalar > > expression, const Eigen::Ref< const Eigen::MatrixX< scalar::Precision< T_Scalar > > > field, const unsigned int offset) const override
    {
      for(unsigned int d = 0; d < field.cols(); ++d) {
        expression.col(d) = Eigen::Map< const Eigen::VectorX< T_Scalar > >((_values[_isConstant ? 0 : offset] * Eigen::Map< const Eigen::Vector3< scalar::Precision< T_Scalar > > >(field.col(d).data(), 3, 1)).eval().data(), 3, 1);
      }
    }

    unsigned int size1() const override
    {
      return 4; // means field size
    }

    unsigned int size2() const override
    {
      return 0; // means number of dofs by element
    }

    void clear() const override
    {
      _values.clear();
      _isConstant = false;
    }

    common::Memory memory() const override
    {
      return common::Memory(_values.size() * sizeof(typename MathObject< T_Scalar, Degree::Degree2 >::Object));
    }
  };

  template< class T_Scalar >
  class EquationEvaluatorCompound< T_Scalar, field::Form::Form3, Degree::Degree2, equation::Product::ScalarProduct, Degree::Empty, equation::Product::Empty, EvaluationOrder::Left > final : public EquationEvaluatorInterface< T_Scalar, field::Form::Form3 >
  {
   private:
    function::Function< T_Scalar, Degree::Degree2 > _left;
    mutable std::vector< typename MathObject< T_Scalar, Degree::Degree2 >::Object, numa::allocator< typename MathObject< T_Scalar, Degree::Degree2 >::Object > > _values;
    mutable bool _isConstant;

   public:
    EquationEvaluatorCompound(const function::Function< T_Scalar, Degree::Degree2 > &left) :
      _left(left), _values(), _isConstant(false)
    {
    }
    ~EquationEvaluatorCompound()
    {
    }

    bool isConstant(const std::pair< int, int > &entity) const override
    {
      return _left.isConstant(entity);
    }

    void initialize(const std::vector< scalar::Precision< T_Scalar >, numa::allocator< scalar::Precision< T_Scalar > > > &points, const std::vector< scalar::Precision< T_Scalar > > &gaussPoints, const int elementType, const std::pair< int, int > &entity) const override
    {
      _isConstant = isConstant(entity);
      if(_isConstant) {
        _values.resize(1);
        _left.evaluate(_values[0], 0., 0., 0., entity);
      }
      else {
#pragma omp parallel num_threads(omp::getMaxThreads())
        {
          _left.evaluate(_values, points, gaussPoints, elementType, entity);
        }
      }
    }

    void operator()(Eigen::Ref< Eigen::MatrixX< T_Scalar > > expression, const Eigen::Ref< const Eigen::MatrixX< scalar::Precision< T_Scalar > > > field, const unsigned int offset) const override
    {
      for(unsigned int d = 0; d < field.cols(); ++d) {
        expression.col(d) = Eigen::Map< const Eigen::VectorX< T_Scalar > >((_values[_isConstant ? 0 : offset] * Eigen::Map< const Eigen::VectorX< scalar::Precision< T_Scalar > > >(field.col(d).data(), 3, 1)).eval().data(), 3, 1);
      }
    }

    unsigned int size1() const override
    {
      return 4; // means field size
    }

    unsigned int size2() const override
    {
      return 0; // means number of dofs by element
    }

    void clear() const override
    {
      _values.clear();
      _isConstant = false;
    }

    common::Memory memory() const override
    {
      return common::Memory(_values.size() * sizeof(typename MathObject< T_Scalar, Degree::Degree2 >::Object));
    }
  };

  template< class T_Scalar, field::Form T_Form >
  class EquationEvaluatorCompound< T_Scalar, T_Form, Degree::Degree4, equation::Product::ScalarProduct, Degree::Empty, equation::Product::Empty, EvaluationOrder::Left > final : public EquationEvaluatorInterface< T_Scalar, T_Form >
  {
   private:
    function::Function< T_Scalar, Degree::Degree4 > _left;
    mutable std::vector< typename MathObject< T_Scalar, Degree::Degree4 >::Object, numa::allocator< typename MathObject< T_Scalar, Degree::Degree4 >::Object > > _values;
    mutable bool _isConstant;

   public:
    EquationEvaluatorCompound(const function::Function< T_Scalar, Degree::Degree4 > &left) :
      _left(left), _values(), _isConstant(false)
    {
    }
    ~EquationEvaluatorCompound()
    {
    }

    bool isConstant(const std::pair< int, int > &entity) const override
    {
      return _left.isConstant(entity);
    }

    void initialize(const std::vector< scalar::Precision< T_Scalar >, numa::allocator< scalar::Precision< T_Scalar > > > &points, const std::vector< scalar::Precision< T_Scalar > > &gaussPoints, const int elementType, const std::pair< int, int > &entity) const override
    {
      _isConstant = isConstant(entity);
      if(_isConstant) {
        _values.resize(1);
        _left.evaluate(_values[0], 0., 0., 0., entity);
      }
      else {
#pragma omp parallel num_threads(omp::getMaxThreads())
        {
          _left.evaluate(_values, points, gaussPoints, elementType, entity);
        }
      }
    }

    void operator()(Eigen::Ref< Eigen::MatrixX< T_Scalar > > expression, const Eigen::Ref< const Eigen::MatrixX< scalar::Precision< T_Scalar > > > field, const unsigned int offset) const override
    {
      // C_{ijkl} : u_{kl} -> v_{ij}
      for(unsigned int d = 0; d < field.cols(); ++d) {
        const Eigen::Matrix3< scalar::Precision< T_Scalar > > mat = Eigen::Map< const Eigen::Matrix3< scalar::Precision< T_Scalar > > >(field.col(d).data(), 3, 3).transpose();
        for(unsigned int j = 0; j < 3; ++j) {
          for(unsigned int i = 0; i < 3; ++i) {
            expression(i + j * 3, d) = (_values[_isConstant ? 0 : offset](i, j) * mat).trace();
          }
        }
      }
    }

    unsigned int size1() const override
    {
      return 4; // means field size
    }

    unsigned int size2() const override
    {
      return 0; // means number of dofs by element
    }

    void clear() const override
    {
      _values.clear();
      _isConstant = false;
    }

    common::Memory memory() const override
    {
      return common::Memory(_values.size() * sizeof(typename MathObject< T_Scalar, Degree::Degree4 >::Object));
    }
  };

  //
  // Right
  //

  template< class T_Scalar, field::Form T_Form, Degree T_DegreeRhs >
  class EquationEvaluatorCompound< T_Scalar, T_Form, Degree::Empty, equation::Product::Empty, T_DegreeRhs, equation::Product::ScalarProduct, EvaluationOrder::Right > final : public EquationEvaluatorInterface< T_Scalar, T_Form >
  {
   private:
    function::Function< T_Scalar, T_DegreeRhs > _right;
    mutable std::vector< typename MathObject< T_Scalar, T_DegreeRhs >::Object, numa::allocator< typename MathObject< T_Scalar, T_DegreeRhs >::Object > > _values;
    mutable bool _isConstant;

   public:
    EquationEvaluatorCompound(const function::Function< T_Scalar, T_DegreeRhs > &right) :
      _right(right), _values(), _isConstant(false)
    {
    }
    ~EquationEvaluatorCompound()
    {
    }

    bool isConstant(const std::pair< int, int > &entity) const override
    {
      return _right.isConstant(entity);
    }

    void initialize(const std::vector< scalar::Precision< T_Scalar >, numa::allocator< scalar::Precision< T_Scalar > > > &points, const std::vector< scalar::Precision< T_Scalar > > &gaussPoints, const int elementType, const std::pair< int, int > &entity) const override
    {
      _isConstant = isConstant(entity);
      if(_isConstant) {
        _values.resize(1);
        _right.evaluate(_values[0], 0., 0., 0., entity);
      }
      else {
#pragma omp parallel num_threads(omp::getMaxThreads())
        {
          _right.evaluate(_values, points, gaussPoints, elementType, entity);
        }
      }
    }

    void operator()(Eigen::Ref< Eigen::MatrixX< T_Scalar > > expression, const Eigen::Ref< const Eigen::MatrixX< scalar::Precision< T_Scalar > > > field, const unsigned int offset) const override
    {
      expression = field.transpose() * _values[_isConstant ? 0 : offset];
    }

    unsigned int size1() const override
    {
      return 0; // means number of dofs by element
    }

    unsigned int size2() const override
    {
      return 4; // means field size
    }

    void clear() const override
    {
      _values.clear();
      _isConstant = false;
    }

    common::Memory memory() const override
    {
      return common::Memory(_values.size() * sizeof(typename MathObject< T_Scalar, T_DegreeRhs >::Object));
    }
  };

  template< class T_Scalar >
  class EquationEvaluatorCompound< T_Scalar, field::Form::Form0, Degree::Empty, equation::Product::Empty, Degree::Degree1, equation::Product::ScalarProduct, EvaluationOrder::Right > final : public EquationEvaluatorInterface< T_Scalar, field::Form::Form0 >
  {
   private:
    function::Function< T_Scalar, Degree::Degree1 > _right;
    mutable std::vector< typename MathObject< T_Scalar, Degree::Degree1 >::Object, numa::allocator< typename MathObject< T_Scalar, Degree::Degree1 >::Object > > _values;
    mutable bool _isConstant;

   public:
    EquationEvaluatorCompound(const function::Function< T_Scalar, Degree::Degree1 > &right) :
      _right(right), _values(), _isConstant(false)
    {
    }
    ~EquationEvaluatorCompound()
    {
    }

    bool isConstant(const std::pair< int, int > &entity) const override
    {
      return _right.isConstant(entity);
    }

    void initialize(const std::vector< scalar::Precision< T_Scalar >, numa::allocator< scalar::Precision< T_Scalar > > > &points, const std::vector< scalar::Precision< T_Scalar > > &gaussPoints, const int elementType, const std::pair< int, int > &entity) const override
    {
      _isConstant = isConstant(entity);
      if(_isConstant) {
        _values.resize(1);
        _right.evaluate(_values[0], 0., 0., 0., entity);
      }
      else {
#pragma omp parallel num_threads(omp::getMaxThreads())
        {
          _right.evaluate(_values, points, gaussPoints, elementType, entity);
        }
      }
    }

    void operator()(Eigen::Ref< Eigen::MatrixX< T_Scalar > > expression, const Eigen::Ref< const Eigen::MatrixX< scalar::Precision< T_Scalar > > > field, const unsigned int offset) const override
    {
      expression = field.transpose() * _values[_isConstant ? 0 : offset];
    }

    unsigned int size1() const override
    {
      return 0; // means number of dofs by element
    }

    unsigned int size2() const override
    {
      return 1;
    }

    void clear() const override
    {
      _values.clear();
      _isConstant = false;
    }

    common::Memory memory() const override
    {
      return common::Memory(_values.size() * sizeof(typename MathObject< T_Scalar, Degree::Degree1 >::Object));
    }
  };

  template< class T_Scalar >
  class EquationEvaluatorCompound< T_Scalar, field::Form::Form3, Degree::Empty, equation::Product::Empty, Degree::Degree1, equation::Product::ScalarProduct, EvaluationOrder::Right > final : public EquationEvaluatorInterface< T_Scalar, field::Form::Form3 >
  {
   private:
    function::Function< T_Scalar, Degree::Degree1 > _right;
    mutable std::vector< typename MathObject< T_Scalar, Degree::Degree1 >::Object, numa::allocator< typename MathObject< T_Scalar, Degree::Degree1 >::Object > > _values;
    mutable bool _isConstant;

   public:
    EquationEvaluatorCompound(const function::Function< T_Scalar, Degree::Degree1 > &right) :
      _right(right), _values(), _isConstant(false)
    {
    }
    ~EquationEvaluatorCompound()
    {
    }

    bool isConstant(const std::pair< int, int > &entity) const override
    {
      return _right.isConstant(entity);
    }

    void initialize(const std::vector< scalar::Precision< T_Scalar >, numa::allocator< scalar::Precision< T_Scalar > > > &points, const std::vector< scalar::Precision< T_Scalar > > &gaussPoints, const int elementType, const std::pair< int, int > &entity) const override
    {
      _isConstant = isConstant(entity);
      if(_isConstant) {
        _values.resize(1);
        _right.evaluate(_values[0], 0., 0., 0., entity);
      }
      else {
#pragma omp parallel num_threads(omp::getMaxThreads())
        {
          _right.evaluate(_values, points, gaussPoints, elementType, entity);
        }
      }
    }

    void operator()(Eigen::Ref< Eigen::MatrixX< T_Scalar > > expression, const Eigen::Ref< const Eigen::MatrixX< scalar::Precision< T_Scalar > > > field, const unsigned int offset) const override
    {
      expression = field.transpose() * _values[_isConstant ? 0 : offset];
    }

    unsigned int size1() const override
    {
      return 0; // means number of dofs by element
    }

    unsigned int size2() const override
    {
      return 1;
    }

    void clear() const override
    {
      _values.clear();
      _isConstant = false;
    }

    common::Memory memory() const override
    {
      return common::Memory(_values.size() * sizeof(typename MathObject< T_Scalar, Degree::Degree1 >::Object));
    }
  };


} // namespace gmshfem::term::evaluator


#endif // H_GMSHFEM_EQUATIONEVALUATOR
