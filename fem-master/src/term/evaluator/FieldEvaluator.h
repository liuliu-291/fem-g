// GmshFEM - Copyright (C) 2019-2022, A. Royer, E. Béchet, C. Geuzaine, Université de Liège
//
// See the LICENSE.txt file for license information. Please report all
// issues on https://gitlab.onelab.info/gmsh/fem/issues

#ifndef H_GMSHFEM_FIELDEVALUATOR
#define H_GMSHFEM_FIELDEVALUATOR

#include "FieldObject.h"
#include "scalar.h"

#include <iostream>


namespace gmshfem::term::evaluator
{


  //
  // Field Evaluator
  //

  template< class T_Scalar, field::Form T_Form >
  class FieldEvaluator
  {
   public:
    FieldEvaluator() {}
    virtual ~FieldEvaluator() {}

    virtual void operator()(Eigen::Ref< Eigen::MatrixX< scalar::Precision< T_Scalar > > > fieldExpression, const scalar::Precision< T_Scalar > *const bf, const scalar::Precision< T_Scalar > *const jacobians, const scalar::Precision< T_Scalar > *const determinants, const unsigned int nbrDofs) const = 0;
    virtual bool needJacobians() const = 0;
    virtual unsigned int size() const = 0;
  };

  // Matrix : 1 x nDof
  template< class T_Scalar >
  class FieldEvaluator< T_Scalar, field::Form::Form0 >
  {
   public:
    FieldEvaluator() {}
    virtual ~FieldEvaluator() {}

    virtual void operator()(Eigen::Ref< Eigen::MatrixX< scalar::Precision< T_Scalar > > > fieldExpression, const scalar::Precision< T_Scalar > *const bf, const scalar::Precision< T_Scalar > *const jacobians, const scalar::Precision< T_Scalar > *const determinants, const unsigned int nbrDofs) const
    {
      fieldExpression = Eigen::Map< const Eigen::MatrixX< scalar::Precision< T_Scalar > > >(bf, 1, nbrDofs);
    }

    virtual bool needJacobians() const
    {
      return false;
    }

    virtual unsigned int size() const
    {
      return 1;
    }
  };

  // Matrix : 3 x nDof
  template< class T_Scalar >
  class FieldEvaluator< T_Scalar, field::Form::Form1 >
  {
   public:
    FieldEvaluator() {}
    virtual ~FieldEvaluator() {}

    virtual void operator()(Eigen::Ref< Eigen::MatrixX< scalar::Precision< T_Scalar > > > fieldExpression, const scalar::Precision< T_Scalar > *const bf, const scalar::Precision< T_Scalar > *const jacobians, const scalar::Precision< T_Scalar > *const determinants, const unsigned int nbrDofs) const
    {
      fieldExpression = Eigen::Map< const Eigen::Matrix3< scalar::Precision< T_Scalar > > >(jacobians).transpose().inverse() * Eigen::Map< const Eigen::MatrixX< scalar::Precision< T_Scalar > > >(bf, 3, nbrDofs);
    }

    virtual bool needJacobians() const
    {
      return true;
    }

    virtual unsigned int size() const
    {
      return 3;
    }
  };

  // Matrix : 3 x nDof
  template< class T_Scalar >
  class FieldEvaluator< T_Scalar, field::Form::Form2 >
  {
   public:
    FieldEvaluator() {}
    virtual ~FieldEvaluator() {}

    virtual void operator()(Eigen::Ref< Eigen::MatrixX< scalar::Precision< T_Scalar > > > fieldExpression, const scalar::Precision< T_Scalar > *const bf, const scalar::Precision< T_Scalar > *const jacobians, const scalar::Precision< T_Scalar > *const determinants, const unsigned int nbrDofs) const
    {
      fieldExpression = 1. / (*determinants) * (Eigen::Map< const Eigen::Matrix3< scalar::Precision< T_Scalar > > >(jacobians) * Eigen::Map< const Eigen::MatrixX< scalar::Precision< T_Scalar > > >(bf, 3, nbrDofs));
    }

    virtual bool needJacobians() const
    {
      return true;
    }

    virtual unsigned int size() const
    {
      return 3;
    }
  };

  // Matrix : 1 x nDof
  template< class T_Scalar >
  class FieldEvaluator< T_Scalar, field::Form::Form3 >
  {
   public:
    FieldEvaluator() {}
    virtual ~FieldEvaluator() {}

    virtual void operator()(Eigen::Ref< Eigen::MatrixX< scalar::Precision< T_Scalar > > > fieldExpression, const scalar::Precision< T_Scalar > *const bf, const scalar::Precision< T_Scalar > *const jacobians, const scalar::Precision< T_Scalar > *const determinants, const unsigned int nbrDofs) const
    {
      fieldExpression = 1. / (*determinants) * Eigen::Map< const Eigen::MatrixX< scalar::Precision< T_Scalar > > >(bf, 1, nbrDofs);
    }

    virtual bool needJacobians() const
    {
      return true;
    }

    virtual unsigned int size() const
    {
      return 1;
    }
  };


  //
  // CompoundField Evaluator
  //

  template< class T_Scalar, field::Form T_Form, unsigned int T_NumFields >
  class CompoundFieldEvaluator : public FieldEvaluator< T_Scalar, T_Form >
  {
  };

  // Matrix : 3 x nDof
  template< class T_Scalar, unsigned int T_NumFields >
  class CompoundFieldEvaluator< T_Scalar, field::Form::Form0, T_NumFields > : public FieldEvaluator< T_Scalar, field::Form::Form0 >
  {
   public:
    CompoundFieldEvaluator() {}
    ~CompoundFieldEvaluator() {}

    virtual void operator()(Eigen::Ref< Eigen::MatrixX< scalar::Precision< T_Scalar > > > fieldExpression, const scalar::Precision< T_Scalar > *const bf, const scalar::Precision< T_Scalar > *const jacobians, const scalar::Precision< T_Scalar > *const determinants, const unsigned int nbrDofs) const override
    {
      const Eigen::Map< const Eigen::MatrixX< scalar::Precision< T_Scalar > > > tmp = Eigen::Map< const Eigen::MatrixX< scalar::Precision< T_Scalar > > >(bf, 1, nbrDofs);

      fieldExpression(0, Eigen::seq(Eigen::fix< 0 >, Eigen::last, Eigen::fix< T_NumFields >)) = tmp;
      fieldExpression(1, Eigen::seq(Eigen::fix< 1 >, Eigen::last, Eigen::fix< T_NumFields >)) = tmp;
      if constexpr(T_NumFields == 3) {
        fieldExpression(2, Eigen::seq(Eigen::fix< 2 >, Eigen::last, Eigen::fix< T_NumFields >)) = tmp;
      }
    }

    virtual bool needJacobians() const override
    {
      return false;
    }

    virtual unsigned int size() const override
    {
      return 3;
    }
  };

  // Matrix : 9 x nDof
  template< class T_Scalar, unsigned int T_NumFields >
  class CompoundFieldEvaluator< T_Scalar, field::Form::Form1, T_NumFields > : public FieldEvaluator< T_Scalar, field::Form::Form1 >
  {
   public:
    CompoundFieldEvaluator() {}
    virtual ~CompoundFieldEvaluator() {}

    void operator()(Eigen::Ref< Eigen::MatrixX< scalar::Precision< T_Scalar > > > fieldExpression, const scalar::Precision< T_Scalar > *const bf, const scalar::Precision< T_Scalar > *const jacobians, const scalar::Precision< T_Scalar > *const determinants, const unsigned int nbrDofs) const override
    {
      const Eigen::MatrixX< scalar::Precision< T_Scalar > > tmp = Eigen::Map< const Eigen::Matrix3< scalar::Precision< T_Scalar > > >(jacobians).transpose().inverse() * Eigen::Map< const Eigen::MatrixX< scalar::Precision< T_Scalar > > >(bf, 3, nbrDofs);

      fieldExpression(Eigen::seq(Eigen::fix< 0 >, Eigen::fix< 2 >), Eigen::seq(Eigen::fix< 0 >, Eigen::last, Eigen::fix< T_NumFields >)) = tmp;
      fieldExpression(Eigen::seq(Eigen::fix< 3 >, Eigen::fix< 5 >), Eigen::seq(Eigen::fix< 1 >, Eigen::last, Eigen::fix< T_NumFields >)) = tmp;
      if constexpr(T_NumFields == 3) {
        fieldExpression(Eigen::seq(Eigen::fix< 6 >, Eigen::fix< 8 >), Eigen::seq(Eigen::fix< 2 >, Eigen::last, Eigen::fix< T_NumFields >)) = tmp;
      }
    }

    bool needJacobians() const override
    {
      return true;
    }

    unsigned int size() const override
    {
      return 9;
    }
  };

  // Matrix : 9 x nDof
  template< class T_Scalar, unsigned int T_NumFields >
  class CompoundFieldEvaluator< T_Scalar, field::Form::Form2, T_NumFields > : public FieldEvaluator< T_Scalar, field::Form::Form2 >
  {
   public:
    CompoundFieldEvaluator() {}
    virtual ~CompoundFieldEvaluator() {}

    void operator()(Eigen::Ref< Eigen::MatrixX< scalar::Precision< T_Scalar > > > fieldExpression, const scalar::Precision< T_Scalar > *const bf, const scalar::Precision< T_Scalar > *const jacobians, const scalar::Precision< T_Scalar > *const determinants, const unsigned int nbrDofs) const override
    {
      const Eigen::MatrixX< scalar::Precision< T_Scalar > > tmp = 1. / (*determinants) * (Eigen::Map< const Eigen::Matrix3< scalar::Precision< T_Scalar > > >(jacobians) * Eigen::Map< const Eigen::MatrixX< scalar::Precision< T_Scalar > > >(bf, 3, nbrDofs));

      fieldExpression(Eigen::seq(Eigen::fix< 0 >, Eigen::fix< 2 >), Eigen::seq(Eigen::fix< 0 >, Eigen::last, Eigen::fix< T_NumFields >)) = tmp;
      fieldExpression(Eigen::seq(Eigen::fix< 3 >, Eigen::fix< 5 >), Eigen::seq(Eigen::fix< 1 >, Eigen::last, Eigen::fix< T_NumFields >)) = tmp;
      if constexpr(T_NumFields == 3) {
        fieldExpression(Eigen::seq(Eigen::fix< 6 >, Eigen::fix< 8 >), Eigen::seq(Eigen::fix< 2 >, Eigen::last, Eigen::fix< T_NumFields >)) = tmp;
      }
    }

    bool needJacobians() const override
    {
      return true;
    }

    unsigned int size() const override
    {
      return 9;
    }
  };

  // Matrix : 3 x nDof
  template< class T_Scalar, unsigned int T_NumFields >
  class CompoundFieldEvaluator< T_Scalar, field::Form::Form3, T_NumFields > : public FieldEvaluator< T_Scalar, field::Form::Form3 >
  {
   public:
    CompoundFieldEvaluator() {}
    virtual ~CompoundFieldEvaluator() {}

    virtual void operator()(Eigen::Ref< Eigen::MatrixX< scalar::Precision< T_Scalar > > > fieldExpression, const scalar::Precision< T_Scalar > *const bf, const scalar::Precision< T_Scalar > *const jacobians, const scalar::Precision< T_Scalar > *const determinants, const unsigned int nbrDofs) const override
    {
      const Eigen::MatrixX< scalar::Precision< T_Scalar > > tmp = 1. / (*determinants) * Eigen::Map< const Eigen::MatrixX< scalar::Precision< T_Scalar > > >(bf, 1, nbrDofs);

      fieldExpression(0, Eigen::seq(Eigen::fix< 0 >, Eigen::last, Eigen::fix< T_NumFields >)) = tmp;
      fieldExpression(1, Eigen::seq(Eigen::fix< 1 >, Eigen::last, Eigen::fix< T_NumFields >)) = tmp;
      if constexpr(T_NumFields == 3) {
        fieldExpression(2, Eigen::seq(Eigen::fix< 2 >, Eigen::last, Eigen::fix< T_NumFields >)) = tmp;
      }
    }

    virtual bool needJacobians() const override
    {
      return true;
    }

    virtual unsigned int size() const override
    {
      return 3;
    }
  };


} // namespace gmshfem::term::evaluator


#endif // H_GMSHFEM_FIELDEVALUATOR
