// GmshFEM - Copyright (C) 2019-2022, A. Royer, E. Béchet, C. Geuzaine, Université de Liège
//
// See the LICENSE.txt file for license information. Please report all
// issues on https://gitlab.onelab.info/gmsh/fem/issues

#ifndef H_GMSHFEM_FIELDNODE
#define H_GMSHFEM_FIELDNODE

#include "ExecutionTree.h"
#include "FieldInterface.h"
#include "IndiceBucket.h"
#include "OperationsInterface.h"

#include <complex>
#include <gmsh.h>
#include <vector>

namespace gmshfem::function
{


  //
  // FieldNode with it operation class
  //

  template< class T_Scalar, Degree T_Degree >
  class FieldOperation : public OperationsInterface< T_Scalar, T_Degree >
  {
   public:
    typedef T_Scalar D_Scalar;
    constexpr static Degree D_Degree = T_Degree;

    FieldOperation() :
      OperationsInterface< T_Scalar, T_Degree >()
    {
    }

    FieldOperation(const FieldOperation &other) :
      OperationsInterface< T_Scalar, T_Degree >(other)
    {
    }
  };

  template< class T_Op >
  class FieldNode final : public ExecutionTreeWithDegreeAndLeaves< typename T_Op::D_Scalar, T_Op::D_Degree >
  {
   protected:
    T_Op _op;
    const field::FieldInterface< typename T_Op::D_Scalar > *_field;
    const unsigned int _fieldTag;
    const std::string _fieldName;

    FieldNode(const FieldNode &other) :
      ExecutionTreeWithDegreeAndLeaves< typename T_Op::D_Scalar, T_Op::D_Degree >(other), _op(other._op), _field(other._field), _fieldTag(other._fieldTag), _fieldName(other._fieldName)
    {
    }

   public:
    FieldNode(const field::FieldInterface< typename T_Op::D_Scalar > *const field) :
      ExecutionTreeWithDegreeAndLeaves< typename T_Op::D_Scalar, T_Op::D_Degree >(), _op(), _field(field), _fieldTag(field->tag()), _fieldName(field->name())
    {
    }

    FieldNode(const T_Op &other, const field::FieldInterface< typename T_Op::D_Scalar > *const field) :
      ExecutionTreeWithDegreeAndLeaves< typename T_Op::D_Scalar, T_Op::D_Degree >(), _op(other), _field(field), _fieldTag(field->tag()), _fieldName(field->name())
    {
    }

    virtual ~FieldNode() {}

    virtual NodeType nodeType() const override
    {
      return NodeType::Field;
    }

    virtual bool isConstant(const std::pair< int, int > &entity) const override
    {
      return false;
    }

    virtual void evaluate(OutputVector< typename T_Op::D_Scalar, T_Op::D_Degree > &values, const std::vector< scalar::Precision< typename T_Op::D_Scalar >, numa::allocator< scalar::Precision< typename T_Op::D_Scalar > > > &points, const std::vector< scalar::Precision< typename T_Op::D_Scalar > > &gaussPoints, const int elementType, const std::pair< int, int > &entity) const override
    {
      if(values.size() < points.size() / 3) {
#pragma omp barrier
        FunctionAllocator< typename MathObject< typename T_Op::D_Scalar, T_Op::D_Degree >::Object >::initalizeNUMAMemory(points.size() / 3);
#pragma omp master
        values.resize(points.size() / 3);
#pragma omp barrier
      }

      _op(_field, values, points, gaussPoints, elementType, entity);
#pragma omp barrier
    }

    virtual FieldNode< T_Op > *copy() const override
    {
      return new FieldNode< T_Op >(*this);
    }

    const T_Op &operation() const
    {
      return _op;
    }

    virtual void exportTree(const domain::Domain &domain, common::HTMLNode *node) const override
    {
      node->name = _op.name() + " " + _field->name();
      node->scalar = scalar::Name< typename T_Op::D_Scalar >::name;
      node->degree = T_Op::D_Degree;
      node->nodeType = this->nodeType();
    }

    virtual common::Memory peakMemoryByGaussPoint(const domain::Domain &domain) const override
    {
      return common::Memory(0);
    }

    virtual bool operator==(const ExecutionTreeInterface &other) const override
    {
      if(this->tag() == other.tag()) {
        return true;
      }

      const std::type_info &ext1 = typeid(*this);
      const std::type_info &ext2 = typeid(other);
      if(ext1.hash_code() == ext2.hash_code()) {
        if(this->_field == static_cast< const FieldNode< T_Op > & >(other)._field) {
          const std::type_info &op1 = typeid(this->_op);
          const std::type_info &op2 = typeid(static_cast< const FieldNode< T_Op > & >(other)._op);
          if(op1.hash_code() == op2.hash_code()) {
            return true;
          }
        }
      }

      return false;
    }
  };

  //
  // FieldScalarTypeNode with it operation class
  //

  template< class T_ScalarOut, class T_ScalarIn, Degree T_Degree >
  class FieldScalarTypeOperation : public OperationsInterface< T_ScalarOut, T_Degree >
  {
   public:
    typedef T_ScalarOut D_ScalarOut;
    typedef T_ScalarIn D_ScalarIn;
    constexpr static Degree D_Degree = T_Degree;

    FieldScalarTypeOperation() :
      OperationsInterface< T_ScalarOut, T_Degree >()
    {
    }

    FieldScalarTypeOperation(const FieldScalarTypeOperation &other) :
      OperationsInterface< T_ScalarOut, T_Degree >(other)
    {
    }
  };

  template< class T_Op >
  class FieldScalarTypeNode final : public ExecutionTreeWithDegreeAndScalarModification< typename T_Op::D_ScalarOut, typename T_Op::D_ScalarIn, T_Op::D_Degree >
  {
   protected:
    T_Op _op;
    const field::FieldInterface< typename T_Op::D_ScalarIn > *_field;

    FieldScalarTypeNode(const FieldScalarTypeNode &other) :
      ExecutionTreeWithDegreeAndScalarModification< typename T_Op::D_ScalarOut, typename T_Op::D_ScalarIn, T_Op::D_Degree >(other), _op(other._op), _field(other._field)
    {
    }

   public:
    FieldScalarTypeNode(const field::FieldInterface< typename T_Op::D_ScalarIn > *const field) :
      ExecutionTreeWithDegreeAndScalarModification< typename T_Op::D_ScalarOut, typename T_Op::D_ScalarIn, T_Op::D_Degree >(), _op(), _field(field)
    {
    }

    FieldScalarTypeNode(const T_Op &other, const field::FieldInterface< typename T_Op::D_ScalarIn > *const field) :
      ExecutionTreeWithDegreeAndScalarModification< typename T_Op::D_ScalarOut, typename T_Op::D_ScalarIn, T_Op::D_Degree >(), _op(other), _field(field)
    {
    }

    virtual ~FieldScalarTypeNode() {}

    virtual NodeType nodeType() const override
    {
      return NodeType::FieldScalarType;
    }

    virtual bool isConstant(const std::pair< int, int > &entity) const override
    {
      return false;
    }

    virtual void evaluate(OutputVector< typename T_Op::D_ScalarOut, T_Op::D_Degree > &values, const std::vector< scalar::Precision< typename T_Op::D_ScalarIn >, numa::allocator< scalar::Precision< typename T_Op::D_ScalarIn > > > &points, const std::vector< scalar::Precision< typename T_Op::D_ScalarIn > > &gaussPoints, const int elementType, const std::pair< int, int > &entity) const override
    {
      if(values.size() < points.size() / 3) {
#pragma omp barrier
        FunctionAllocator< typename MathObject< typename T_Op::D_ScalarOut, T_Op::D_Degree >::Object >::initalizeNUMAMemory(points.size() / 3);
#pragma omp master
        values.resize(points.size() / 3);
#pragma omp barrier
      }

      _op(_field, values, points, gaussPoints, elementType, entity);
#pragma omp barrier
    }

    virtual FieldScalarTypeNode< T_Op > *copy() const override
    {
      return new FieldScalarTypeNode< T_Op >(*this);
    }

    const T_Op &operation() const
    {
      return _op;
    }

    virtual void exportTree(const domain::Domain &domain, common::HTMLNode *node) const override
    {
      node->name = _op.name() + " " + _field->name();
      node->scalar = scalar::Name< typename T_Op::D_ScalarOut >::name;
      node->degree = T_Op::D_Degree;
      node->nodeType = this->nodeType();
      node->constant = false;
    }

    virtual common::Memory peakMemoryByGaussPoint(const domain::Domain &domain) const override
    {
      return common::Memory(0);
    }

    virtual bool operator==(const ExecutionTreeInterface &other) const override
    {
      if(this->tag() == other.tag()) {
        return true;
      }

      const std::type_info &ext1 = typeid(*this);
      const std::type_info &ext2 = typeid(other);
      if(ext1.hash_code() == ext2.hash_code()) {
        if(this->_field == static_cast< const FieldScalarTypeNode< T_Op > & >(other)._field) {
          return true;
        }
      }

      return false;
    }
  };


} // namespace gmshfem::function


#endif // H_GMSHFEM_FIELDNODE
