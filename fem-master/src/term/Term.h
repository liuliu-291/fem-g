// GmshFEM - Copyright (C) 2019-2022, A. Royer, E. Béchet, C. Geuzaine, Université de Liège
//
// See the LICENSE.txt file for license information. Please report all
// issues on https://gitlab.onelab.info/gmsh/fem/issues

#ifndef H_GMSHFEM_TERM
#define H_GMSHFEM_TERM

#include "Domain.h"
#include "EquationEvaluator.h"
#include "Function.h"
#include "FunctionSpaceInterface.h"
#include "MathObject.h"
#include "Memory.h"
#include "numa.h"
#include "scalar.h"

#include <gmsh.h>
#include <set>
#include <string>
#include <vector>

namespace gmshfem::problem
{
  template< class T_PScalar >
  class FunctionSpaceBucket;
}

namespace gmshfem::system
{
  template< class T_Scalar >
  class MatrixFactory;

  template< class T_Scalar >
  class VectorFactory;
} // namespace gmshfem::system

namespace gmshfem::field
{
  template< class T_Scalar >
  class FieldInterface;
}

namespace gmshfem::term
{


  template< class T_Scalar >
  class Term
  {
   protected:
    const unsigned int _tag;
    domain::GeometricObject _domain;
    std::string _integrationType;
    ProductType _productType;
    unsigned int _numberOfGaussPoints;
    bool _activated;

   public:
    Term(const domain::GeometricObject &domain, const std::string &integrationType, const ProductType productType);
    virtual ~Term();

    virtual std::string termName() const = 0;
    unsigned int tag() const;
    std::string integrationType() const;
    unsigned int nbrGaussPoints() const;
    virtual common::Memory memory() const = 0;
    bool isActivated() const;
    void activate();
    void deactivate();

    virtual unsigned int nbrDofsByElement(unsigned int field = 0) const = 0;

    virtual bool isLinear() const = 0;
    virtual bool isBilinear() const = 0;

    virtual bool needJacobians() const = 0;
    virtual bool needGaussCoordinates(const std::pair< int, int > &entity) const = 0;

    const domain::GeometricObject &domain() const;

    virtual void assemblyInitialization(problem::FunctionSpaceBucket< scalar::Precision< T_Scalar > > &functionSpaces, const std::vector< scalar::Precision< T_Scalar >, numa::allocator< scalar::Precision< T_Scalar > > > &points, const int elementType, const std::pair< int, int > &entity) = 0;
    virtual Term< T_Scalar > *switchTestFunctionField(const field::FieldInterface< T_Scalar > *primal, field::FieldInterface< T_Scalar > *dual) const = 0;
  };


} // namespace gmshfem::term

#include "BilinearTermInterface.h"
#include "LinearTermInterface.h"

#endif // H_GMSHFEM_TERM
