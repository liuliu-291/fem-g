// GmshFEM - Copyright (C) 2019-2022, A. Royer, E. Béchet, C. Geuzaine, Université de Liège
//
// See the LICENSE.txt file for license information. Please report all
// issues on https://gitlab.onelab.info/gmsh/fem/issues

#ifndef H_GMSHFEM_LINEARTERMINTERFACE
#define H_GMSHFEM_LINEARTERMINTERFACE

#include "Term.h"

namespace gmshfem::term
{


  template< class T_Scalar >
  class LinearTermInterface : public Term< T_Scalar >
  {
   protected:
    unsigned int _nbrOfDofsByElement;
    const field::FieldInterface< T_Scalar > *const _fieldRhs;
    unsigned _rhsIdx;

   public:
    LinearTermInterface(const domain::GeometricObject &domain, const std::string &integrationType, const field::FieldInterface< T_Scalar > *const fieldRhs, const ProductType productType, unsigned rhsIdx = 0);
    virtual ~LinearTermInterface();

    virtual std::string termName() const override;
    virtual unsigned int nbrDofsByElement(unsigned int field = 0) const final;

    virtual bool isLinear() const final;
    virtual bool isBilinear() const final;

    unsigned int rhsIdx() const {return _rhsIdx;}

    const field::FieldInterface< T_Scalar > *field() const;

    virtual void evaluate(Eigen::VectorX< T_Scalar > &b_e, const scalar::Precision< T_Scalar > *const determinants, const scalar::Precision< T_Scalar > *const jacobians, const unsigned int elementIndex) = 0;
  };


} // namespace gmshfem::term

#include "LinearTerm.h"

#endif // H_GMSHFEM_LINEARTERMINTERFACE
