// GmshFEM - Copyright (C) 2019-2022, A. Royer, E. Béchet, C. Geuzaine, Université de Liège
//
// See the LICENSE.txt file for license information. Please report all
// issues on https://gitlab.onelab.info/gmsh/fem/issues

#ifndef H_GMSHFEM_BILINEARTERMINTERFACE
#define H_GMSHFEM_BILINEARTERMINTERFACE

#include "Term.h"
#include "UnknownField.h"

namespace gmshfem::term
{


  template< class T_Scalar >
  class BilinearTermInterface : public Term< T_Scalar >
  {
   protected:
    unsigned int _nbrOfDofsByElement[2];
    const field::FieldInterface< T_Scalar > *const _field[2];

   public:
    BilinearTermInterface(const domain::GeometricObject &domain, const std::string &integrationType, const field::FieldInterface< T_Scalar > *const fieldLhs, const field::FieldInterface< T_Scalar > *const fieldRhs, const ProductType productType);
    virtual ~BilinearTermInterface();

    virtual std::string termName() const override;
    virtual unsigned int nbrDofsByElement(unsigned int field = 0) const final;

    virtual bool isLinear() const final;
    virtual bool isBilinear() const final;

    virtual equation::UnknownFieldType unknownFieldType() const = 0;

    const field::FieldInterface< T_Scalar > *field(const unsigned int field) const;

    virtual void evaluate(Eigen::MatrixX< T_Scalar > &A_e, const scalar::Precision< T_Scalar > *const determinants, const scalar::Precision< T_Scalar > *const jacobians, const unsigned int elementIndex) = 0;
  };


} // namespace gmshfem::term

#include "BilinearTerm.h"

#endif // H_GMSHFEM_BILINEARTERMINTERFACE
