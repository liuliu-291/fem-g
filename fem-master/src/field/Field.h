// GmshFEM - Copyright (C) 2019-2022, A. Royer, E. Béchet, C. Geuzaine, Université de Liège
//
// See the LICENSE.txt file for license information. Please report all
// issues on https://gitlab.onelab.info/gmsh/fem/issues

#ifndef H_GMSHFEM_FIELD
#define H_GMSHFEM_FIELD

#include "Domain.h"
#include "FieldInterface.h"
#include "GeneralEvaluableObject.h"

#include <string>

namespace gmshfem::field
{


  template< class T_Scalar, field::Form T_Form >
  class Field : public FieldInterface< T_Scalar >, public function::GeneralEvaluableObject< T_Scalar, field::DegreeOfForm< T_Form >::value >
  {
   protected:
    std::vector< Constraint< T_Scalar, DegreeOfForm< T_Form >::value > > _constraints;
    function::Function< T_Scalar, field::DegreeOfForm< T_Form >::value > _myself;

    virtual bool _checkDomain() const override;

   public:
    Field();
    Field(const std::string &name, const domain::Domain &domain, const field::FunctionSpaceOfForm< T_Form > &type, const unsigned int order = 0, const std::string &model = "");
    Field(const std::string &name, const domain::Domain &domain, const FunctionSpace< scalar::Precision< T_Scalar >, T_Form > &functionSpace, const std::string &model = "");
    Field(const Field< T_Scalar, T_Form > &other);
    virtual ~Field();

    virtual void clear() override;

    virtual field::Form form() const override;
    virtual unsigned int multiplicity() const override;

    virtual field::FunctionSpace< scalar::Precision< T_Scalar >, T_Form > *getFunctionSpace() const override;

    Field< T_Scalar, T_Form > &operator=(const Field< T_Scalar, T_Form > &other);

    void addConstraint(const domain::Domain &domain, const function::Function< T_Scalar, DegreeOfForm< T_Form >::value > &function);
    std::vector< Constraint< T_Scalar, DegreeOfForm< T_Form >::value > > getConstraints() const;
    void removeConstraints();

    // GeneralEvaluableObject
    virtual const function::Function< T_Scalar, field::DegreeOfForm< T_Form >::value > &getEvaluableFunction() const override;
  };


} // namespace gmshfem::field


#endif // H_GMSHFEM_FIELD
