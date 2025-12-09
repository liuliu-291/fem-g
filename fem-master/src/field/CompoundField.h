// GmshFEM - Copyright (C) 2019-2022, A. Royer, E. Béchet, C. Geuzaine, Université de Liège
//
// See the LICENSE.txt file for license information. Please report all
// issues on https://gitlab.onelab.info/gmsh/fem/issues

#ifndef H_GMSHFEM_COMPOUNDFIELD
#define H_GMSHFEM_COMPOUNDFIELD

#include "Domain.h"
#include "FieldInterface.h"

#include <string>

namespace gmshfem::field
{


  template< class T_Scalar, field::Form T_Form, unsigned int T_NumFields >
  class CompoundField : public FieldInterface< T_Scalar >, public function::GeneralEvaluableObject< T_Scalar, field::DegreeOfCompoundForm< T_Form >::value >
  {
   protected:
    std::vector< Constraint< T_Scalar, DegreeOfCompoundForm< T_Form >::value > > _constraints;
    std::array< std::vector< Constraint< T_Scalar, DegreeOfForm< T_Form >::value > >, T_NumFields > _constraintsOnComponent;
    function::Function< T_Scalar, field::DegreeOfCompoundForm< T_Form >::value > _myself;

    virtual bool _checkDomain() const override;

   public:
    CompoundField();
    CompoundField(const std::string &name, const domain::Domain &domain, const field::FunctionSpaceOfForm< T_Form > &type, const unsigned int order = 0, const std::string &model = "");
    CompoundField(const CompoundField< T_Scalar, T_Form, T_NumFields > &other);
    virtual ~CompoundField();

    virtual void clear() override;

    virtual field::Form form() const override;
    virtual unsigned int multiplicity() const override;

    virtual field::FunctionSpace< scalar::Precision< T_Scalar >, T_Form > *getFunctionSpace() const override;

    CompoundField< T_Scalar, T_Form, T_NumFields > &operator=(const CompoundField< T_Scalar, T_Form, T_NumFields > &other);

    void addConstraint(const domain::Domain &domain, const function::Function< T_Scalar, DegreeOfCompoundForm< T_Form >::value > &function);
    std::vector< Constraint< T_Scalar, DegreeOfCompoundForm< T_Form >::value > > getConstraints() const;

    template< unsigned int T_Component, class = std::enable_if< (T_Component < T_NumFields) > >
    void addConstraintOnComponent(const domain::Domain &domain, const function::Function< T_Scalar, DegreeOfForm< T_Form >::value > &function);
    template< unsigned int T_Component, class = std::enable_if< (T_Component < T_NumFields) > >
    std::vector< Constraint< T_Scalar, DegreeOfForm< T_Form >::value > > getConstraintsOnComponent() const;
    
    void removeConstraints();

    // GeneralEvaluableObject
    virtual const function::Function< T_Scalar, field::DegreeOfCompoundForm< T_Form >::value > &getEvaluableFunction() const override;
  };


} // namespace gmshfem::field


#endif // H_GMSHFEM_COMPOUNDFIELD
