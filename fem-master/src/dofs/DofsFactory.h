// GmshFEM - Copyright (C) 2019-2022, A. Royer, E. Béchet, C. Geuzaine, Université de Liège
//
// See the LICENSE.txt file for license information. Please report all
// issues on https://gitlab.onelab.info/gmsh/fem/issues

#ifndef H_GMSHFEM_DOFSFACTORY
#define H_GMSHFEM_DOFSFACTORY

#include "FieldInterface.h"
#include "scalar.h"

#include <string>
#include <vector>

namespace gmshfem::dofs
{
  template< class T_Scalar >
  class DofsManager;
}

namespace gmshfem::term
{
  template< class T_Scalar >
  class Term;
}

namespace gmshfem::domain
{
  class Domain;
}

namespace gmshfem::dofs
{


  template< class T_Scalar >
  class DofsFactory
  {
   private:
    DofsManager< T_Scalar > *_dofsM;
    unsigned long long _nbrFixedDofs;
    unsigned long long _nbrFixedGlobalDofs;
    unsigned long long _nbrUnknownDofs;
    unsigned long long _nbrBubbleUnknownDofs;
    unsigned long long _nbrUnknownGlobalDofs;
    unsigned long long _nbrLinkedDofs;
    unsigned long long _nbrBubbleLinkedDofs;

    // FixedDofs
    template< field::Form T_Form >
    void _generateFixedDofs(field::Field< T_Scalar, T_Form > *const field, const domain::Domain &domain, const function::Function< T_Scalar, field::DegreeOfForm< T_Form >::value > &function);
    template< field::Form T_Form, unsigned int T_NumFields >
    void _generateCompoundFixedDofs(field::CompoundField< T_Scalar, T_Form, T_NumFields > *const field, const domain::Domain &domain, const function::Function< T_Scalar, field::DegreeOfCompoundForm< T_Form >::value > &function);
    template< field::Form T_Form, unsigned int T_NumFields, unsigned int T_Component, class = std::enable_if< (T_Component < T_NumFields) > >
    void _generateCompoundFixedDofsOnComponent(field::CompoundField< T_Scalar, T_Form, T_NumFields > *const field, const domain::Domain &domain, const function::Function< T_Scalar, field::DegreeOfForm< T_Form >::value > &function);

    void _generateFixedGlobalDofs(field::GlobalQuantity< T_Scalar > *const globalQuantity, const domain::Domain &domain, const bool primal);

    // UnknownDofs
    void _generateUnknownDofs(field::FieldInterface< T_Scalar > *const field, const domain::Domain &domain, const unsigned int maxDim);
    void _generateUnknownGlobalDofs(field::GlobalQuantity< T_Scalar > *const globalQuantity, const domain::Domain &domain, const bool primal);

    // LinkedDofs
    void _generateLinkedDofs(field::FieldInterface< T_Scalar > *const field, const domain::Domain &master, const domain::Domain &slave, const T_Scalar &coefficient, const unsigned int maxDim);

   public:
    explicit DofsFactory(DofsManager< T_Scalar > *const dofsM);
    void generate(const std::vector< term::Term< T_Scalar > * > &terms, const std::vector< field::FieldInterface< T_Scalar > * > &fields);

    unsigned long long nbrFixedDofsCreated() const;
    unsigned long long nbrFixedGlobalDofsCreated() const;
    unsigned long long nbrUnknownDofsCreated() const;
    unsigned long long nbrBubbleUnknownDofsCreated() const;
    unsigned long long nbrUnknownGlobalDofsCreated() const;
    unsigned long long nbrLinkedDofsCreated() const;
    unsigned long long nbrBubbleLinkedDofsCreated() const;

    template< field::Form T_Form >
    struct GenerateFixedDofs {
      static void run(DofsFactory< T_Scalar > &enclose, field::FieldInterface< T_Scalar > *field);
    };
  };


} // namespace gmshfem::dofs

#endif // H_GMSHFEM_DOFSFACTORY
