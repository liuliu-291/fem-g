// GmshFEM - Copyright (C) 2019-2022, A. Royer, E. Béchet, C. Geuzaine, Université de Liège
//
// See the LICENSE.txt file for license information. Please report all
// issues on https://gitlab.onelab.info/gmsh/fem/issues

#ifndef H_GMSHFEM_GLOBALQUANTITY
#define H_GMSHFEM_GLOBALQUANTITY

#include "Dof.h"
#include "Domain.h"
#include "FieldObject.h"

#include <string>

namespace gmshfem::problem
{
  template< class T_Scalar >
  class Formulation;
}

namespace gmshfem::dofs
{
  template< class T_Scalar >
  class DofsFactory;
}

namespace gmshfem::field
{
  template< class T_Scalar >
  class FieldInterface;
}

namespace gmshfem::field
{


  enum class FixedComponent {
    Primal,
    Dual,
    None
  };


  template< class T_Scalar >
  class GlobalQuantity
  {
   private:
    std::string _name;
    std::string _model;
    const unsigned int _tag;
    domain::Domain _domain;
    T_Scalar _primal;
    T_Scalar _dual;
    FixedComponent _fixedComponent;
    FieldInterface< T_Scalar > *_associatedPrimalField;
    FieldInterface< T_Scalar > *_associatedDualField;
    bool _isActivated;
    const dofs::Dof *_associatedPrimalDof;
    const dofs::Dof *_associatedDualDof;

    void _setPrimalValue(const T_Scalar &value);
    void _setDualValue(const T_Scalar &value);
    void _setFixedComponent(const FixedComponent fixedComponent);

   public:
    GlobalQuantity();
    GlobalQuantity(const std::string &name, const domain::Domain &domain, const std::string &model = "");
    GlobalQuantity(const GlobalQuantity &other);
    ~GlobalQuantity();

    GlobalQuantity< T_Scalar > &operator=(const GlobalQuantity &other);

    unsigned int tag() const;
    std::string name() const;
    std::string model() const;
    void name(const std::string &name);
    domain::Domain domain() const;
    void domain(const domain::Domain &domain);

    T_Scalar getPrimalValue() const;
    T_Scalar getDualValue() const;

    FixedComponent fixedComponent() const;
    bool isActivated() const;

    void setAssociatedPrimalDof(const dofs::Dof *dof);
    unsigned int getNumPrimalDof() const;
    void setAssociatedDualDof(const dofs::Dof *dof);
    unsigned int getNumDualDof() const;
    void setAssociatedDualField(FieldInterface< T_Scalar > *const field);
    FieldInterface< T_Scalar > *getAssociatedDualField() const;
    void setAssociatedPrimalField(FieldInterface< T_Scalar > *const field);
    FieldInterface< T_Scalar > *getAssociatedPrimalField() const;

    friend class problem::Formulation< T_Scalar >;
  };

  bool globalValueIsStillValid(const unsigned int tag);


} // namespace gmshfem::field

#endif // H_GMSHFEM_GLOBALQUANTITY
