// GmshFEM - Copyright (C) 2019-2022, A. Royer, E. Béchet, C. Geuzaine, Université de Liège
//
// See the LICENSE.txt file for license information. Please report all
// issues on https://gitlab.onelab.info/gmsh/fem/issues

#ifndef H_GMSHFEM_DOFSMANAGER
#define H_GMSHFEM_DOFSMANAGER

#include "scalar.h"

#include <string>
#include <unordered_map>
#include <unordered_set>
#include <vector>

namespace gmshfem::dofs
{
  class Dof;
}

namespace gmshfem::term
{
  template< class T_Scalar >
  class Term;
}

namespace gmshfem::field
{
  template< class T_Scalar >
  class FieldInterface;
}

namespace gmshfem::system
{
  template< class T_Scalar >
  class MatrixFactory;
}

namespace gmshfem::dofs
{


  template< class T_Scalar >
  class DofsManager
  {
   private:
    std::unordered_set< const field::FieldInterface< T_Scalar > * > _fields;

    unsigned long long _nbrDofs;
    unsigned long long _nbrUnknownDofs;
    unsigned long long _nbrFixedDofs;
    unsigned long long _nbrBubbleDofs;
    unsigned long long _nbrUnknownGlobalDofs;
    unsigned long long _nbrFixedGlobalDofs;
    unsigned long long _nbrLinkedDofs;
    unsigned long long _nbrBubbleLinkedDofs;

   public:
    DofsManager();
    ~DofsManager();

    DofsManager(DofsManager< T_Scalar > &&other);

    void clear();

    void addField(const field::FieldInterface< T_Scalar > *field);

    unsigned long long nbrDofs() const;
    unsigned long long nbrUnknownDofs() const;
    unsigned long long nbrFixedDofs() const;
    unsigned long long nbrBubbleDofs() const;
    unsigned long long nbrUnknownGlobalDofs() const;
    unsigned long long nbrFixedGlobalDofs() const;
    unsigned long long nbrLinkedDofs() const;
    unsigned long long nbrBubbleLinkedDofs() const;

    void reorderWithHilbert(const bool bubble = false);
    void reorderWithRCM(system::MatrixFactory< T_Scalar > *const A);

    void built(const std::vector< term::Term< T_Scalar > * > &terms, const std::vector< std::pair< unsigned int, field::FieldInterface< T_Scalar > * > > &fields);
  };


} // namespace gmshfem::dofs

#endif // H_GMSHFEM_DOFSMANAGER
