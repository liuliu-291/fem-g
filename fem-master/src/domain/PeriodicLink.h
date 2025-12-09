// GmshFEM - Copyright (C) 2019-2022, A. Royer, E. Béchet, C. Geuzaine, Université de Liège
//
// See the LICENSE.txt file for license information. Please report all
// issues on https://gitlab.onelab.info/gmsh/fem/issues

#ifndef H_GMSHFEM_PERIODICLINK
#define H_GMSHFEM_PERIODICLINK

#include "Domain.h"

#include <string>

namespace gmshfem::domain
{


  class PeriodicLink : protected Domain
  {
   protected:
    Domain _master;

    void _buildLink();

   public:
    PeriodicLink();
    PeriodicLink(const int dim, const int tag);
    PeriodicLink(const std::pair< int, int > &dimTag);
    PeriodicLink(const int dim, const int tag, const JacobiansModificator &jacMod);
    PeriodicLink(const std::pair< int, int > &dimTag, const JacobiansModificator &jacMod);
    PeriodicLink(const std::string &name);
    PeriodicLink(const std::string &name, const JacobiansModificator &jacMod);
    PeriodicLink(const PeriodicLink &other);
    PeriodicLink(PeriodicLink &&other);
    virtual ~PeriodicLink();

    PeriodicLink &operator=(const PeriodicLink &other);
    PeriodicLink &operator=(PeriodicLink &&other);
    bool operator==(const PeriodicLink &other) const;
    bool operator!=(const PeriodicLink &other) const;
    PeriodicLink &operator|=(const PeriodicLink &other); // union
    PeriodicLink operator|(const PeriodicLink &other) const;
    PeriodicLink &operator&=(const PeriodicLink &other); // intersection
    PeriodicLink operator&(const PeriodicLink &other) const;
    PeriodicLink &operator^=(const PeriodicLink &other); // symmetric difference
    PeriodicLink operator^(const PeriodicLink &other) const;
    PeriodicLink operator~() const; // complement

    void printDebug() const;
    void saveDebug(const std::string &name = "") const;

    Domain slave() const;
    Domain master() const;

    bool isChainedTo(const PeriodicLink &other) const;
  };


} // namespace gmshfem::domain

#endif // H_GMSHFEM_PERIODICLINK
