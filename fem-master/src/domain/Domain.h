// GmshFEM - Copyright (C) 2019-2022, A. Royer, E. Béchet, C. Geuzaine, Université de Liège
//
// See the LICENSE.txt file for license information. Please report all
// issues on https://gitlab.onelab.info/gmsh/fem/issues

#ifndef H_GMSHFEM_DOMAIN
#define H_GMSHFEM_DOMAIN

#include "GeometricObject.h"
#include "SkinLayer.h"

#include <string>

namespace gmshfem::domain
{


  class Domain : public GeometricObject
  {
   public:
    Domain();
    Domain(const int dim, const int tag);
    Domain(const std::pair< int, int > &dimTag);
    Domain(const int dim, const int tag, const JacobiansModificator &jacMod);
    Domain(const std::pair< int, int > &dimTag, const JacobiansModificator &jacMod);
    Domain(const std::string &name);
    Domain(const std::string &name, const JacobiansModificator &jacMod);
    Domain(const Domain &other);
    Domain(Domain &&other);
    virtual ~Domain();

    void addEntity(const int dimEntity, const int tagEntity);

    Domain &operator=(const Domain &other);
    Domain &operator=(Domain &&other);
    bool operator==(const Domain &other) const;
    bool operator!=(const Domain &other) const;
    Domain &operator|=(const Domain &other); // union
    Domain operator|(const Domain &other) const;
    Domain &operator&=(const Domain &other); // intersection
    Domain operator&(const Domain &other) const;
    Domain &operator^=(const Domain &other); // symmetric difference
    Domain operator^(const Domain &other) const;
    Domain operator~() const; // complement

    Domain getBoundary(const bool combined = true) const;
    SkinLayer getSkinLayer(const Domain &tool) const;

    bool isIncludedInto(const Domain &other) const;
    bool included(const Domain &other) const;

    void printDebug() const;
    void saveDebug(const std::string &name = "") const;
  };


} // namespace gmshfem::domain

#include "PeriodicLink.h"

#endif // H_GMSHFEM_DOMAIN
