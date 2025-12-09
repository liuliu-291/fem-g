// GmshFEM - Copyright (C) 2019-2022, A. Royer, E. Béchet, C. Geuzaine, Université de Liège
//
// See the LICENSE.txt file for license information. Please report all
// issues on https://gitlab.onelab.info/gmsh/fem/issues

#ifndef H_GMSHFEM_SKINLAYER
#define H_GMSHFEM_SKINLAYER

#include "GeometricObject.h"

#include <string>

namespace gmshfem::domain
{
  class Domain;
}

namespace gmshfem::domain
{


  class SkinLayer : public GeometricObject
  {
   protected:
    const int _entity;
    const Domain *_object;
    const Domain *_tool;

    SkinLayer(const int dim, const int entity, const Domain &object, const Domain &tool, const std::map< std::pair< int, int >, JacobiansModificator * > &jacobiansModificators);
    friend class Domain;

   public:
    SkinLayer();
    ~SkinLayer();

    void printDebug() const;
    void saveDebug(const std::string &name = "") const;
  };


} // namespace gmshfem::domain

#endif // H_GMSHFEM_SKINLAYER
