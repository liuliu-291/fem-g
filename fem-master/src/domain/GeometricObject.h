// GmshFEM - Copyright (C) 2019-2022, A. Royer, E. Béchet, C. Geuzaine, Université de Liège
//
// See the LICENSE.txt file for license information. Please report all
// issues on https://gitlab.onelab.info/gmsh/fem/issues

#ifndef H_GMSHFEM_GEOMETRICOBJECT
#define H_GMSHFEM_GEOMETRICOBJECT

#include "JacobiansModificator.h"
#include "numa.h"

#include <cstddef>
#include <map>
#include <set>

namespace gmshfem::domain
{


  class GeometricObject
  {
   protected:
    std::set< std::pair< int, int > > _entities;
    std::map< std::pair< int, int >, JacobiansModificator * > _jacobiansModificators;

    void _addEntity(const int dim, const int tag);
    void _addEntity(const int dim, const int tag, JacobiansModificator *jacMod);

   public:
    GeometricObject();
    GeometricObject(const GeometricObject &other);
    GeometricObject(GeometricObject &&other);
    virtual ~GeometricObject();

    virtual bool isEmpty() const;
    virtual void clear();
    std::set< std::pair< int, int > >::const_iterator cbegin() const;
    std::set< std::pair< int, int > >::const_iterator cend() const;
    unsigned int numberOfEntities() const;

    bool have(const std::pair< int, int > &entity) const;
    unsigned int maxDim() const;
    unsigned int minDim() const;

    size_t hash() const;

    bool haveJacobiansModificators(const std::pair< int, int > &dimTag) const;
    bool needGaussCoordinates() const;

    void applyJacobiansModificator(const std::vector< double, numa::allocator< double > > &points, std::vector< double, numa::allocator< double > > &determinants, std::vector< double, numa::allocator< double > > &jacobians, const std::pair< int, int > &dimTag) const;
    void applyJacobiansModificator(const std::vector< float, numa::allocator< float > > &points, std::vector< float, numa::allocator< float > > &determinants, std::vector< float, numa::allocator< float > > &jacobians, const std::pair< int, int > &dimTag) const;
  };


} // namespace gmshfem::domain

#endif // H_GMSHFEM_GEOMETRICOBJECT
