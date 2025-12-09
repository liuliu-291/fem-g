// GmshFEM - Copyright (C) 2019-2022, A. Royer, E. Béchet, C. Geuzaine, Université de Liège
//
// See the LICENSE.txt file for license information. Please report all
// issues on https://gitlab.onelab.info/gmsh/fem/issues

#include "DofsManager.h"

#include "CSVio.h"
#include "Dof.h"
#include "DofsFactory.h"
#include "Exception.h"
#include "FieldInterface.h"
#include "Hilbert.h"
#include "MatrixFactory.h"
#include "OmpInterface.h"
#include "RCM.h"
#include "Term.h"
#include "instantiate.h"

#include <cstring>
#include <gmsh.h>
#include <string>
#include <unordered_map>
#include <utility>

namespace gmshfem::dofs
{


  template< class T_Scalar >
  DofsManager< T_Scalar >::DofsManager() :
    _fields(), _nbrDofs(0), _nbrUnknownDofs(0), _nbrFixedDofs(0), _nbrBubbleDofs(0), _nbrUnknownGlobalDofs(0), _nbrFixedGlobalDofs(0), _nbrLinkedDofs(0), _nbrBubbleLinkedDofs(0)
  {
  }

  template< class T_Scalar >
  DofsManager< T_Scalar >::~DofsManager()
  {
    clear();
  }

  template< class T_Scalar >
  DofsManager< T_Scalar >::DofsManager(DofsManager< T_Scalar > &&other) :
    _fields(std::move(other._fields)), _nbrDofs(std::move(other._nbrDofs)), _nbrUnknownDofs(std::move(other._nbrUnknownDofs)), _nbrFixedDofs(std::move(other._nbrFixedDofs)), _nbrBubbleDofs(std::move(other._nbrBubbleDofs)), _nbrUnknownGlobalDofs(std::move(other._nbrUnknownGlobalDofs)), _nbrFixedGlobalDofs(std::move(other._nbrFixedGlobalDofs)), _nbrLinkedDofs(std::move(other._nbrLinkedDofs)), _nbrBubbleLinkedDofs(std::move(other._nbrBubbleLinkedDofs))
  {
    other.clear();
  }

  template< class T_Scalar >
  void DofsManager< T_Scalar >::clear()
  {
    _nbrDofs = 0;
    _nbrUnknownDofs = 0;
    _nbrFixedDofs = 0;
    _nbrBubbleDofs = 0;
    _nbrUnknownGlobalDofs = 0;
    _nbrFixedGlobalDofs = 0;
    _nbrLinkedDofs = 0;
    _nbrBubbleLinkedDofs = 0;

    _fields.clear();
  }

  template< class T_Scalar >
  void DofsManager< T_Scalar >::addField(const field::FieldInterface< T_Scalar > *field)
  {
    _fields.insert(field);
  }

  template< class T_Scalar >
  unsigned long long DofsManager< T_Scalar >::nbrDofs() const
  {
    return _nbrDofs;
  }

  template< class T_Scalar >
  unsigned long long DofsManager< T_Scalar >::nbrUnknownDofs() const
  {
    return _nbrUnknownDofs;
  }

  template< class T_Scalar >
  unsigned long long DofsManager< T_Scalar >::nbrFixedDofs() const
  {
    return _nbrFixedDofs;
  }

  template< class T_Scalar >
  unsigned long long DofsManager< T_Scalar >::nbrBubbleDofs() const
  {
    return _nbrBubbleDofs;
  }

  template< class T_Scalar >
  unsigned long long DofsManager< T_Scalar >::nbrUnknownGlobalDofs() const
  {
    return _nbrUnknownGlobalDofs;
  }

  template< class T_Scalar >
  unsigned long long DofsManager< T_Scalar >::nbrFixedGlobalDofs() const
  {
    return _nbrFixedGlobalDofs;
  }

  template< class T_Scalar >
  unsigned long long DofsManager< T_Scalar >::nbrLinkedDofs() const
  {
    return _nbrLinkedDofs;
  }

  template< class T_Scalar >
  unsigned long long DofsManager< T_Scalar >::nbrBubbleLinkedDofs() const
  {
    return _nbrBubbleLinkedDofs;
  }

  template< class T_Scalar >
  void DofsManager< T_Scalar >::reorderWithHilbert(const bool bubble)
  {
    if(_nbrUnknownDofs == 0)
      return;

    const int dim = gmsh::model::getDimension();
    double xmin = 0., ymin = 0., zmin = 0., xmax = 0., ymax = 0., zmax = 0.;

    gmsh::model::getBoundingBox(-1, -1, xmin, ymin, zmin, xmax, ymax, zmax);

    if(bubble) {
      if(dim == 2) {
        std::vector< reorder::SortedEntity< 2 > > sortedEntitiesUnknownDofs;
        std::vector< reorder::SortedEntity< 2 > > sortedEntitiesUnknownBubbleDofs;
        sortedEntitiesUnknownDofs.reserve(_nbrUnknownDofs - _nbrBubbleDofs);
        sortedEntitiesUnknownBubbleDofs.reserve(_nbrBubbleDofs);

        reorder::SortedEntity< 2 > se;
        for(auto itField = _fields.begin(); itField != _fields.end(); ++itField) {
          for(auto it = (*itField)->begin(); it != (*itField)->end(); ++it) {
            if(it->first->type() == Type::Unknown) {
              se.f[0] = it->first->x();
              se.f[1] = it->first->y();
              se.ptrA = it->first;
              if(static_cast< const UnknownDof * >(it->first)->isBubble()) {
                sortedEntitiesUnknownBubbleDofs.push_back(se);
              }
              else {
                sortedEntitiesUnknownDofs.push_back(se);
              }
            }
          }
        }

        float min[2] = {static_cast< float >(xmin - 0.000001 * std::abs(xmin)),
                        static_cast< float >(ymin - 0.000001 * std::abs(ymin))};
        float max[2] = {static_cast< float >(xmax + 0.000001 * std::abs(xmax)),
                        static_cast< float >(ymax + 0.000001 * std::abs(ymax))};
        reorder::Hilbert< 2, 4 > hilbert;
        hilbert.apply(sortedEntitiesUnknownDofs, min, max);
        hilbert.apply(sortedEntitiesUnknownBubbleDofs, min, max);

        std::unordered_map< unsigned long long, unsigned long long > globalDofs;
        unsigned long long globalDofShift = 0;
        for(auto i = 0ULL; i < sortedEntitiesUnknownDofs.size(); ++i) {
          if(static_cast< UnknownDof * >(sortedEntitiesUnknownDofs[i].ptrB)->assembleType() == 4) {
            auto itFind = globalDofs.find(static_cast< UnknownDof * >(sortedEntitiesUnknownDofs[i].ptrB)->numDof());
            if(itFind == globalDofs.end()) {
              globalDofs.insert(std::make_pair(static_cast< UnknownDof * >(sortedEntitiesUnknownDofs[i].ptrB)->numDof(), i + 1 - globalDofShift));
              static_cast< UnknownDof * >(sortedEntitiesUnknownDofs[i].ptrB)->numDof(i + 1 - globalDofShift);
            }
            else {
              static_cast< UnknownDof * >(sortedEntitiesUnknownDofs[i].ptrB)->numDof(itFind->second);
              ++globalDofShift;
            }
          }
          else {
            static_cast< UnknownDof * >(sortedEntitiesUnknownDofs[i].ptrB)->numDof(i + 1 - globalDofShift);
          }
        }
#pragma omp parallel for num_threads(omp::getMaxThreads())
        for(auto i = 0ULL; i < sortedEntitiesUnknownBubbleDofs.size(); ++i) {
          static_cast< UnknownDof * >(sortedEntitiesUnknownBubbleDofs[i].ptrB)->numDof(_nbrUnknownDofs - _nbrBubbleDofs - globalDofShift + i + 1);
        }
      }
      else if(dim == 3) {
        std::vector< reorder::SortedEntity< 3 > > sortedEntitiesUnknownDofs;
        std::vector< reorder::SortedEntity< 3 > > sortedEntitiesUnknownBubbleDofs;
        sortedEntitiesUnknownDofs.reserve(_nbrUnknownDofs - _nbrBubbleDofs);
        sortedEntitiesUnknownBubbleDofs.reserve(_nbrBubbleDofs);

        reorder::SortedEntity< 3 > se;
        for(auto itField = _fields.begin(); itField != _fields.end(); ++itField) {
          for(auto it = (*itField)->begin(); it != (*itField)->end(); ++it) {
            if(it->first->type() == Type::Unknown) {
              se.f[0] = it->first->x();
              se.f[1] = it->first->y();
              se.f[2] = it->first->z();
              se.ptrA = it->first;
              if(static_cast< const UnknownDof * >(it->first)->isBubble()) {
                sortedEntitiesUnknownBubbleDofs.push_back(se);
              }
              else {
                sortedEntitiesUnknownDofs.push_back(se);
              }
            }
          }
        }

        float min[3] = {static_cast< float >(xmin - 0.000001 * std::abs(xmin)),
                        static_cast< float >(ymin - 0.000001 * std::abs(ymin)),
                        static_cast< float >(zmin - 0.000001 * std::abs(zmin))};
        float max[3] = {static_cast< float >(xmax + 0.000001 * std::abs(xmax)),
                        static_cast< float >(ymax + 0.000001 * std::abs(ymax)),
                        static_cast< float >(zmax + 0.000001 * std::abs(zmax))};
        reorder::Hilbert< 3, 8 > hilbert;
        hilbert.apply(sortedEntitiesUnknownDofs, min, max);
        hilbert.apply(sortedEntitiesUnknownBubbleDofs, min, max);

        std::unordered_map< unsigned long long, unsigned long long > globalDofs;
        unsigned long long globalDofShift = 0;
        for(auto i = 0ULL; i < sortedEntitiesUnknownDofs.size(); ++i) {
          if(static_cast< UnknownDof * >(sortedEntitiesUnknownDofs[i].ptrB)->assembleType() == 4) {
            auto itFind = globalDofs.find(static_cast< UnknownDof * >(sortedEntitiesUnknownDofs[i].ptrB)->numDof());
            if(itFind == globalDofs.end()) {
              globalDofs.insert(std::make_pair(static_cast< UnknownDof * >(sortedEntitiesUnknownDofs[i].ptrB)->numDof(), i + 1 - globalDofShift));
              static_cast< UnknownDof * >(sortedEntitiesUnknownDofs[i].ptrB)->numDof(i + 1 - globalDofShift);
            }
            else {
              static_cast< UnknownDof * >(sortedEntitiesUnknownDofs[i].ptrB)->numDof(itFind->second);
              ++globalDofShift;
            }
          }
          else {
            static_cast< UnknownDof * >(sortedEntitiesUnknownDofs[i].ptrB)->numDof(i + 1 - globalDofShift);
          }
        }
#pragma omp parallel for num_threads(omp::getMaxThreads())
        for(auto i = 0ULL; i < sortedEntitiesUnknownBubbleDofs.size(); ++i) {
          static_cast< UnknownDof * >(sortedEntitiesUnknownBubbleDofs[i].ptrB)->numDof(_nbrUnknownDofs - _nbrBubbleDofs - globalDofShift + i + 1);
        }
      }
    }
    else {
      if(dim == 2) {
        std::vector< reorder::SortedEntity< 2 > > sortedEntitiesUnknownDofs;
        sortedEntitiesUnknownDofs.reserve(_nbrUnknownDofs);

        reorder::SortedEntity< 2 > se;
        for(auto itField = _fields.begin(); itField != _fields.end(); ++itField) {
          for(auto it = (*itField)->begin(); it != (*itField)->end(); ++it) {
            if(it->first->type() == Type::Unknown) {
              se.f[0] = it->first->x();
              se.f[1] = it->first->y();
              se.ptrA = it->first;
              sortedEntitiesUnknownDofs.push_back(se);
            }
          }
        }

        float min[2] = {static_cast< float >(xmin - 0.000001 * std::abs(xmin)),
                        static_cast< float >(ymin - 0.000001 * std::abs(ymin))};
        float max[2] = {static_cast< float >(xmax + 0.000001 * std::abs(xmax)),
                        static_cast< float >(ymax + 0.000001 * std::abs(ymax))};
        reorder::Hilbert< 2, 4 > hilbert;
        hilbert.apply(sortedEntitiesUnknownDofs, min, max);

        std::unordered_map< unsigned long long, unsigned long long > globalDofs;
        unsigned long long globalDofShift = 0;
        for(auto i = 0ULL; i < sortedEntitiesUnknownDofs.size(); ++i) {
          if(static_cast< UnknownDof * >(sortedEntitiesUnknownDofs[i].ptrB)->assembleType() == 4) {
            auto itFind = globalDofs.find(static_cast< UnknownDof * >(sortedEntitiesUnknownDofs[i].ptrB)->numDof());
            if(itFind == globalDofs.end()) {
              globalDofs.insert(std::make_pair(static_cast< UnknownDof * >(sortedEntitiesUnknownDofs[i].ptrB)->numDof(), i + 1 - globalDofShift));
              static_cast< UnknownDof * >(sortedEntitiesUnknownDofs[i].ptrB)->numDof(i + 1 - globalDofShift);
            }
            else {
              static_cast< UnknownDof * >(sortedEntitiesUnknownDofs[i].ptrB)->numDof(itFind->second);
              ++globalDofShift;
            }
          }
          else {
            static_cast< UnknownDof * >(sortedEntitiesUnknownDofs[i].ptrB)->numDof(i + 1 - globalDofShift);
          }
        }
      }
      else if(dim == 3) {
        std::vector< reorder::SortedEntity< 3 > > sortedEntitiesUnknownDofs;
        sortedEntitiesUnknownDofs.reserve(_nbrUnknownDofs);

        reorder::SortedEntity< 3 > se;
        for(auto itField = _fields.begin(); itField != _fields.end(); ++itField) {
          for(auto it = (*itField)->begin(); it != (*itField)->end(); ++it) {
            if(it->first->type() == Type::Unknown) {
              se.f[0] = it->first->x();
              se.f[1] = it->first->y();
              se.f[2] = it->first->z();
              se.ptrA = it->first;
              sortedEntitiesUnknownDofs.push_back(se);
            }
          }
        }

        float min[3] = {static_cast< float >(xmin - 0.000001 * std::abs(xmin)),
                        static_cast< float >(ymin - 0.000001 * std::abs(ymin)),
                        static_cast< float >(zmin - 0.000001 * std::abs(zmin))};
        float max[3] = {static_cast< float >(xmax + 0.000001 * std::abs(xmax)),
                        static_cast< float >(ymax + 0.000001 * std::abs(ymax)),
                        static_cast< float >(zmax + 0.000001 * std::abs(zmax))};
        reorder::Hilbert< 3, 8 > hilbert;
        hilbert.apply(sortedEntitiesUnknownDofs, min, max);

        std::unordered_map< unsigned long long, unsigned long long > globalDofs;
        unsigned long long globalDofShift = 0;
        for(auto i = 0ULL; i < sortedEntitiesUnknownDofs.size(); ++i) {
          if(static_cast< UnknownDof * >(sortedEntitiesUnknownDofs[i].ptrB)->assembleType() == 4) {
            auto itFind = globalDofs.find(static_cast< UnknownDof * >(sortedEntitiesUnknownDofs[i].ptrB)->numDof());
            if(itFind == globalDofs.end()) {
              globalDofs.insert(std::make_pair(static_cast< UnknownDof * >(sortedEntitiesUnknownDofs[i].ptrB)->numDof(), i + 1 - globalDofShift));
              static_cast< UnknownDof * >(sortedEntitiesUnknownDofs[i].ptrB)->numDof(i + 1 - globalDofShift);
            }
            else {
              static_cast< UnknownDof * >(sortedEntitiesUnknownDofs[i].ptrB)->numDof(itFind->second);
              ++globalDofShift;
            }
          }
          else {
            static_cast< UnknownDof * >(sortedEntitiesUnknownDofs[i].ptrB)->numDof(i + 1 - globalDofShift);
          }
        }
      }
    }
  }

  template< class T_Scalar >
  void DofsManager< T_Scalar >::reorderWithRCM(system::MatrixFactory< T_Scalar > *const A)
  {
    // dofs::DofsManager is a friend of system::Matrix
    unsigned long long size = A->_size;
    unsigned long long *row = A->_ai;
    unsigned long long *indices = A->_aj;

    std::vector< unsigned long long > sorted(size, 0);
    reorder::RCM rcm;
    rcm.apply(sorted, row, indices);

    for(auto itField = _fields.begin(); itField != _fields.end(); ++itField) {
      for(auto it = (*itField)->begin(); it != (*itField)->end(); ++it) {
        if(it->first->type() == Type::Unknown) {
          static_cast< UnknownDof * >(it->first)->numDof(sorted[it->first->numDof() - 1]);
        }
      }
    }

    A->_reorder(sorted);
  }

  template< class T_Scalar >
  void DofsManager< T_Scalar >::built(const std::vector< term::Term< T_Scalar > * > &terms, const std::vector< std::pair< unsigned int, field::FieldInterface< T_Scalar > * > > &fields)
  {
    std::vector< field::FieldInterface< T_Scalar > * > fieldsCp(fields.size());
    for(auto i = 0ULL; i < fields.size(); ++i) {
      fieldsCp[i] = fields[i].second;
    }
    DofsFactory< T_Scalar > factory(this);
    factory.generate(terms, fieldsCp);

    _nbrFixedDofs = factory.nbrFixedDofsCreated();
    _nbrUnknownDofs = factory.nbrUnknownDofsCreated();
    _nbrBubbleDofs = factory.nbrBubbleUnknownDofsCreated();
    _nbrUnknownGlobalDofs = factory.nbrUnknownGlobalDofsCreated();
    _nbrFixedGlobalDofs = factory.nbrFixedGlobalDofsCreated();
    _nbrLinkedDofs = factory.nbrLinkedDofsCreated();
    _nbrBubbleLinkedDofs = factory.nbrBubbleLinkedDofsCreated();
    _nbrDofs = _nbrFixedDofs + _nbrUnknownDofs + _nbrLinkedDofs;
  }

  INSTANTIATE_CLASS(DofsManager, 4, TEMPLATE_ARGS(std::complex< double >, std::complex< float >, double, float))


} // namespace gmshfem::dofs
