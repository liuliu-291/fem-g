// GmshFEM - Copyright (C) 2019-2022, A. Royer, E. Béchet, C. Geuzaine, Université de Liège
//
// See the LICENSE.txt file for license information. Please report all
// issues on https://gitlab.onelab.info/gmsh/fem/issues

#include "PeriodicLink.h"

#include "Exception.h"
#include "Message.h"
#include "Post.h"
#include "ScalarFunction.h"

#include <algorithm>
#include <gmsh.h>
#include <unordered_set>

namespace gmshfem::domain
{


  PeriodicLink::PeriodicLink() :
    Domain(), _master()
  {
    _buildLink();
  }

  PeriodicLink::PeriodicLink(const int dim, const int tag) :
    Domain(dim, tag), _master()
  {
    _buildLink();
  }

  PeriodicLink::PeriodicLink(const std::pair< int, int > &dimTag) :
    Domain(dimTag), _master()
  {
    _buildLink();
  }

  PeriodicLink::PeriodicLink(const int dim, const int tag, const JacobiansModificator &jacMod) :
    Domain(dim, tag, jacMod), _master()
  {
    _buildLink();
  }

  PeriodicLink::PeriodicLink(const std::pair< int, int > &dimTag, const JacobiansModificator &jacMod) :
    Domain(dimTag, jacMod), _master()
  {
    _buildLink();
  }

  PeriodicLink::PeriodicLink(const std::string &name) :
    Domain(name), _master()
  {
    _buildLink();
  }

  PeriodicLink::PeriodicLink(const std::string &name, const JacobiansModificator &jacMod) :
    Domain(name, jacMod), _master()
  {
    _buildLink();
  }

  PeriodicLink::PeriodicLink(const PeriodicLink &other) :
    Domain(other), _master(other._master)
  {
  }

  PeriodicLink::PeriodicLink(PeriodicLink &&other) :
    Domain(other), _master(std::move(other._master))
  {
  }

  PeriodicLink::~PeriodicLink()
  {
  }

  void PeriodicLink::_buildLink()
  {
    std::set< std::pair< int, int > > newEntities;
    for(auto it = _entities.begin(); it != _entities.end(); ++it) {
      std::vector< int > tagMaster;
      gmsh::model::mesh::getPeriodic(it->first, {it->second}, tagMaster);
      if(tagMaster[0] == it->second) {
        msg::warning << "Entity (" << it->first << ", " << it->second << ") does not have master entity" << msg::endl;
      }
      else {
        newEntities.insert(*it);
        auto itJacMod = _jacobiansModificators.find(*it);
        if(itJacMod != _jacobiansModificators.end()) {
          reinterpret_cast< PeriodicLink * >(&_master)->_addEntity(it->first, tagMaster[0], itJacMod->second->copy()); // Hack because protected _addEntity can only by access through a PeriodicLink class
        }
        else {
          reinterpret_cast< PeriodicLink * >(&_master)->_addEntity(it->first, tagMaster[0]); // Hack because protected _addEntity can only by access through a PeriodicLink class
        }
      }
    }
    _entities = std::move(newEntities);
  }

  PeriodicLink &PeriodicLink::operator=(const PeriodicLink &other)
  {
    Domain::operator=(other);
    _master = other._master;
    return *this;
  }

  PeriodicLink &PeriodicLink::operator=(PeriodicLink &&other)
  {
    Domain::operator=(other);
    _master = other._master;
    return *this;
  }

  bool PeriodicLink::operator==(const PeriodicLink &other) const
  {
    return Domain::operator==(other) && _master == other._master;
  }

  bool PeriodicLink::operator!=(const PeriodicLink &other) const
  {
    return !((*this) == other);
  }

  PeriodicLink &PeriodicLink::operator|=(const PeriodicLink &other)
  {
    Domain::operator|=(other);
    _master |= other._master;
    return *this;
  }

  PeriodicLink PeriodicLink::operator|(const PeriodicLink &other) const
  {
    PeriodicLink link(*this);
    link |= other;

    return link;
  }

  PeriodicLink &PeriodicLink::operator&=(const PeriodicLink &other)
  {
    Domain::operator&=(other);
    _master.clear();
    _buildLink();
    return *this;
  }

  PeriodicLink PeriodicLink::operator&(const PeriodicLink &other) const
  {
    PeriodicLink link(*this);
    link &= other;

    return link;
  }

  PeriodicLink &PeriodicLink::operator^=(const PeriodicLink &other)
  {
    Domain::operator^=(other);
    _master.clear();
    _buildLink();
    return *this;
  }

  PeriodicLink PeriodicLink::operator^(const PeriodicLink &other) const
  {
    PeriodicLink link(*this);
    link ^= other;

    return link;
  }

  PeriodicLink PeriodicLink::operator~() const
  {
    Domain slave = Domain::operator~();
    PeriodicLink link;
    static_cast< Domain >(link) = slave;
    link._buildLink();
    return link;
  }

  void PeriodicLink::printDebug() const
  {
    msg::debug << "PeriodicLink has:" << msg::endl;
    msg::debug << " * slave elementary entities (" << _entities.size() << "):" << msg::endl;
    for(auto it = _entities.begin(); it != _entities.end(); ++it) {
      msg::debug << "  - " << it->first << " / " << it->second << msg::endl;
    }
    msg::debug << " * master elementary entities (" << _master.numberOfEntities() << "):" << msg::endl;
    for(auto it = _master.cbegin(); it != _master.cend(); ++it) {
      msg::debug << "  - " << it->first << " / " << it->second << msg::endl;
    }
  }

  void PeriodicLink::saveDebug(const std::string &name) const
  {
    Domain::saveDebug("slave " + name);
    _master.saveDebug("master " + name);
  }

  Domain PeriodicLink::slave() const
  {
    return *this;
  }

  Domain PeriodicLink::master() const
  {
    return _master;
  }

  bool PeriodicLink::isChainedTo(const PeriodicLink &other) const
  {
    return (_master & Domain(other)).isEmpty();
  }


} // namespace gmshfem::domain
