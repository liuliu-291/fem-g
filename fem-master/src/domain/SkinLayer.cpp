// GmshFEM - Copyright (C) 2019-2022, A. Royer, E. Béchet, C. Geuzaine, Université de Liège
//
// See the LICENSE.txt file for license information. Please report all
// issues on https://gitlab.onelab.info/gmsh/fem/issues

#include "SkinLayer.h"

#include "Domain.h"
#include "Message.h"
#include "Post.h"
#include "ScalarFunction.h"

namespace gmshfem::domain
{


  SkinLayer::SkinLayer(const int dim, const int entity, const Domain &object, const Domain &tool, const std::map< std::pair< int, int >, JacobiansModificator * > &jacobiansModificators) :
    GeometricObject(), _entity(entity), _object(new Domain(object)), _tool(new Domain(tool))
  {
    _entities.insert(std::make_pair(dim, entity));
    for(auto it = jacobiansModificators.begin(); it != jacobiansModificators.end(); ++it) {
      _jacobiansModificators.insert(std::pair< std::pair< int, int >, JacobiansModificator * >(it->first, it->second->copy()));
    }
  }

  SkinLayer::SkinLayer() :
    GeometricObject(), _entity(-1), _object(nullptr), _tool(nullptr)
  {
  }

  SkinLayer::~SkinLayer()
  {
    if(_object != nullptr) {
      delete _object;
    }
    if(_tool != nullptr) {
      delete _tool;
    }
  }

  void SkinLayer::printDebug() const
  {
    msg::debug << "Skin layer has:" << msg::endl;
    msg::debug << " * elementary entities (" << _entities.size() << "):" << msg::endl;
    for(auto it = _entities.begin(); it != _entities.end(); ++it) {
      msg::debug << "  - " << it->first << " / " << it->second << msg::endl;
    }
  }

  void SkinLayer::saveDebug(const std::string &name) const
  {
    post::save(function::ScalarFunction< double >(1.), *this, "Skin layer debug" + (name == "" ? "" : " : " + name), "pos");
  }


} // namespace gmshfem::domain
