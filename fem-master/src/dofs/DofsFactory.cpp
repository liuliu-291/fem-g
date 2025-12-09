// GmshFEM - Copyright (C) 2019-2022, A. Royer, E. Béchet, C. Geuzaine, Université de Liège
//
// See the LICENSE.txt file for license information. Please report all
// issues on https://gitlab.onelab.info/gmsh/fem/issues

#include "DofsFactory.h"

#include "Dof.h"
#include "DofsManager.h"
#include "Domain.h"
#include "Exception.h"
#include "FieldInterface.h"
#include "Formulation.h"
#include "MathObject.h"
#include "Message.h"
#include "OmpInterface.h"
#include "Term.h"
#include "instantiate.h"
#include "numa.h"
#include "scalar.h"

#include <gmsh.h>
#include <unordered_map>

namespace gmshfem::dofs
{


  template< class T_Scalar >
  DofsFactory< T_Scalar >::DofsFactory(DofsManager< T_Scalar > *const dofsM) :
    _dofsM(dofsM), _nbrFixedDofs(0), _nbrFixedGlobalDofs(0), _nbrUnknownDofs(0), _nbrBubbleUnknownDofs(0), _nbrUnknownGlobalDofs(0), _nbrLinkedDofs(0), _nbrBubbleLinkedDofs(0)
  {
  }

  // FixedDofs

  template< class T_Scalar >
  template< field::Form T_Form >
  void DofsFactory< T_Scalar >::_generateFixedDofs(field::Field< T_Scalar, T_Form > *const field, const domain::Domain &domain, const function::Function< T_Scalar, field::DegreeOfForm< T_Form >::value > &function)
  {
    msg::debug << "Generate " << field::NameOfForm< T_Form >::value << " fixed dofs of field '" << field->name() << "' on domain:" << msg::endl;
    domain.printDebug();
    
    if constexpr(T_Form == field::Form::Form0) {
      if(field->getFunctionSpace()->getGmshFemOrientationName() == "None") {
        for(auto itEntity = domain.cbegin(); itEntity != domain.cend(); ++itEntity) {
          if(function.isConstant(*itEntity)) {
            typename MathObject< T_Scalar, Degree::Degree0 >::Object value;
            function.evaluate(value, 0., 0., 0., *itEntity);

            std::vector< int > elementTypes;
            gmsh::model::mesh::getElementTypes(elementTypes, itEntity->first, itEntity->second);
            for(auto typeIndex = 0ULL; typeIndex < elementTypes.size(); ++typeIndex) {
              std::vector< int > typeKeys;
              std::vector< unsigned long long > entityKeys;
              std::vector< scalar::Precision< T_Scalar >, numa::allocator< scalar::Precision< T_Scalar > > > coord;
              std::vector< scalar::Precision< T_Scalar > > coordNoNuma;
              field->getFunctionSpace()->getKeys(true, typeKeys, entityKeys, coordNoNuma, elementTypes[typeIndex], *itEntity);
              numa::copy(coord, coordNoNuma);

              field->reserve(field->numberOfDofs() + typeKeys.size());

              for(auto j = 0ULL; j < typeKeys.size(); ++j) {
                FixedDof *dof = new(field->getNextFixedDofMemoryPlace()) FixedDof(typeKeys[j] + GMSHFEM_DOF_FIELD_OFFSET * field->tag(), entityKeys[j]);
                if(field->setValue(dof, value)) {
                  dof->coordinates(coord[3 * j + 0], coord[3 * j + 1], coord[3 * j + 2]);
                  dof->numDof(++_nbrFixedDofs);
                }
              }
            }
          }
          else {
            std::vector< typename MathObject< T_Scalar, Degree::Degree0 >::Object, numa::allocator< typename MathObject< T_Scalar, Degree::Degree0 >::Object > > values;
            std::vector< int > elementTypes;
            gmsh::model::mesh::getElementTypes(elementTypes, itEntity->first, itEntity->second);
            for(auto typeIndex = 0ULL; typeIndex < elementTypes.size(); ++typeIndex) {
              std::vector< int > typeKeys;
              std::vector< unsigned long long > entityKeys;
              std::vector< scalar::Precision< T_Scalar > > gaussPoints;
              std::vector< scalar::Precision< T_Scalar >, numa::allocator< scalar::Precision< T_Scalar > > > coord;
              std::vector< scalar::Precision< T_Scalar > > coordNoNuma;
              field->getFunctionSpace()->getKeys(true, typeKeys, entityKeys, coordNoNuma, elementTypes[typeIndex], *itEntity);
              numa::copy(coord, coordNoNuma);
#pragma omp parallel num_threads(omp::getMaxThreads())
              {
                function.evaluate(values, coord, gaussPoints, elementTypes[typeIndex], *itEntity);
              }

              field->reserve(field->numberOfDofs() + typeKeys.size());

              for(auto j = 0ULL; j < typeKeys.size(); ++j) {
                FixedDof *dof = new(field->getNextFixedDofMemoryPlace()) FixedDof(typeKeys[j] + GMSHFEM_DOF_FIELD_OFFSET * field->tag(), entityKeys[j]);
                if(field->setValue(dof, values[j])) {
                  dof->coordinates(coord[3 * j + 0], coord[3 * j + 1], coord[3 * j + 2]);
                  dof->numDof(++_nbrFixedDofs);
                }
              }
            }
          }
        }
        _dofsM->addField(field);
        return;
      }
    }
 
    field::Field< T_Scalar, T_Form > fieldCp("fieldCp", domain, field->getFunctionSpace()->type(), field->getFunctionSpace()->order());
    const std::string gauss = "Gauss" + std::to_string(2 * (field->getFunctionSpace()->order() + 1));
    problem::Formulation< T_Scalar > formulation("boundary");

    formulation.integral(equation::dof(fieldCp), equation::tf(fieldCp), domain, gauss);
    formulation.integral(-function, equation::tf(fieldCp), domain, gauss);

    const int pastVerbose = common::Options::instance()->verbose;
    common::Options::instance()->verbose = 2;
    formulation.pre(problem::DofsSort::Algorithm::None);
    formulation.assemble();
    formulation.solve();
    common::Options::instance()->verbose = pastVerbose;

    for(auto it = fieldCp.begin(); it != fieldCp.end(); ++it) {
      const unsigned long long keyFirst = it->first->numType() - GMSHFEM_DOF_FIELD_OFFSET * fieldCp.tag();
      const unsigned long long keySecond = it->first->entity();
      const T_Scalar value = it->second;
      const scalar::Precision< T_Scalar > coord[3]{it->first->x(), it->first->y(), it->first->z()};

      FixedDof *dof = new(field->getNextFixedDofMemoryPlace()) FixedDof(keyFirst + GMSHFEM_DOF_FIELD_OFFSET * field->tag(), keySecond);
      if(field->setValue(dof, value)) {
        dof->coordinates(coord[0], coord[1], coord[2]);
        dof->numDof(++_nbrFixedDofs);
      }
    }
    _dofsM->addField(field);
  }

  template< class T_Scalar >
  template< field::Form T_Form, unsigned int T_NumFields >
  void DofsFactory< T_Scalar >::_generateCompoundFixedDofs(field::CompoundField< T_Scalar, T_Form, T_NumFields > *const field, const domain::Domain &domain, const function::Function< T_Scalar, field::DegreeOfCompoundForm< T_Form >::value > &function)
  {
    msg::debug << "Generate " << field::NameOfForm< T_Form >::value << " fixed dofs of field '" << field->name() << "' on domain:" << msg::endl;
    domain.printDebug();
    
    if constexpr(T_Form == field::Form::Form0) {
      if(field->getFunctionSpace()->getGmshFemOrientationName() == "None") {
        for(auto itEntity = domain.cbegin(); itEntity != domain.cend(); ++itEntity) {
          if(function.isConstant(*itEntity)) {
            typename MathObject< T_Scalar, Degree::Degree1 >::Object value;
            function.evaluate(value, 0., 0., 0., *itEntity);

            std::vector< int > elementTypes;
            gmsh::model::mesh::getElementTypes(elementTypes, itEntity->first, itEntity->second);
            for(auto typeIndex = 0ULL; typeIndex < elementTypes.size(); ++typeIndex) {
              std::vector< int > typeKeys;
              std::vector< unsigned long long > entityKeys;
              std::vector< scalar::Precision< T_Scalar >, numa::allocator< scalar::Precision< T_Scalar > > > coord;
              std::vector< scalar::Precision< T_Scalar > > coordNoNuma;
              field->getFunctionSpace()->getKeys(true, typeKeys, entityKeys, coordNoNuma, elementTypes[typeIndex], *itEntity);
              numa::copy(coord, coordNoNuma);

              field->reserve(field->numberOfDofs() + typeKeys.size());

              for(auto j = 0ULL; j < typeKeys.size(); ++j) {
                for(unsigned int k = 0; k < field->multiplicity(); ++k) {
                  FixedDof *dof = new(field->getNextFixedDofMemoryPlace()) FixedDof(typeKeys[j] + GMSHFEM_DOF_FIELD_OFFSET * (field->tag() + k), entityKeys[j]);
                  if(field->setValue(dof, value(k))) {
                    dof->coordinates(coord[3 * j + 0], coord[3 * j + 1], coord[3 * j + 2]);
                    dof->numDof(++_nbrFixedDofs);
                  }
                }
              }
            }
          }
          else {
            std::vector< typename MathObject< T_Scalar, Degree::Degree1 >::Object, numa::allocator< typename MathObject< T_Scalar, Degree::Degree1 >::Object > > values;
            std::vector< int > elementTypes;
            gmsh::model::mesh::getElementTypes(elementTypes, itEntity->first, itEntity->second);
            for(auto typeIndex = 0ULL; typeIndex < elementTypes.size(); ++typeIndex) {
              std::vector< int > typeKeys;
              std::vector< unsigned long long > entityKeys;
              std::vector< scalar::Precision< T_Scalar > > gaussPoints;
              std::vector< double > gmshNodesCoord, gmshNodesParametricCoord;
              std::vector< std::size_t > gmshNodesTags;
              gmsh::model::mesh::getNodesByElementType(elementTypes[typeIndex], gmshNodesTags, gmshNodesCoord, gmshNodesParametricCoord, itEntity->second, false);
              scalar::move(gaussPoints, gmshNodesCoord);

              std::vector< scalar::Precision< T_Scalar >, numa::allocator< scalar::Precision< T_Scalar > > > coord;
              std::vector< scalar::Precision< T_Scalar > > coordNoNuma;
              field->getFunctionSpace()->getKeys(true, typeKeys, entityKeys, coordNoNuma, elementTypes[typeIndex], *itEntity);
              numa::copy(coord, coordNoNuma);
  #pragma omp parallel num_threads(omp::getMaxThreads())
              {
                function.evaluate(values, coord, gaussPoints, elementTypes[typeIndex], *itEntity);
              }

              field->reserve(field->numberOfDofs() + typeKeys.size());

              for(auto j = 0ULL; j < typeKeys.size(); ++j) {
                for(unsigned int k = 0; k < field->multiplicity(); ++k) {
                  FixedDof *dof = new(field->getNextFixedDofMemoryPlace()) FixedDof(typeKeys[j] + GMSHFEM_DOF_FIELD_OFFSET * (field->tag() + k), entityKeys[j]);
                  if(field->setValue(dof, values[j](k))) {
                    dof->coordinates(coord[3 * j + 0], coord[3 * j + 1], coord[3 * j + 2]);
                    dof->numDof(++_nbrFixedDofs);
                  }
                }
              }
            }
          }
        }
        _dofsM->addField(field);
        return;
      }
    }
    
    field::CompoundField< T_Scalar, T_Form, T_NumFields > fieldCp("fieldCp", domain, field->getFunctionSpace()->type(), field->getFunctionSpace()->order());
    const std::string gauss = "Gauss" + std::to_string(2 * field->getFunctionSpace()->order());
    problem::Formulation< T_Scalar > formulation("boundary");

    formulation.integral(equation::dof(fieldCp), equation::tf(fieldCp), domain, gauss);
    formulation.integral(-function, equation::tf(fieldCp), domain, gauss);

    const int pastVerbose = common::Options::instance()->verbose;
    common::Options::instance()->verbose = 2;
    formulation.pre(problem::DofsSort::Algorithm::None);
    formulation.assemble();
    formulation.solve();
    common::Options::instance()->verbose = pastVerbose;

    for(auto it = fieldCp.begin(); it != fieldCp.end(); ++it) {
      const unsigned long long keyFirst = it->first->numType() - GMSHFEM_DOF_FIELD_OFFSET * fieldCp.tag();
      const unsigned long long keySecond = it->first->entity();
      const T_Scalar value = it->second;
      const scalar::Precision< T_Scalar > coord[3]{it->first->x(), it->first->y(), it->first->z()};

      FixedDof *dof = new(field->getNextFixedDofMemoryPlace()) FixedDof(keyFirst + GMSHFEM_DOF_FIELD_OFFSET * field->tag(), keySecond);
      if(field->setValue(dof, value)) {
        dof->coordinates(coord[0], coord[1], coord[2]);
        dof->numDof(++_nbrFixedDofs);
      }
    }
    _dofsM->addField(field);
  }
  
  template< class T_Scalar >
  template< field::Form T_Form, unsigned int T_NumFields, unsigned int T_Component, class >
  void DofsFactory< T_Scalar >::_generateCompoundFixedDofsOnComponent(field::CompoundField< T_Scalar, T_Form, T_NumFields > *const field, const domain::Domain &domain, const function::Function< T_Scalar, field::DegreeOfForm< T_Form >::value > &function)
  {
    msg::debug << "Generate " << field::NameOfForm< T_Form >::value << " fixed dofs of the " << T_Component << "-component of field '" << field->name() << "' on domain:" << msg::endl;
    domain.printDebug();
    
    if constexpr(T_Form == field::Form::Form0) {
      if(field->getFunctionSpace()->getGmshFemOrientationName() == "None") {
        for(auto itEntity = domain.cbegin(); itEntity != domain.cend(); ++itEntity) {
          if(function.isConstant(*itEntity)) {
            typename MathObject< T_Scalar, Degree::Degree0 >::Object value;
            function.evaluate(value, 0., 0., 0., *itEntity);

            std::vector< int > elementTypes;
            gmsh::model::mesh::getElementTypes(elementTypes, itEntity->first, itEntity->second);
            for(auto typeIndex = 0ULL; typeIndex < elementTypes.size(); ++typeIndex) {
              std::vector< int > typeKeys;
              std::vector< unsigned long long > entityKeys;
              std::vector< scalar::Precision< T_Scalar >, numa::allocator< scalar::Precision< T_Scalar > > > coord;
              std::vector< scalar::Precision< T_Scalar > > coordNoNuma;
              field->getFunctionSpace()->getKeys(true, typeKeys, entityKeys, coordNoNuma, elementTypes[typeIndex], *itEntity);
              numa::copy(coord, coordNoNuma);

              field->reserve(field->numberOfDofs() + typeKeys.size());

              for(auto j = 0ULL; j < typeKeys.size(); ++j) {
                FixedDof *dof = new(field->getNextFixedDofMemoryPlace()) FixedDof(typeKeys[j] + GMSHFEM_DOF_FIELD_OFFSET * (field->tag() + T_Component), entityKeys[j]);
                if(field->setValue(dof, value)) {
                  dof->coordinates(coord[3 * j + 0], coord[3 * j + 1], coord[3 * j + 2]);
                  dof->numDof(++_nbrFixedDofs);
                }
              }
            }
          }
          else {
            std::vector< typename MathObject< T_Scalar, Degree::Degree0 >::Object, numa::allocator< typename MathObject< T_Scalar, Degree::Degree0 >::Object > > values;
            std::vector< int > elementTypes;
            gmsh::model::mesh::getElementTypes(elementTypes, itEntity->first, itEntity->second);
            for(auto typeIndex = 0ULL; typeIndex < elementTypes.size(); ++typeIndex) {
              std::vector< int > typeKeys;
              std::vector< unsigned long long > entityKeys;
              std::vector< scalar::Precision< T_Scalar > > gaussPoints;
              std::vector< double > gmshNodesCoord, gmshNodesParametricCoord;
              std::vector< std::size_t > gmshNodesTags;
              gmsh::model::mesh::getNodesByElementType(elementTypes[typeIndex], gmshNodesTags, gmshNodesCoord, gmshNodesParametricCoord, itEntity->second, false);
              scalar::move(gaussPoints, gmshNodesCoord);

              std::vector< scalar::Precision< T_Scalar >, numa::allocator< scalar::Precision< T_Scalar > > > coord;
              std::vector< scalar::Precision< T_Scalar > > coordNoNuma;
              field->getFunctionSpace()->getKeys(true, typeKeys, entityKeys, coordNoNuma, elementTypes[typeIndex], *itEntity);
              numa::copy(coord, coordNoNuma);
  #pragma omp parallel num_threads(omp::getMaxThreads())
              {
                function.evaluate(values, coord, gaussPoints, elementTypes[typeIndex], *itEntity);
              }

              field->reserve(field->numberOfDofs() + typeKeys.size());

              for(auto j = 0ULL; j < typeKeys.size(); ++j) {
                FixedDof *dof = new(field->getNextFixedDofMemoryPlace()) FixedDof(typeKeys[j] + GMSHFEM_DOF_FIELD_OFFSET * (field->tag() + T_Component), entityKeys[j]);
                if(field->setValue(dof, values[j])) {
                  dof->coordinates(coord[3 * j + 0], coord[3 * j + 1], coord[3 * j + 2]);
                  dof->numDof(++_nbrFixedDofs);
                }
              }
            }
          }
        }
        _dofsM->addField(field);
        return;
      }
    }
    
    field::Field< T_Scalar, T_Form > fieldCp("fieldCp", domain, field->getFunctionSpace()->type(), field->getFunctionSpace()->order());
    const std::string gauss = "Gauss" + std::to_string(2 * field->getFunctionSpace()->order());
    problem::Formulation< T_Scalar > formulation("boundary");

    formulation.integral(equation::dof(fieldCp), equation::tf(fieldCp), domain, gauss);
    formulation.integral(-function, equation::tf(fieldCp), domain, gauss);

    const int pastVerbose = common::Options::instance()->verbose;
    common::Options::instance()->verbose = 2;
    formulation.pre(problem::DofsSort::Algorithm::None);
    formulation.assemble();
    formulation.solve();
    common::Options::instance()->verbose = pastVerbose;

    for(auto it = fieldCp.begin(); it != fieldCp.end(); ++it) {
      const unsigned long long keyFirst = it->first->numType() - GMSHFEM_DOF_FIELD_OFFSET * fieldCp.tag();
      const unsigned long long keySecond = it->first->entity();
      const T_Scalar value = it->second;
      const scalar::Precision< T_Scalar > coord[3]{it->first->x(), it->first->y(), it->first->z()};

      FixedDof *dof = new(field->getNextFixedDofMemoryPlace()) FixedDof(keyFirst + GMSHFEM_DOF_FIELD_OFFSET * (field->tag() + T_Component), keySecond);
      if(field->setValue(dof, value)) {
        dof->coordinates(coord[0], coord[1], coord[2]);
        dof->numDof(++_nbrFixedDofs);
      }
    }
    _dofsM->addField(field);
  }

  template< class T_Scalar >
  void DofsFactory< T_Scalar >::_generateFixedGlobalDofs(field::GlobalQuantity< T_Scalar > *const globalQuantity, const domain::Domain &domain, const bool primal)
  {
    msg::debug << "Generate fixed global dofs of gloabl quantity '" << globalQuantity->name() << "' on domain:" << msg::endl;
    domain.printDebug();
    std::vector< int > typeKeys;
    std::vector< unsigned long long > entityKeys;
    std::vector< scalar::Precision< T_Scalar > > coord;
    gmsh::vectorpair infoKey;
    bool hasZero = false;
    bool hasCreatedDof = false;
    field::FieldInterface< T_Scalar > *const field = (primal ? globalQuantity->getAssociatedPrimalField() : globalQuantity->getAssociatedDualField());
    for(auto itEntity = domain.cbegin(); itEntity != domain.cend(); ++itEntity) {
      std::vector< int > elementTypes;
      gmsh::model::mesh::getElementTypes(elementTypes, itEntity->first, itEntity->second);
      for(auto typeIndex = 0ULL; typeIndex < elementTypes.size(); ++typeIndex) {
        typeKeys.clear();
        entityKeys.clear();
        coord.clear();
        field->getFunctionSpace()->getKeys(true, typeKeys, entityKeys, coord, elementTypes[typeIndex], *itEntity);
        field->getFunctionSpace()->getKeyInformation(typeKeys, entityKeys, elementTypes[typeIndex], infoKey);

        field->reserve(field->numberOfDofs() + typeKeys.size());

        for(auto j = 0ULL; j < typeKeys.size(); ++j) {
          if(infoKey[j].second == 1 || infoKey[j].second == -1) {
            FixedDof *dof = new(field->getNextFixedDofMemoryPlace()) FixedDof(typeKeys[j] + GMSHFEM_DOF_FIELD_OFFSET * field->tag(), entityKeys[j]);
            if(field->setValue(dof, (primal ? globalQuantity->getPrimalValue() : globalQuantity->getDualValue()))) {
              hasCreatedDof = true;
              dof->coordinates(coord[3 * j + 0], coord[3 * j + 1], coord[3 * j + 2]);
              dof->numDof(_nbrFixedDofs + 1);
              if(primal) {
                globalQuantity->setAssociatedPrimalDof(dof);
              }
              else {
                globalQuantity->setAssociatedDualDof(dof);
              }
            }
          }
          else {
            FixedDof *dof = new(field->getNextFixedDofMemoryPlace()) FixedDof(typeKeys[j] + GMSHFEM_DOF_FIELD_OFFSET * field->tag(), entityKeys[j]);
            if(field->setValue(dof, 0.)) {
              hasZero = true;
              hasCreatedDof = true;
              dof->coordinates(coord[3 * j + 0], coord[3 * j + 1], coord[3 * j + 2]);
              dof->numDof(_nbrFixedDofs + 2);
            }
          }
        }
      }
    }
    _dofsM->addField(field);
    if(hasCreatedDof) {
      _nbrFixedDofs += (hasZero ? 2 : 1);
      ++_nbrFixedGlobalDofs;
    }
  }

  // UnknownDofs

  template< class T_Scalar >
  void DofsFactory< T_Scalar >::_generateUnknownDofs(field::FieldInterface< T_Scalar > *const field, const domain::Domain &domain, const unsigned int maxDim)
  {
    msg::debug << "Generate unknown dofs of field '" << field->name() << "' on domain:" << msg::endl;
    domain.printDebug();
    std::vector< int > typeKeys;
    std::vector< unsigned long long > entityKeys;
    gmsh::vectorpair infoKey;
    std::vector< scalar::Precision< T_Scalar > > coord;
    for(auto itEntity = domain.cbegin(); itEntity != domain.cend(); ++itEntity) {
      std::vector< int > elementTypes;
      gmsh::model::mesh::getElementTypes(elementTypes, itEntity->first, itEntity->second);
      for(auto typeIndex = 0ULL; typeIndex < elementTypes.size(); ++typeIndex) {
        typeKeys.clear();
        entityKeys.clear();
        infoKey.clear();
        coord.clear();
        field->getFunctionSpace()->getKeys(true, typeKeys, entityKeys, coord, elementTypes[typeIndex], *itEntity);
        if(static_cast< unsigned int >(itEntity->first) == maxDim) {
          field->getFunctionSpace()->getKeyInformation(typeKeys, entityKeys, elementTypes[typeIndex], infoKey);
        }

        field->reserve(field->numberOfDofs() + typeKeys.size() * field->multiplicity());

        if(static_cast< unsigned int >(itEntity->first) == maxDim) {
          for(auto j = 0ULL; j < typeKeys.size(); ++j) {
            for(unsigned int k = 0; k < field->multiplicity(); ++k) {
              UnknownDof *dof = nullptr;
              const bool isBubble = (static_cast< unsigned int >(infoKey[j].first) == maxDim);
              if(isBubble) {
                dof = new(field->getNextUnknownDofMemoryPlace()) UnknownBubbleDof(typeKeys[j] + GMSHFEM_DOF_FIELD_OFFSET * (field->tag() + k), entityKeys[j]);
              }
              else {
                dof = new(field->getNextUnknownDofMemoryPlace()) UnknownDof(typeKeys[j] + GMSHFEM_DOF_FIELD_OFFSET * (field->tag() + k), entityKeys[j]);
              }
              if(field->setDof(dof)) {
                dof->coordinates(coord[3 * j + 0], coord[3 * j + 1], coord[3 * j + 2]);
                dof->numDof(++_nbrUnknownDofs);
                if(isBubble) ++_nbrBubbleUnknownDofs;
              }
            }
          }
        }
        else {
          for(auto j = 0ULL; j < typeKeys.size(); ++j) {
            for(unsigned int k = 0; k < field->multiplicity(); ++k) {
              UnknownDof *dof = new(field->getNextUnknownDofMemoryPlace()) UnknownDof(typeKeys[j] + GMSHFEM_DOF_FIELD_OFFSET * (field->tag() + k), entityKeys[j]);
              if(field->setDof(dof)) {
                dof->coordinates(coord[3 * j + 0], coord[3 * j + 1], coord[3 * j + 2]);
                dof->numDof(++_nbrUnknownDofs);
              }
            }
          }
        }
      }
    }
    _dofsM->addField(field);
  }

  template< class T_Scalar >
  void DofsFactory< T_Scalar >::_generateUnknownGlobalDofs(field::GlobalQuantity< T_Scalar > *const globalQuantity, const domain::Domain &domain, const bool primal)
  {
    msg::debug << "Generate unknown global dofs of global quantity '" << globalQuantity->name() << "' on domain:" << msg::endl;
    domain.printDebug();
    std::vector< int > typeKeys;
    std::vector< unsigned long long > entityKeys;
    std::vector< scalar::Precision< T_Scalar > > coord;
    gmsh::vectorpair infoKey;
    bool hasZero = false;
    bool hasCreatedDof = false;
    field::FieldInterface< T_Scalar > *const field = (primal ? globalQuantity->getAssociatedPrimalField() : globalQuantity->getAssociatedDualField());
    for(auto itEntity = domain.cbegin(); itEntity != domain.cend(); ++itEntity) {
      std::vector< int > elementTypes;
      gmsh::model::mesh::getElementTypes(elementTypes, itEntity->first, itEntity->second);
      for(auto typeIndex = 0ULL; typeIndex < elementTypes.size(); ++typeIndex) {
        typeKeys.clear();
        entityKeys.clear();
        coord.clear();
        field->getFunctionSpace()->getKeys(true, typeKeys, entityKeys, coord, elementTypes[typeIndex], *itEntity);
        field->getFunctionSpace()->getKeyInformation(typeKeys, entityKeys, elementTypes[typeIndex], infoKey);

        field->reserve(field->numberOfDofs() + typeKeys.size());

        for(auto j = 0ULL; j < typeKeys.size(); ++j) {
          if(infoKey[j].second == 1 || infoKey[j].second == -1) {
            UnknownGlobalDof *dof = new(field->getNextUnknownDofMemoryPlace()) UnknownGlobalDof(typeKeys[j] + GMSHFEM_DOF_FIELD_OFFSET * field->tag(), entityKeys[j]);
            if(field->setDof(dof)) {
              hasCreatedDof = true;
              dof->coordinates(coord[3 * j + 0], coord[3 * j + 1], coord[3 * j + 2]);
              dof->numDof(_nbrUnknownDofs + 1);
              if(primal) {
                globalQuantity->setAssociatedPrimalDof(dof);
              }
              else {
                globalQuantity->setAssociatedDualDof(dof);
              }
            }
          }
          else {
            FixedDof *dof = new(field->getNextFixedDofMemoryPlace()) FixedDof(typeKeys[j] + GMSHFEM_DOF_FIELD_OFFSET * field->tag(), entityKeys[j]);
            if(field->setValue(dof, 0.)) {
              hasZero = true;
              hasCreatedDof = true;
              dof->coordinates(coord[3 * j + 0], coord[3 * j + 1], coord[3 * j + 2]);
              dof->numDof(_nbrFixedDofs + 1);
            }
          }
        }
      }
    }
    _dofsM->addField(field);
    if(hasCreatedDof) {
      _nbrUnknownDofs++;
      _nbrUnknownGlobalDofs++;
      if(hasZero) {
        _nbrFixedDofs++;
        _nbrFixedGlobalDofs++;
      }
    }
  }

  // LinkedDofs

  template< class T_Scalar >
  void DofsFactory< T_Scalar >::_generateLinkedDofs(field::FieldInterface< T_Scalar > *const field, const domain::Domain &master, const domain::Domain &slave, const T_Scalar &coefficient, const unsigned int maxDim)
  {
    msg::debug << "Generate linked dofs of field '" << field->name() << "' on domain:" << msg::endl;
    slave.printDebug();
    std::vector< int > slaveTypeKeys, masterTypeKeys;
    std::vector< unsigned long long > slaveEntityKeys, masterEntityKeys;
    gmsh::vectorpair infoKey;
    std::vector< scalar::Precision< T_Scalar > > slaveCoord, masterCoord;
    for(auto itEntity = slave.cbegin(); itEntity != slave.cend(); ++itEntity) {
      std::vector< int > elementTypes;
      gmsh::model::mesh::getElementTypes(elementTypes, itEntity->first, itEntity->second);
      for(auto typeIndex = 0ULL; typeIndex < elementTypes.size(); ++typeIndex) {
        slaveTypeKeys.clear();
        masterTypeKeys.clear();
        slaveEntityKeys.clear();
        masterEntityKeys.clear();
        infoKey.clear();
        slaveCoord.clear();
        masterCoord.clear();
        bool isPeriodic = field->getFunctionSpace()->getPeriodicKeys(true, slaveTypeKeys, slaveEntityKeys, slaveCoord, masterTypeKeys, masterEntityKeys, masterCoord, elementTypes[typeIndex], *itEntity);
        if(!isPeriodic) {
          msg::error << "Entity (" << itEntity->first << ", " << itEntity->second << ") is not periodic" << msg::endl;
          typeIndex = elementTypes.size();
          continue;
        }
        if(static_cast< unsigned int >(itEntity->first) == maxDim) {
          field->getFunctionSpace()->getKeyInformation(slaveTypeKeys, slaveEntityKeys, elementTypes[typeIndex], infoKey);
        }

        field->reserve(field->numberOfDofs() + slaveEntityKeys.size() + masterEntityKeys.size());

        if(static_cast< unsigned int >(itEntity->first) == maxDim) {
          for(auto j = 0ULL; j < masterTypeKeys.size(); ++j) {
            UnknownDof *dof = nullptr;
            const bool isBubble = (static_cast< unsigned int >(infoKey[j].first) == maxDim);
            if(isBubble) {
              dof = new(field->getNextUnknownDofMemoryPlace()) UnknownBubbleDof(masterTypeKeys[j] + GMSHFEM_DOF_FIELD_OFFSET * field->tag(), masterEntityKeys[j]);
            }
            else {
              dof = new(field->getNextUnknownDofMemoryPlace()) UnknownDof(masterTypeKeys[j] + GMSHFEM_DOF_FIELD_OFFSET * field->tag(), masterEntityKeys[j]);
            }
            if(field->setDof(dof)) {
              dof->coordinates(masterCoord[3 * j + 0], masterCoord[3 * j + 1], masterCoord[3 * j + 2]);
              dof->numDof(++_nbrUnknownDofs);
              if(isBubble) ++_nbrBubbleUnknownDofs;

              LinkedDof *ldof = nullptr;
              if(isBubble) {
                ldof = new(field->getNextLinkedDofMemoryPlace()) LinkedBubbleDof(slaveTypeKeys[j] + GMSHFEM_DOF_FIELD_OFFSET * field->tag(), slaveEntityKeys[j], dof);
              }
              else {
                ldof = new(field->getNextLinkedDofMemoryPlace()) LinkedDof(slaveTypeKeys[j] + GMSHFEM_DOF_FIELD_OFFSET * field->tag(), slaveEntityKeys[j], dof);
              }
              if(field->setDof(ldof, coefficient)) {
                ldof->coordinates(slaveCoord[3 * j + 0], slaveCoord[3 * j + 1], slaveCoord[3 * j + 2]);
                ldof->numDof(++_nbrLinkedDofs);
                if(isBubble) ++_nbrBubbleLinkedDofs;
              }
            }
            else {
              const dofs::Dof *masterDof = field->searchDof(masterTypeKeys[j], masterEntityKeys[j]);
              if(masterDof == nullptr) {
                continue;
              }

              if(masterDof->type() == dofs::Type::Fixed) {
                FixedDof *dof = new(field->getNextFixedDofMemoryPlace()) FixedDof(slaveTypeKeys[j] + GMSHFEM_DOF_FIELD_OFFSET * field->tag(), slaveEntityKeys[j]);
                const T_Scalar value = coefficient * field->getValue(masterDof);
                if(field->setValue(dof, value)) {
                  dof->coordinates(slaveCoord[0], slaveCoord[1], slaveCoord[2]);
                  dof->numDof(++_nbrFixedDofs);
                }
              }
            }
          }
        }
        else {
          for(auto j = 0ULL; j < masterTypeKeys.size(); ++j) {
            UnknownDof *dof = new(field->getNextUnknownDofMemoryPlace()) UnknownDof(masterTypeKeys[j] + GMSHFEM_DOF_FIELD_OFFSET * field->tag(), masterEntityKeys[j]);
            if(field->setDof(dof)) {
              dof->coordinates(masterCoord[3 * j + 0], masterCoord[3 * j + 1], masterCoord[3 * j + 2]);
              dof->numDof(++_nbrUnknownDofs);

              LinkedDof *ldof = new(field->getNextLinkedDofMemoryPlace()) LinkedDof(slaveTypeKeys[j] + GMSHFEM_DOF_FIELD_OFFSET * field->tag(), slaveEntityKeys[j], dof);
              if(field->setDof(ldof, coefficient)) {
                ldof->coordinates(slaveCoord[3 * j + 0], slaveCoord[3 * j + 1], slaveCoord[3 * j + 2]);
                ldof->numDof(++_nbrLinkedDofs);
              }
            }
            else {
              const dofs::Dof *masterDof = field->searchDof(masterTypeKeys[j], masterEntityKeys[j]);
              if(masterDof == nullptr) {
                continue;
              }

              if(masterDof->type() == dofs::Type::Fixed) {
                FixedDof *dof = new(field->getNextFixedDofMemoryPlace()) FixedDof(slaveTypeKeys[j] + GMSHFEM_DOF_FIELD_OFFSET * field->tag(), slaveEntityKeys[j]);
                const T_Scalar value = coefficient * field->getValue(masterDof);
                if(field->setValue(dof, value)) {
                  dof->coordinates(slaveCoord[0], slaveCoord[1], slaveCoord[2]);
                  dof->numDof(++_nbrFixedDofs);
                }
              }
            }
          }
        }
      }
    }
  }

  template< class T_Scalar >
  struct FlattenPeriodicConstraint {
    domain::Domain _master;
    domain::Domain _slave;
    T_Scalar _bondageCoefficient;

    FlattenPeriodicConstraint(const domain::Domain &master, const domain::Domain &slave, const T_Scalar &bondageCoefficient) :
      _master(master), _slave(slave), _bondageCoefficient(bondageCoefficient)
    {
    }
  };

  template< class T_Scalar >
  static std::vector< FlattenPeriodicConstraint< T_Scalar > > s_flattenRecursivePeriodicConstraints(std::vector< field::PeriodicConstraint< T_Scalar > > &constraints)
  {
    std::vector< FlattenPeriodicConstraint< T_Scalar > > flattenConstraints;
    flattenConstraints.reserve(constraints.size());
    for(auto i = 0ULL; i < constraints.size(); ++i) {
      flattenConstraints.push_back(FlattenPeriodicConstraint< T_Scalar >(constraints[i]._link.master(), constraints[i]._link.slave(), constraints[i]._bondageCoefficient));
    }

    for(auto i = 0ULL; i < flattenConstraints.size(); ++i) {
      for(auto j = 0ULL; j < flattenConstraints.size(); ++j) {
        if(i == j) continue;
        if(!(flattenConstraints[i]._master & flattenConstraints[j]._slave).isEmpty()) {
          const domain::Domain overlap = flattenConstraints[i]._master & flattenConstraints[j]._slave;
          // I condition
          domain::Domain freeI, linkedI;
          for(auto it = flattenConstraints[i]._slave.cbegin(); it != flattenConstraints[i]._slave.cend(); ++it) {
            int tagMaster;
            std::vector< std::size_t > nodeTags;
            std::vector< std::size_t > nodeTagsMaster;
            std::vector< double > affineTransform;
            gmsh::model::mesh::getPeriodicNodes(it->first, it->second, tagMaster, nodeTags, nodeTagsMaster, affineTransform, false);
            const domain::Domain parent = domain::Domain(it->first, tagMaster);
            if((parent & overlap) == domain::Domain()) {
              freeI |= parent;
            }
            else {
              linkedI |= parent;
            }
          }
          // J condition
          domain::Domain freeJ, linkedJ;
          for(auto it = flattenConstraints[j]._slave.cbegin(); it != flattenConstraints[j]._slave.cend(); ++it) {
            int tagMaster;
            std::vector< std::size_t > nodeTags;
            std::vector< std::size_t > nodeTagsMaster;
            std::vector< double > affineTransform;
            gmsh::model::mesh::getPeriodicNodes(it->first, it->second, tagMaster, nodeTags, nodeTagsMaster, affineTransform, false);
            const domain::Domain subset = domain::Domain(it->first, it->second);
            const domain::Domain parent = domain::Domain(it->first, tagMaster);
            if((subset & overlap) == domain::Domain()) {
              linkedJ |= parent;
            }
            else {
              freeJ |= parent;
            }
          }

          {
            // free I link side
            FlattenPeriodicConstraint< T_Scalar > constraint(flattenConstraints[i]._master ^ overlap, freeI, flattenConstraints[i]._bondageCoefficient);
            flattenConstraints.push_back(constraint);
          }
          {
            // chained
            FlattenPeriodicConstraint< T_Scalar > constraint(linkedJ, linkedI, flattenConstraints[i]._bondageCoefficient * flattenConstraints[j]._bondageCoefficient);
            flattenConstraints.push_back(constraint);
          }
          {
            // free J link side
            FlattenPeriodicConstraint< T_Scalar > constraint(freeJ, flattenConstraints[j]._slave ^ overlap, flattenConstraints[j]._bondageCoefficient);
            flattenConstraints.push_back(constraint);
          }
          flattenConstraints.erase(flattenConstraints.begin() + i);
          flattenConstraints.erase(flattenConstraints.begin() + j + (j > i ? -1 : 0));
          --i;
          --j;
        }
      }
    }

    return flattenConstraints;
  }

  template< class T_Scalar >
  void DofsFactory< T_Scalar >::generate(const std::vector< term::Term< T_Scalar > * > &terms, const std::vector< field::FieldInterface< T_Scalar > * > &fields)
  {
    std::vector< std::pair< int, int > > entities;
    gmsh::model::getEntities(entities);

    unsigned int maxDim = 0;
    for(auto it = fields.begin(); it != fields.end(); ++it) {
      if(maxDim < (*it)->domain().maxDim()) {
        maxDim = (*it)->domain().maxDim();
      }
    }

    // Fixed Dofs
    try {
      // global quantities
      for(auto it = fields.begin(); it != fields.end(); ++it) {
        for(auto itG = (*it)->firstGlobalQuantity(); itG != (*it)->lastGlobalQuantity(); ++itG) {
          if((*it)->multiplicity() > 1) {
            throw common::Exception("Global quantities are not implemented for 'CompoundField'");
          }
          if(field::globalValueIsStillValid(itG->first)) {
            if(itG->second->isActivated()) {
              if(itG->second->fixedComponent() == field::FixedComponent::Primal) {
                _generateFixedGlobalDofs(itG->second, itG->second->domain(), true);
              }
              else if(itG->second->fixedComponent() == field::FixedComponent::Dual) {
                _generateFixedGlobalDofs(itG->second, itG->second->domain(), false);
              }
            }
          }
          else {
            msg::warning << "An unvalid global quantity is skipped" << msg::endl;
          }
        }
      }

      // normal fixed dofs
      for(auto it = fields.begin(); it != fields.end(); ++it) {
        switch((*it)->form()) {
        case field::Form::Form0:
            GenerateFixedDofs< field::Form::Form0 >::run(*this, *it);
          break;
        case field::Form::Form1:
            GenerateFixedDofs< field::Form::Form1 >::run(*this, *it);
          break;
        case field::Form::Form2:
            GenerateFixedDofs< field::Form::Form2 >::run(*this, *it);
          break;
        case field::Form::Form3:
            GenerateFixedDofs< field::Form::Form3 >::run(*this, *it);
          break;
        default:
          break;
        }
      }
    }
    catch(const std::exception &exc) {
      msg::error << "Unable to create fixed dofs" << msg::endl;
      throw;
    }

    // Linked Dofs
    try {
      for(auto it = fields.begin(); it != fields.end(); ++it) {
        std::vector< field::PeriodicConstraint< T_Scalar > > constraints = (*it)->getPeriodicConstraints();
        std::vector< FlattenPeriodicConstraint< T_Scalar > > flattenConstraints = s_flattenRecursivePeriodicConstraints(constraints);
        for(auto itP = flattenConstraints.begin(); itP != flattenConstraints.end(); ++itP) {
          if((*it)->multiplicity() > 1) {
            throw common::Exception("Periodic conditions are not implemented for 'CompoundField'");
          }
          _generateLinkedDofs(*it, itP->_master, itP->_slave, itP->_bondageCoefficient, maxDim);
        }
      }
    }
    catch(const std::exception &exc) {
      msg::error << "Unable to create linked dofs" << msg::endl;
      throw;
    }

    // Unknown Dofs
    try {
      // global quantities
      for(auto it = fields.begin(); it != fields.end(); ++it) {
        for(auto itG = (*it)->firstGlobalQuantity(); itG != (*it)->lastGlobalQuantity(); ++itG) {
          if((*it)->multiplicity() > 1) {
            throw common::Exception("Global quantities are not implemented for 'CompoundField'");
          }
          if(field::globalValueIsStillValid(itG->first)) {
            if(itG->second->isActivated()) {
              if(itG->second->fixedComponent() == field::FixedComponent::Primal) {
                _generateUnknownGlobalDofs(itG->second, itG->second->domain(), false);
              }
              else if(itG->second->fixedComponent() == field::FixedComponent::Dual) {
                _generateUnknownGlobalDofs(itG->second, itG->second->domain(), true);
              }
            }
          }
          else {
            msg::warning << "An unvalid global quantity is skipped" << msg::endl;
          }
        }
      }

      // normal unknown
      for(auto it = fields.begin(); it != fields.end(); ++it) {
        domain::Domain domain = (*it)->domain();
        _generateUnknownDofs(*it, domain, maxDim);
      }
    }
    catch(const std::exception &exc) {
      msg::error << "Unable to create unknown dofs" << msg::endl;
      throw;
    }
  }

  template< class T_Scalar >
  unsigned long long DofsFactory< T_Scalar >::nbrFixedDofsCreated() const
  {
    return _nbrFixedDofs;
  }

  template< class T_Scalar >
  unsigned long long DofsFactory< T_Scalar >::nbrFixedGlobalDofsCreated() const
  {
    return _nbrFixedGlobalDofs;
  }

  template< class T_Scalar >
  unsigned long long DofsFactory< T_Scalar >::nbrUnknownDofsCreated() const
  {
    return _nbrUnknownDofs;
  }

  template< class T_Scalar >
  unsigned long long DofsFactory< T_Scalar >::nbrBubbleUnknownDofsCreated() const
  {
    return _nbrBubbleUnknownDofs;
  }

  template< class T_Scalar >
  unsigned long long DofsFactory< T_Scalar >::nbrUnknownGlobalDofsCreated() const
  {
    return _nbrUnknownGlobalDofs;
  }

  template< class T_Scalar >
  unsigned long long DofsFactory< T_Scalar >::nbrLinkedDofsCreated() const
  {
    return _nbrLinkedDofs;
  }

  template< class T_Scalar >
  unsigned long long DofsFactory< T_Scalar >::nbrBubbleLinkedDofsCreated() const
  {
    return _nbrBubbleLinkedDofs;
  }
  
  template< class T_Scalar >
  template< field::Form T_Form >
  void DofsFactory< T_Scalar >::GenerateFixedDofs< T_Form >::run(DofsFactory< T_Scalar > &enclose, field::FieldInterface< T_Scalar > *field)
  {
    if((field)->multiplicity() == 1) {
      field::Field< T_Scalar, T_Form > *f = static_cast< field::Field< T_Scalar, T_Form > * >(field);
      auto constraints = f->getConstraints();
      for(auto i = 0ULL; i < constraints.size(); ++i) {
        enclose._generateFixedDofs< T_Form >(f, constraints[i]._domain, constraints[i]._function);
      }
    }
    else if((field)->multiplicity() == 2) {
      field::CompoundField< T_Scalar, T_Form, 2 > *f = static_cast< field::CompoundField< T_Scalar, T_Form, 2 > * >(field);
      // constraints on the 0 component
      auto constraintsOn0 = f->template getConstraintsOnComponent< 0 >();
      for(auto i = 0ULL; i < constraintsOn0.size(); ++i) {
        enclose._generateCompoundFixedDofsOnComponent< T_Form, 2, 0 >(f, constraintsOn0[i]._domain, constraintsOn0[i]._function);
      }
      // constraints on the 1 component
      auto constraintsOn1 = f->template getConstraintsOnComponent< 1 >();
      for(auto i = 0ULL; i < constraintsOn1.size(); ++i) {
        enclose._generateCompoundFixedDofsOnComponent< T_Form, 2, 1 >(f, constraintsOn1[i]._domain, constraintsOn1[i]._function);
      }
      // constraints on the whole comonent
      auto constraints = f->getConstraints();
      for(auto i = 0ULL; i < constraints.size(); ++i) {
        enclose._generateCompoundFixedDofs< T_Form, 2 >(f, constraints[i]._domain, constraints[i]._function);
      }
    }
    else if((field)->multiplicity() == 3) {
      field::CompoundField< T_Scalar, T_Form, 3 > *f = static_cast< field::CompoundField< T_Scalar, T_Form, 3 > * >(field);
      // constraints on the 0 component
      auto constraintsOn0 = f->template getConstraintsOnComponent< 0 >();
      for(auto i = 0ULL; i < constraintsOn0.size(); ++i) {
        enclose._generateCompoundFixedDofsOnComponent< T_Form, 3, 0 >(f, constraintsOn0[i]._domain, constraintsOn0[i]._function);
      }
      // constraints on the 1 component
      auto constraintsOn1 = f->template getConstraintsOnComponent< 1 >();
      for(auto i = 0ULL; i < constraintsOn1.size(); ++i) {
        enclose._generateCompoundFixedDofsOnComponent< T_Form, 3, 1 >(f, constraintsOn1[i]._domain, constraintsOn1[i]._function);
      }
      // constraints on the 2 component
      auto constraintsOn2 = f->template getConstraintsOnComponent< 2 >();
      for(auto i = 0ULL; i < constraintsOn2.size(); ++i) {
        enclose._generateCompoundFixedDofsOnComponent< T_Form, 3, 2 >(f, constraintsOn2[i]._domain, constraintsOn2[i]._function);
      }
      // constraints on the whole comonent
      auto constraints = f->getConstraints();
      for(auto i = 0ULL; i < constraints.size(); ++i) {
        enclose._generateCompoundFixedDofs< T_Form, 3 >(f, constraints[i]._domain, constraints[i]._function);
      }
    }
  }

  INSTANTIATE_CLASS(DofsFactory, 4, TEMPLATE_ARGS(std::complex< double >, std::complex< float >, double, float))


} // namespace gmshfem::dofs
