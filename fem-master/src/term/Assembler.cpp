// GmshFEM - Copyright (C) 2019-2022, A. Royer, E. Béchet, C. Geuzaine, Université de Liège
//
// See the LICENSE.txt file for license information. Please report all
// issues on https://gitlab.onelab.info/gmsh/fem/issues

#include "Assembler.h"

#include "ElementBucket.h"
#include "Exception.h"
#include "FieldInterface.h"
#include "Function.h"
#include "FunctionSpaceInterface.h"
#include "IndiceBucket.h"
#include "MatrixFactory.h"
#include "Message.h"
#include "VectorFactory.h"
#include "instantiate.h"
#include "numa.h"

#include <unordered_map>
#include <vector>

namespace gmshfem::term
{


  template< class T_Scalar >
  Assembler< T_Scalar >::Assembler() :
    _terms(), _bilinearTerms(), _linearTerms(), _sortedBilinearTerms(), _sortedLinearTerms()
  {
  }

  template< class T_Scalar >
  Assembler< T_Scalar >::~Assembler()
  {
  }

  template< class T_Scalar >
  void Assembler< T_Scalar >::addBilinearTerm(BilinearTermInterface< T_Scalar > *const term)
  {
    _bilinearTerms.push_back(term);
    _terms.push_back(term);
  }

  template< class T_Scalar >
  void Assembler< T_Scalar >::addLinearTerm(LinearTermInterface< T_Scalar > *const term)
  {
    _linearTerms.push_back(term);
    _terms.push_back(term);
  }

  template< class T_Scalar >
  void Assembler< T_Scalar >::sort()
  {
    // Sort bilinear terms by fields and by integrationType
    for(auto i = 0ULL; i < _bilinearTerms.size(); ++i) {
      const std::string integrationType = _bilinearTerms[i]->integrationType();
      _sortedBilinearTerms[std::make_pair(_bilinearTerms[i]->field(0), _bilinearTerms[i]->field(1))][integrationType].push_back(_bilinearTerms[i]);
    }

    // Sort linear terms by field and by integrationType
    for(auto i = 0ULL; i < _linearTerms.size(); ++i) {
      const std::string integrationType = _linearTerms[i]->integrationType();
      _sortedLinearTerms[_linearTerms[i]->field()][integrationType].push_back(_linearTerms[i]);
    }
  }

  template< class T_Scalar >
  void Assembler< T_Scalar >::assemblyInitialization(problem::IndiceBucket &indices, problem::ElementBucket< scalar::Precision< T_Scalar > > &elements, problem::FunctionSpaceBucket< scalar::Precision< T_Scalar > > &functionSpaces, const int elementType, const std::pair< int, int > &entity)
  {
    // get all fields and integrationTypes
    std::vector< const field::FieldInterface< T_Scalar > * > fields;
    std::vector< std::string > integrationTypes;
    std::unordered_map< std::string, std::vector< Term< T_Scalar > * > > termsByIntegrationType;

    for(auto i = 0ULL; i < _bilinearTerms.size(); ++i) {
      for(auto j = 0; j < 2; ++j) {
        auto it = std::find(fields.begin(), fields.end(), _bilinearTerms[i]->field(j));
        if(it == fields.end()) {
          fields.push_back(_bilinearTerms[i]->field(j));
        }
      }
      auto it = std::find(integrationTypes.begin(), integrationTypes.end(), _bilinearTerms[i]->integrationType());
      if(it == integrationTypes.end()) {
        integrationTypes.push_back(_bilinearTerms[i]->integrationType());
      }

      termsByIntegrationType[_bilinearTerms[i]->integrationType()].push_back(_bilinearTerms[i]);
    }

    for(auto i = 0ULL; i < _linearTerms.size(); ++i) {
      {
        auto it = std::find(fields.begin(), fields.end(), _linearTerms[i]->field());
        if(it == fields.end()) {
          fields.push_back(_linearTerms[i]->field());
        }
      }
      auto it = std::find(integrationTypes.begin(), integrationTypes.end(), _linearTerms[i]->integrationType());
      if(it == integrationTypes.end()) {
        integrationTypes.push_back(_linearTerms[i]->integrationType());
      }

      termsByIntegrationType[_linearTerms[i]->integrationType()].push_back(_linearTerms[i]);
    }

    // get jacobians info
    std::unordered_map< std::string, bool > needJacobians;
    std::unordered_map< std::string, bool > needPoints;

    for(auto i = 0ULL; i < integrationTypes.size(); ++i) {
      needJacobians.insert(std::make_pair(integrationTypes[i], false));
      needPoints.insert(std::make_pair(integrationTypes[i], false));
    }

    for(auto i = 0ULL; i < _terms.size(); ++i) {
      needJacobians[_terms[i]->integrationType()] = needJacobians[_terms[i]->integrationType()] || _terms[i]->needJacobians();
      needPoints[_terms[i]->integrationType()] = needPoints[_terms[i]->integrationType()] || _terms[i]->needGaussCoordinates(entity);
    }

    std::vector< scalar::Precision< T_Scalar >, numa::allocator< scalar::Precision< T_Scalar > > > jacobians;
    std::vector< double > gmshJacobians;
    std::vector< scalar::Precision< T_Scalar >, numa::allocator< scalar::Precision< T_Scalar > > > determinants;
    std::vector< double > gmshDeterminants;
    std::vector< scalar::Precision< T_Scalar >, numa::allocator< scalar::Precision< T_Scalar > > > points;
    std::vector< double > gmshPoints;

    std::vector< double > integrationWeights;
    std::vector< double > integrationPoints;

    for(auto i = 0ULL; i < integrationTypes.size(); ++i) {
      unsigned int nbrOfGaussPoints = field::FunctionSpaceInterface< double >::GetGaussInfo(integrationTypes[i], elementType, integrationWeights, integrationPoints);
      gmsh::model::mesh::preallocateJacobians(elementType, nbrOfGaussPoints, needJacobians[integrationTypes[i]], true, needPoints[integrationTypes[i]], gmshJacobians, gmshDeterminants, gmshPoints, entity.second);
#pragma omp parallel num_threads(omp::getMaxThreads())
      {
        const unsigned int numThreads = omp::getNumThreads();
        const unsigned int myThreadID = omp::getThreadNum();

        gmsh::model::mesh::getJacobians(elementType, integrationPoints, gmshJacobians, gmshDeterminants, gmshPoints, entity.second, myThreadID, numThreads);
      }

      numa::copy(determinants, gmshDeterminants);
      elements.determinants(integrationTypes[i], determinants);
      numa::copy(jacobians, gmshJacobians);
      elements.jacobians(integrationTypes[i], jacobians);
      numa::copy(points, gmshPoints);

      for(auto j = 0ULL; j < termsByIntegrationType[integrationTypes[i]].size(); ++j) {
        termsByIntegrationType[integrationTypes[i]][j]->assemblyInitialization(functionSpaces, points, elementType, entity);
      }
    }

    for(auto i = 0ULL; i < fields.size(); ++i) {
      fields[i]->fillIndices(indices, entity, elementType);
    }
  }

  template< class T_Scalar >
  void Assembler< T_Scalar >::_assembleBilinear(const problem::ElementBucket< scalar::Precision< T_Scalar > > &elements, const problem::IndiceBucket &indices, std::vector<system::VectorFactory< T_Scalar >>& b, system::MatrixFactory< T_Scalar > *const A, const int elementType, const std::pair< int, int > &entity) const
  {
    for(auto itField = _sortedBilinearTerms.begin(); itField != _sortedBilinearTerms.end(); ++itField) {
      const field::FieldInterface< T_Scalar > *fieldLhs = itField->first.first;
      const field::FieldInterface< T_Scalar > *fieldRhs = itField->first.second;

      const std::vector< std::pair< unsigned long long, int > > *const indicesLhs = indices.indices(fieldLhs->tag());
      const std::vector< std::pair< unsigned long long, int > > *const indicesRhs = indices.indices(fieldRhs->tag());

      std::map< unsigned int, const std::vector< scalar::Precision< T_Scalar >, numa::allocator< scalar::Precision< T_Scalar > > > * > determinants;
      std::map< unsigned int, const std::vector< scalar::Precision< T_Scalar >, numa::allocator< scalar::Precision< T_Scalar > > > * > jacobians;

      std::vector< BilinearTermInterface< T_Scalar > * > terms;

      // pre-assembly part
      for(auto itIntegration = itField->second.begin(); itIntegration != itField->second.end(); ++itIntegration) {
        std::vector< scalar::Precision< T_Scalar > > gaussWeights;
        std::vector< scalar::Precision< T_Scalar > > gaussPoints;
        unsigned int nbrOfGaussPoints = field::FunctionSpaceInterface< scalar::Precision< T_Scalar > >::GetGaussInfo(itIntegration->first, elementType, gaussWeights, gaussPoints);

        bool needPoints = false;
        for(auto termIndex = 0ULL; termIndex < itIntegration->second.size(); ++termIndex) {
          if(itIntegration->second[termIndex]->domain().haveJacobiansModificators(entity)) {
            needPoints = true;
            break;
          }
        }

        std::vector< scalar::Precision< T_Scalar >, numa::allocator< scalar::Precision< T_Scalar > > > points;

        if(needPoints) {
          std::vector< double > gmshJacobians;
          std::vector< double > gmshDeterminants;
          std::vector< double > gmshPoints;

          std::vector< double > gmshGaussPoints;
          scalar::copy(gmshGaussPoints, gaussPoints);

          gmsh::model::mesh::preallocateJacobians(elementType, nbrOfGaussPoints, false, true, true, gmshJacobians, gmshDeterminants, gmshPoints, entity.second);
#pragma omp parallel num_threads(omp::getMaxThreads())
          {
            const unsigned int numThreads = omp::getNumThreads();
            const unsigned int myThreadID = omp::getThreadNum();

            gmsh::model::mesh::getJacobians(elementType, gmshGaussPoints, gmshJacobians, gmshDeterminants, gmshPoints, entity.second, myThreadID, numThreads);
          }
          gmshDeterminants.clear();
          gmshDeterminants.shrink_to_fit();
          numa::copy(points, gmshPoints);
        }

        for(auto termIndex = 0ULL; termIndex < itIntegration->second.size(); ++termIndex) {
          if(itIntegration->second[termIndex]->domain().haveJacobiansModificators(entity)) {
            std::vector< scalar::Precision< T_Scalar >, numa::allocator< scalar::Precision< T_Scalar > > > *determinantMod = new std::vector< scalar::Precision< T_Scalar >, numa::allocator< scalar::Precision< T_Scalar > > >();
            std::vector< scalar::Precision< T_Scalar >, numa::allocator< scalar::Precision< T_Scalar > > > *jacobianJac = new std::vector< scalar::Precision< T_Scalar >, numa::allocator< scalar::Precision< T_Scalar > > >();
            if(elements.determinants(itIntegration->first) != nullptr) {
              determinantMod->insert(determinantMod->begin(), elements.determinants(itIntegration->first)->begin(), elements.determinants(itIntegration->first)->end());
            }
            if(elements.jacobians(itIntegration->first) != nullptr) {
              jacobianJac->insert(jacobianJac->begin(), elements.jacobians(itIntegration->first)->begin(), elements.jacobians(itIntegration->first)->end());
            }
#pragma omp parallel num_threads(omp::getMaxThreads())
            itIntegration->second[termIndex]->domain().applyJacobiansModificator(points, *determinantMod, *jacobianJac, entity);
            determinants[itIntegration->second[termIndex]->tag()] = determinantMod;
            jacobians[itIntegration->second[termIndex]->tag()] = jacobianJac;
          }
          else {
            determinants[itIntegration->second[termIndex]->tag()] = elements.determinants(itIntegration->first);
            jacobians[itIntegration->second[termIndex]->tag()] = elements.jacobians(itIntegration->first);
          }

          terms.push_back(itIntegration->second[termIndex]);
        }
      }

      // assembly part
      if(A->getModule() == nullptr) {
        throw common::Exception("There is no module in the matrix factory");
      }
      const std::string moduleName = A->getModule()->name();
      if(moduleName == "A") {
#pragma omp parallel num_threads(omp::getMaxThreads())
        {
          const unsigned int nbrDofsByElementLhs = terms[0]->nbrDofsByElement(0);
          const unsigned int nbrDofsByElementRhs = terms[0]->nbrDofsByElement(1);
          const unsigned long long nbrElements = determinants[terms[0]->tag()]->size() / terms[0]->nbrGaussPoints();
          Eigen::MatrixX< T_Scalar > A_e(nbrDofsByElementLhs, nbrDofsByElementRhs);

#pragma omp for
          for(auto elementIndex = 0ULL; elementIndex < nbrElements; ++elementIndex) {
            A_e.setZero();
            for(auto termIndex = 0ULL; termIndex < terms.size(); ++termIndex) {
              const unsigned long long offset = elementIndex * terms[termIndex]->nbrGaussPoints();
              const unsigned int tag = terms[termIndex]->tag();

              terms[termIndex]->evaluate(A_e, &((*determinants[tag])[offset]), &(*jacobians[tag])[9 * offset], elementIndex);
            }

            if(A->addValues(nbrDofsByElementRhs, nbrDofsByElementLhs, &(*indicesRhs)[elementIndex * nbrDofsByElementRhs], &(*indicesLhs)[elementIndex * nbrDofsByElementLhs], A_e.data())) {
              for(auto &bi : b)
                bi.addValuesDC(nbrDofsByElementRhs, nbrDofsByElementLhs, &(*indicesRhs)[elementIndex * nbrDofsByElementRhs], &(*indicesLhs)[elementIndex * nbrDofsByElementLhs], A_e.data());
            }
          }
        }
      }
      else if(moduleName == "AFrequency") {
        if constexpr(scalar::IsComplex< T_Scalar >::value) {
          const scalar::ComplexPrecision< T_Scalar > im(0., 1.);
          const scalar::Precision< T_Scalar > frequency = static_cast< const system::AFrequencyModule< T_Scalar > * >(A->getModule())->getFrequency();
#pragma omp parallel num_threads(omp::getMaxThreads())
          {
            const unsigned int nbrDofsByElementLhs = terms[0]->nbrDofsByElement(0);
            const unsigned int nbrDofsByElementRhs = terms[0]->nbrDofsByElement(1);
            const unsigned long long nbrElements = determinants[terms[0]->tag()]->size() / terms[0]->nbrGaussPoints();
            Eigen::MatrixX< T_Scalar > A_e(nbrDofsByElementLhs, nbrDofsByElementRhs);
            Eigen::MatrixX< T_Scalar > A_e_tmp(nbrDofsByElementLhs, nbrDofsByElementRhs);

#pragma omp for
            for(auto elementIndex = 0ULL; elementIndex < nbrElements; ++elementIndex) {
              A_e.setZero();
              A_e_tmp.setZero();
              for(auto termIndex = 0ULL; termIndex < terms.size(); ++termIndex) {
                const unsigned long long offset = elementIndex * terms[termIndex]->nbrGaussPoints();
                const unsigned int tag = terms[termIndex]->tag();

                if(terms[termIndex]->unknownFieldType() == equation::UnknownFieldType::DtDt) {
                  terms[termIndex]->evaluate(A_e_tmp, &((*determinants[tag])[offset]), &(*jacobians[tag])[9 * offset], elementIndex);
                  A_e += -frequency * frequency * A_e_tmp;
                  A_e_tmp.setZero();
                }
                else if(terms[termIndex]->unknownFieldType() == equation::UnknownFieldType::Dt) {
                  terms[termIndex]->evaluate(A_e_tmp, &((*determinants[tag])[offset]), &(*jacobians[tag])[9 * offset], elementIndex);
                  A_e += -im * frequency * A_e_tmp;
                  A_e_tmp.setZero();
                }
                else if(terms[termIndex]->unknownFieldType() == equation::UnknownFieldType::NotDt) {
                  terms[termIndex]->evaluate(A_e, &((*determinants[tag])[offset]), &(*jacobians[tag])[9 * offset], elementIndex);
                }
              }

              if(A->addValues(nbrDofsByElementRhs, nbrDofsByElementLhs, &(*indicesRhs)[elementIndex * nbrDofsByElementRhs], &(*indicesLhs)[elementIndex * nbrDofsByElementLhs], A_e.data())) {
                for(auto &bi : b)
                  bi.addValuesDC(nbrDofsByElementRhs, nbrDofsByElementLhs, &(*indicesRhs)[elementIndex * nbrDofsByElementRhs], &(*indicesLhs)[elementIndex * nbrDofsByElementLhs], A_e.data());
              }
            }
          }
        }
      }
      else {
#pragma omp parallel num_threads(omp::getMaxThreads())
        for(auto idMatrix = 0ULL; idMatrix < moduleName.size(); ++idMatrix) {
          const char matrixType = moduleName[idMatrix];
#pragma omp single
          A->getModule()->activate(matrixType);
          const unsigned int nbrDofsByElementLhs = terms[0]->nbrDofsByElement(0);
          const unsigned int nbrDofsByElementRhs = terms[0]->nbrDofsByElement(1);
          const unsigned long long nbrElements = determinants[terms[0]->tag()]->size() / terms[0]->nbrGaussPoints();
          Eigen::MatrixX< T_Scalar > A_e(nbrDofsByElementLhs, nbrDofsByElementRhs);

#pragma omp for
          for(auto elementIndex = 0ULL; elementIndex < nbrElements; ++elementIndex) {
            A_e.setZero();
            for(auto termIndex = 0ULL; termIndex < terms.size(); ++termIndex) {
              const unsigned long long offset = elementIndex * terms[termIndex]->nbrGaussPoints();
              const unsigned int tag = terms[termIndex]->tag();

              if(matrixType == 'M' && terms[termIndex]->unknownFieldType() == equation::UnknownFieldType::DtDt) {
                terms[termIndex]->evaluate(A_e, &((*determinants[tag])[offset]), &(*jacobians[tag])[9 * offset], elementIndex);
              }
              else if(matrixType == 'C' && terms[termIndex]->unknownFieldType() == equation::UnknownFieldType::Dt) {
                terms[termIndex]->evaluate(A_e, &((*determinants[tag])[offset]), &(*jacobians[tag])[9 * offset], elementIndex);
              }
              else if(matrixType == 'K' && terms[termIndex]->unknownFieldType() == equation::UnknownFieldType::NotDt) {
                terms[termIndex]->evaluate(A_e, &((*determinants[tag])[offset]), &(*jacobians[tag])[9 * offset], elementIndex);
              }
            }

            if(matrixType == 'K') {
              if(A->addValues(nbrDofsByElementRhs, nbrDofsByElementLhs, &(*indicesRhs)[elementIndex * nbrDofsByElementRhs], &(*indicesLhs)[elementIndex * nbrDofsByElementLhs], A_e.data())) {
                for(auto &bi : b)

                  bi.addValuesDC(nbrDofsByElementRhs, nbrDofsByElementLhs, &(*indicesRhs)[elementIndex * nbrDofsByElementRhs], &(*indicesLhs)[elementIndex * nbrDofsByElementLhs], A_e.data());
              }
            }
            else {
              A->addValues(nbrDofsByElementRhs, nbrDofsByElementLhs, &(*indicesRhs)[elementIndex * nbrDofsByElementRhs], &(*indicesLhs)[elementIndex * nbrDofsByElementLhs], A_e.data());
            }
          }
        }
      }

      // free jacobians that are modified
      for(auto itIntegration = itField->second.begin(); itIntegration != itField->second.end(); ++itIntegration) {
        for(auto termIndex = 0ULL; termIndex < itIntegration->second.size(); ++termIndex) {
          if(itIntegration->second[termIndex]->domain().haveJacobiansModificators(entity)) {
            delete determinants[itIntegration->second[termIndex]->tag()];
            delete jacobians[itIntegration->second[termIndex]->tag()];
          }
        }
      }
    }
  }

  template< class T_Scalar >
  void Assembler< T_Scalar >::_assembleLinear(const problem::ElementBucket< scalar::Precision< T_Scalar > > &elements, const problem::IndiceBucket &indices, std::vector<system::VectorFactory< T_Scalar >>& b, system::MatrixFactory< T_Scalar > *const A, const int elementType, const std::pair< int, int > &entity) const
  {
    for(auto fieldAndTerms: _sortedLinearTerms) {
      // fieldAndTerms is a key-value pair where the key is the field (whose test function is used) and value is another map
      // that maps an integrationType (e.g. "Gauss10") to the list of terms needing those Gauss points
      auto& integrationTermsDict = fieldAndTerms.second; // All terms, regrouped by integration
      const field::FieldInterface< T_Scalar > *field = fieldAndTerms.first; // Field whose TF is used

      const std::vector< std::pair< unsigned long long, int > > *const indicesRhs = indices.indices(field->tag());

      std::map< unsigned int, const std::vector< scalar::Precision< T_Scalar >, numa::allocator< scalar::Precision< T_Scalar > > > * > determinants;
      std::map< unsigned int, const std::vector< scalar::Precision< T_Scalar >, numa::allocator< scalar::Precision< T_Scalar > > > * > jacobians;

      std::vector< LinearTermInterface< T_Scalar > * > terms; // Why duplicate ?

      // pre-assembly part
      // Iterate on pairs <string, TermInterface*>
      for(auto integrationTermsPair: integrationTermsDict) {
        std::string integrationType = integrationTermsPair.first;
        auto& termsList = integrationTermsPair.second;
        std::vector< scalar::Precision< T_Scalar > > gaussWeights;
        std::vector< scalar::Precision< T_Scalar > > gaussPoints;
        unsigned int nbrOfGaussPoints = field::FunctionSpaceInterface< scalar::Precision< T_Scalar > >::GetGaussInfo(integrationType, elementType, gaussWeights, gaussPoints);

        bool needPoints = false;
        for(auto termPtr: termsList) {
          if(termPtr->domain().haveJacobiansModificators(entity)) {
            needPoints = true;
            break;
          }
        }

        std::vector< scalar::Precision< T_Scalar >, numa::allocator< scalar::Precision< T_Scalar > > > points;

        if(needPoints) {
          std::vector< double > gmshJacobians;
          std::vector< double > gmshDeterminants;
          std::vector< double > gmshPoints;

          std::vector< double > gmshGaussPoints;
          scalar::copy(gmshGaussPoints, gaussPoints);

          gmsh::model::mesh::preallocateJacobians(elementType, nbrOfGaussPoints, false, true, true, gmshJacobians, gmshDeterminants, gmshPoints, entity.second);
#pragma omp parallel num_threads(omp::getMaxThreads())
          {
            const unsigned int numThreads = omp::getNumThreads();
            const unsigned int myThreadID = omp::getThreadNum();

            gmsh::model::mesh::getJacobians(elementType, gmshGaussPoints, gmshJacobians, gmshDeterminants, gmshPoints, entity.second, myThreadID, numThreads);
          }
          gmshDeterminants.clear();
          gmshDeterminants.shrink_to_fit();
          numa::copy(points, gmshPoints);
        }

        for(auto termPtr: termsList) {
          if(termPtr->domain().haveJacobiansModificators(entity)) {
            std::vector< scalar::Precision< T_Scalar >, numa::allocator< scalar::Precision< T_Scalar > > > *determinantMod = new std::vector< scalar::Precision< T_Scalar >, numa::allocator< scalar::Precision< T_Scalar > > >();
            std::vector< scalar::Precision< T_Scalar >, numa::allocator< scalar::Precision< T_Scalar > > > *jacobianJac = new std::vector< scalar::Precision< T_Scalar >, numa::allocator< scalar::Precision< T_Scalar > > >();
            if(elements.determinants(integrationType) != nullptr) {
              determinantMod->insert(determinantMod->begin(), elements.determinants(integrationType)->begin(), elements.determinants(integrationType)->end());
            }
            if(elements.jacobians(integrationType) != nullptr) {
              jacobianJac->insert(jacobianJac->begin(), elements.jacobians(integrationType)->begin(), elements.jacobians(integrationType)->end());
            }
#pragma omp parallel num_threads(omp::getMaxThreads())
            termPtr->domain().applyJacobiansModificator(points, *determinantMod, *jacobianJac, entity);
            determinants[termPtr->tag()] = determinantMod;
            jacobians[termPtr->tag()] = jacobianJac;
          }
          else {
            determinants[termPtr->tag()] = elements.determinants(integrationType);
            jacobians[termPtr->tag()] = elements.jacobians(integrationType);
          }

          terms.push_back(termPtr);
        }
      }

      // assembly part
#pragma omp parallel num_threads(omp::getMaxThreads())
      {
        const unsigned int nbrDofsByElement = terms[0]->nbrDofsByElement();
        const unsigned long long nbrElements = determinants[terms[0]->tag()]->size() / terms[0]->nbrGaussPoints();
        //Eigen::VectorX< T_Scalar > b_e(nbrDofsByElement);
        std::vector< Eigen::VectorX< T_Scalar > > b_es(b.size(), Eigen::VectorX< T_Scalar >(nbrDofsByElement));
#pragma omp for
        for(auto elementIndex = 0ULL; elementIndex < nbrElements; ++elementIndex) {
          for(auto &b_e : b_es)
            b_e.setZero();

          for(auto termIndex = 0ULL; termIndex < terms.size(); ++termIndex) {
            unsigned iRHS = terms[termIndex]->rhsIdx();
            const unsigned long long offset = elementIndex * terms[termIndex]->nbrGaussPoints();
            const unsigned int tag = terms[termIndex]->tag();

            terms[termIndex]->evaluate(b_es.at(iRHS), &((*determinants[tag])[offset]), &(*jacobians[tag])[9 * offset], elementIndex);
          }

          for(unsigned iRhs = 0; iRhs < b.size(); ++iRhs) {
            b[iRhs].addValues(nbrDofsByElement, &(*indicesRhs)[elementIndex * nbrDofsByElement], b_es[iRhs].data());
          }
        }
      }

      // free jacobians that are modified
      for(auto itIntegration = integrationTermsDict.begin(); itIntegration != integrationTermsDict.end(); ++itIntegration) {
        auto &termsList = itIntegration->second;

        for(auto termIndex = 0ULL; termIndex < termsList.size(); ++termIndex) {
          if(termsList[termIndex]->domain().haveJacobiansModificators(entity)) {
            delete determinants[termsList[termIndex]->tag()];
            delete jacobians[termsList[termIndex]->tag()];
          }
        }
      }
    }
  }

  template< class T_Scalar >
  void Assembler< T_Scalar >::assemble(const problem::ElementBucket< scalar::Precision< T_Scalar > > &elements, const problem::IndiceBucket &indices, std::vector<system::VectorFactory< T_Scalar >>& b, system::MatrixFactory< T_Scalar > *const A, const int elementType, const std::pair< int, int > &entity) const
  {
    // Could be done:
    // Loop over N memory buckets (n= 1 ... N) and then called terms to evaluate: points[points(1), points(2), ..., points(n)], elms[elms(1), elms(2), ..., elms(n)]
    //_equationEvaluatorLhs->initialize(points(n), gaussPoints, elementType, entity);
    //_equationEvaluatorRhs->initialize(points(n), gaussPoints, elementType, entity);
    //
    //_assembleBilinear(elements, indices, b, A, elementType, entity, n);
    //_assembleLinear(elements, indices, b, A, elementType, entity, n);

    _assembleBilinear(elements, indices, b, A, elementType, entity);
    _assembleLinear(elements, indices, b, A, elementType, entity);
  }

  INSTANTIATE_CLASS(Assembler, 4, TEMPLATE_ARGS(std::complex< double >, std::complex< float >, double, float))


} // namespace gmshfem::term
