// GmshFEM - Copyright (C) 2019-2022, A. Royer, E. Béchet, C. Geuzaine, Université de Liège
//
// See the LICENSE.txt file for license information. Please report all
// issues on https://gitlab.onelab.info/gmsh/fem/issues

#include "Integrate.h"

#include "FunctionSpaceInterface.h"
#include "KahanSum.h"
#include "Message.h"
#include "OmpInterface.h"
#include "instantiate.h"

#include <gmsh.h>

namespace gmshfem::post
{


  template< class T_Scalar, Degree T_Degree >
  Integrate< T_Scalar, T_Degree >::Integrate(typename MathObject< T_Scalar, T_Degree >::Object *const value, const function::Function< T_Scalar, T_Degree > &function, const domain::GeometricObject &domain, const std::string &integrationType) :
    PostInterface(), _value(value), _function(function), _domain(domain), _integrationType(integrationType)
  {
  }

  template< class T_Scalar, Degree T_Degree >
  Integrate< T_Scalar, T_Degree >::~Integrate()
  {
  }

  template< class T_Scalar, Degree T_Degree >
  int Integrate< T_Scalar, T_Degree >::run()
  {
    typename MathObject< T_Scalar, T_Degree >::Object total;
    MathObject< T_Scalar, T_Degree >::zero(total);

    for(auto itEntity = _domain.cbegin(); itEntity != _domain.cend(); ++itEntity) {
      const bool needPoints = !_function.isConstant(*itEntity);
      std::vector< int > elementTypes;
      gmsh::model::mesh::getElementTypes(elementTypes, itEntity->first, itEntity->second);

      for(auto typeIndex = 0ULL; typeIndex < elementTypes.size(); ++typeIndex) {
        std::vector< common::KahanSum< typename MathObject< T_Scalar, T_Degree >::Object > > integral(omp::getMaxThreads(), total);

        std::vector< scalar::Precision< T_Scalar >, numa::allocator< scalar::Precision< T_Scalar > > > jacobians;
        std::vector< double > gmshJacobians;
        std::vector< scalar::Precision< T_Scalar >, numa::allocator< scalar::Precision< T_Scalar > > > determinants;
        std::vector< double > gmshDeterminants;
        std::vector< scalar::Precision< T_Scalar >, numa::allocator< scalar::Precision< T_Scalar > > > points;
        std::vector< double > gmshPoints;

        std::vector< scalar::Precision< T_Scalar > > gaussWeights;
        std::vector< scalar::Precision< T_Scalar > > gaussPoints;
        std::vector< double > gmshGaussPoints;

        unsigned int nbrOfGaussPoints = field::FunctionSpaceInterface< scalar::Precision< T_Scalar > >::GetGaussInfo(_integrationType, elementTypes[typeIndex], gaussWeights, gaussPoints);

        scalar::copy(gmshGaussPoints, gaussPoints);

        gmsh::model::mesh::preallocateJacobians(elementTypes[typeIndex], nbrOfGaussPoints, false, true, needPoints, gmshJacobians, gmshDeterminants, gmshPoints, itEntity->second);

#pragma omp parallel num_threads(omp::getMaxThreads())
        {
          const unsigned int numThreads = omp::getNumThreads();
          const unsigned int myThreadID = omp::getThreadNum();

          gmsh::model::mesh::getJacobians(elementTypes[typeIndex], gmshGaussPoints, gmshJacobians, gmshDeterminants, gmshPoints, itEntity->second, myThreadID, numThreads);
        }
        numa::copy(jacobians, gmshJacobians);
        numa::copy(determinants, gmshDeterminants);
        numa::copy(points, gmshPoints);

        if(_domain.haveJacobiansModificators(*itEntity)) {
          _domain.applyJacobiansModificator(points, determinants, jacobians, *itEntity);
        }

        typename MathObject< T_Scalar, T_Degree >::Object constant; // Used only if _function is constant
        if(_function.isConstant(*itEntity)) {
          _function.evaluate(constant, 0., 0., 0., *itEntity);
        }
        std::vector< typename MathObject< T_Scalar, T_Degree >::Object, numa::allocator< typename MathObject< T_Scalar, T_Degree >::Object > > values; // Used only if _function isn't constant
#pragma omp parallel num_threads(omp::getMaxThreads())
        {
          const unsigned int myThreadID = omp::getThreadNum();
          if(_function.isConstant(*itEntity)) {
            const unsigned int nbrElements = determinants.size() / nbrOfGaussPoints;
#pragma omp for
            for(auto j = 0U; j < nbrElements; ++j) {
              for(auto k = 0U; k < nbrOfGaussPoints; ++k) {
                integral[myThreadID] += constant * gaussWeights[k] * determinants[j * nbrOfGaussPoints + k];
              }
            }
          }
          else {
            const unsigned int nbrElements = determinants.size() / nbrOfGaussPoints;
            _function.evaluate(values, points, gaussPoints, elementTypes[typeIndex], *itEntity);
#pragma omp for
            for(auto j = 0U; j < nbrElements; ++j) {
              for(auto k = 0U; k < nbrOfGaussPoints; ++k) {
                integral[myThreadID] += values[j * nbrOfGaussPoints + k] * gaussWeights[k] * determinants[j * nbrOfGaussPoints + k];
              }
            }
          }
        }

        for(auto threadIndex = 0U; threadIndex < omp::getMaxThreads(); ++threadIndex) {
          total += integral[threadIndex].sum();
        }
      }
    }

    *_value = std::move(total);
    return 0;
  }

  INSTANTIATE_CLASS_2(Integrate, 4, 3, TEMPLATE_ARGS(std::complex< double >, std::complex< float >, double, float), TEMPLATE_ARGS(Degree::Degree0, Degree::Degree1, Degree::Degree2))


} // namespace gmshfem::post
