// GmshFEM - Copyright (C) 2019-2022, A. Royer, E. Béchet, C. Geuzaine, Université de Liège
//
// See the LICENSE.txt file for license information. Please report all
// issues on https://gitlab.onelab.info/gmsh/fem/issues

#include "Pattern.h"

#include "FieldInterface.h"
#include "IndiceBucket.h"
#include "MatrixFactory.h"
#include "OmpInterface.h"
#include "VectorFactory.h"
#include "instantiate.h"

namespace gmshfem::term
{


  template< class T_Scalar >
  Pattern< T_Scalar >::Pattern() :
    _bilinearTerm(), _linearTerm(), _selectedField()
  {
  }

  template< class T_Scalar >
  Pattern< T_Scalar >::~Pattern()
  {
  }

  template< class T_Scalar >
  void Pattern< T_Scalar >::addBilinearTerm(BilinearTermInterface< T_Scalar > *const term)
  {
    _bilinearTerm.push_back(term);
  }

  template< class T_Scalar >
  void Pattern< T_Scalar >::addLinearTerm(LinearTermInterface< T_Scalar > *const term)
  {
    _linearTerm.push_back(term);
  }

  template< class T_Scalar >
  void Pattern< T_Scalar >::sort(const std::pair< int, int > &entity)
  {
    // The pattern is only computed on bilinear term
    for(auto i = 0ULL; i < _bilinearTerm.size(); ++i) {
      if(_bilinearTerm[i]->domain().have(entity)) {
        _selectedField.insert(std::make_pair(_bilinearTerm[i]->field(0), _bilinearTerm[i]->field(1)));
      }
    }
  }

  template< class T_Scalar >
  void Pattern< T_Scalar >::build(system::VectorFactory< T_Scalar > *const b, system::MatrixFactory< T_Scalar > *const A, const int elementType, const std::pair< int, int > &entity)
  {
    problem::IndiceBucket indices;
    for(auto it = _selectedField.begin(); it != _selectedField.end(); ++it) {
      unsigned int nbrDofsByElement[2];
      for(auto i = 0; i < 2; ++i) {
        const field::FieldInterface< T_Scalar > *field = (i == 0 ? it->first : it->second);
        nbrDofsByElement[i] = field->multiplicity() * field->fillIndices(indices, entity, elementType);
      }

      const std::vector< std::pair< unsigned long long, int > > *const indicesLhs = indices.indices(it->first->tag());
      const std::vector< std::pair< unsigned long long, int > > *const indicesRhs = indices.indices(it->second->tag());

#pragma omp parallel for num_threads(omp::getMaxThreads())
      for(auto elementIndex = 0ULL; elementIndex < indicesLhs->size() / nbrDofsByElement[0]; ++elementIndex) {
        A->addPattern(nbrDofsByElement[1], nbrDofsByElement[0], &(*indicesRhs)[elementIndex * nbrDofsByElement[1]], &(*indicesLhs)[elementIndex * nbrDofsByElement[0]]);
      }
    }
  }

  INSTANTIATE_CLASS(Pattern, 4, TEMPLATE_ARGS(std::complex< double >, std::complex< float >, double, float))


} // namespace gmshfem::term
