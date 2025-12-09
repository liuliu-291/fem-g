// GmshFEM - Copyright (C) 2019-2022, A. Royer, E. Béchet, C. Geuzaine, Université de Liège
//
// See the LICENSE.txt file for license information. Please report all
// issues on https://gitlab.onelab.info/gmsh/fem/issues

#ifndef H_GMSHFEM_ASSEMBLER
#define H_GMSHFEM_ASSEMBLER

#include "BilinearTermInterface.h"
#include "LinearTermInterface.h"
#include "MathObject.h"

#include <map>
#include <vector>

namespace gmshfem::problem
{
  template< class T_PScalar >
  class ElementBucket;

  template< class T_PScalar >
  class FunctionSpaceBucket;

  class IndiceBucket;
} // namespace gmshfem::problem

namespace gmshfem::term
{


  template< class T_Scalar >
  class Assembler
  {
   private:
    std::vector< Term< T_Scalar > * > _terms;
    std::vector< BilinearTermInterface< T_Scalar > * > _bilinearTerms;
    std::vector< LinearTermInterface< T_Scalar > * > _linearTerms;

    std::map< std::pair< const field::FieldInterface< T_Scalar > *, const field::FieldInterface< T_Scalar > * >, std::map< std::string, std::vector< BilinearTermInterface< T_Scalar > * > > > _sortedBilinearTerms;
    std::map< const field::FieldInterface< T_Scalar > *, std::map< std::string, std::vector< LinearTermInterface< T_Scalar > * > > > _sortedLinearTerms;

    void _assembleBilinear(const problem::ElementBucket< scalar::Precision< T_Scalar > > &elements, const problem::IndiceBucket &indices, std::vector<system::VectorFactory< T_Scalar >>& b, system::MatrixFactory< T_Scalar > *const A, const int elementType, const std::pair< int, int > &entity) const;
    void _assembleLinear(const problem::ElementBucket< scalar::Precision< T_Scalar > > &elements, const problem::IndiceBucket &indices, std::vector<system::VectorFactory< T_Scalar >>& b, system::MatrixFactory< T_Scalar > *const A, const int elementType, const std::pair< int, int > &entity) const;

   public:
    Assembler();
    ~Assembler();

    void addBilinearTerm(BilinearTermInterface< T_Scalar > *const term);
    void addLinearTerm(LinearTermInterface< T_Scalar > *const term);
    void sort();

    void assemblyInitialization(problem::IndiceBucket &indices, problem::ElementBucket< scalar::Precision< T_Scalar > > &elements, problem::FunctionSpaceBucket< scalar::Precision< T_Scalar > > &functionSpaces, const int elementType, const std::pair< int, int > &entity);

    void assemble(const problem::ElementBucket< scalar::Precision< T_Scalar > > &elements, const problem::IndiceBucket &indices, std::vector<system::VectorFactory< T_Scalar >>& b, system::MatrixFactory< T_Scalar > *const A, const int elementType, const std::pair< int, int > &entity) const;
  };


} // namespace gmshfem::term

#endif // H_GMSHFEM_ASSEMBLER
