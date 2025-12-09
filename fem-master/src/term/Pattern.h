// GmshFEM - Copyright (C) 2019-2022, A. Royer, E. Béchet, C. Geuzaine, Université de Liège
//
// See the LICENSE.txt file for license information. Please report all
// issues on https://gitlab.onelab.info/gmsh/fem/issues

#ifndef H_GMSHFEM_PATTERN
#define H_GMSHFEM_PATTERN

#include "BilinearTermInterface.h"
#include "LinearTermInterface.h"
#include "MathObject.h"

#include <set>
#include <vector>

namespace gmshfem::term
{


  template< class T_Scalar >
  class Pattern
  {
   private:
    std::vector< BilinearTermInterface< T_Scalar > * > _bilinearTerm;
    std::vector< LinearTermInterface< T_Scalar > * > _linearTerm;
    std::set< std::pair< const field::FieldInterface< T_Scalar > *, const field::FieldInterface< T_Scalar > * > > _selectedField;

   public:
    Pattern();
    ~Pattern();

    void addBilinearTerm(BilinearTermInterface< T_Scalar > *const term);
    void addLinearTerm(LinearTermInterface< T_Scalar > *const term);
    void sort(const std::pair< int, int > &entity);

    void build(system::VectorFactory< T_Scalar > *const b, system::MatrixFactory< T_Scalar > *const A, const int elementType, const std::pair< int, int > &entity);
  };


} // namespace gmshfem::term

#endif // H_GMSHFEM_PATTERN
