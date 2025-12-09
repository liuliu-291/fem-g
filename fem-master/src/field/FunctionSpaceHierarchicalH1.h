// GmshFEM - Copyright (C) 2019-2022, A. Royer, E. Béchet, C. Geuzaine, Université de Liège
//
// See the LICENSE.txt file for license information. Please report all
// issues on https://gitlab.onelab.info/gmsh/fem/issues

#ifndef H_GMSHFEM_FUNCTIONSPACEHIERARCHICALH1
#define H_GMSHFEM_FUNCTIONSPACEHIERARCHICALH1

#include "FunctionSpace.h"

#include <vector>

namespace gmshfem::field
{


  template< class T_PScalar >
  class FunctionSpaceHierarchicalH1 final : public FunctionSpace< T_PScalar, Form::Form0 >
  {
   private:
    const unsigned int _order;

   public:
    FunctionSpaceHierarchicalH1(const unsigned int order);
    ~FunctionSpaceHierarchicalH1();

    virtual std::string _getGmshName(bool derivative) const override;
    virtual std::string getGmshFemName(bool derivative) const override;
    virtual std::string getGmshFemOrientationName() const override;
    virtual FunctionSpaceTypeForm0 type() const override;
    virtual unsigned int order() const override;
    virtual bool isConstantByElements() const override;
    virtual bool isFormP() const override;

    virtual FunctionSpaceHierarchicalH1 *copy() const override;

    virtual unsigned int getBasisFunctions(const bool derivative, std::vector< T_PScalar > &functions, std::vector< int > &orientations, const std::vector< T_PScalar > &integrationPoints, const int elementType, const std::pair< int, int > &entity) const override;
    virtual unsigned int getBasisFunction_noexcept(const bool derivative, std::vector< T_PScalar > &functions, const std::vector< T_PScalar > &integrationPoints, const int elementType, const unsigned int orientation = 0) const noexcept override;
  };


} // namespace gmshfem::field

#endif // H_GMSHFEM_FUNCTIONSPACEHIERARCHICALH1
