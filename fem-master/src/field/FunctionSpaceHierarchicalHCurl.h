// GmshFEM - Copyright (C) 2019-2022, A. Royer, E. Béchet, C. Geuzaine, Université de Liège
//
// See the LICENSE.txt file for license information. Please report all
// issues on https://gitlab.onelab.info/gmsh/fem/issues

#ifndef H_GMSHFEM_FUNCTIONSPACEHIERARCHICALHCURL
#define H_GMSHFEM_FUNCTIONSPACEHIERARCHICALHCURL

#include "FunctionSpace.h"

#include <vector>

namespace gmshfem::field
{


  template< class T_PScalar >
  class FunctionSpaceHierarchicalHCurl final : public FunctionSpace< T_PScalar, Form::Form1 >
  {
   private:
    const unsigned int _order;

   public:
    FunctionSpaceHierarchicalHCurl(const unsigned int order);
    ~FunctionSpaceHierarchicalHCurl();

    virtual std::string _getGmshName(bool derivative) const override;
    virtual std::string getGmshFemName(bool derivative) const override;
    virtual std::string getGmshFemOrientationName() const override;
    virtual FunctionSpaceTypeForm1 type() const override;
    virtual unsigned int order() const override;
    virtual bool isConstantByElements() const override;
    virtual bool isFormP() const override;

    virtual FunctionSpaceHierarchicalHCurl *copy() const override;

    virtual unsigned int getBasisFunctions(const bool derivative, std::vector< T_PScalar > &functions, std::vector< int > &orientations, const std::vector< T_PScalar > &integrationPoints, const int elementType, const std::pair< int, int > &entity) const override;
    virtual unsigned int getBasisFunction_noexcept(const bool derivative, std::vector< T_PScalar > &functions, const std::vector< T_PScalar > &integrationPoints, const int elementType, const unsigned int orientation = 0) const noexcept override;
  };


} // namespace gmshfem::field

#endif // H_GMSHFEM_FUNCTIONSPACEHIERARCHICALHCURL
