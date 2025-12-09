// GmshFEM - Copyright (C) 2019-2022, A. Royer, E. Béchet, C. Geuzaine, Université de Liège
//
// See the LICENSE.txt file for license information. Please report all
// issues on https://gitlab.onelab.info/gmsh/fem/issues

#ifndef H_GMSHFEM_FUNCTIONSPACE
#define H_GMSHFEM_FUNCTIONSPACE

#include "FunctionSpaceInterface.h"

namespace gmshfem::field
{


  enum class FunctionSpaceTypeForm0 {
    Lagrange,
    HierarchicalH1
  };

  enum class FunctionSpaceTypeForm1 {
    HierarchicalHCurl,
    // Perpendicular
    P_Lagrange,
    P_HierarchicalH1,
    // Derivative
    D_Lagrange,
    D_HierarchicalH1
  };

  enum class FunctionSpaceTypeForm2 {
    P_HierarchicalHCurl,
    // Derivative
    D_HierarchicalHCurl,
    // Derivative and Perpendicular
    DP_Lagrange,
    DP_HierarchicalH1,
  };

  enum class FunctionSpaceTypeForm3 {
    Constant,
    P_Constant,
    // Derivative and Perpendicular
    DP_HierarchicalHCurl
  };

  template< class T_PScalar, field::Form T_Form >
  class FunctionSpace : public FunctionSpaceInterface< T_PScalar >
  {
   public:
    FunctionSpace();
    virtual ~FunctionSpace();

    virtual field::Form form() const override;
    virtual field::FunctionSpaceOfForm< T_Form > type() const = 0;
    virtual unsigned int order() const = 0;

    virtual FunctionSpace *copy() const override = 0;

    //  FunctionSpace - form0:
    //     Values return by gauss points:
    //      - BF : matrix(1, nbrDofsByElements)
    //      - Grad(BF) : matrix(3, nbrDofsByElements)
    //  FunctionSpace - form1:
    //     Values return by gauss points:
    //      - BF : matrix(3, nbrDofsByElements)
    //      - Curl(BF) : matrix(3, nbrDofsByElements)
    //  FunctionSpace - form2:
    //     Values return by gauss points:
    //      - BF : matrix(3, nbrDofsByElements)
    //      - Div(BF) : matrix(1, nbrDofsByElements)
    //  FunctionSpace - form3:
    //     Values return by gauss points:
    //      - BF : vector(1, nbrDofsByElements)

    virtual unsigned int numberOfComponents(const bool derivative) const override;
    virtual unsigned int getBasisFunctions(const bool derivative, std::vector< T_PScalar > &functions, std::vector< int > &orientations, const std::vector< T_PScalar > &integrationPoints, const int elementType, const std::pair< int, int > &entity) const override = 0;
    virtual unsigned int getBasisFunction_noexcept(const bool derivative, std::vector< T_PScalar > &functions, const std::vector< T_PScalar > &integrationPoints, const int elementType, const unsigned int orientations = 0) const noexcept override = 0;
  };


} // namespace gmshfem::field

namespace gmshfem
{
  typedef field::FunctionSpaceTypeForm0 functionSpaceH1;
  typedef field::FunctionSpaceTypeForm1 functionSpaceHCurl;
  typedef field::FunctionSpaceTypeForm2 functionSpaceHDiv;
  typedef field::FunctionSpaceTypeForm3 functionSpaceH0;
} // namespace gmshfem

#include "FunctionSpaceConstant.h"
#include "FunctionSpaceCurlHCurl.h"
#include "FunctionSpaceDivHDiv.h"
#include "FunctionSpaceGradH1.h"
#include "FunctionSpaceHierarchicalH1.h"
#include "FunctionSpaceHierarchicalHCurl.h"
#include "FunctionSpaceLagrange.h"
#include "FunctionSpaceParallel.h"
#include "FunctionSpacePerpendicular.h"

#endif // H_GMSHFEM_FUNCTIONSPACE
