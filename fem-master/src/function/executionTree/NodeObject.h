// GmshFEM - Copyright (C) 2019-2022, A. Royer, E. Béchet, C. Geuzaine, Université de Liège
//
// See the LICENSE.txt file for license information. Please report all
// issues on https://gitlab.onelab.info/gmsh/fem/issues

#ifndef H_GMSHFEM_NODEOBJECT
#define H_GMSHFEM_NODEOBJECT

namespace gmshfem::function
{


  enum class Arity {
    Nullary = 0,
    Unary = 1,
    Binary = 2,
    Ternary = 3,
    Novenary = 9,
    Multary
  };

  enum class NodeType {
    Analytical,
    Binary,
    ChangeOfCoordinates,
    Field,
    FieldScalarType,
    Multary,
    Novenary,
    Nullary,
    ScalarType,
    Ternary,
    Unary,
    UserDefined
  };


} // namespace gmshfem::function


#endif // H_GMSHFEM_NODEOBJECT
