// GmshFEM - Copyright (C) 2019-2021 Université de Liège
//
// See the LICENSE.txt file for license information. Please report all
// issues on https://gitlab.onelab.info/gmsh/fem/issues

#ifndef H_GMSHFEM_TEST_FIELDS
#define H_GMSHFEM_TEST_FIELDS

#include <GmshFem.h>

void fields(gmshfem::common::GmshFem &gmshfem, unsigned int &numTest);

template< class T_Scalar >
void trigoField();
template< class T_Scalar >
void hyperbolicField();
template< class T_Scalar >
void logarithmField();
template< class T_Scalar >
void rootField();
template< class T_Scalar >
void normField();
template< class T_Scalar >
void complexFieldFunction();
template< class T_Scalar >
void derivativeFieldFunction();
template< class T_Scalar >
void compoundField();


#endif // H_GMSHFEM_TEST_FIELDS
