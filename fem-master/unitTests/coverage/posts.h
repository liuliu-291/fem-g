// GmshFEM - Copyright (C) 2019-2021 Université de Liège
//
// See the LICENSE.txt file for license information. Please report all
// issues on https://gitlab.onelab.info/gmsh/fem/issues

#ifndef H_GMSHFEM_TEST_POSTS
#define H_GMSHFEM_TEST_POSTS

#include <GmshFem.h>

void posts(gmshfem::common::GmshFem &gmshfem, unsigned int &numTest);

template< class T_Scalar >
void integrate();
template< class T_Scalar >
void save();
template< class T_Scalar >
void pointEvaluationAndFieldCopy();
template< class T_Scalar >
void algebra();
template< class T_Scalar >
void basisFunctionsDraw();


#endif // H_GMSHFEM_TEST_POSTS
