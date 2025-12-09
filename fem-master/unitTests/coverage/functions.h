// GmshFEM - Copyright (C) 2019-2021 Université de Liège
//
// See the LICENSE.txt file for license information. Please report all
// issues on https://gitlab.onelab.info/gmsh/fem/issues

#ifndef H_GMSHFEM_TEST_FUNCTIONS
#define H_GMSHFEM_TEST_FUNCTIONS

#include <GmshFem.h>

void functions(gmshfem::common::GmshFem &gmshfem, unsigned int &numTest);

template< class T_Scalar >
void trigo();
template< class T_Scalar >
void hyperbolic();
template< class T_Scalar >
void logarithm();
template< class T_Scalar >
void root();
template< class T_Scalar >
void bessel();
template< class T_Scalar >
void scalarFunction();
template< class T_Scalar >
void vectorFunction();
template< class T_Scalar >
void tensorFunction();
template< class T_Scalar >
void piecewiseFunction();
template< class T_Scalar >
void normalAndTangentFunction();
template< class T_Scalar >
void interpolation();
template< class T_Scalar >
void norm();
template< class T_Scalar >
void complexFunction();
template< class T_Scalar >
void coordinates();
template< class T_Scalar >
void others();


#endif // H_GMSHFEM_TEST_FUNCTIONS
