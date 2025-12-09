// GmshFEM - Copyright (C) 2019-2021 Université de Liège
//
// See the LICENSE.txt file for license information. Please report all
// issues on https://gitlab.onelab.info/gmsh/fem/issues

#ifndef H_GMSHFEM_TEST_FORMULATIONS
#define H_GMSHFEM_TEST_FORMULATIONS

#include <GmshFem.h>

void formulations(gmshfem::common::GmshFem &gmshfem, unsigned int &numTest);

template< class T_Scalar >
void globalQuantity();
template< class T_Scalar >
void eigensolve();
template< class T_Scalar >
void nonsymmetricMatrix();
template< class T_Scalar >
void periodic();


#endif // H_GMSHFEM_TEST_FORMULATIONS
