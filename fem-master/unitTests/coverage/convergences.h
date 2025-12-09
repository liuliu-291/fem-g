// GmshFEM - Copyright (C) 2019-2021 Université de Liège
//
// See the LICENSE.txt file for license information. Please report all
// issues on https://gitlab.onelab.info/gmsh/fem/issues

#ifndef H_GMSHFEM_TEST_CONVERGENCES
#define H_GMSHFEM_TEST_CONVERGENCES

#include <GmshFem.h>

void convergences(gmshfem::common::GmshFem &gmshfem, unsigned int &numTest);

void convergencesH1();

#endif // H_GMSHFEM_TEST_CONVERGENCES
