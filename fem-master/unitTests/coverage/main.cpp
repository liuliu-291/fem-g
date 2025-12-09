// GmshFEM - Copyright (C) 2019-2021 Université de Liège
//
// See the LICENSE.txt file for license information. Please report all
// issues on https://gitlab.onelab.info/gmsh/fem/issues

#include <GmshFem.h>
#include <Message.h>
#include <Timer.h>

#include "convergences.h"
#include "fields.h"
#include "formulations.h"
#include "functions.h"
#include "posts.h"

int main(int argc, char **argv)
{
  gmshfem::common::GmshFem gmshFem(argc, argv);
  
  gmshfem::msg::info << "Running continuous integration test..." << gmshfem::msg::endl;
  gmshfem::common::Timer timer;
  timer.tick();
  
  unsigned int numTest = 0;

  functions(gmshFem, numTest);
  formulations(gmshFem, numTest);
  posts(gmshFem, numTest);
  convergences(gmshFem, numTest);
  fields(gmshFem, numTest);
  
  timer.tock();
  gmshfem::msg::info << "Running continuous done in " << timer << "s" << gmshfem::msg::endl;
  return 0;
}
