#ifndef INIT
#define INIT

#include "petsc.h"

struct TopTenInit {

  TopTenInit(int& argc, char**& argv, char help[])
  {
    PetscInitialize(&argc,&argv,(char*) 0,help);
  }

  ~TopTenInit()
  {
    PetscFinalize();
  }

};

#endif
