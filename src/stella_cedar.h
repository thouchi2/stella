#ifndef CEDAR_STELLA_H
#define CEDAR_STELLA_H

#define chkerr(cedar_err) do{if (PetscUnlikely(cedar_err)) { SETERRQ1(PETSC_COMM_WORLD, PETSC_ERR_LIB, "Cedar error code: %d", cedar_err); }} while(0)

#ifdef WITH_CEDAR
#include <cedar/capi.h>
void cedar_copyto(const double *src, cedar_vec dst);
void cedar_copyfrom(cedar_vec src, double *dst);
#endif

#endif
