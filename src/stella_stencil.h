#ifndef STELLA_STENCIL_H
#define STELLA_STENCIL_H

#include <petscsys.h>


typedef struct {
	double *v;
	int len;
	int off;
} stella_stencil;


PetscErrorCode stella_stencil_create(stella_stencil **sten, int len, int off);


PetscErrorCode stella_stencil_destroy(stella_stencil *sten);


#endif
