#include <stdlib.h>

#include "stella_stencil.h"


PetscErrorCode stella_stencil_create(stella_stencil **sten, int len, int off)
{
	stella_stencil *st;

	st = (stella_stencil*) malloc(sizeof(stella_stencil));
	st->v = (double*) malloc(len*sizeof(double));
	st->v += off;

	st->len = len;
	st->off = off;

	*sten = st;
	return 0;
}


PetscErrorCode stella_stencil_destroy(stella_stencil *sten)
{
	double *dummy;

	dummy = sten->v - sten->off;
	free(dummy);

	return 0;
}
