#ifndef STELLA_MAT_H
#define STELLA_MAT_H

#include "petscmat.h"

#ifdef WITH_BOXMG
#include <boxmg/capi.h>

typedef struct
{
	bmg2_operator op;
} stella_bmg2_mat;
#endif

PetscErrorCode stella_bmg2_SetValuesStencil(Mat mat, PetscInt m, const MatStencil idxm[], PetscInt n,
                                        const MatStencil idxn[], const PetscScalar v[],
                                        InsertMode addv);


PetscErrorCode stella_bmg2_mult(Mat mat, Vec, Vec);

#endif
