#ifndef STELLA_MAT_H
#define STELLA_MAT_H

#include "petscmat.h"

#ifdef WITH_BOXMG
#include <cedar/capi.h>

typedef struct
{
	int nd;
	bmg2_operator op2;
	bmg3_operator op3;
} stella_bmg2_mat;
#endif

PetscErrorCode stella_bmg_SetValuesStencil(Mat mat, PetscInt m, const MatStencil idxm[], PetscInt n,
                                        const MatStencil idxn[], const PetscScalar v[],
                                        InsertMode addv);


PetscErrorCode stella_bmg_mult(Mat mat, Vec, Vec);

#endif
