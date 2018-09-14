#ifndef STELLA_PC_H
#define STELLA_PC_H

#include "petscksp.h"

#ifdef WITH_CEDAR
#include <cedar/capi.h>
#endif

/**
 * @file: stella_pc.h
 *
 * Preconditioner object
 */

typedef struct {
	int nd;
	#ifdef WITH_CEDAR
	cedar_solver solver;
	cedar_vec solvec;
	cedar_vec rhsvec;
	#endif
} stella_pc;

PetscErrorCode stella_pc_create(stella_pc**);

PetscErrorCode stella_pc_setup(PC);

PetscErrorCode stella_pc_apply(PC pc, Vec x, Vec y);

PetscErrorCode stella_pc_destroy(PC);

#endif
