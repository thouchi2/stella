#ifndef STELLA_PC_H
#define STELLA_PC_H

#include "petscksp.h"

#ifdef WITH_BOXMG
#include <boxmg/capi.h>
#endif

/**
 * @file: stella_pc.h
 *
 * Preconditioner object
 */

typedef struct {
	#ifdef WITH_BOXMG
	bmg2_solver solver;
	#endif
} stella_pc;

PetscErrorCode stella_pc_create(stella_pc**);

PetscErrorCode stella_pc_setup(PC);

PetscErrorCode stella_pc_apply(PC pc, Vec x, Vec y);

PetscErrorCode stella_pc_destroy(PC);

#endif
