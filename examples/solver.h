#ifndef SOLVER_H
#define SOLVER_H

#include "stella_interface.h"

#include "boundary.h"
#include "state.h"
#include "grid.h"
#include "problem.h"

#ifdef __cplusplus
extern "C" {
#endif

typedef struct {
	blist *bnd;
	state *state;
	void  *ptr;
	int axisymmetric;
} solver;

solver *solver_create(grid*, problem*);

PetscErrorCode solver_init(solver*, grid*);

PetscErrorCode solver_run(solver*);

PetscErrorCode solver_destroy(solver*);

#ifdef __cplusplus
}
#endif

#endif
