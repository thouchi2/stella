#ifndef BOUNDARY_H
#define BOUNDARY_H

#include "grid.h"
#include "problem.h"

#ifdef __cplusplus
extern "C" {
#endif

typedef enum {
	NORTH = 0,
	SOUTH = 1,
	WEST = 2,
	EAST = 3,
	BACK = 4,
	FRONT = 5
} bdir;

typedef enum {
	DIRICHLET=30, NEUMANN=31
} btype;

typedef struct {
	char *classify;  /** Byte array that classifies grid points across domain */
	char *norm_dir;  /** Specifies outward normal of bc that may be applied at a point */
	double *values;  /** Value of boundary condition to be applied.  For example, x=value at a point for Dirichlet.*/
} boundary;

boundary *boundary_create(grid*, problem*);

void boundary_destroy(boundary*);

#ifdef __cplusplus
}
#endif

#endif
