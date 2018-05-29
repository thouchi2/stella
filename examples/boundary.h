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
	DIRICHLET=600, NEUMANN=602
} btype;

typedef struct {
	int is[3];
	int ie[3];
	btype type;
	double *dirichlet;
	int norm_dir;
} boundary;

typedef struct {
	boundary *v;
	int len;
} blist;

blist *boundary_create(grid*, problem*);

void boundary_destroy(blist*);

#ifdef __cplusplus
}
#endif

#endif
