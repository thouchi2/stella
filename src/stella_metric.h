#ifndef STELLA_METRIC_H
#define STELLA_METRIC_H

#include <petscdmda.h>

#include "stella_fd.h"

typedef struct {
	Vec coef[10];
	Vec lcoef[10];
	int t2map[2][2];
	int t3map[3][3];
	Vec jac_v[9];
	Vec ljac_v[9];
	Vec jac;
	int nd;
} stella_metric;

typedef enum {
	MET_C=0,
	MET_S=1,
	MET_W=2,
	MET_SW=3,
	MET_SE=4,
	MET_B=5,
	MET_WB=6,
	MET_EB=7,
	MET_SB=8,
	MET_NB=9
} stella_met_dir;

PetscErrorCode stella_metric_create(stella_metric**,DM,stella_fd*);


PetscErrorCode stella_metric_destroy(stella_metric*);


#endif
