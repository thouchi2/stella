#ifndef STELLA_LEVEL_H
#define STELLA_LEVEL_H

#include <petscvec.h>

#include "stella_metric.h"
#include "stella_classify.h"

/**
 * Data structure containing data
 * needed on each MG level (previously an array because PDE was discretized on each level)
 */
typedef struct {
	Vec dcoef;           /** Permittivity */
	Vec ldcoef;          /** Local vector for dcoef */
	Vec nscale;        /** Neumann scaling for symmetry */
	Vec add_cont;      /** Additive contribution to rhs */
	Vec ladd_cont;     /** Local vector for add_cont */
	Vec pm;            /** Pointwise multiplication for add_cont */
	stella_metric *metric;
	Vec bcoef;
	Vec lbcoef;
	DM dm;

	stella_classify *classify;
} stella_level;


#endif
