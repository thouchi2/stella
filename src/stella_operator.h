#ifndef STELLA_OPERATOR_H
#define STELLA_OPERATOR_H

#include <petscdmda.h>

#include "stella_metric.h"
#include "stella_level.h"
#include "stella_fd.h"

/**
 * @file: stella_operator.h
 *
 * This is used to create an operator object that is
 * responsible for assembling the matrix for all rows
 * that are not owned by boundary conditions.
 */

static inline double have(double v1, double v2)
{
	return 2.0/((1.0/v1) + (1.0/v2));
}


typedef struct stella_operator_{
	PetscErrorCode (*assemble)(struct stella_operator_ *op, Mat A, DM da);
	stella_level *level;
	stella_fd *fd;
	int axisymmetric;
} stella_operator;


/**
 * Constructor for the operator object
 *
 * @param[out] efop operator object to create
 * @param level level object
 * @param fd finite difference coefficient object
 * @param nd number of dimensions
 */
PetscErrorCode stella_operator_create(stella_operator **efop, stella_level *level, stella_fd *fd, int nd);


/**
 * Used to assemble the interior of the matrix
 * (everywhere not owned by a boundary condition)
 *
 * @param op operator object
 * @param A PETSc matrix
 * @param da PETSc DM object
 */
PetscErrorCode stella_operator_assemble(stella_operator *op, Mat A, DM da);


/**
 * Destroys data structures owned by the operator
 */
PetscErrorCode stella_operator_destroy(stella_operator*);


#endif
