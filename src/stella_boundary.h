#ifndef STELLA_BOUNDARY_H
#define STELLA_BOUNDARY_H

#include "stella_bc.h"


/**
 * Contains all boundary conditions for a grid.
 */
typedef struct stella_boundary {
	stella_bc **bc;
	int num_bc;       /**< number of boundary conditions */
	stella_level *level; /**< data needed on each level */
	stella_dmap *slv_dmap;/**< efield solver's datamapping object */
	stella_state *state;  /**< state data for efield solver */
	int axisymmetric; /**< flag to indicate axisymmetric operator */
	stella_fd *fd;        /**< finite difference coefficients */
} stella_boundary;


/**
 * Creates the boundary object for a grid.
 *
 * @param stella_bnd location to put newly created object
 * @param level data needed on each MG level
 * @param dmap efield solver's data mapping object
 * @param state state data for efield solver
 * @param fd finite difference coefficients
 */
PetscErrorCode stella_boundary_create(stella_boundary **stella_bnd, stella_level *level,
                                      stella_dmap *dmap, stella_state *state, stella_fd *fd);


/**
 * Set boundary conditions for elliptic solve.
 *
 * @param bnd boundary object
 * @param norm_dir grid array specifying the outward normal of a boundary condition (if applicable)
 * @param values grid array of boundary values
 */
PetscErrorCode stella_boundary_set(stella_boundary *bnd, char *norm_dir, double *values);


/**
 * Adds boundary conditions to matrix.
 *
 * @param bnd boundary object
 * @param A matrix to add boundary conditions
 * @param da data structure for structured grid
 */
PetscErrorCode stella_boundary_apply(stella_boundary *bnd, Mat A, DM da);


/**
 * Adds boundary conditions to the right hand side.
 *
 * Depending on the boundary condition and choice of implementation,
 * the matrix and right hand side may need to be modified.
 * @param bnd boundary object
 * @param da data structure for structured grid
 * @param rhs right hand side vector
 */
PetscErrorCode stella_boundary_apply_rhs(stella_boundary *bnd, DM da, Vec rhs);


/**
 * Destructor for boundary.
 */
PetscErrorCode stella_boundary_destroy(stella_boundary *bnd);


#endif
