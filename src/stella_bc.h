#ifndef STELLA_BC_H
#define STELLA_BC_H

#include <petscdmda.h>
#include <petscmat.h>
#include <petscvec.h>

#include "stella_level.h"
#include "stella_dmap.h"
#include "stella_state.h"
#include "stella_grid.h"

/**
 * Base class for boundary conditions
 *
 * This module is used for defining common functions among boundary conditions.
 */
typedef enum {
	STELLA_DIRICHLET = 600,
	STELLA_NEUMANN = 602,
	STELLA_SCHWARZ = 698
} stella_bctype;


/**
 * Base class for boundary condition objects
 */
typedef struct stella_bc_ {
	PetscErrorCode (*apply)(struct stella_bc_ *bc, Mat A, DM da);
	PetscErrorCode (*apply_rhs)(struct stella_bc_ *bc, DM da, Vec rhs);
	PetscErrorCode (*symmetric)(struct stella_bc_ *bc, Mat A, DM da);
	PetscErrorCode (*destroy)(struct stella_bc_ *bc);
	stella_level *level;
	int axisymmetric;
	stella_dmap *slv_dmap;
	stella_state *state;
	stella_fd *fd;
	char *norm_dirs;
	double *values;
	stella_bctype btype;
	void *sub; /** pointer for child specific data */
} stella_bc;


/**
 * Creates a boundary condition object
 *
 * @param[out] efbc boundary condition object that is created
 * @param btype flag to specify type of boundary condition (e.g. Dirichlet or Neumann)
 * @param norm_dir direction of the normal for the boundary condition
 * @param is global indicies of local grid starts by dimension
 * @param ie global indicies of local grid ends by dimension
 * @param dirichlet grid function of dirichlet values (not used if btype is not dirichlet)
 * @param level stella_level object
 * @param dmap data mapping object for solver
 * @param state state variable object
 * @param fd object that contains finite difference coefficients
 */
PetscErrorCode stella_bc_create(stella_bc **efbc, stella_bctype btype, char *norm_dir, double *values,
                                stella_level *level, stella_dmap *dmap,
                                stella_state *state, stella_fd *fd);


/**
 * Applies boundary condition to matrix
 *
 * @param bc boundary condition object
 * @param A matrix to apply bc to
 * @param da PETSc DM object
 */
PetscErrorCode stella_bc_apply(stella_bc *bc, Mat A, DM da);


/**
 * Applies boundary condition to the rhs
 *
 * @param bc boundary condition object
 * @param rhs right hand side for solve
 */
PetscErrorCode stella_bc_apply_rhs(stella_bc *bc, DM da, Vec rhs);


/**
 * Callback used to make the operator symmetric in the presence of this bc
 *
 * @param bc boundary condition object
 * @param A PETSc matrix
 * @param da PETSC DM object
 */
PetscErrorCode stella_bc_symmetric(stella_bc *bc, Mat A, DM da);


/**
 * Destroys data structures owned by the bc
 *
 * @param bc boundary condition object
 */
PetscErrorCode stella_bc_destroy(stella_bc *bc);


#endif
