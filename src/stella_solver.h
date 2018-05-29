#ifndef STELLA_SOLVER_H
#define STELLA_SOLVER_H

#include <petscdmda.h>
#include <petscksp.h>

#include "stella_dmap.h"
#include "stella_grid.h"
#include "stella_metric.h"
#include "stella_operator.h"
#include "stella_boundary.h"
#include "stella_level.h"
#include "stella_state.h"
#include "stella_fd.h"

/**
 * Efield solver log level.
 *
 * Sets what the efield solver will log.
 */
typedef enum {
	EFS_LOG_STATUS   = 1,      /**< log status messages */
	EFS_LOG_PROBLEM  = (1<<1), /**< Write the problem (matrix and rhs) to disk */
	EFS_LOG_RESIDUAL = (1<<2), /**< Write residual to stdout */
	EFS_LOG_VTK      = (1<<3), /**< Write vtk files of computed solution */
	EFS_LOG_EIGS     = (1<<4), /**< Output approximate eigenvalues from Krylov solver */
	EFS_LOG_ALL      = (1<<5)-1
} efs_log_level;


/**
 * Options for efield solver.
 */
typedef struct {
	int matfree;             /**< Sets whether the geometric solver uses matrix free operations */
	int algebraic;           /**< Use AMG vs structured */
	int axisymmetric;        /**< Use axisymmetric operator */
	PetscBool galerkin;      /**< Use Galerking coarsening */
	efs_log_level log_level; /**< Specify loglevel */
	PetscInt levels;         /**< Depth of multilevel solve */
} efs_option;


/**
 * Main datastructure for efield solver
 */
typedef struct {
	Mat A;   /**< Matrix for linear solve */
	DM  dm;  /**< Datastructure for structured grid */
	KSP ksp;
	PC pc;
	Vec x;
	Vec rhs;
	Vec sol;
	Vec dcoef;
	Vec bcoef;
	MPI_Comm comm;
	stella_state   state;
	stella_grid    grid;
	efs_option options;
	stella_level level;
	stella_dmap  *dmap;
	stella_operator *op;
	stella_boundary *boundary;
	stella_fd *fd;
	int num_patches;
	int ts; /**< Timestep (only used for logging). */
} efs;


PetscErrorCode efs_create(efs**, MPI_Comm);


PetscErrorCode efs_setup(efs*, int offset[], int stride[]);


PetscErrorCode efs_solve(efs*);


/**
 * Sets efield solver state.
 *
 * @param slv solver object
 * @param phi array where solver will put the solution
 * @param dcoef permittivity of electric field
 * @param jump  jump condition
 */
PetscErrorCode efs_set_state(efs *slv, double phi[], double dcoef[],
                             double bcoef[], double jump[]);


/**
 * Sets efield solver rhs.
 *
 * @param slv solver object
 * @param rhs right hand side for linear solve
 */
PetscErrorCode efs_set_rhs(efs *slv, double rhs[]);


/**
 * Sets anaylitical solution for linear solve.
 *
 * This is used in testing to set the analytical solution
 * in order to compute the error of the solve.
 */
PetscErrorCode efs_set_sol(efs *slv, double sol[]);


/**
 * Sets grid corners and Cartesian coordinates.
 *
 * @param slv solver object
 * @param is global index of first grid point on this processor in each dimension
 * @param ie global index of last grid point on this processor in each dimension
 * @param num_cells number of grid points on local processor
 * @param xyz array of Cartesian coordinates for each grid point on this processor
 */
PetscErrorCode efs_set_grid(efs *slv, int is[], int ie[], int num_cells, double xyz[]);


/**
 * Checks if log level is active
 */
int efs_log(efs*, efs_log_level);


/**
 * Set log level for efield solver
 */
void efs_set_log(efs*, efs_log_level);


/**
 * Calls all callbacks needed for setting up the matrix
 */
PetscErrorCode efs_setup_op(efs*);


/**
 * Calls callbacks required when the rhs is changed
 */
PetscErrorCode efs_setup_rhs(efs*);


/**
 * Destroys data structures owned by the efield solver
 */
PetscErrorCode efs_destroy(efs*);

#endif
