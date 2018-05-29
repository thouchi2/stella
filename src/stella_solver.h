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
 * Stella solver log level.
 *
 * Sets what the stella solver will log.
 */
typedef enum {
	STELLA_LOG_STATUS   = 1,      /**< log status messages */
	STELLA_LOG_PROBLEM  = (1<<1), /**< Write the problem (matrix and rhs) to disk */
	STELLA_LOG_RESIDUAL = (1<<2), /**< Write residual to stdout */
	STELLA_LOG_VTK      = (1<<3), /**< Write vtk files of computed solution */
	STELLA_LOG_EIGS     = (1<<4), /**< Output approximate eigenvalues from Krylov solver */
	STELLA_LOG_ALL      = (1<<5)-1
} stella_log_level;


/**
 * Options for stella solver.
 */
typedef struct {
	int algebraic;              /**< Use AMG vs structured */
	int axisymmetric;           /**< Use axisymmetric operator */
	stella_log_level log_level; /**< Specify loglevel */
} stella_option;


/**
 * Main datastructure for stella
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
	stella_option options;
	stella_level level;
	stella_dmap  *dmap;
	stella_operator *op;
	stella_boundary *boundary;
	stella_fd *fd;
	int num_patches;
	int ts; /**< Timestep (only used for logging). */
} stella;


/**
 * Initializes structured elliptic solver
 *
 * @param[out] solver_ctx pointer to the solver object allocated in this call
 * @param[in] comm MPI Communicator to be used for the efield solver
 * @param[in] nGlobal number of global grid points in each dimension
 * @param[in] nProcs number of processors in each dimension
 * @param[in] nLocal number of local grid points in each dimension
 * @param[in] offset offsets of data inside Fortran arrays (e.g. if ghosts are included)
 * @param[in] stride strides of data inside Fortran arrays
 * @param[in] cartCoord Cartesian coordinate of my rank
 * @param[in] periodic flag for setting periodicity of a dimension
 * @param[in] periodic_storage flag for duplicating grid points if periodic
 * @param[in] nd number of dimensions
 * @param[in] ng grid number
 * @param[in] axisymmetric flag for setting the solver to be axisymmetric
 */
PetscErrorCode stella_init(stella **solver_ctx, MPI_Comm comm,
                           int nGlobal[], int nProcs[], int nLocal[],
                           int offset[], int stride[],
                           int cartCoord[], int periodic[], int periodic_storage,
                           int nd, int axisymmetric);


PetscErrorCode stella_solve(stella*);


/**
 * Sets solver state.
 *
 * @param slv solver object
 * @param phi array where solver will put the solution
 * @param dcoef permittivity of electric field
 * @param jump  jump condition
 */
PetscErrorCode stella_set_state(stella *slv, double phi[], double dcoef[],
                                double bcoef[], double jump[]);


/**
 * Sets solver rhs.
 *
 * @param slv solver object
 * @param rhs right hand side for linear solve
 */
PetscErrorCode stella_set_rhs(stella *slv, double rhs[]);


/**
 * Sets anaylitical solution for linear solve.
 *
 * This is used in testing to set the analytical solution
 * in order to compute the error of the solve.
 */
PetscErrorCode stella_set_sol(stella *slv, double sol[]);


/**
 * Sets grid corners and Cartesian coordinates.
 *
 * @param slv solver object
 * @param is global index of first grid point on this processor in each dimension
 * @param ie global index of last grid point on this processor in each dimension
 * @param num_cells number of grid points on local processor
 * @param xyz array of Cartesian coordinates for each grid point on this processor
 */
PetscErrorCode stella_set_grid(stella *slv, int is[], int ie[], int num_cells, double xyz[]);


/**
 * Checks if log level is active
 */
int stella_log(stella*, stella_log_level);


/**
 * Set log level for stella
 */
void stella_set_log(stella*, stella_log_level);


/**
 * Calls all callbacks needed for setting up the matrix
 */
PetscErrorCode stella_setup_op(stella*);


/**
 * Destroys data structures owned by the elliptic discretization
 */
PetscErrorCode stella_cleanup(stella*);

#endif
