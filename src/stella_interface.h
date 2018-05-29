#ifndef STELLA_INTERFACE_H
#define STELLA_INTERFACE_H

#include "stella_solver.h"
#include "stella_bc.h"
/**
 * file: stella_interface.h
 *
 * This file contains routines to interface with PlasComCM.
 * These routines are meant to be called from Fortran.
 */


/**
 * Initializes efield solver
 *
 * @param[out] solverCtx pointer to the solver object allocated in this call
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
 * @param[out] ierr error code
 */
void stella_init(void **solverCtx, MPI_Comm *comm,
             int nGlobal[], int nProcs[], int nLocal[],
             int offset[], int stride[],
             int cartCoord[], int periodic[], int periodic_storage,
             int *nd, int *ng, int axisymmetric, PetscErrorCode *ierr);


/**
 * Sets state variables for efield solver
 *    -div dcoef grad phi + bcoef phi = rhs
 *     [n . dcoef grad phi] = jump
 *                      phi = g_{dirichlet} on \Gamma_D
 *               n grad phi = 0 on \Gamma_
 *
 * @param[in] solver_ctx handle for the solver object
 * @param[in] phi where we should store the solution
 * @param[in] dcoef grid function (see above)
 * @param[in] bcoef grid function (see above)
 * @param[in] jump grid function (see above)
 * @param[out] ierr error code
 */
void stella_set_state(void **solver_ctx, double **phi, double **dcoef,
                  double **bcoef, double **jump, PetscErrorCode *ierr);


/**
 * interface routine to set the right hand side
 *
 * @param[in] solver_ctx handle for the solver object
 * @param[in] rhs grid function to set the rhs
 * @param[out] ierr error code
 */
void stella_set_rhs(void **solver_ctx, double **rhs, PetscErrorCode *ierr);


/**
 * interface routine to set the analytic solution (for testing)
 *
 * @param[in] solver_ctx handle for the solver object
 * @param[in] sol grid function to set the solution
 * @param[out] ierr error code
 */
void stella_set_sol(void **solver_ctx, double sol[], PetscErrorCode *ierr);


/**
 * interface routine to set the analytic solution (for testing)
 *
 * @param[in] solver_ctx handle for the solver object
 * @param[in] is global indices of local grid starts by dimension
 * @param[in] ie global indices of local grid ends by dimension
 * @param[in] nCells number of d.o.f.
 * @param[in] xyz grid coordinates
 * @param[out] ierr error code
 */
void stella_set_grid(void **solver_ctx, int is[], int ie[], int *nCells,
                 double **xyz, PetscErrorCode *ierr);


/**
 * Used to tell the solver to setup the operator
 *
 * @param[in] solver_ctx handle for the solver object
 * @param[out] ierr error code
 */
void stella_setup_op(void **solver_ctx, PetscErrorCode *ierr);

/**
 * Runs the efield solver (see ModElectricSolve in ModElectric.fpp
 *
 * @param[in] solver_ctx handle for the solver object
 * @param[out] ierr error code
 */
void stella_solve(void **solver_ctx, PetscErrorCode *ierr);


/**
 * Used to add a boundary condition (called patch in PlasComCM
 *
 * @param[in] solver_ctx handle for the solver object
 * @param[out] ierr error code
 */
void stella_set_patch(void **solver_ctx, int is[], int ie[], int *bc_type,
                  int *norm_dir, double **dirichlet);


/**
 * Destroys efield data structures
 *
 * @param[in] solver_ctx handle for the solver object
 * @param[out] ierr error code
 */
void stella_cleanup(void **solver_ctx, PetscErrorCode *ierr);

#endif
