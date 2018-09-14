#ifndef STELLA_GRID_H
#define STELLA_GRID_H

#include <petscdmda.h>
#include "mpi.h"

#ifdef WITH_CEDAR
#include <cedar/capi.h>
#endif

/**
 * @file: stella_grid.h
 *
 * Class used to store grid information from PlasComCM
 */


/**
 * Holds grid information from PlasComCM
 */
typedef struct {
	int nd;        /**< number of dimensions */
	int id;        /**< grid id */
	int is[3];     /**< global indicies of local grid starts by dimension */
	int ie[3];     /**< global indicies of local grid ends by dimension */
	int num_cells; /**< cells number of d.o.f. in grid */
	int overlap_periodic;  /**< are there duplicate grid points if periodic */
	double *xyz;   /**< grid coordinates */
	#ifdef WITH_CEDAR
	cedar_topo topo;/**< holds cedar topology handle */
	#endif
} stella_grid;


/**
 * Sets up the grid data structure
 * The primary responsibilities are:
 *   1. Renumber ranks to match PlasComCM rank to Cartesian grid ordering
 *   2. Create PETSc DMDA object with parallel grid topology info from PlasComCM
 *   3. Initialize DMDA coordinates (overwritten by efs_set_grid in stella_solver.c)
 *
 * @param grid grid object
 * @param dm PETSc DM object to setup
 * @param[in/out] comm MPI Communicator from PlasComCM
 * @param nGlobal number of global grid points by dimension
 * @param nProcs number of processors by dimension
 * @param nLocal number of local grid points by dimension
 * @param cartCoord Cartesian coordinates of grid
 * @param periodic flag to set periodicity by dimension
 * @param periodic_storage flag for duplicate grid points if periodic
 * @param nd number of dimensions
 */
PetscErrorCode stella_grid_setup(stella_grid *grid, DM *dm, MPI_Comm *comm, int nGlobal[], int nProcs[], int nLocal[],
                             int cartCoord[], int periodic[], int periodic_storage, int nd);


#endif
