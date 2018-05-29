#ifndef STELLA_IO_H
#define STELLA_IO_H

#include <petscdmda.h>

/**
 * file: stella_io.h
 *
 * General input/output utility functions for efield solver
 */


/**
 * Writes a function over a structured grid to a vtk file
 *
 * @param da DM object that holds structured grid information
 * @param v function to output
 * @param name filename for vtkfile
 * @param grid_id unique id for multiple grid cases
 * @param ts time step
 */
PetscErrorCode stella_io_vtkwrite(DM da, Vec v, char *name, int grid_id, int ts);


/**
 * Writes estimation of the operators eigenvalues
 */
PetscErrorCode stella_io_eigwrite(PetscReal r[], PetscReal c[], PetscInt neig);


/**
 * Prints message to console
 */
PetscErrorCode stella_io_print(MPI_Comm comm, char *msg);

#endif
