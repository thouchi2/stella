#ifndef STELLA_UTIL_H
#define STELLA_UTIL_H

#include "stella_solver.h"

/**
 * Copy external array to PETSc vec
 */
PetscErrorCode stella_store_external_array(stella *slv, double src[], Vec dest);

#endif
