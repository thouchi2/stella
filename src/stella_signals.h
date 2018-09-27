#ifndef STELLA_SIGNALS_H
#define STELLA_SIGNALS_H

#include "stella_solver.h"

#ifdef __cplusplus
extern "C" {
#endif

/**
 * Signal boundary conditions have changed
 */
PetscErrorCode stella_changed_bc(stella * slv);


/**
 * Signal rhs has changed
 */
PetscErrorCode stella_changed_rhs(stella *slv);


/**
 * Signal dcoef has changed
 *
 * Changing dcoef is expensive, requiring the following:
 * 1. Reassembly of the linear operator (rerun MG setup phase)
 * 2. Reprocess boundary conditions
 */
PetscErrorCode stella_changed_dcoef(stella *slv);


#ifdef __cplusplus
}
#endif

#endif
