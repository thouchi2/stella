#ifndef STELLA_CLASSIFY_H
#define STELLA_CLASSIFY_H

#include <petscdmda.h>

#include "stella_dmap.h"


/**
 * Stella grid point types.
 */
typedef struct {
	char interior;
	char dirichlet;
	char neumann;
} stella_ptypes;


typedef struct {
	int owns_data;
	stella_ptypes ptypes;
	stella_dmap *dmap;
	char *classify;
} stella_classify;

/**
 * Helper function to use consume classify array from external library.
 *
 * If the classify array given from an external library has ghosts, use that directly without copying.
 * Otherwise, create a new classify array and fill ghost region (currently expensive).
 * @param[in] da DMDA object
 * @param[in] dmap data mapping object for mapping user application data
 * @param[out] classify_dmap data mapping object for using classify data
 * @param[out] classify array to be for classifying grid points
 * @params[in] ptypes key indicating values for classify array
 */
PetscErrorCode stella_classify_create(stella_classify **cls, DM da,
                                      stella_dmap *dmap_ext, char *classify_ext,
                                      stella_ptypes ptypes);


PetscErrorCode stella_classify_get(stella_classify *cls, void *array);
PetscErrorCode stella_classify_restore(stella_classify *cls, void *array);


/**
 * Clean up stella_classify data
 */
PetscErrorCode stella_classify_destroy(stella_classify *cls);

#endif
