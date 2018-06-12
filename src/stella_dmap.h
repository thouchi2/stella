#ifndef STELLA_DMAP_H
#define STELLA_DMAP_H

#include "petscsys.h"


/**
 * Data mapping object.
 *
 * This object manages mapping data from a given source.  For example,
 * arrays on 2 or 3 dimensional grids may include ghost points or
 * have different strides.  This object provides a general interface from
 * which the efield solver can map contiguous data arrays from different
 * sources to two or three dimensional arrays on a structured grid.
 */
typedef struct stella_dmap {
	int xs, ys, zs, xm, ym, zm;
	int nd;        /**< number of dimensions of structured grid */
	int stride[3]; /**< strides of array that will be mapped */
} stella_dmap;


/**
 * Constructor for a data mapping object on a three dimensional grid.
 *
 * @param xs index of first grid point in the first dimension
 * @param ys index of first grid point in the second dimension
 * @param zs index of first grid point in the third dimension
 * @param xm number of grid points in the first dimension
 * @param ym number of grid points in the second dimension
 * @param zm number of grid points in the third dimension
 * @param stride stride of array in each dimension
 */
stella_dmap *stella_dmap_create_3d(int xs, int ys, int zs, int xm, int ym, int zm, int stride[]);


/**
 * Constructor for a data mapping object on a two dimensional grid.
 *
 * @param xs index of first grid point in the first dimension
 * @param ys index of first grid point in the second dimension
 * @param xm number of grid points in the first dimension
 * @param ym number of grid points in the second dimension
 * @param stride stride of array in each dimension
 */
stella_dmap *stella_dmap_create_2d(int xs, int ys, int xm, int ym, int stride[]);


/**
 * Gives a multidimensional array of mapped data indexed globally
 * on the two or three dimensional structured grid.
 *
 * @param mp data mapping object
 * @param src data mapping source array
 * @param array multidimensional array of mapped data
 */
PetscErrorCode stella_dmap_get(stella_dmap *mp, double *src, void *array);


/**
 * Frees memory created when using stella_dmap_get
 */
PetscErrorCode stella_dmap_restore(stella_dmap *mp, void *array);


/**
 * Gives a multidimensional array of mapped data indexed globally
 * on the two or three dimensional structured grid.
 *
 * @param mp data mapping object
 * @param src data mapping source array
 * @param array multidimensional array of mapped data
 */
PetscErrorCode stella_dmap_get_char(stella_dmap *mp, char *src, void *array);


/**
 * Frees memory created when using stella_dmap_get
 */
PetscErrorCode stella_dmap_restore_char(stella_dmap *mp, void *array);


PetscErrorCode stella_dmap_destroy(stella_dmap *mp);


#endif
