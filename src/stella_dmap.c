#include <stdlib.h>

#include "stella_dmap.h"


static PetscErrorCode get_3d_double(stella_dmap *idx, double *src, double ****dest)
{
	PetscErrorCode ierr;
	double **b;
	int i, j;
	ierr = PetscMalloc1(idx->stride[2]*sizeof(double**) + idx->stride[2]*idx->stride[1], dest);CHKERRQ(ierr);

	b = (double**)((*dest) + idx->stride[2]);
	for (i = 0; i < idx->stride[2]; i++) (*dest)[i] = b + i*idx->stride[1] - idx->ys;
	for (i = 0; i < idx->stride[2]; i++) {
		for (j = 0; j < idx->stride[1]; j++) {
			b[i*idx->stride[1] + j] = src + i*idx->stride[0]*idx->stride[1] + j*idx->stride[0] - idx->xs;
		}
	}

	*dest -= idx->zs;

	return 0;
}


static PetscErrorCode get_3d_char(stella_dmap *idx, char *src, char ****dest)
{
	PetscErrorCode ierr;
	char **b;
	int i, j;
	ierr = PetscMalloc1(idx->stride[2]*sizeof(char**) + idx->stride[2]*idx->stride[1], dest);CHKERRQ(ierr);

	b = (char**)((*dest) + idx->stride[2]);
	for (i = 0; i < idx->stride[2]; i++) (*dest)[i] = b + i*idx->stride[1] - idx->ys;
	for (i = 0; i < idx->stride[2]; i++) {
		for (j = 0; j < idx->stride[1]; j++) {
			b[i*idx->stride[1] + j] = src + i*idx->stride[0]*idx->stride[1] + j*idx->stride[0] - idx->xs;
		}
	}

	*dest -= idx->zs;

	return 0;
}


static PetscErrorCode restore_3d_double(stella_dmap *idx, double ****arr)
{
	PetscErrorCode ierr;
	void *dummy;

	dummy = (void*)(*arr + idx->zs);
	ierr = PetscFree(dummy);CHKERRQ(ierr);

	return 0;
}


static PetscErrorCode restore_3d_char(stella_dmap *idx, char ****arr)
{
	PetscErrorCode ierr;
	void *dummy;

	dummy = (void*)(*arr + idx->zs);
	ierr = PetscFree(dummy);CHKERRQ(ierr);

	return 0;
}


static PetscErrorCode get_2d_double(stella_dmap *idx, double *src, double ***dest)
{
	PetscErrorCode ierr;
	int i;

	ierr = PetscMalloc1(idx->stride[1], dest);CHKERRQ(ierr);

	for (i = 0; i < idx->stride[1]; i++) {
		(*dest)[i] = src + i*idx->stride[0] - idx->xs;
	}
	*dest -= idx->ys;

	return 0;
}


static PetscErrorCode get_2d_char(stella_dmap *idx, char *src, char ***dest)
{
	PetscErrorCode ierr;
	int i;

	ierr = PetscMalloc1(idx->stride[1], dest);CHKERRQ(ierr);

	for (i = 0; i < idx->stride[1]; i++) {
		(*dest)[i] = src + i*idx->stride[0] - idx->xs;
	}
	*dest -= idx->ys;

	return 0;
}


static PetscErrorCode restore_2d_double(stella_dmap *idx, double ***arr)
{
	PetscErrorCode ierr;
	void *dummy;

	dummy = (void*)(*arr + idx->ys);
	ierr = PetscFree(dummy);CHKERRQ(ierr);

	return 0;
}


static PetscErrorCode restore_2d_char(stella_dmap *idx, char ***arr)
{
	PetscErrorCode ierr;
	void *dummy;

	dummy = (void*)(*arr + idx->ys);
	ierr = PetscFree(dummy);CHKERRQ(ierr);

	return 0;
}


stella_dmap *stella_dmap_create_3d(int xs, int ys, int zs, int xm, int ym, int zm, int stride[])
{
	stella_dmap *mp = (stella_dmap*) malloc(sizeof(stella_dmap));

	mp->xs = xs; mp->ys = ys; mp->zs = zs;
	mp->xm = xm; mp->ym = ym; mp->zm = zm;
	mp->nd = 3;
	mp->stride[0] = stride[0];
	mp->stride[1] = stride[1];
	mp->stride[2] = stride[2];

	return mp;
}


stella_dmap *stella_dmap_create_2d(int xs, int ys, int xm, int ym, int stride[])
{
	stella_dmap *mp = (stella_dmap*) malloc(sizeof(stella_dmap));
	mp->xs = xs; mp->ys = ys;
	mp->xm = xm; mp->ym = ym;
	mp->nd = 2;
	mp->stride[0] = stride[0];
	mp->stride[1] = stride[1];

	return mp;
}


PetscErrorCode stella_dmap_get(stella_dmap *mp, double *src, void *array)
{
	PetscErrorCode ierr;

	if (mp->nd == 3) {
		ierr = get_3d_double(mp, src, (double ****) array);CHKERRQ(ierr);
	} else {
		ierr = get_2d_double(mp, src, (double ***) array);CHKERRQ(ierr);
	}

	return 0;
}


PetscErrorCode stella_dmap_get_char(stella_dmap *mp, char *src, void *array)
{
	PetscErrorCode ierr;

	if (mp->nd == 3) {
		ierr = get_3d_char(mp, src, (char ****) array);CHKERRQ(ierr);
	} else {
		ierr = get_2d_char(mp, src, (char ***) array);CHKERRQ(ierr);
	}

	return 0;
}


PetscErrorCode stella_dmap_restore(stella_dmap *mp, void *array)
{
	PetscErrorCode ierr;

	if (mp->nd == 3) {
		ierr = restore_3d_double(mp, (double ****) array);CHKERRQ(ierr);
	} else {
		ierr = restore_2d_double(mp, (double ***) array);CHKERRQ(ierr);
	}

	return 0;
}


PetscErrorCode stella_dmap_restore_char(stella_dmap *mp, void *array)
{
	PetscErrorCode ierr;

	if (mp->nd == 3) {
		ierr = restore_3d_char(mp, (char ****) array);CHKERRQ(ierr);
	} else {
		ierr = restore_2d_char(mp, (char ***) array);CHKERRQ(ierr);
	}

	return 0;
}


PetscErrorCode stella_dmap_destroy(stella_dmap *mp)
{
	return 0;
}

