#include "stella_classify.h"


static int classify_has_ghosts(stella_dmap *dmap)
{
	return !(stride[0] == xm);
}


/**
 * Create data mapper and classify array with room for halo region.
 */
static PetscErrorCode create_ghosted_classify(stella_classify *cls, DM da)
{
	PetscErrorCode ierr;
	PetscInt dim;

	ierr = DMGetDimension(da, &dim);CHKERRQ(ierr);

	if (dim == 2) {
		PetscInt xs, ys, xm, ym;
		ierr = DMDAGetCorners(da, &xs, &ys, 0, &xm, &ym, 0);CHKERRQ(ierr);

		int stride[2];
		stride[0] = xm + 2;
		stride[1] = ym + 2;
		cls->dmap = stella_dmap_create_2d(xs - 1, ys - 1, xm, ym, stride);CHKERRQ(ierr);
		cls->classify = (char*) malloc((xm+2)*(ym+2)*sizeof(char));
	} else {
		PetscInt xs, ys, zs, xm, ym, zm;
		ierr = DMDAGetCorners(da, &xs, &ys, &zs, &xm, &ym, &zm);CHKERRQ(ierr);

		int stride[3];
		stride[0] = xm + 2;
		stride[1] = ym + 2;
		stride[2] = zm + 2;
		cls->dmap = stella_dmap_create_3d(xs - 1, ys - 1, zs - 1, xm, ym, zm, stride);CHKERRQ(ierr);
		cls->classify = (char*) malloc((xm+2)*(ym+2)*(zm+2)*sizeof(char));
	}

	return 0;
}


/**
 * This uses DMDA to transmit halo region for newly created DMDA.
 * This involves an extra data copy and a vector of type PetscScalar (expensive).
 */
static PetscErrorCode setup_ghosted_classify(stella_classify *cls, DM da,
                                             stella_dmap *dmap_ext, char *classify_ext)
{
	PetscErrorCode ierr;

	// Use DMDA to transmit halo
	Vec g, l;
	PetscInt dim;

	ierr = DMGetDimension(da, &dim);CHKERRQ(ierr);
	ierr = DMCreateGlobalVector(da, &g);CHKERRQ(ierr);
	ierr = DMCreateLocalVector(da, &l);CHKERRQ(ierr);

	if (dim == 2) {
		PetscScalar **arr;
		char **classify;
		PetscInt i, j;
		PetscInt xs, ys, xm, ym;

		ierr = DMDAGetCorners(da, &xs, &ys, 0, &xm, &ym, 0);CHKERRQ(ierr);
		ierr = DMDAVecGetArray(da, g, &arr);CHKERRQ(ierr);
		ierr = stella_dmap_get_char(dmap_ext, classify_ext, &classify);CHKERRQ(ierr);
		for (j = ys; j < ys + ym; j++) {
			for (i = xs; i < xs + xm; i++) {
				arr[j][i] = classify[j][i];
			}
		}
		ierr = DMDAVecRestoreArray(da, g, &arr);CHKERRQ(ierr);
		ierr = stella_dmap_restore_char(dmap_ext, &classify);CHKERRQ(ierr);
	} else {
		PetscScalar ***arr;
		char ***classify;
		PetscInt i, j, k;
		PetscInt xs, ys, zs, xm, ym, zm;

		ierr = DMDAGetCorners(da, &xs, &ys, &zs, &xm, &ym, &zm);CHKERRQ(ierr);
		ierr = DMDAVecGetArray(da, g, &arr);CHKERRQ(ierr);
		ierr = stella_dmap_get_char(dmap_ext, classify_ext, &classify);CHKERRQ(ierr);
		for (k = zs; k < zs + zm; k++) {
			for (j = ys; j < ys + ym; j++) {
				for (i = xs; i < xs + xm; i++) {
					arr[k][j][i] = classify[k][j][i];
				}
			}
		}
		ierr = DMDAVecRestoreArray(da, g, &arr);CHKERRQ(ierr);
		ierr = stella_dmap_restore_char(dmap_ext, &classify);CHKERRQ(ierr);
	}

	// Transmit halo
	ierr = DMGlobalToLocalBegin(da, g, INSERT_VALUES, l);CHKERRQ(ierr);
	ierr = DMGlobalToLocalEnd(da, g, INSERT_VALUES, l);CHKERRQ(ierr);

	if (dim == 2) {
		PetscScalar **arr;
		char **classify;
		PetscInt i, j;
		PetscInt xs, ys, xm, ym;
		ierr = DMDAGetCorners(da, &xs, &ys, 0, &xm, &ym, 0);CHKERRQ(ierr);

		ierr = DMDAVecGetArray(da, l, &arr);CHKERRQ(ierr);
		ierr = stella_classify_get(cls, &classify);CHKERRQ(ierr);
		for (j = ys - 1; j < ys + ym + 1; j++) {
			for (i = xs - 1; i < xs + xm + 1; i++) {
				classify[j][i] = arr[j][i];
			}
		}
		ierr = DMDAVecRestoreArray(da, l, &arr);CHKERRQ(ierr);
		ierr = stella_classify_restore(cls, &classify);CHKERRQ(ierr);
	} else {
		PetscScalar ***arr;
		char ***classify;
		PetscInt i, j, k;
		PetscInt xs, ys, zs, xm, ym, zm;
		ierr = DMDAGetCorners(da, &xs, &ys, &zs, &xm, &ym, &zm);CHKERRQ(ierr);

		ierr = DMDAVecGetArray(da, l, &arr);CHKERRQ(ierr);
		ierr = stella_classify_get(cls, &classify);CHKERRQ(ierr);
		for (k = zs - 1; k < zs + zm + 1; k++) {
			for (j = ys - 1; j < ys + ym + 1; j++) {
				for (i = xs - 1; i < xs + xm + 1; i++) {
					classify[k][j][i] = arr[k][j][i];
				}
			}
		}
		ierr = DMDAVecRestoreArray(da, l, &arr);CHKERRQ(ierr);
		ierr = stella_classify_restore(cls, &classify);CHKERRQ(ierr);
	}

	return 0;
}


PetscErrorCode stella_classify_create(stella_classify **cls_ret, DM da,
                                      stella_dmap *dmap_ext, char *classify_ext,
                                      stella_ptypes ptypes)
{
	PetscErrorCode ierr;

	stella_classify *cls = (stella_classify*) malloc(sizeof(stella_classify));
	cls->ptypes = ptypes;

	if (classify_has_ghosts(dmap_ext)) {
		// use application array for classify
		cls->owns_data = 0;
		cls->dmap = dmap_ext;
		cls->classify = classify_ext;
	} else {
		cls->owns_data = 1;
		ierr = create_ghosted_classify(cls, da);CHKERRQ(ierr);
		ierr = setup_ghosted_classify(cls, da, dmap_ext, classify_ext);CHKERRQ(ierr);
	}

	*cls_ret = cls;

	return 0;
}


PetscErrorCode stella_classify_get(stella_classify *cls, void *array)
{
	PetscErrorCode ierr;

	ierr = stella_dmap_get_char(cls->dmap, cls->classify, &array);CHKERRQ(ierr);

	return 0;
}


PetscErrorCode stella_classify_restore(stella_classify *cls, void *array)
{
	PetscErrorCode ierr;

	ierr = stella_dmap_restore_char(cls->dmap, &array);CHKERRQ(ierr);

	return 0;
}


PetscErrorCode stella_classify_destroy(stella_classify *cls)
{
	PetscErrorCode ierr;

	if (cls->owns_data) {
		free(cls->classify);
		ierr = stella_dmap_destroy(cls->dmap);CHKERRQ(ierr);
		free(cls->dmap);
	}

	return 0;
}
