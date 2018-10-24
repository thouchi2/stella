#include "stella_util.h"


PetscErrorCode stella_store_external_array(stella *slv, double src[], Vec dest)
{
	PetscErrorCode ierr;

	if (slv->grid.nd == 2) {
		PetscScalar **dest_loc;
		double **src_loc;
		PetscInt i, j;
		PetscInt xs, ys, xm, ym;

		ierr = DMDAGetCorners(slv->dm, &xs, &ys, 0, &xm, &ym, 0);CHKERRQ(ierr);
		ierr = stella_dmap_get(slv->dmap, src, &src_loc);CHKERRQ(ierr);
		ierr = DMDAVecGetArray(slv->dm, dest, &dest_loc);CHKERRQ(ierr);
		for (j = ys; j < ys + ym; j++) {
			for (i = xs; i < xs + xm; i++) {
				dest_loc[j][i] = src_loc[j][i];
			}
		}
		ierr = stella_dmap_restore(slv->dmap, &src_loc);CHKERRQ(ierr);
		ierr = DMDAVecRestoreArray(slv->dm, dest, &dest_loc);CHKERRQ(ierr);
	} else { // slv->grid.nd == 3
		PetscScalar ***dest_loc;
		double ***src_loc;
		PetscInt i, j, k;
		PetscInt xs, ys, zs, xm, ym, zm;
		ierr = DMDAGetCorners(slv->dm, &xs, &ys, &zs, &xm, &ym, &zm);CHKERRQ(ierr);
		ierr = stella_dmap_get(slv->dmap, src, &src_loc);CHKERRQ(ierr);
		ierr = DMDAVecGetArray(slv->dm, dest, &dest_loc);CHKERRQ(ierr);
		for (k = zs; k < zs + zm; k++) {
			for (j = ys; j < ys + ym; j++) {
				for (i = xs; i < xs + xm; i++) {
					dest_loc[k][j][i] = src_loc[k][j][i];
				}
			}
		}
		ierr = stella_dmap_restore(slv->dmap, &src_loc);CHKERRQ(ierr);
		ierr = DMDAVecRestoreArray(slv->dm, dest, &dest_loc);CHKERRQ(ierr);
	}

	return 0;
}
