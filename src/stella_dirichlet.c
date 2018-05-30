#include "stella_dirichlet.h"
#include "stella_mat.h"
#include "stella_operator.h"


/**
 * Deletes connections from interior points that neighbor a Dirichlet boundary.
 * This step is needed to get a symmetric operator.
 */
static PetscErrorCode update_dirichlet_sym(stella_bc *bc, Mat A, DM da)
{
	PetscErrorCode ierr;
	PetscInt i,j,cnt;
	PetscInt ngx, ngy;
	PetscScalar v[9];
	MatStencil row[9], col;
	PetscInt level;
	stella_patch *patch = bc->patch;

	ierr = DMGetCoarsenLevel(da, &level);CHKERRQ(ierr);

	ierr = DMDAGetInfo(da, 0, &ngx, &ngy, 0,0,0,0,0,0,0,0,0,0);CHKERRQ(ierr);

	for (i = 0; i < 8; i++) v[i] = 0;
	for (j = patch->corners.is[1]; j <= patch->corners.ie[1]; j++) {
		for (i = patch->corners.is[0]; i <= patch->corners.ie[0]; i++) {
			col.i = i; col.j = j;
			cnt = 0;

			if (i-1 < patch->corners.is[0] &&
				i-1 >= 0) {
				row[cnt].i = i-1; row[cnt].j = j;
				cnt++;
				if (j-1 >= 0) {
					row[cnt].i = i-1; row[cnt].j = j-1;
					cnt++;
				}
				if (j+1 <= ngy-1) {
					row[cnt].i = i-1; row[cnt].j = j+1;
					cnt++;
				}
			}
			if (i+1 > patch->corners.ie[0] &&
				i+1 <= ngx-1) {
				row[cnt].i = i+1; row[cnt].j = j;
				cnt++;
				if (j-1 >= 0) {
					row[cnt].i = i+1; row[cnt].j = j-1;
					cnt++;
				}
				if (j+1 <= ngy-1) {
					row[cnt].i = i+1; row[cnt].j = j+1;
					cnt++;
				}
			}
			if (j-1 < patch->corners.is[1] &&
			    j-1 >= 0) {
				row[cnt].i = i; row[cnt].j = j-1;
				cnt++;
				if (i-1 >= 0) {
					row[cnt].i = i-1; row[cnt].j = j-1;
					cnt++;
				}
				if (i+1 <= ngx - 1) {
					row[cnt].i = i+1; row[cnt].j = j-1;
					cnt++;
				}
			}
			if (j+1 > patch->corners.ie[1] &&
				j+1 <= ngy-1) {
				row[cnt].i = i; row[cnt].j = j+1;
				cnt++;
				if (i-1 >= 0) {
					row[cnt].i = i-1; row[cnt].j = j+1;
					cnt++;
				}
				if (i+1 <= ngx - 1) {
					row[cnt].i = i+1; row[cnt].j = j+1;
					cnt++;
				}
			}

			ierr = MatSetValuesStencil(A, cnt, row, 1, &col, v, INSERT_VALUES);CHKERRQ(ierr);
		}
	}

	return 0;
}


static PetscErrorCode update_dirichlet_sym_3d(stella_bc *bc, Mat A, DM da)
{
	PetscErrorCode ierr;
	PetscInt i, j, k, cnt;
	PetscInt ngx, ngy, ngz;
	PetscScalar v[64];
	MatStencil row[64], col;
	stella_patch *patch = bc->patch;

	ierr = DMDAGetInfo(da, 0, &ngx, &ngy, &ngz,0,0,0,0,0,0,0,0,0);CHKERRQ(ierr);

	for (i = 0; i < 19; i++) v[i] = 0;
	for (k = patch->corners.is[2]; k <= patch->corners.ie[2]; k++) {
		for (j = patch->corners.is[1]; j <= patch->corners.ie[1]; j++) {
			for (i = patch->corners.is[0]; i <= patch->corners.ie[0]; i++) {
				col.i = i; col.j = j; col.k = k;
				cnt = 0;

				if (i-1 < patch->corners.is[0] && i-1 >= 0) {
					row[cnt].i = i-1; row[cnt].j = j; row[cnt].k = k;
					cnt++;
					if (j-1 >= 0) {
						row[cnt].i = i-1; row[cnt].j = j-1; row[cnt].k = k;
						cnt++;
					}
					if (j+1 <= ngy-1) {
						row[cnt].i = i-1; row[cnt].j = j+1; row[cnt].k = k;
						cnt++;
					}
					if (k-1 >= 0) {
						row[cnt].i = i-1; row[cnt].j = j; row[cnt].k = k-1;
						cnt++;
					}
					if (k+1 <= ngz - 1) {
						row[cnt].i = i-1; row[cnt].j = j; row[cnt].k = k+1;
						cnt++;
					}
				}
				if (i+1 > patch->corners.ie[0] && i+1 <= ngx - 1) {
					row[cnt].i = i+1; row[cnt].j = j; row[cnt].k = k;
					cnt++;
					if (j-1 >= 0) {
						row[cnt].i = i+1; row[cnt].j = j-1; row[cnt].k = k;
						cnt++;
					}
					if (j+1 <= ngy-1) {
						row[cnt].i = i+1; row[cnt].j = j+1; row[cnt].k = k;
						cnt++;
					}
					if (k-1 >= 0) {
						row[cnt].i = i+1; row[cnt].j = j; row[cnt].k = k-1;
						cnt++;
					}
					if (k+1 <= ngz - 1) {
						row[cnt].i = i+1; row[cnt].j = j; row[cnt].k = k+1;
						cnt++;
					}
				}
				if (j-1 < patch->corners.is[1] && j-1 >= 0) {
					row[cnt].i = i; row[cnt].j = j-1; row[cnt].k = k;
					cnt++;
					if (i-1 >= 0) {
						row[cnt].i = i-1; row[cnt].j = j-1; row[cnt].k = k;
						cnt++;
					}
					if (i+1 <= ngx - 1) {
						row[cnt].i = i+1; row[cnt].j = j-1; row[cnt].k = k;
						cnt++;
					}
					if (k-1 >= 0) {
						row[cnt].i = i; row[cnt].j = j-1; row[cnt].k = k-1;
						cnt++;
					}
					if (k+1 <= ngz - 1) {
						row[cnt].i = i; row[cnt].j = j-1; row[cnt].k = k+1;
						cnt++;
					}
				}
				if (j+1 > patch->corners.ie[1] && j+1 <= ngy - 1) {
					row[cnt].i = i; row[cnt].j = j+1; row[cnt].k = k;
					cnt++;
					if (i-1 >= 0) {
						row[cnt].i = i-1; row[cnt].j = j+1; row[cnt].k = k;
						cnt++;
					}
					if (i+1 <= ngx - 1) {
						row[cnt].i = i+1; row[cnt].j = j+1; row[cnt].k = k;
						cnt++;
					}
					if (k-1 >= 0) {
						row[cnt].i = i; row[cnt].j = j+1; row[cnt].k = k-1;
						cnt++;
					}
					if (k+1 <= ngz - 1) {
						row[cnt].i = i; row[cnt].j = j+1; row[cnt].k = k+1;
						cnt++;
					}
				}
				if (k-1 < patch->corners.is[2] && k-1 >= 0) {
					row[cnt].i = i; row[cnt].j = j; row[cnt].k = k-1;
					cnt++;
					if (i-1 >= 0) {
						row[cnt].i = i-1; row[cnt].j = j; row[cnt].k = k-1;
						cnt++;
					}
					if (i+1 <= ngx - 1) {
						row[cnt].i = i+1; row[cnt].j = j; row[cnt].k = k-1;
						cnt++;
					}
					if (j-1 >= 0) {
						row[cnt].i = i; row[cnt].j = j-1; row[cnt].k = k-1;
						cnt++;
					}
					if (j+1 <= ngy - 1) {
						row[cnt].i = i; row[cnt].j = j+1; row[cnt].k = k-1;
						cnt++;
					}
				}
				if (k+1 > patch->corners.ie[2] && k+1 <= ngz - 1) {
					row[cnt].i = i; row[cnt].j = j; row[cnt].k = k+1;
					cnt++;
					if (i-1 >= 0) {
						row[cnt].i = i-1; row[cnt].j = j; row[cnt].k = k+1;
						cnt++;
					}
					if (i+1 <= ngx - 1) {
						row[cnt].i = i+1; row[cnt].j = j; row[cnt].k = k+1;
						cnt++;
					}
					if (j-1 >= 0) {
						row[cnt].i = i; row[cnt].j = j-1; row[cnt].k = k+1;
						cnt++;
					}
					if (j+1 <= ngy - 1) {
						row[cnt].i = i; row[cnt].j = j+1; row[cnt].k = k+1;
						cnt++;
					}
				}

				ierr = MatSetValuesStencil(A, cnt, row, 1, &col, v, INSERT_VALUES);CHKERRQ(ierr);
			}
		}
	}

	return 0;
}



/**
 * Eliminates connections from boundary to neighboring points.
 * Saves the weight of this connection as additive contribution to RHS.
 */
static PetscErrorCode apply_dirichlet(stella_bc *bc, Mat A, DM da)
{
	PetscErrorCode ierr;
	DM cda;
	DMDACoor2d **coors;
	Vec lc;
	PetscInt i,j,cnt;
	PetscInt ngx, ngy;
	PetscScalar v[9];
	MatStencil row, col[9];
	stella_patch *patch = bc->patch;
	MatType mtype;
	PetscBool is_boxmg;
	PetscScalar **acont;
	stella_metric *met;
	PetscScalar **dcoef;
	PetscScalar **coef[5];
	PetscInt xs,ys,xe,ye,ym,xm;
	int stride[2];
	int offset[2];
	double **dirichlet;
	stella_dmap *dmap;

	met = bc->level->metric;

	ys = patch->corners.is[1]; xs = patch->corners.is[0];
	ye = patch->corners.ie[1]+1; xe = patch->corners.ie[0]+1;
	ym = patch->corners.ie[1] - patch->corners.is[1] + 1;
	xm = patch->corners.ie[0] - patch->corners.is[0] + 1;

	stride[0] = patch->stride[0];
	stride[1] = patch->stride[1];
	offset[0] = patch->offset[0];
	offset[1] = patch->offset[1];
	dmap = stella_dmap_create_2d(xs - offset[0], ys - offset[1], xm, ym, stride);
	ierr = stella_dmap_get(dmap, patch->dirichlet, &dirichlet);CHKERRQ(ierr);

	ierr = MatGetType(A, &mtype);CHKERRQ(ierr);
	ierr = DMDAGetInfo(da, 0, &ngx, &ngy, 0,0,0,0,0,0,0,0,0,0);CHKERRQ(ierr);
	ierr = DMDAVecGetArray(da, bc->level->ladd_cont, &acont);CHKERRQ(ierr);
	ierr = DMDAVecGetArray(da, bc->level->ldcoef, &dcoef);CHKERRQ(ierr);

	if (bc->axisymmetric) {
		ierr = DMGetCoordinateDM(da, &cda);CHKERRQ(ierr);
		ierr = DMGetCoordinatesLocal(da, &lc);CHKERRQ(ierr);
		ierr = DMDAVecGetArray(cda, lc, &coors);CHKERRQ(ierr);
	}

	for (i = 0; i < 5; i++) {
		ierr = DMDAVecGetArray(da, met->lcoef[i], &coef[i]);CHKERRQ(ierr);
	}

	int ibeg = patch->corners.is[0];
	int jbeg = patch->corners.is[1];
	int iend = patch->corners.ie[0];
	int jend = patch->corners.ie[1];

	v[0] = 1.0;
	for (j = jbeg; j <= jend; j++) {
		for (i = ibeg; i <= iend; i++) {
			row.i = i; row.j = j;
			col[0].i = i; col[0].j = j;
			cnt = 1;

			if (i != 0) {
				col[cnt].i = i-1; col[cnt].j = j;
				v[cnt] = 0; cnt++;
				if ((i-1) < patch->corners.is[0]) {
					double weight = coef[MET_W][j][i]*have(dcoef[j][i],dcoef[j][i-1])*dirichlet[j][i];
					if (bc->axisymmetric) weight *= (coors[j][i].x + coors[j][i-1].x) * .5;
					acont[j][i-1] += weight;
				}
			}
			if (j != 0) {
				col[cnt].i = i; col[cnt].j = j-1;
				v[cnt] = 0; cnt++;
				if ((j-1) < patch->corners.is[1]) {
					double weight = coef[MET_S][j][i]*have(dcoef[j][i],dcoef[j-1][i])*dirichlet[j][i];
					if (bc->axisymmetric) weight *= coors[j][i].x;
					acont[j-1][i] += weight;
				}
			}
			// SW
			if (i != 0 && j != 0) {
				col[cnt].i = i-1; col[cnt].j = j-1;
				v[cnt] = 0; cnt++;
				if ((i-1) < patch->corners.is[0] || (j-1) < patch->corners.is[1]) {
					double weight = coef[MET_SW][j][i]*have(have(dcoef[j][i], dcoef[j][i-1]),
					                                  have(dcoef[j-1][i], dcoef[j-1][i-1])) * dirichlet[j][i];
					acont[j-1][i-1] += weight;
				}
			}
			if (i != ngx - 1) {
				col[cnt].i = i+1; col[cnt].j = j;
				v[cnt] = 0; cnt++;
				if ((i+1) > patch->corners.ie[0]) {
					double weight = coef[MET_W][j][i+1]*have(dcoef[j][i],dcoef[j][i+1])*dirichlet[j][i];
					if (bc->axisymmetric) weight *= (coors[j][i+1].x + coors[j][i].x) * .5;
					acont[j][i+1] += weight;
				}
			}
			// SE
			if (i != ngx - 1 && j != 0) {
				col[cnt].i = i+1; col[cnt].j = j-1;
				v[cnt] = 0; cnt++;
				if ((i+1) > patch->corners.ie[0] || (j-1) < patch->corners.is[1]) {
					double weight = coef[MET_SE][j][i]*have(have(dcoef[j][i], dcoef[j][i+1]),
					                                        have(dcoef[j-1][i], dcoef[j-1][i+1]))*dirichlet[j][i];
					acont[j-1][i+1] += weight;
				}
			}
			// NW
			if (i !=0 && j != ngy - 1) {
				col[cnt].i = i-1; col[cnt].j = j+1;
				v[cnt] = 0; cnt++;
				if  ((i-1) < patch->corners.is[0] || (j+1) > patch->corners.ie[1]) {
					double weight = coef[MET_SE][j+1][i-1]*have(have(dcoef[j][i], dcoef[j][i-1]),
					                                            have(dcoef[j+1][i], dcoef[j+1][i-1]))*dirichlet[j][i];
					acont[j+1][i-1] += weight;
				}
			}
			if (j != ngy - 1) {
				col[cnt].i = i; col[cnt].j = j+1;
				v[cnt] = 0; cnt++;
				if ((j+1) > patch->corners.ie[1]) {
					double weight = coef[MET_S][j+1][i]*have(dcoef[j][i],dcoef[j+1][i])*dirichlet[j][i];
					if (bc->axisymmetric) weight *= coors[j][i].x;
					acont[j+1][i] += weight;
				}
			}
			// NE
			if ((i != ngx - 1) && (j != ngy - 1)) {
				col[cnt].i = i+1; col[cnt].j = j+1;
				v[cnt] = 0; cnt++;
				if ((i+1) > patch->corners.ie[0] || (j+1) > patch->corners.ie[1]) {
					double weight = coef[MET_SW][j+1][i+1]*have(have(dcoef[j][i], dcoef[j][i+1]),
					                                            have(dcoef[j+1][i], dcoef[j+1][i+1]))*dirichlet[j][i];
					acont[j+1][i+1] += weight;
				}
			}

			ierr = PetscStrcmp(mtype, MATSHELL, &is_boxmg);CHKERRQ(ierr);
			if (is_boxmg) {
				ierr = stella_bmg_SetValuesStencil(A, 1, &row, cnt, col, v, INSERT_VALUES);CHKERRQ(ierr);
			} else {
				ierr = MatSetValuesStencil(A, 1, &row, cnt, col, v, INSERT_VALUES);CHKERRQ(ierr);
			}
		}
	}

	ierr = DMDAVecRestoreArray(da, bc->level->ladd_cont, &acont);CHKERRQ(ierr);
	ierr = DMDAVecRestoreArray(da, bc->level->ldcoef, &dcoef);CHKERRQ(ierr);
	for (i = 0; i < 5; i++) {
		ierr = DMDAVecRestoreArray(da, met->lcoef[i], &coef[i]);CHKERRQ(ierr);
	}

	if (bc->axisymmetric) {
		ierr = DMDAVecRestoreArray(cda, lc, &coors);CHKERRQ(ierr);
	}

	ierr = stella_dmap_restore(dmap, &dirichlet);CHKERRQ(ierr);
	stella_dmap_destroy(dmap); free(dmap);

	return 0;
}


static PetscErrorCode apply_dirichlet_3d(stella_bc *bc, Mat A, DM da)
{
	PetscErrorCode ierr;
	PetscInt i,j,k,cnt;
	PetscInt ngx, ngy, ngz;
	PetscScalar v[19];
	MatStencil row, col[19];
	MatType mtype;
	PetscBool is_boxmg;
	double dcoefh[19];
	stella_patch *patch = bc->patch;
	PetscScalar ***acont;
	stella_metric *met;
	PetscScalar ***dcoef;
	PetscScalar ***coef[10];
	PetscInt xs, ys, zs, xe, ye, ze, xm, ym, zm;
	int stride[3];
	int offset[3];
	double ***dirichlet;
	stella_dmap *dmap;

	ierr = MatGetType(A, &mtype);
	PetscErrorCode (*set_stencil)(Mat mat, PetscInt m, const MatStencil idxm[], PetscInt n,
	                              const MatStencil idxn[], const PetscScalar v[],
	                              InsertMode addv);
	ierr = PetscStrcmp(mtype, MATSHELL, &is_boxmg);CHKERRQ(ierr);
	if (is_boxmg) {
		set_stencil = &stella_bmg_SetValuesStencil;
	} else {
		set_stencil = &MatSetValuesStencil;
	}

	met = bc->level->metric;
	ys = patch->corners.is[1]; xs = patch->corners.is[0]; zs = patch->corners.is[2];
	ye = patch->corners.ie[1]+1; xe = patch->corners.ie[0]+1; ze = patch->corners.ie[2]+1;
	ym = patch->corners.ie[1] - patch->corners.is[1] + 1;
	xm = patch->corners.ie[0] - patch->corners.is[0] + 1;
	zm = patch->corners.ie[2] - patch->corners.is[2] + 1;
	stride[0] = patch->stride[0]; stride[1] = patch->stride[1]; stride[2] = patch->stride[2];
	offset[0] = patch->offset[0]; offset[1] = patch->offset[1]; offset[2] = patch->offset[2];

	dmap = stella_dmap_create_3d(xs - offset[0], ys - offset[1], zs - offset[2], xm, ym, zm, stride);
	ierr = stella_dmap_get(dmap, patch->dirichlet, &dirichlet);CHKERRQ(ierr);
	ierr = DMDAGetInfo(da, 0, &ngx, &ngy, &ngz,0,0,0,0,0,0,0,0,0);CHKERRQ(ierr);
	ierr = DMDAVecGetArray(da, bc->level->ladd_cont, &acont);CHKERRQ(ierr);
	ierr = DMDAVecGetArray(da, bc->level->ldcoef, &dcoef);CHKERRQ(ierr);

	enum dir {
		O = 0, N, S, E, W, F, B, NW, NE, SW, SE, NF, NB, SF, SB, FE, FW, BE, BW
	};

	for (i = 0; i < 10; i++) {
		ierr = DMDAVecGetArray(da, met->lcoef[i], &coef[i]);CHKERRQ(ierr);
	}

	int ibeg = patch->corners.is[0];
	int jbeg = patch->corners.is[1];
	int kbeg = patch->corners.is[2];
	int iend = patch->corners.ie[0];
	int jend = patch->corners.ie[1];
	int kend = patch->corners.ie[2];
        int il, jl, ir, jr, kl, kr;

	v[0] = 1.0;
	for (k = kbeg; k <= kend; k++) {
		for (j = jbeg; j <= jend; j++) {
			for (i = ibeg; i <= iend; i++) {

                          il = i>0     ? i-1 : 0;
                          jl = j>0     ? j-1 : 0;
                          kl = k>0     ? k-1 : 0;
                          ir = i<ngx-1 ? i+1 : ngx-1;
                          jr = j<ngy-1 ? j+1 : ngy-1;
                          kr = k<ngz-1 ? k+1 : ngz-1;

				dcoefh[N] = have(dcoef[k][j][i], dcoef[k][jr][i]);
				dcoefh[S] = have(dcoef[k][j][i], dcoef[k][jl][i]);
				dcoefh[W] = have(dcoef[k][j][i], dcoef[k][j][il]);
				dcoefh[E] = have(dcoef[k][j][i], dcoef[k][j][ir]);
				dcoefh[F] = have(dcoef[k][j][i], dcoef[kr][j][i]);
				dcoefh[B] = have(dcoef[k][j][i], dcoef[kl][j][i]);
				dcoefh[SW] = have(dcoefh[W],
				                have(dcoef[k][jl][i], dcoef[k][jl][il]));
				dcoefh[NW] = have(dcoefh[W],
				                have(dcoef[k][jr][i], dcoef[k][jr][il]));
				dcoefh[SE] = have(dcoefh[E],
				                have(dcoef[k][jl][i], dcoef[k][jl][ir]));
				dcoefh[NE] = have(dcoefh[E],
				                have(dcoef[k][jr][i], dcoef[k][jr][ir]));
				dcoefh[FE] = have(dcoefh[E],
				                have(dcoef[kr][j][i], dcoef[kr][j][ir]));
				dcoefh[BE] = have(dcoefh[E],
				                have(dcoef[kl][j][i], dcoef[kl][j][ir]));
				dcoefh[FW] = have(dcoefh[W],
				                have(dcoef[kr][j][i], dcoef[kr][j][il]));
				dcoefh[BW] = have(dcoefh[W],
				                have(dcoef[kl][j][i], dcoef[kl][j][il]));
				dcoefh[NF] = have(dcoefh[F],
				                have(dcoef[k][jr][i], dcoef[kr][jr][i]));
				dcoefh[SF] = have(dcoefh[F],
				                have(dcoef[k][jl][i], dcoef[kr][jl][i]));
				dcoefh[NB] = have(dcoefh[B],
				                have(dcoef[k][jr][i], dcoef[kl][jr][i]));
				dcoefh[SB] = have(dcoefh[B],
				                have(dcoef[k][jl][i], dcoef[kl][jl][i]));
				row.i = i; row.j = j; row.k = k;
				col[0].i = i; col[0].j = j; col[0].k = k;
				cnt = 1;

				if (i != 0) {
					col[cnt].i = i-1; col[cnt].j = j; col[cnt].k = k;
					v[cnt] = 0; cnt++;
					if ((i-1) < patch->corners.is[0]) {
						double weight = coef[MET_W][k][j][i] * dcoefh[W] * dirichlet[k][j][i];
						acont[k][j][i-1] += weight;
					}
				}
				if (j != 0) {
					col[cnt].i = i; col[cnt].j = j-1; col[cnt].k = k;
					v[cnt] = 0; cnt++;
					if ((j-1) < patch->corners.is[1]) {
						double weight = coef[MET_S][k][j][i] * dcoefh[S] * dirichlet[k][j][i];
						acont[k][j-1][i] += weight;
					}
				}
				if (k != 0) {
					col[cnt].i = i; col[cnt].j = j; col[cnt].k = k-1;
					v[cnt] = 0; cnt++;
					if ((k-1) < patch->corners.is[2]) {
						double weight = coef[MET_B][k][j][i] * dcoefh[B] * dirichlet[k][j][i];
						acont[k-1][j][i] += weight;
					}
				}
				if (i != ngx - 1) {
					col[cnt].i = i+1; col[cnt].j = j; col[cnt].k = k;
					v[cnt] = 0; cnt++;
					if ((i+1) > patch->corners.ie[0]) {
						double weight = coef[MET_W][k][j][i+1] * dcoefh[E] * dirichlet[k][j][i];
						acont[k][j][i+1] += weight;
					}
				}
				if (j != ngy - 1) {
					col[cnt].i = i; col[cnt].j = j+1; col[cnt].k = k;
					v[cnt] = 0; cnt++;
					if ((j+1) > patch->corners.ie[1]) {
						double weight = coef[MET_S][k][j+1][i] * dcoefh[N] * dirichlet[k][j][i];
						acont[k][j+1][i] += weight;
					}
				}
				if (k != ngz - 1) {
					col[cnt].i = i; col[cnt].j = j; col[cnt].k = k+1;
					v[cnt] = 0; cnt++;
					if ((k+1) > patch->corners.ie[2]) {
						double weight = coef[MET_B][k+1][j][i] * dcoefh[F] * dirichlet[k][j][i];
						acont[k+1][j][i] += weight;
					}
				}
				if (i != 0 && j != 0) {
					col[cnt].i = i-1; col[cnt].j = j-1; col[cnt].k = k;
					v[cnt] = 0; cnt++;
					if ((i-1) < patch->corners.is[0] || (j-1) < patch->corners.is[1]) {
						double weight = coef[MET_SW][k][j][i] * dcoefh[SW] * dirichlet[k][j][i];
						acont[k][j-1][i-1] += weight;
					}
				}
				if (i != ngx - 1 && j != 0) {
					col[cnt].i = i+1; col[cnt].j = j-1; col[cnt].k = k;
					v[cnt] = 0; cnt++;
					if ((i+1) > patch->corners.ie[0] || (j-1) < patch->corners.is[1]) {
						double weight = coef[MET_SE][k][j][i] * dcoefh[SE] * dirichlet[k][j][i];
						acont[k][j-1][i+1] += weight;
					}
				}
				if (i != 0 && j != ngy-1) {
					col[cnt].i = i-1; col[cnt].j = j+1; col[cnt].k = k;
					v[cnt] = 0; cnt++;
					if ((i-1) < patch->corners.is[0] || (j+1) > patch->corners.ie[1]) {
						double weight = coef[MET_SE][k][j+1][i-1] * dcoefh[NW] * dirichlet[k][j][i];
						acont[k][j+1][i-1] += weight;
					}
				}
				if ((i != ngx -1) && (j != ngy - 1)) {
					col[cnt].i = i+1; col[cnt].j = j+1; col[cnt].k = k;
					v[cnt] = 0; cnt++;
					if ((i+1) > patch->corners.ie[0] || (j+1) > patch->corners.ie[1]) {
						double weight = coef[MET_SW][k][j+1][i+1] * dcoefh[NW] * dirichlet[k][j][i];
						acont[k][j+1][i+1] += weight;
					}
				}
				if (i != 0 && k != 0) {
					col[cnt].i = i-1; col[cnt].j = j; col[cnt].k = k-1;
					v[cnt] = 0; cnt++;
					if ((i-1) < patch->corners.is[0] || (k-1) < patch->corners.is[2]) {
						double weight = coef[MET_WB][k][j][i] * dcoefh[BW] * dirichlet[k][j][i];
						acont[k-1][j][i-1] += weight;
					}
				}
				if (i != ngx - 1 && k != 0) {
					col[cnt].i = i+1; col[cnt].j = j; col[cnt].k = k-1;
					v[cnt] = 0; cnt++;
					if ((i+1) > patch->corners.ie[0] || (k-1) < patch->corners.is[2]) {
						double weight = coef[MET_EB][k][j][i] * dcoefh[BE] * dirichlet[k][j][i];
						acont[k-1][j][i+1] += weight;
					}
				}
				if (i != 0 && k != ngz-1) {
					col[cnt].i = i-1; col[cnt].j = j; col[cnt].k = k+1;
					v[cnt] = 0; cnt++;
					if ((i-1) < patch->corners.is[0] || (k+1) > patch->corners.ie[2]) {
						double weight = coef[MET_EB][k+1][j][i-1]  * dcoefh[BE] * dirichlet[k][j][i];
						acont[k+1][j][i-1] += weight;
					}
				}
				if ((i != ngx-1) && (k != ngz -1)) {
					col[cnt].i = i+1; col[cnt].j = j; col[cnt].k = k+1;
					v[cnt] = 0; cnt++;
					if ((i+1) > patch->corners.ie[0] || (k+1) > patch->corners.ie[2]) {
						double weight = coef[MET_WB][k+1][j][i+1] * dcoefh[FE] * dirichlet[k][j][i];
						acont[k+1][j][i+1] += weight;
					}
				}
				if (j != 0 && k != 0) {
					col[cnt].i = i; col[cnt].j = j-1; col[cnt].k = k-1;
					v[cnt] = 0; cnt++;
					if ((j-1) < patch->corners.is[1] || (k-1) < patch->corners.is[2]) {
						double weight = coef[MET_SB][k][j][i] * dcoefh[SB] * dirichlet[k][j][i];
						acont[k-1][j-1][i] += weight;
					}
				}
				if (j != ngy - 1 && k != 0) {
					col[cnt].i = i; col[cnt].j = j+1; col[cnt].k = k-1;
					v[cnt] = 0; cnt++;
					if ((j+1) > patch->corners.ie[1] || (k-1) < patch->corners.is[2]) {
						double weight = coef[MET_NB][k][j][i] * dcoefh[NB] * dirichlet[k][j][i];
						acont[k-1][j+1][i] += weight;
					}
				}
				if (j != 0 && k != ngz-1) {
					col[cnt].i = i; col[cnt].j = j-1; col[cnt].k = k+1;
					v[cnt] = 0; cnt++;
					if ((j-1) < patch->corners.is[1] || (k+1) > patch->corners.ie[2]) {
						double weight = coef[MET_NB][k+1][j-1][i]  * dcoefh[SF] * dirichlet[k][j][i];
						acont[k+1][j-1][i] += weight;
					}
				}
				if ((j != ngy-1) && (k != ngz -1)) {
					col[cnt].i = i; col[cnt].j = j+1; col[cnt].k = k+1;
					v[cnt] = 0; cnt++;
					if ((j+1) > patch->corners.ie[1] || (k+1) > patch->corners.ie[2]) {
						double weight = coef[MET_SB][k+1][j+1][i] * dcoefh[NF] * dirichlet[k][j][i];
						acont[k+1][j+1][i] += weight;
					}
				}
				ierr = set_stencil(A, 1, &row, cnt, col, v, INSERT_VALUES);CHKERRQ(ierr);
			}
		}
	}

	for (i = 0; i < 10; i++) {
		ierr = DMDAVecRestoreArray(da, met->lcoef[i], &coef[i]);CHKERRQ(ierr);
	}

	ierr = stella_dmap_restore(dmap, &dirichlet);CHKERRQ(ierr);
	ierr = DMDAVecRestoreArray(da, bc->level->ladd_cont, &acont);CHKERRQ(ierr);
	ierr = DMDAVecRestoreArray(da, bc->level->ldcoef, &dcoef);CHKERRQ(ierr);

	stella_dmap_destroy(dmap); free(dmap);

	return 0;
}




static PetscErrorCode apply_dirichlet_rhs(stella_bc *bc, DM da, Vec rhs)
{
	PetscErrorCode ierr;
	PetscInt i,j;
	stella_patch *patch = bc->patch;
	PetscScalar **rhs_vec, **pm;
	PetscInt xs,ys,xe,ye,ym,xm;
	int stride[2];
	int offset[2];
	double **dirichlet;
	stella_dmap *dmap;

	ys = patch->corners.is[1]; xs = patch->corners.is[0];
	ye = patch->corners.ie[1]+1; xe = patch->corners.ie[0]+1;
	ym = patch->corners.ie[1] - patch->corners.is[1] + 1;
	xm = patch->corners.ie[0] - patch->corners.is[0] + 1;

	stride[0] = patch->stride[0];
	stride[1] = patch->stride[1];
  offset[0] = patch->offset[0];
  offset[1] = patch->offset[1];
	dmap = stella_dmap_create_2d(xs - offset[0], ys - offset[1], xm, ym, stride);
	ierr = stella_dmap_get(dmap, patch->dirichlet, &dirichlet);CHKERRQ(ierr);

	ierr = DMDAVecGetArray(da, rhs, &rhs_vec);CHKERRQ(ierr);
	ierr = DMDAVecGetArray(da, bc->level->pm, &pm);CHKERRQ(ierr);
	for (j = ys; j < ye; j++) {
		for (i = xs; i < xe; i++) {
			rhs_vec[j][i] = dirichlet[j][i];
			pm[j][i] = 0.0;
		}
	}
	ierr = DMDAVecRestoreArray(da, rhs, &rhs_vec);CHKERRQ(ierr);
	ierr = DMDAVecRestoreArray(da, bc->level->pm, &pm);CHKERRQ(ierr);
	ierr = stella_dmap_restore(dmap, &dirichlet);CHKERRQ(ierr);
	stella_dmap_destroy(dmap); free(dmap);

	return 0;
}


static PetscErrorCode apply_dirichlet_rhs_3d(stella_bc *bc, DM da, Vec rhs)
{
	PetscErrorCode ierr;
	PetscInt i,j,k;
	stella_patch *patch = bc->patch;
	PetscScalar ***rhs_vec, ***pm;
	PetscInt xs, ys, zs, xe, ye, ze, xm, ym, zm;
	int stride[3];
	int offset[3];
	double ***dirichlet;
	stella_dmap *dmap;

	ys = patch->corners.is[1]; xs = patch->corners.is[0]; zs = patch->corners.is[2];
	ye = patch->corners.ie[1]+1; xe = patch->corners.ie[0]+1; ze = patch->corners.ie[2]+1;
	ym = patch->corners.ie[1] - patch->corners.is[1] + 1;
	xm = patch->corners.ie[0] - patch->corners.is[0] + 1;
	zm = patch->corners.ie[2] - patch->corners.is[2] + 1;
	stride[0] = patch->stride[0]; stride[1] = patch->stride[1]; stride[2] = patch->stride[2];
	offset[0] = patch->offset[0]; offset[1] = patch->offset[1]; offset[2] = patch->offset[2];

	dmap = stella_dmap_create_3d(xs - offset[0], ys - offset[1], zs - offset[2], xm, ym, zm, stride);
	ierr = stella_dmap_get(dmap, patch->dirichlet, &dirichlet);CHKERRQ(ierr);

	ierr = DMDAVecGetArray(da, rhs, &rhs_vec);CHKERRQ(ierr);
	ierr = DMDAVecGetArray(da, bc->level->pm, &pm);CHKERRQ(ierr);
	for (k = zs; k < ze; k++) {
		for (j = ys; j < ye; j++) {
			for (i = xs; i < xe; i++) {
				rhs_vec[k][j][i] = dirichlet[k][j][i];
				pm[k][j][i] = 0.0;
			}
		}
	}
	ierr = DMDAVecRestoreArray(da, rhs, &rhs_vec);CHKERRQ(ierr);
	ierr = DMDAVecRestoreArray(da, bc->level->pm, &pm);CHKERRQ(ierr);

	ierr = stella_dmap_restore(dmap, &dirichlet);CHKERRQ(ierr);
	stella_dmap_destroy(dmap); free(dmap);

	return 0;
}


static PetscErrorCode apply_op(stella_bc *bc, Mat A, DM da)
{
	PetscErrorCode ierr;
	PetscInt dim;

	ierr = DMDAGetInfo(da, &dim, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL);CHKERRQ(ierr);

	if (dim == 2) {
		ierr = apply_dirichlet(bc, A, da);CHKERRQ(ierr);
	} else if (dim == 3) {
		ierr = apply_dirichlet_3d(bc, A, da);CHKERRQ(ierr);
	}

	return 0;
}


static PetscErrorCode symmetric(stella_bc *bc, Mat A, DM da)
{
	PetscErrorCode ierr;
	PetscInt dim;

	ierr = DMDAGetInfo(da, &dim, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL);CHKERRQ(ierr);

	if (dim == 2) {
		ierr = update_dirichlet_sym(bc, A, da);CHKERRQ(ierr);
	} else if (dim == 3) {
		ierr = update_dirichlet_sym_3d(bc, A, da);CHKERRQ(ierr);
	}

	return 0;
}


static PetscErrorCode apply_rhs(stella_bc *bc, DM da, Vec rhs)
{
	PetscErrorCode ierr;
	PetscInt dim;

	ierr = DMDAGetInfo(da, &dim, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL);CHKERRQ(ierr);

	if (dim == 2) {
		ierr = apply_dirichlet_rhs(bc, da, rhs);CHKERRQ(ierr);
	} else if (dim == 3) {
		ierr = apply_dirichlet_rhs_3d(bc, da, rhs);CHKERRQ(ierr);
	}

	return 0;
}


PetscErrorCode stella_dirichlet_create(stella_bc *bc)
{
	bc->apply = &apply_op;
	bc->apply_rhs = &apply_rhs;
	bc->symmetric = &symmetric;

	return 0;
}
