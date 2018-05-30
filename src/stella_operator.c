#include <stdlib.h>

#include "stella_operator.h"
#include "stella_mat.h"
#include "stella_io.h"


static PetscErrorCode centered_2d(stella_operator *op, Mat A, DM da)
{
	PetscErrorCode ierr;
	DMBoundaryType bx, by;
	DM cda;
	DMDACoor2d **coors;
	Vec lc;
	PetscInt i,j,k,xs,ys,xm,ym,ngx,ngy;
	PetscInt ncol;
	PetscScalar v[9];
	MatStencil row, col[9];
	double dcoefh[9];
	PetscScalar **dcoef;
	PetscScalar **bcoef;
	PetscScalar **coef[5];
	PetscScalar **jac[4];
    PetscScalar **rootg;
	double acc;
	MatType mtype;
	PetscBool is_boxmg;

	stella_metric *met;

	enum dir {
		O = 0, N=1, S=2, E=3, W=4, NW=5, NE=6, SW=7, SE=8
	};

    // alias 
	met = op->level->metric;

	ierr = stella_io_print(MPI_COMM_WORLD, "start assembling");CHKERRQ(ierr);
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

	ierr = DMDAGetInfo(da, 0, &ngx, &ngy, 0,0,0,0,0,0, &bx, &by,0,0);CHKERRQ(ierr);

	ierr = DMDAGetCorners(da, &xs, &ys, 0, &xm, &ym, 0);CHKERRQ(ierr);

	ierr = DMDAVecGetArray(da, op->level->ldcoef, &dcoef);CHKERRQ(ierr);
	ierr = DMDAVecGetArray(da, op->level->lbcoef, &bcoef);CHKERRQ(ierr);
	for (i = 0; i < 5; i++) {
		ierr = DMDAVecGetArray(da, met->lcoef[i], &coef[i]);CHKERRQ(ierr);
	}

    ierr = DMDAVecGetArray(da, met->jac, &rootg);CHKERRQ(ierr);

	if (op->axisymmetric) {
		ierr = DMGetCoordinateDM(da, &cda);CHKERRQ(ierr);
		ierr = DMGetCoordinatesLocal(da, &lc);CHKERRQ(ierr);
		ierr = DMDAVecGetArray(cda, lc, &coors);CHKERRQ(ierr);
	}

	PetscInt ibeg, iend, jbeg, jend;
	ibeg = xs; jbeg = ys; iend = xs + xm; jend = ys + ym;
	if (xs == 0 && bx != DM_BOUNDARY_PERIODIC)
		ibeg++;
	if (ys == 0 && by != DM_BOUNDARY_PERIODIC)
		jbeg++;
	if (xs + xm == ngx && bx != DM_BOUNDARY_PERIODIC)
		iend--;
	if (ys + ym == ngy && by != DM_BOUNDARY_PERIODIC)
		jend--;

	for (j = jbeg; j < jend; j++) {
		for (i = ibeg; i < iend; i++) {
			dcoefh[N] = have(dcoef[j][i], dcoef[j+1][i]);
			dcoefh[S] = have(dcoef[j][i], dcoef[j-1][i]);
			dcoefh[W] = have(dcoef[j][i], dcoef[j][i-1]);
			dcoefh[E] = have(dcoef[j][i], dcoef[j][i+1]);
			dcoefh[SW] = have(have(dcoef[j][i], dcoef[j][i-1]),
			                have(dcoef[j-1][i], dcoef[j-1][i-1]));
			dcoefh[NW] = have(have(dcoef[j][i], dcoef[j][i-1]),
			                have(dcoef[j+1][i], dcoef[j+1][i-1]));
			dcoefh[SE] = have(have(dcoef[j][i], dcoef[j][i+1]),
			                have(dcoef[j-1][i], dcoef[j-1][i+1]));
			dcoefh[NE] = have(have(dcoef[j][i], dcoef[j][i+1]),
			                have(dcoef[j+1][i], dcoef[j+1][i+1]));

			if (op->axisymmetric) {
				dcoefh[E] = dcoefh[E] * (coors[j][i+1].x + coors[j][i].x) * .5;
				dcoefh[W] = dcoefh[W] * (coors[j][i].x + coors[j][i-1].x) * .5;
				dcoefh[N] = dcoefh[N] * coors[j][i].x;
				dcoefh[S] = dcoefh[S] * coors[j][i].x;
			}

			row.i = i; row.j = j;

			col[SE].i = i+1; col[SE].j = j-1;
			col[S].i = i; col[S].j = j-1;
			col[SW].i = i-1; col[SW].j = j-1;
			col[W].i = i-1; col[W].j = j;
			col[O].i = i; col[O].j = j;
			col[E].i = i+1; col[E].j = j;
			col[NE].i = i+1; col[NE].j = j+1;
			col[N].i = i; col[N].j = j+1;
			col[NW].i = i-1; col[NW].j = j+1;

			ncol = 9;

			v[S] = -1*coef[MET_S][j][i]*dcoefh[S];
			v[N] = -1*coef[MET_S][j+1][i]*dcoefh[N];
			v[W] = -1*coef[MET_W][j][i]*dcoefh[W];
			v[E] = -1*coef[MET_W][j][i+1]*dcoefh[E];
			v[SW] = -1*coef[MET_SW][j][i]*dcoefh[SW];
			v[SE] = -1*coef[MET_SE][j][i]*dcoefh[SE];
			v[NW] = -1*coef[MET_SE][j+1][i-1]*dcoefh[NW];
			v[NE] = -1*coef[MET_SW][j+1][i+1]*dcoefh[NE];
			v[O] = -1 * (v[S] + v[N] + v[E] + v[W]);

            // v[O] = v[O] + (PetscScalar)(100);
            // we are solving
            // -laplacian (phi) = f ;

			PetscScalar diag_cont;
			diag_cont = rootg[j][i]*bcoef[j][i];
			if (op->axisymmetric) {
			  diag_cont *= coors[j][i].x;
			}
			v[O] += PetscSign(v[O]) * PetscAbsScalar(diag_cont);

			ierr = PetscStrcmp(mtype, MATSHELL, &is_boxmg);CHKERRQ(ierr);
			ierr = set_stencil(A, 1, &row, ncol, col, v, INSERT_VALUES);CHKERRQ(ierr);
		}
	}

	ierr = DMDAVecRestoreArray(da, op->level->ldcoef, &dcoef);CHKERRQ(ierr);
	ierr = DMDAVecRestoreArray(da, op->level->lbcoef, &bcoef);CHKERRQ(ierr);
	for (i = 0; i < 5; i++) {
		ierr = DMDAVecRestoreArray(da, met->lcoef[i], &coef[i]);CHKERRQ(ierr);
	}
        ierr = DMDAVecRestoreArray(da, met->jac, &rootg);CHKERRQ(ierr);

	if (op->axisymmetric) {
		ierr = DMDAVecRestoreArray(cda, lc, &coors);CHKERRQ(ierr);
	}

	return 0;
}


static PetscErrorCode centered_3d(stella_operator *op, Mat A, DM da)
{
	PetscErrorCode ierr;
	DMBoundaryType bx, by, bz;
	PetscInt i,j,k,l,xs,ys,zs,xm,ym,zm,ngx,ngy,ngz;
	PetscScalar v[19];
	PetscInt level;
	PetscInt ncol;
	MatStencil row, col[19];
	double dcoefh[19];
	Vec ldcoef;
	PetscScalar ***dcoef;
	PetscScalar ***bcoef;
	PetscScalar ***rootg;
	PetscScalar ***coef[19];
	PetscScalar ***jac[9];
	MatType mtype;
	PetscBool is_boxmg;

	stella_metric *met;


	met = op->level->metric;

	ierr = stella_io_print(MPI_COMM_WORLD, "start assembling");CHKERRQ(ierr);
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


	ierr = DMDAGetInfo(da, 0, &ngx, &ngy, &ngz,0,0,0,0,0, &bx, &by, &bz,0);CHKERRQ(ierr);

	ierr = DMDAGetCorners(da, &xs, &ys, &zs, &xm, &ym, &zm);CHKERRQ(ierr);

	ierr = DMDAVecGetArray(da, op->level->ldcoef, &dcoef);CHKERRQ(ierr);
	ierr = DMDAVecGetArray(da, op->level->lbcoef, &bcoef);CHKERRQ(ierr);
	for (i = 0; i < 10; i++) {
		ierr = DMDAVecGetArray(da, met->lcoef[i], &coef[i]);CHKERRQ(ierr);
	}

	ierr = DMDAVecGetArray(da, met->jac, &rootg);CHKERRQ(ierr);

	enum dir {
		O = 0, N, S, E, W, F, B, NW, NE, SW, SE, NF, NB, SF, SB, FE, FW, BE, BW
	};

	PetscInt ibeg, iend, jbeg, jend, kbeg, kend;
	ibeg = xs; jbeg = ys; kbeg = zs; iend = xs + xm; jend = ys + ym; kend = zs + zm;
	if (xs == 0 && bx != DM_BOUNDARY_PERIODIC)
		ibeg++;
	if (ys == 0 && by != DM_BOUNDARY_PERIODIC)
		jbeg++;
	if (zs == 0 && bz != DM_BOUNDARY_PERIODIC)
		kbeg++;
	if (xs + xm == ngx && bx != DM_BOUNDARY_PERIODIC)
		iend--;
	if (ys + ym == ngy && by != DM_BOUNDARY_PERIODIC)
		jend--;
	if (zs + zm == ngz && bz != DM_BOUNDARY_PERIODIC)
		kend--;

	ncol = 19;

	for (k = kbeg; k < kend; k++) {
		for (j = jbeg; j < jend; j++) {
			for (i = ibeg; i < iend; i++) {
					dcoefh[N] = have(dcoef[k][j][i], dcoef[k][j+1][i]);
					dcoefh[S] = have(dcoef[k][j][i], dcoef[k][j-1][i]);
					dcoefh[W] = have(dcoef[k][j][i], dcoef[k][j][i-1]);
					dcoefh[E] = have(dcoef[k][j][i], dcoef[k][j][i+1]);
					dcoefh[F] = have(dcoef[k][j][i], dcoef[k+1][j][i]);
					dcoefh[B] = have(dcoef[k][j][i], dcoef[k-1][j][i]);
					dcoefh[SW] = have(dcoefh[W],
					                have(dcoef[k][j-1][i], dcoef[k][j-1][i-1]));
					dcoefh[NW] = have(dcoefh[W],
					                have(dcoef[k][j+1][i], dcoef[k][j+1][i-1]));
					dcoefh[SE] = have(dcoefh[E],
					                have(dcoef[k][j-1][i], dcoef[k][j-1][i+1]));
					dcoefh[NE] = have(dcoefh[E],
					                have(dcoef[k][j+1][i], dcoef[k][j+1][i+1]));
					dcoefh[FE] = have(dcoefh[E],
					                have(dcoef[k+1][j][i], dcoef[k+1][j][i+1]));
					dcoefh[BE] = have(dcoefh[E],
					                have(dcoef[k-1][j][i], dcoef[k-1][j][i+1]));
					dcoefh[FW] = have(dcoefh[W],
					                have(dcoef[k+1][j][i], dcoef[k+1][j][i-1]));
					dcoefh[BW] = have(dcoefh[W],
					                have(dcoef[k-1][j][i], dcoef[k-1][j][i-1]));
					dcoefh[NF] = have(dcoefh[F],
					                have(dcoef[k][j+1][i], dcoef[k+1][j+1][i]));
					dcoefh[SF] = have(dcoefh[F],
					                have(dcoef[k][j-1][i], dcoef[k+1][j-1][i]));
					dcoefh[NB] = have(dcoefh[B],
					                have(dcoef[k][j+1][i], dcoef[k-1][j+1][i]));
					dcoefh[SB] = have(dcoefh[B],
					                have(dcoef[k][j-1][i], dcoef[k-1][j-1][i]));

					row.i = i; row.j = j; row.k = k;
					for (l = 0; l < 19; l++) {
						col[l].i = i;
						col[l].j = j;
						col[l].k = k;
					}
					col[N].j++;
					col[S].j--;
					col[E].i++;
					col[W].i--;
					col[F].k++;
					col[B].k--;
					col[NW].j++; col[NW].i--;
					col[NE].j++; col[NE].i++;
					col[SW].j--; col[SW].i--;
					col[SE].j--; col[SE].i++;
					col[NF].j++; col[NF].k++;
					col[NB].j++; col[NB].k--;
					col[SF].j--; col[SF].k++;
					col[SB].j--; col[SB].k--;
					col[FE].k++; col[FE].i++;
					col[FW].k++; col[FW].i--;
					col[BE].k--; col[BE].i++;
					col[BW].k--; col[BW].i--;

					v[S] = -1*coef[MET_S][k][j][i]*dcoefh[S];
					v[N] = -1*coef[MET_S][k][j+1][i]*dcoefh[N];
					v[W] = -1*coef[MET_W][k][j][i]*dcoefh[W];
					v[E] = -1*coef[MET_W][k][j][i+1]*dcoefh[E];
					v[SW] = -1*coef[MET_SW][k][j][i]*dcoefh[SW];
					v[SE] = -1*coef[MET_SE][k][j][i]*dcoefh[SE];
					v[NW] = -1*coef[MET_SE][k][j+1][i-1]*dcoefh[NW];
					v[NE] = -1*coef[MET_SW][k][j+1][i+1]*dcoefh[NE];
					v[F] = -1*coef[MET_B][k+1][j][i]*dcoefh[F];
					v[B] = -1*coef[MET_B][k][j][i]*dcoefh[B];
					v[NF] = -1*coef[MET_SB][k+1][j+1][i]*dcoefh[NF];
					v[NB] = -1*coef[MET_NB][k][j][i]*dcoefh[NB];
					v[SF] = -1*coef[MET_NB][k+1][j-1][i]*dcoefh[SF];
					v[SB] = -1*coef[MET_SB][k][j][i]*dcoefh[SB];
					v[FE] = -1*coef[MET_WB][k+1][j][i+1]*dcoefh[FE];
					v[FW] = -1*coef[MET_EB][k+1][j][i-1]*dcoefh[FW];
					v[BE] = -1*coef[MET_EB][k][j][i]*dcoefh[BE];
					v[BW] = -1*coef[MET_WB][k][j][i]*dcoefh[BW];

					v[O] = 0.0;
					for (l = 1; l < ncol; l++) {
						v[O] += v[l];
					}

					v[O] = -1*v[O];

					PetscScalar diag_cont;
					diag_cont = rootg[k][j][i]*bcoef[k][j][i];
					v[O] += PetscSign(v[O]) * PetscAbsScalar(diag_cont);

					/* printf("(%d %d %d) %g %g %g %g %g %g %g %g %g %g %g %g\n", i,j,k, v[SW], v[SE], v[NW], v[NE], v[NF], v[NB], v[SF], v[SB],v[FE], v[FW], v[BE], v[BW]); */
					/* printf("(%d %d %d) %g\n", i,j,k,v[SB]); */

					/* printf("(%d %d %d) [%f] %f %f %f %f %f %f\n", i,j,k, v[O], v[S], v[N], v[W], v[E], v[F], v[B]); */
					/* printf("(%d %d %d) %f %f %f %f %f %f %f %f %f %f %f %f\n", i, j, k, */
					/*        v[SW], v[SE], v[NW], v[NE], v[NF], v[NB], v[SF], */
					/*        v[SB], v[FE], v[FW], v[BE], v[BW]); */

					ierr = set_stencil(A, 1, &row, ncol, col, v, INSERT_VALUES);CHKERRQ(ierr);
			}
		}
	}

	ierr = DMDAVecRestoreArray(da, op->level->ldcoef, &dcoef);CHKERRQ(ierr);
	ierr = DMDAVecRestoreArray(da, op->level->lbcoef, &bcoef);CHKERRQ(ierr);
	for (i = 0; i < 10; i++) {
		ierr = DMDAVecRestoreArray(da, met->lcoef[i], &coef[i]);CHKERRQ(ierr);
	}
	ierr = DMDAVecRestoreArray(da, met->jac, &rootg);CHKERRQ(ierr);

	return 0;
}


PetscErrorCode stella_operator_create(stella_operator **efop, stella_level *level, stella_fd *fd, int nd)
{
	stella_operator *op = (stella_operator*) malloc(sizeof(stella_operator));

	op->axisymmetric = 0;

	if (nd == 3)
		op->assemble = &centered_3d;
	else{
		op->assemble = &centered_2d;
}
	op->level = level;
	op->fd = fd;

	*efop = op;
	return 0;
}


PetscErrorCode stella_operator_assemble(stella_operator *op, Mat A, DM da)
{
	PetscErrorCode ierr;
	ierr = op->assemble(op, A, da);CHKERRQ(ierr);
	return 0;
}


PetscErrorCode stella_operator_destroy(stella_operator *efop)
{
	return 0;
}
