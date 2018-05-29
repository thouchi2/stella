#include <stdlib.h>
#include <math.h>

#include "stella_gen.h"
#include "stella_metric.h"


static PetscErrorCode compute_3d_interior(stella_metric *met, DM da, stella_fd *fd)
{
	PetscErrorCode ierr;
	DM cda;
	DMBoundaryType bx, by, bz;
        int per_x, per_y, per_z;
	DMDACoor3d ***coors;
	PetscInt i,j,k,xs,ys,zs,xm,ym,zm,ngx,ngy,ngz;
	PetscInt ibeg, jbeg, kbeg, iend, jend, kend;
	PetscScalar ***coef[10], ***jac;
	PetscScalar ***jac_v[9];
	Vec cv;
	double *si, *sp, *sm;
	double feps = 1e-13;


	ierr = DMGetCoordinateDM(da, &cda);CHKERRQ(ierr);
	ierr = DMGetCoordinatesLocal(da, &cv);CHKERRQ(ierr);
	ierr = DMDAVecGetArray(cda, cv, &coors);CHKERRQ(ierr);

	ierr = DMDAGetCorners(cda, &xs, &ys, &zs, &xm, &ym, &zm);CHKERRQ(ierr);
	ierr = DMDAGetInfo(da, 0, &ngx, &ngy, &ngz,0,0,0,0,0, &bx, &by,&bz,0);CHKERRQ(ierr);
        per_x = bx==DM_BOUNDARY_PERIODIC;
        per_y = by==DM_BOUNDARY_PERIODIC;
        per_z = bz==DM_BOUNDARY_PERIODIC;

	ierr = DMDAVecGetArray(da, met->jac, &jac);CHKERRQ(ierr);

        ibeg = xs; iend = xs + xm; jbeg = ys; jend = ys + ym; kbeg = zs; kend = zs + zm;
        if(!per_x) {
          if (xs == 0) ibeg++;
          if (xs + xm == ngx) iend--;
        }
        if(!per_y) {
          if (ys == 0) jbeg++;
          if (ys + ym == ngy) jend--;
        }
        if(!per_z) {
          if (zs == 0) kbeg++;
          if (zs + zm == ngz) kend--;
        }

	for (i = 0; i < 10; i++) {
		ierr = DMDAVecGetArray(da, met->lcoef[i], &coef[i]);CHKERRQ(ierr);
	}

	for (i = 0; i < 9; i++) {
		ierr = DMDAVecGetArray(da, met->jac_v[i], &jac_v[i]);CHKERRQ(ierr);
	}

	si = fd->first[0]->v;
	sp = fd->first[1]->v;
	sm = fd->first[-1]->v;

	INTERIOR_METRICS_3D;

	for (i = 0; i < 10; i++) {
		ierr = DMDAVecRestoreArray(da, met->lcoef[i], &coef[i]);CHKERRQ(ierr);
	}

	for (i = 0; i < 9; i++) {
		ierr = DMDAVecRestoreArray(da, met->jac_v[i], &jac_v[i]);CHKERRQ(ierr);
	}

	ierr = DMDAVecRestoreArray(da, met->jac, &jac);CHKERRQ(ierr);

	ierr = DMDAVecRestoreArray(cda, cv, &coors);CHKERRQ(ierr);

	return 0;
}


static PetscErrorCode compute_2d_interior(stella_metric *met, DM da, stella_fd *fd)
{
	PetscErrorCode ierr;
	DM cda;
	DMBoundaryType bx, by;
        int per_x, per_y;
	DMDACoor2d **coors;
	PetscInt i,j,xs,ys,xm,ym,ngx,ngy;
	int k;
	PetscScalar **coef[5], **jac;
	PetscScalar **jac_v[4];
	PetscScalar **vec__x_s, **vec__x_t, **vec__y_s, **vec__y_t;
	Vec cv;
	double *si, *sp, *sm;
	double feps = 1e-13;

	ierr = DMGetCoordinateDM(da, &cda);CHKERRQ(ierr);
	ierr = DMGetCoordinatesLocal(da, &cv);CHKERRQ(ierr);
	ierr = DMDAVecGetArray(cda, cv, &coors);CHKERRQ(ierr);

	ierr = DMDAGetCorners(cda, &xs, &ys, 0, &xm, &ym, 0);CHKERRQ(ierr);
	ierr = DMDAGetInfo(da, 0, &ngx, &ngy, 0,0,0,0,0,0, &bx, &by,0,0);CHKERRQ(ierr);
        per_x = bx==DM_BOUNDARY_PERIODIC;
        per_y = by==DM_BOUNDARY_PERIODIC;

	for (i = 0; i < 5; i++) {
		ierr = DMDAVecGetArray(da, met->lcoef[i], &coef[i]);CHKERRQ(ierr);
	}

	for (i = 0; i < 4; i++) {
		ierr = DMDAVecGetArray(da, met->jac_v[i], &jac_v[i]);CHKERRQ(ierr);
	}

	vec__x_s = jac_v[met->t2map[0][0]];
	vec__x_t = jac_v[met->t2map[0][1]];
	vec__y_s = jac_v[met->t2map[1][0]];
	vec__y_t = jac_v[met->t2map[1][1]];

	ierr = DMDAVecGetArray(da, met->jac, &jac);CHKERRQ(ierr);

	si = fd->first[0]->v;
	sp = fd->first[1]->v;
	sm = fd->first[-1]->v;

	enum {c=0,x_h=1,y_h=2};
	double x_s[3], x_t[3], y_s[3], y_t[3], J[3];
	DMDACoor2d coors_xh_data[3], coors_yh_data[3];
	DMDACoor2d *coors_xh = &coors_xh_data[1];
	DMDACoor2d *coors_yh = &coors_yh_data[1];
	double g_22, g_11, g_12, g11, g22, g12;

	PetscInt ibeg, jbeg, iend, jend;

	ibeg = xs; iend = xs + xm; jbeg = ys; jend = ys + ym;
        if(!per_x) {
          if (xs == 0) ibeg++;
          if (xs + xm == ngx) iend--;
        }
        if(!per_y) {
          if (ys == 0) jbeg++;
          if (ys + ym == ngy) jend--;
        }

	for (j = jbeg; j < jend; j++) {
		for (i = ibeg; i < iend; i++) {
			x_s[c] = si[1]*coors[j][i+1].x + si[-1]*coors[j][i-1].x;
			y_s[c] = si[1]*coors[j][i+1].y + si[-1]*coors[j][i-1].y;
			x_t[c] = si[1]*coors[j+1][i].x + si[-1]*coors[j-1][i].x;
			y_t[c] = si[1]*coors[j+1][i].y + si[-1]*coors[j-1][i].y;

			J[c] = x_s[c]*y_t[c] - x_t[c]*y_s[c];
			jac[j][i] = J[c];
			vec__x_s[j][i] = x_s[c];
			vec__y_s[j][i] = y_s[c];
			vec__x_t[j][i] = x_t[c];
			vec__y_t[j][i] = y_t[c];

			g_12 = x_s[c]*x_t[c] + y_s[c]*y_t[c];
			if (fabs(J[c]) < feps) g12 = 0.0;
			else g12 = -1./(J[c]*J[c]) * g_12;

			coef[MET_SW][j][i+1] += .25 * J[c] * g12;
			coef[MET_SW][j+1][i] += .25 * J[c] * g12;
			coef[MET_SE][j+1][i] += -.25 * J[c] * g12;
			coef[MET_SE][j][i-1] += -.25 * J[c] * g12;
		}
	}


	for (j = jbeg; j < ys + ym; j++) {
		for (i = ibeg; i < iend; i++) {
			for (k = -1; k < 2; k = k + 2) {
				coors_yh[k].x = .5 * (coors[j-1][i+k].x + coors[j][i+k].x);
				coors_yh[k].y = .5 * (coors[j-1][i+k].y + coors[j][i+k].y);
			}

			x_s[y_h] = si[-1]*coors_yh[-1].x + si[1]*coors_yh[1].x;
			y_s[y_h] = si[-1]*coors_yh[-1].y + si[1]*coors_yh[1].y;
			x_t[y_h] = coors[j][i].x - coors[j-1][i].x;
			y_t[y_h] = coors[j][i].y - coors[j-1][i].y;

			J[y_h] = x_s[y_h]*y_t[y_h] - x_t[y_h]*y_s[y_h];

			g_11 = x_s[y_h]*x_s[y_h] + y_s[y_h]*y_s[y_h];
			if (fabs(J[y_h]) < feps) g22 = 0.0;
			else g22 = 1./(J[y_h]*J[y_h]) * g_11;
			coef[MET_S][j][i] = J[y_h] * g22;
		}
	}


	for (j = jbeg; j < jend; j++) {
		for (i = ibeg; i < xs + xm; i++) {
			for (k = -1; k < 2; k = k + 2) {
				coors_xh[k].x = .5 * (coors[j+k][i-1].x + coors[j+k][i].x);
				coors_xh[k].y = .5 * (coors[j+k][i-1].y + coors[j+k][i].y);
			}

			x_s[x_h] = coors[j][i].x - coors[j][i-1].x;
			y_s[x_h] = coors[j][i].y - coors[j][i-1].y;
			x_t[x_h] = si[-1]*coors_xh[-1].x + si[1]*coors_xh[1].x;
			y_t[x_h] = si[-1]*coors_xh[-1].y + si[1]*coors_xh[1].y;

			J[x_h] = x_s[x_h]*y_t[x_h] - x_t[x_h]*y_s[x_h];

			g_22 = x_t[x_h]*x_t[x_h] + y_t[x_h]*y_t[x_h];
			if (fabs(J[x_h]) < feps) g11 = 0.0;
			else g11 = 1./(J[x_h]*J[x_h]) * g_22;
			coef[MET_W][j][i] = J[x_h] * g11;
		}
	}


	for (i = 0; i < 5; i++) {
		ierr = DMDAVecRestoreArray(da, met->lcoef[i], &coef[i]);CHKERRQ(ierr);
	}

	for (i = 0; i < 4; i++) {
		ierr = DMDAVecRestoreArray(da, met->jac_v[i], &jac_v[i]);CHKERRQ(ierr);
	}

	ierr = DMDAVecRestoreArray(da, met->jac, &jac);CHKERRQ(ierr);

	ierr = DMDAVecRestoreArray(cda, cv, &coors);CHKERRQ(ierr);

	return 0;
}


static PetscErrorCode compute_3d_boundary(stella_metric *met, DM da, stella_fd *fd)
{
	PetscErrorCode ierr;
	DM cda;
	DMBoundaryType bx, by, bz;
	DMDACoor3d ***coors;
	PetscInt i,j,k,xs,ys,zs,xm,ym,zm,ngx,ngy,ngz;
	PetscInt ibeg, jbeg, kbeg, iend, jend, kend;
	PetscScalar ***coef[10], ***jac;
	PetscScalar ***jac_v[9];
	Vec cv;
	double *si, *sp, *sm;
	double feps = 1e-13;
	int per_x, per_y, per_z;


	ierr = DMGetCoordinateDM(da, &cda);CHKERRQ(ierr);
	ierr = DMGetCoordinatesLocal(da, &cv);CHKERRQ(ierr);
	ierr = DMDAVecGetArray(cda, cv, &coors);CHKERRQ(ierr);

	ierr = DMDAGetCorners(cda, &xs, &ys, &zs, &xm, &ym, &zm);CHKERRQ(ierr);
	ierr = DMDAGetInfo(da, 0, &ngx, &ngy, &ngz,0,0,0,0,0, &bx, &by,&bz,0);CHKERRQ(ierr);
        per_x = bx==DM_BOUNDARY_PERIODIC;
        per_y = by==DM_BOUNDARY_PERIODIC;
        per_z = bz==DM_BOUNDARY_PERIODIC;

	ierr = DMDAVecGetArray(da, met->jac, &jac);CHKERRQ(ierr);
	for (i = 0; i < 10; i++) {
		ierr = DMDAVecGetArray(da, met->lcoef[i], &coef[i]);CHKERRQ(ierr);
	}

	for (i = 0; i < 9; i++) {
		ierr = DMDAVecGetArray(da, met->jac_v[i], &jac_v[i]);CHKERRQ(ierr);
	}

	si = fd->first[0]->v;
	sp = fd->first[1]->v;
	sm = fd->first[-1]->v;

	BOUNDARY_METRICS_3D;

	for (i = 0; i < 10; i++) {
		ierr = DMDAVecRestoreArray(da, met->lcoef[i], &coef[i]);CHKERRQ(ierr);
		ierr = VecSet(met->coef[i], 0.0);CHKERRQ(ierr);
		ierr = DMLocalToGlobalBegin(da, met->lcoef[i], ADD_VALUES, met->coef[i]);CHKERRQ(ierr);
		ierr = DMLocalToGlobalEnd(da, met->lcoef[i], ADD_VALUES, met->coef[i]);CHKERRQ(ierr);
		ierr = DMGlobalToLocalBegin(da, met->coef[i], INSERT_VALUES, met->lcoef[i]);CHKERRQ(ierr);
		ierr = DMGlobalToLocalEnd(da, met->coef[i], INSERT_VALUES, met->lcoef[i]);CHKERRQ(ierr);
	}

	for (i = 0; i < 9; i++) {
		ierr = DMDAVecRestoreArray(da, met->jac_v[i], &jac_v[i]);CHKERRQ(ierr);
	}

	ierr = DMDAVecRestoreArray(da, met->jac, &jac);CHKERRQ(ierr);

	ierr = DMDAVecRestoreArray(cda, cv, &coors);CHKERRQ(ierr);

	return 0;
}


static PetscErrorCode compute_2d_boundary(stella_metric *met, DM da, stella_fd *fd)
{
	PetscErrorCode ierr;
	DM cda;
	DMBoundaryType bx, by;
	DMDACoor2d **coors;
	PetscInt i,j,xs,ys,xm,ym,ngx,ngy;
	int k;
	PetscScalar **coef[5], **jac;
	PetscScalar **jac_v[4];
	PetscScalar **vec__x_s, **vec__x_t, **vec__y_s, **vec__y_t;
	Vec cv;
	double *si, *sp, *sm;
	double feps = 1e-13;
	int per_x, per_y;

	ierr = DMGetCoordinateDM(da, &cda);CHKERRQ(ierr);
	ierr = DMGetCoordinatesLocal(da, &cv);CHKERRQ(ierr);
	ierr = DMDAVecGetArray(cda, cv, &coors);CHKERRQ(ierr);

	ierr = DMDAGetCorners(cda, &xs, &ys, 0, &xm, &ym, 0);CHKERRQ(ierr);
	ierr = DMDAGetInfo(da, 0, &ngx, &ngy, 0,0,0,0,0,0, &bx, &by,0,0);CHKERRQ(ierr);

	per_x = bx==DM_BOUNDARY_PERIODIC;
	per_y = by==DM_BOUNDARY_PERIODIC;

	for (i = 0; i < 5; i++) {
		ierr = DMDAVecGetArray(da, met->lcoef[i], &coef[i]);CHKERRQ(ierr);
	}

	for (i = 0; i < 4; i++) {
		ierr = DMDAVecGetArray(da, met->jac_v[i], &jac_v[i]);CHKERRQ(ierr);
	}

	vec__x_s = jac_v[met->t2map[0][0]];
	vec__x_t = jac_v[met->t2map[0][1]];
	vec__y_s = jac_v[met->t2map[1][0]];
	vec__y_t = jac_v[met->t2map[1][1]];

	ierr = DMDAVecGetArray(da, met->jac, &jac);CHKERRQ(ierr);

	si = fd->first[0]->v;
	sp = fd->first[1]->v;
	sm = fd->first[-1]->v;

	enum {c=0,x_h=1,y_h=2};
	double x_s[3], x_t[3], y_s[3], y_t[3], J[3];
	DMDACoor2d coors_xh_data[5], coors_yh_data[5];
	DMDACoor2d *coors_xh = &coors_xh_data[1];
	DMDACoor2d *coors_yh = &coors_yh_data[1];
	double g_22, g_11, g_12, g12, g11, g22;

	PetscInt ibeg, jbeg, iend, jend;

	ibeg = xs; iend = xs + xm; jbeg = ys; jend = ys + ym;
	if (xs == 0) ibeg++;
	if (xs + xm == ngx) iend--;
	if (ys == 0) jbeg++;
	if (ys + ym == ngy) jend--;


	// left boundary
	if (xs == 0 && !per_x) {
		for (j = jbeg; j < jend; j++) {
			x_s[c] = sm[0]*coors[j][0].x + sm[1]*coors[j][1].x + sm[2]*coors[j][2].x;
			y_s[c] = sm[0]*coors[j][0].y + sm[1]*coors[j][1].y + sm[2]*coors[j][2].y;
			x_t[c] = si[-1]*coors[j-1][0].x + si[1]*coors[j+1][0].x;
			y_t[c] = si[-1]*coors[j-1][0].y + si[1]*coors[j+1][0].y;

			J[c] = x_s[c]*y_t[c] - x_t[c]*y_s[c];
			jac[j][0] = J[c];
			vec__x_s[j][0] = x_s[c];
			vec__y_s[j][0] = y_s[c];
			vec__x_t[j][0] = x_t[c];
			vec__y_t[j][0] = y_t[c];

			g_12 = x_s[c]*x_t[c] + y_s[c]*y_t[c];
			if (fabs(J[c] < feps)) g12 = 0.0;
			else g12 = -1./(J[c]*J[c]) * g_12;

			coef[MET_SW][j][1] += .25 * J[c] * g12;
			coef[MET_SW][j+1][0] += .25 * J[c] * g12;
			coef[MET_SE][j+1][0] += -.25 * J[c] * g12;
		}

		for (j = jbeg; j < ys + ym; j++) {
			for (k = 0; k < 3; k++) {
				coors_yh[k].x = .5 * (coors[j-1][k].x + coors[j][k].x);
				coors_yh[k].y = .5 * (coors[j-1][k].y + coors[j][k].y);
			}

			x_s[y_h] = sm[0]*coors_yh[0].x + sm[1]*coors_yh[1].x + sm[2]*coors_yh[2].x;
			y_s[y_h] = sm[0]*coors_yh[0].y + sm[1]*coors_yh[1].y + sm[2]*coors_yh[2].y;
			x_t[y_h] = coors[j][0].x - coors[j-1][0].x;
			y_t[y_h] = coors[j][0].y - coors[j-1][0].y;

			J[y_h] = x_s[y_h]*y_t[y_h] - x_t[y_h]*y_s[y_h];

			g_11 = x_s[y_h]*x_s[y_h] + y_s[y_h]*y_s[y_h];
			if (fabs(J[y_h]) < feps) g22 = 0.0;
			else g22 = 1./(J[y_h]*J[y_h]) * g_11;
			coef[MET_S][j][0] = J[y_h] * g22;
		}
	}

	// right boundary
	if (xs + xm == ngx) {
		int idx = ngx - 1;
		for (j = jbeg; j < jend; j++) {
			x_s[c] = sp[0]*coors[j][idx].x + sp[-1]*coors[j][idx-1].x + sp[-2]*coors[j][idx-2].x;
			y_s[c] = sp[0]*coors[j][idx].y + sp[-1]*coors[j][idx-1].y + sp[-2]*coors[j][idx-2].y;
			x_t[c] = si[-1]*coors[j-1][idx].x + si[1]*coors[j+1][idx].x;
			y_t[c] = si[-1]*coors[j-1][idx].y + si[1]*coors[j+1][idx].y;

			J[c] = x_s[c]*y_t[c] - x_t[c]*y_s[c];
			jac[j][idx] = J[c];
			vec__x_s[j][idx] = x_s[c];
			vec__y_s[j][idx] = y_s[c];
			vec__x_t[j][idx] = x_t[c];
			vec__y_t[j][idx] = y_t[c];

			g_12 = x_s[c]*x_t[c] + y_s[c]*y_t[c];
			if (fabs(J[c]) < feps) g12 = 0.0;
			else g12 = -1./(J[c]*J[c]) * g_12;

			coef[MET_SW][j+1][idx] += .25 * J[c] * g12;
			coef[MET_SE][j+1][idx] += -.25 * J[c] * g12;
			coef[MET_SE][j][idx-1] += -.25 * J[c] * g12;
		}

		for (j = jbeg; j < ys + ym; j++) {
			for (k = 0; k < 3; k++) {
				coors_yh[k].x = .5 * (coors[j-1][idx-k].x + coors[j][idx-k].x);
				coors_yh[k].y = .5 * (coors[j-1][idx-k].y + coors[j][idx-k].y);
			}

			x_s[y_h] = sp[0]*coors_yh[0].x + sp[-1]*coors_yh[1].x + sp[-2]*coors_yh[2].x;
			y_s[y_h] = sp[0]*coors_yh[0].y + sp[-1]*coors_yh[1].y + sp[-2]*coors_yh[2].y;
			x_t[y_h] = coors[j][idx].x - coors[j-1][idx].x;
			y_t[y_h] = coors[j][idx].y - coors[j-1][idx].y;

			J[y_h] = x_s[y_h]*y_t[y_h] - x_t[y_h]*y_s[y_h];

			g_11 = x_s[y_h]*x_s[y_h] + y_s[y_h]*y_s[y_h];
			if (fabs(J[y_h]) < feps) g22 = 0.0;
			else g22 = 1./(J[y_h]*J[y_h]) * g_11;
			coef[MET_S][j][idx] = J[y_h] * g22;
		}
	}

	// bottom boundary
	if (ys == 0 && !per_y) {
		for (i = ibeg; i < iend; i++) {
			x_s[c] = si[-1]*coors[0][i-1].x + si[1]*coors[0][i+1].x;
			y_s[c] = si[-1]*coors[0][i-1].y + si[1]*coors[0][i+1].y;
			x_t[c] = sm[0]*coors[0][i].x + sm[1]*coors[1][i].x + sm[2]*coors[2][i].x;
			y_t[c] = sm[0]*coors[0][i].y + sm[1]*coors[1][i].y + sm[2]*coors[2][i].y;

			J[c] = x_s[c]*y_t[c] - x_t[c]*y_s[c];
			jac[0][i] = J[c];
			vec__x_s[0][i] = x_s[c];
			vec__y_s[0][i] = y_s[c];
			vec__x_t[0][i] = x_t[c];
			vec__y_t[0][i] = y_t[c];

			g_12 = x_s[c]*x_t[c] + y_s[c]*y_t[c];
			if (fabs(J[c]) < feps) g12 = 0.0;
			else g12 = -1./(J[c]*J[c]) * g_12;

			coef[MET_SW][0][i+1] += .25 * J[c] * g12;
			coef[MET_SW][1][i] += .25 * J[c] * g12;
			coef[MET_SE][1][i] += -.25 * J[c] * g12;
			coef[MET_SE][0][i-1] += -.25 * J[c] * g12;

		}

		for (i = ibeg; i < xs + xm; i++) {
			for (k = 0; k < 3; k++) {
				coors_xh[k].x = .5 * (coors[k][i-1].x + coors[k][i].x);
				coors_xh[k].y = .5 * (coors[k][i-1].y + coors[k][i].y);
			}

//Andrew can you check the four lines below???
			x_s[x_h] = coors[0][i].x - coors[0][i-1].x;
			y_s[x_h] = coors[0][i].y - coors[0][i-1].y;
			x_t[x_h] = sm[0]*coors_xh[0].x + sm[1]*coors_xh[1].x + sm[2]*coors_xh[2].x;
			y_t[x_h] = sm[0]*coors_xh[0].y + sm[1]*coors_xh[1].y + sm[2]*coors_xh[2].y;

			J[x_h] = x_s[x_h]*y_t[x_h] - x_t[x_h]*y_s[x_h];

			g_22 = x_t[x_h]*x_t[x_h] + y_t[x_h]*y_t[x_h];
			if (fabs(J[x_h]) > feps)
				g11 = 1./(J[x_h]*J[x_h]) * g_22;
			else
				g11 = 0;
			coef[MET_W][0][i] = J[x_h] * g11;
		}
	}

	// top boundary
	if (ys + ym == ngy && !per_y) {
		int jdx = ngy - 1;
		for (i = ibeg; i < iend; i++) {
			x_s[c] = si[-1]*coors[jdx][i-1].x + si[1]*coors[jdx][i+1].x;
			y_s[c] = si[-1]*coors[jdx][i-1].y + si[1]*coors[jdx][i+1].y;
			x_t[c] = sp[0]*coors[jdx][i].x + sp[-1]*coors[jdx-1][i].x + sp[-2]*coors[jdx-2][i].x;
			y_t[c] = sp[0]*coors[jdx][i].y + sp[-1]*coors[jdx-1][i].y + sp[-2]*coors[jdx-2][i].y;

			J[c] = x_s[c]*y_t[c] - x_t[c]*y_s[c];
			jac[jdx][i] = J[c];
			vec__x_s[jdx][i] = x_s[c];
			vec__y_s[jdx][i] = y_s[c];
			vec__x_t[jdx][i] = x_t[c];
			vec__y_t[jdx][i] = y_t[c];

			g_12 = x_s[c]*x_t[c] + y_s[c]*y_t[c];
			if (fabs(J[c]) < feps) g12 = 0.0;
			else g12 = -1./(J[c]*J[c]) * g_12;

			coef[MET_SW][jdx][i+1] += .25 * J[c] * g12;
			coef[MET_SE][jdx][i-1] += -.25 * J[c] * g12;
		}

		for (i = ibeg; i < xs + xm; i++) {
			for (k = 0; k < 3; k++) {
				coors_xh[k].x = .5 * (coors[jdx-k][i-1].x + coors[jdx-k][i].x);
				coors_xh[k].y = .5 * (coors[jdx-k][i-1].y + coors[jdx-k][i].y);
			}

			x_s[x_h] = coors[jdx][i].x - coors[jdx][i-1].x;
			y_s[x_h] = coors[jdx][i].y - coors[jdx][i-1].y;
			x_t[x_h] = sp[0]*coors_xh[0].x + sp[-1]*coors_xh[1].x + sp[-2]*coors_xh[2].x;
			y_t[x_h] = sp[0]*coors_xh[0].y + sp[-1]*coors_xh[1].y + sp[-2]*coors_xh[2].y;

			J[x_h] = x_s[x_h]*y_t[x_h] - x_t[x_h]*y_s[x_h];

			g_22 = x_t[x_h]*x_t[x_h] + y_t[x_h]*y_t[x_h];
			if (fabs(J[x_h]) < feps) g11 = 0.0;
			else g11 = 1./(J[x_h]*J[x_h]) * g_22;
			coef[MET_W][jdx][i] = J[x_h] * g11;
		}
	}


	// South West corner
	if (xs == 0 && ys == 0) {
		x_s[c] = sm[0]*coors[0][0].x + sm[1]*coors[0][1].x + sm[2]*coors[0][2].x;
		y_s[c] = sm[0]*coors[0][0].y + sm[1]*coors[0][1].y + sm[2]*coors[0][2].y;
		x_t[c] = sm[0]*coors[0][0].x + sm[1]*coors[1][0].x + sm[2]*coors[2][0].x;
		y_t[c] = sm[0]*coors[0][0].y + sm[1]*coors[1][0].y + sm[2]*coors[2][0].y;

		J[c] = x_s[c]*y_t[c] - x_t[c]*y_s[c];
		jac[0][0] = J[c];
		vec__x_s[0][0] = x_s[c];
		vec__y_s[0][0] = y_s[c];
		vec__x_t[0][0] = x_t[c];
		vec__y_t[0][0] = y_t[c];

		g_12 = x_s[c]*x_t[c] + y_s[c]*y_t[c];
		if (fabs(J[c]) < feps) g12 = 0.0;
		else g12 = -1./(J[c]*J[c]) * g_12;

		coef[MET_SW][0][1] += .25 * J[c] * g12;
		coef[MET_SW][1][0] += .25 * J[c] * g12;
		coef[MET_SE][1][0] += -.25 * J[c] * g12;
	}

	// North West corner
	if (xs == 0 && ys + ym == ngy) {
		int jdx = ngy - 1;
		x_s[c] = sm[0]*coors[jdx][0].x + sm[1]*coors[jdx][1].x + sm[2]*coors[jdx][2].x;
		y_s[c] = sm[0]*coors[jdx][0].y + sm[1]*coors[jdx][1].y + sm[2]*coors[jdx][2].y;
		x_t[c] = sp[0]*coors[jdx][0].x + sp[-1]*coors[jdx-1][0].x + sp[-2]*coors[jdx-2][0].x;
		y_t[c] = sp[0]*coors[jdx][0].y + sp[-1]*coors[jdx-1][0].y + sp[-2]*coors[jdx-2][0].y;

		J[c] = x_s[c]*y_t[c] - x_t[c]*y_s[c];
		jac[jdx][0] = J[c];
		vec__x_s[jdx][0] = x_s[c];
		vec__y_s[jdx][0] = y_s[c];
		vec__x_t[jdx][0] = x_t[c];
		vec__y_t[jdx][0] = y_t[c];

		g_12 = x_s[c]*x_t[c] + y_s[c]*y_t[c];
		if (fabs(J[c]) < feps) g12 = 0.0;
		else g12 = -1./(J[c]*J[c]) * g_12;

		coef[MET_SW][jdx][1] += .25 * J[c] * g12;
	}

	// South East corner
	if (xs + xm == ngx && ys == 0) {
		int idx = ngx - 1;
		x_s[c] = sp[0]*coors[0][idx].x + sp[-1]*coors[0][idx-1].x + sp[-2]*coors[0][idx-2].x;
		y_s[c] = sp[0]*coors[0][idx].y + sp[-1]*coors[0][idx-1].y + sp[-2]*coors[0][idx-2].y;
		x_t[c] = sm[0]*coors[0][idx].x + sm[1]*coors[1][idx].x + sm[2]*coors[2][idx].x;
		y_t[c] = sm[0]*coors[0][idx].y + sm[1]*coors[1][idx].y + sm[2]*coors[2][idx].y;

		J[c] = x_s[c]*y_t[c] - x_t[c]*y_s[c];
		jac[0][idx] = J[c];
		vec__x_s[0][idx] = x_s[c];
		vec__y_s[0][idx] = y_s[c];
		vec__x_t[0][idx] = x_t[c];
		vec__y_t[0][idx] = y_t[c];


		g_12 = x_s[c]*x_t[c] + y_s[c]*y_t[c];
		if (fabs(J[c]) < feps) g12 = 0.0;
		else g12 = -1./(J[c]*J[c]) * g_12;

		coef[MET_SW][1][idx] += .25 * J[c] * g12;
		coef[MET_SE][1][idx] += -.25 * J[c] * g12;
		coef[MET_SE][0][idx-1] += -.25 * J[c] * g12;
	}

	// North East corner
	if (xs + xm == ngx && ys + ym == ngy) {
		int idx = ngx - 1;
		int jdx = ngy - 1;

		x_s[c] = sp[0]*coors[jdx][idx].x + sp[-1]*coors[jdx][idx-1].x + sp[-2]*coors[jdx][idx-2].x;
		y_s[c] = sp[0]*coors[jdx][idx].y + sp[-1]*coors[jdx][idx-1].y + sp[-2]*coors[jdx][idx-2].y;
		x_t[c] = sp[0]*coors[jdx][idx].x + sp[-1]*coors[jdx-1][idx].x + sp[-2]*coors[jdx-2][idx].x;
		y_t[c] = sp[0]*coors[jdx][idx].y + sp[-1]*coors[jdx-1][idx].y + sp[-2]*coors[jdx-2][idx].y;

		J[c] = x_s[c]*y_t[c] - x_t[c]*y_s[c];
		jac[jdx][idx] = J[c];
		vec__x_s[jdx][idx] = x_s[c];
		vec__y_s[jdx][idx] = y_s[c];
		vec__x_t[jdx][idx] = x_t[c];
		vec__y_t[jdx][idx] = y_t[c];

		g_12 = x_s[c]*x_t[c] + y_s[c]*y_t[c];
		if (fabs(J[c]) < feps) g12 = 0.0;
		else g12 = -1./(J[c]*J[c]) * g_12;

		coef[MET_SE][jdx][idx-1] += -.25 * J[c] * g12;
	}


	for (i = 0; i < 5; i++) {
		ierr = DMDAVecRestoreArray(da, met->lcoef[i], &coef[i]);CHKERRQ(ierr);
		ierr = VecSet(met->coef[i], 0.0);CHKERRQ(ierr);
		ierr = DMLocalToGlobalBegin(da, met->lcoef[i], ADD_VALUES, met->coef[i]);CHKERRQ(ierr);
		ierr = DMLocalToGlobalEnd(da, met->lcoef[i], ADD_VALUES, met->coef[i]);CHKERRQ(ierr);
		ierr = DMGlobalToLocalBegin(da, met->coef[i], INSERT_VALUES, met->lcoef[i]);CHKERRQ(ierr);
		ierr = DMGlobalToLocalEnd(da, met->coef[i], INSERT_VALUES, met->lcoef[i]);CHKERRQ(ierr);
	}

	for (i = 0; i < 4; i++) {
		ierr = DMDAVecRestoreArray(da, met->jac_v[i], &jac_v[i]);CHKERRQ(ierr);
	}

	ierr = DMDAVecRestoreArray(da, met->jac, &jac);CHKERRQ(ierr);

	ierr = DMDAVecRestoreArray(cda, cv, &coors);CHKERRQ(ierr);

	return 0;
}


static PetscErrorCode compute_2d(stella_metric *met, DM da, stella_fd *fd)
{
	PetscErrorCode ierr;

	ierr = compute_2d_interior(met, da, fd);CHKERRQ(ierr);
	ierr = compute_2d_boundary(met, da, fd);CHKERRQ(ierr);

	return 0;
}


static PetscErrorCode compute_3d(stella_metric *met, DM da, stella_fd *fd)
{
	PetscErrorCode ierr;

	ierr = compute_3d_interior(met, da, fd);CHKERRQ(ierr);
	ierr = compute_3d_boundary(met, da, fd);CHKERRQ(ierr);

	return 0;
}


PetscErrorCode stella_metric_create(stella_metric **metptr, DM da, stella_fd *fd)
{
	PetscErrorCode ierr;
	int i;
	PetscInt dim;
	stella_metric *met = (stella_metric*) malloc(sizeof(stella_metric));

	ierr = DMDAGetInfo(da, &dim, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL);CHKERRQ(ierr);

	met->t2map[0][0] = 0;
	met->t2map[0][1] = 1;
	met->t2map[1][0] = 2;
	met->t2map[1][1] = 3;

	met->t3map[0][0] = 0;
	met->t3map[0][1] = 1;
	met->t3map[0][2] = 2;
	met->t3map[1][0] = 3;
	met->t3map[2][0] = 4;
	met->t3map[1][1] = 5;
	met->t3map[1][2] = 6;
	met->t3map[2][1] = 7;
	met->t3map[2][2] = 8;

	ierr = DMCreateGlobalVector(da, &met->coef[0]);CHKERRQ(ierr);
	ierr = DMCreateLocalVector(da, &met->lcoef[0]);CHKERRQ(ierr);
	ierr = VecDuplicate(met->coef[0], &met->jac);CHKERRQ(ierr);

	met->nd = dim;

	if (dim == 2) {
		for (i=0; i < 4; i++) {
			ierr = VecDuplicate(met->coef[0], &met->jac_v[i]);CHKERRQ(ierr);
		}

		for (i = 1; i < 5; i++) {
			ierr = VecDuplicate(met->coef[0], &met->coef[i]);CHKERRQ(ierr);
			ierr = VecDuplicate(met->lcoef[0], &met->lcoef[i]);CHKERRQ(ierr);
		}

		ierr = compute_2d(met, da, fd);CHKERRQ(ierr);
	} else if (dim == 3) {
		for (i=0; i < 9; i++) {
			ierr = VecDuplicate(met->coef[0], &met->jac_v[i]);CHKERRQ(ierr);
		}

		for (i = 1; i < 10; i++) {
			ierr = VecDuplicate(met->coef[0], &met->coef[i]);CHKERRQ(ierr);
			ierr = VecDuplicate(met->lcoef[0], &met->lcoef[i]);CHKERRQ(ierr);
		}

		ierr = compute_3d(met, da, fd);CHKERRQ(ierr);
	}

	*metptr = met;

	return 0;
}


PetscErrorCode stella_metric_destroy(stella_metric *metric)
{
	PetscErrorCode ierr;
	int i;

	int ncoef = 0;
	int njac = 0;

	if (metric->nd == 2) {
		ncoef = 5;
		njac = 4;
	} else if (metric->nd == 3) {
		ncoef = 10;
		njac = 9;
	}

	for (i = 0; i < ncoef; i++) {
		ierr = VecDestroy(&metric->coef[i]);CHKERRQ(ierr);
	}

	for (i = 0; i < njac; i++) {
		ierr = VecDestroy(&metric->jac_v[i]);CHKERRQ(ierr);
	}

	for (i = 0; i < ncoef; i++) {
		ierr = VecDestroy(&metric->lcoef[i]);CHKERRQ(ierr);
	}

	ierr = VecDestroy(&metric->jac);CHKERRQ(ierr);

	return 0;
}
