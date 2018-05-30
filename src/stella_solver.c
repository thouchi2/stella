#include <stdlib.h>

#include "stella_solver.h"
#include "stella_io.h"

#include "stella_mat.h"
#include "stella_pc.h"

#ifdef WITH_BOXMG
#include <boxmg/capi.h>
#endif

#include "petscmat.h"

static inline int smallerValue(int i1, int i2, double v1, double v2)
{
	return   (v1==v2) ? i1 : ( v1<v2 ? i1 : i2 );
}

static inline PetscErrorCode petsc_options_has_name(const char name[], PetscBool *set) {

	PetscErrorCode ierr;

	#ifdef PETSC_OPTIONS_TAKES_DATABASE
	ierr = PetscOptionsHasName(PETSC_NULL, PETSC_NULL, name, set);CHKERRQ(ierr);
	#else
	ierr = PetscOptionsHasName(PETSC_NULL, name, set);CHKERRQ(ierr);
	#endif

	return 0;
}

static inline PetscErrorCode petsc_options_set_value(const char name[], const char value[]) {

	PetscErrorCode ierr;

	#ifdef PETSC_OPTIONS_TAKES_DATABASE
	ierr = PetscOptionsSetValue(PETSC_NULL, name, value);CHKERRQ(ierr);
	#else
	ierr = PetscOptionsSetValue(name, value);CHKERRQ(ierr);
	#endif

	return 0;
}


/**
 * Calls callbacks required when the rhs is changed
 */
static PetscErrorCode stella_setup_rhs(stella *slv)
{
	PetscErrorCode ierr;
	PetscInt i, j, k, xs, ys, zs, xm, ym, zm, ngx, ngy, ngz;
	Vec b;

	b = slv->rhs;

	ierr = DMDAGetCorners(slv->dm, &xs, &ys, &zs, &xm, &ym, &zm);CHKERRQ(ierr);
	ierr = DMDAGetInfo(slv->dm, 0, &ngx, &ngy, &ngz,0,0,0,0,0,0,0,0,0);CHKERRQ(ierr);

	if (slv->grid.nd == 2) {
		double **rhs, **jump;
		PetscScalar **bvec, **dcoef;

		ierr = stella_dmap_get(slv->dmap, slv->state.rhs, &rhs);CHKERRQ(ierr);
		ierr = stella_dmap_get(slv->dmap, slv->state.jump, &jump);CHKERRQ(ierr);
		ierr = DMDAVecGetArray(slv->dm, slv->level.ldcoef, &dcoef);CHKERRQ(ierr);
		ierr = DMDAVecGetArray(slv->dm, b, &bvec);CHKERRQ(ierr);
		if (slv->op->axisymmetric) {
			DM cda;
			DMDACoor2d **coors;
			Vec gc;

			ierr = DMGetCoordinateDM(slv->dm, &cda);CHKERRQ(ierr);
			ierr = DMGetCoordinates(slv->dm, &gc);CHKERRQ(ierr);
			ierr = DMDAVecGetArray(cda, gc, &coors);CHKERRQ(ierr);

			for (j = ys; j < ys + ym; j++) {
				for (i = xs; i < xs + xm; i++) {
					bvec[j][i] = rhs[j][i] * coors[j][i].x;
				}
			}

			ierr = DMDAVecRestoreArray(cda, gc, &coors);CHKERRQ(ierr);
		} else {
			double fxp, fxm, fyp, fym;
			PetscScalar **jac[4];
			PetscScalar **x_s, **x_t, **y_s, **y_t;
			stella_metric *met;
			met = slv->level.metric;

			for (i = 0; i < 4; i++) {
				ierr = DMDAVecGetArray(slv->dm, met->jac_v[i], &jac[i]);CHKERRQ(ierr);
			}

			x_s = jac[met->t2map[0][0]];
			x_t = jac[met->t2map[0][1]];
			y_s = jac[met->t2map[1][0]];
			y_t = jac[met->t2map[1][1]];

			for (j = ys; j < ys + ym; j++) {
				for (i = xs; i < xs + xm; i++) {

					int ip = smallerValue(i, i+1, dcoef[j][i], dcoef[j][i+1]);
					int im = smallerValue(i, i-1, dcoef[j][i], dcoef[j][i-1]);
					int jp = smallerValue(j, j+1, dcoef[j][i], dcoef[j+1][i]);
					int jm = smallerValue(j, j-1, dcoef[j][i], dcoef[j-1][i]);

					if ((i != ngx-1) && (dcoef[j][i] != dcoef[j][i+1]))
						fxp = 2.0*dcoef[j][i]*jump[j][ip] / (dcoef[j][i] + dcoef[j][i+1]) / (x_s[j][i]+x_s[j][i+1]);
					else fxp = 0;
					if ((i != 0) && (dcoef[j][i] != dcoef[j][i-1]))
						fxm = 2.0*dcoef[j][i]*jump[j][im] / (dcoef[j][i] + dcoef[j][i-1]) / (x_s[j][i]+x_s[j][i-1]);
					else fxm = 0;
					if ((j != ngy-1) && (dcoef[j][i] != dcoef[j+1][i]))
						fyp = 2.0*dcoef[j][i]*jump[jp][i] / (dcoef[j][i] + dcoef[j+1][i]) / (y_t[j][i]+y_t[j+1][i]);
					else fyp = 0;
					if ((j != 0) && (dcoef[j][i] != dcoef[j-1][i]))
						fym = 2.0*dcoef[j][i]*jump[jm][i] / (dcoef[j][i] + dcoef[j-1][i]) / (y_t[j][i]+y_t[j-1][i]);
					else fym = 0;
					bvec[j][i] = rhs[j][i] - fxp - fxm - fyp - fym;
				}
			}

			for (i = 0; i < 4; i++) {
				ierr = DMDAVecRestoreArray(slv->dm, met->jac_v[i], &jac[i]);CHKERRQ(ierr);

			}
		}

		ierr = stella_dmap_restore(slv->dmap, &jump);CHKERRQ(ierr);
		ierr = stella_dmap_restore(slv->dmap, &rhs);CHKERRQ(ierr);
		ierr = DMDAVecRestoreArray(slv->dm, b, &bvec);CHKERRQ(ierr);
	} else {
		double ***rhs, ***jump;
		PetscScalar ***bvec, ***dcoef;

		double fxp, fxm, fyp, fym, fzp, fzm;
		PetscScalar ***jac[6];
		PetscScalar ***x_r, ***y_s, ***z_t;
		stella_metric *met;
		met = slv->level.metric;

		for (i = 0; i < 6; i++) {
			ierr = DMDAVecGetArray(slv->dm, met->jac_v[i], &jac[i]);CHKERRQ(ierr);
		}

		x_r = jac[met->t3map[0][0]];
		y_s = jac[met->t3map[1][1]];
		z_t = jac[met->t3map[2][2]];

		ierr = DMDAVecGetArray(slv->dm, slv->level.ldcoef, &dcoef);CHKERRQ(ierr);
		ierr = stella_dmap_get(slv->dmap, slv->state.rhs, &rhs);CHKERRQ(ierr);
		ierr = stella_dmap_get(slv->dmap, slv->state.jump, &jump);CHKERRQ(ierr);
		ierr = DMDAVecGetArray(slv->dm, b, &bvec);CHKERRQ(ierr);
		for (k = zs; k < zs + zm; k++) {
			for (j = ys; j < ys + ym; j++) {
				for (i = xs; i < xs + xm; i++) {
					int ip = smallerValue(i, i+1, dcoef[k][j][i], dcoef[k][j][i+1]);
					int im = smallerValue(i, i-1, dcoef[k][j][i], dcoef[k][j][i-1]);
					int jp = smallerValue(j, j+1, dcoef[k][j][i], dcoef[k][j+1][i]);
					int jm = smallerValue(j, j-1, dcoef[k][j][i], dcoef[k][j-1][i]);
					int kp = smallerValue(k, k+1, dcoef[k][j][i], dcoef[k+1][j][i]);
					int km = smallerValue(k, k-1, dcoef[k][j][i], dcoef[k-1][j][i]);

					if ((i != ngx-1) && (dcoef[k][j][i] != dcoef[k][j][i+1]))
						fxp = 2.0*dcoef[k][j][i]*jump[k][j][ip] / (dcoef[k][j][i] + dcoef[k][j][i+1]) / (x_r[k][j][i]+x_r[k][j][i+1]);
					else fxp = 0;
					if ((i != 0) && (dcoef[k][j][i] != dcoef[k][j][i-1]))
						fxm = 2.0*dcoef[k][j][i]*jump[k][j][im] / (dcoef[k][j][i] + dcoef[k][j][i-1]) / (x_r[k][j][i]+x_r[k][j][i-1]);
					else fxm = 0;
					if ((j != ngy-1) && (dcoef[k][j][i] != dcoef[k][j+1][i]))
						fyp = 2.0*dcoef[k][j][i]*jump[k][jp][i] / (dcoef[k][j][i] + dcoef[k][j+1][i]) / (y_s[k][j][i]+y_s[k][j+1][i]);
					else fyp = 0;
					if ((j != 0) && (dcoef[k][j][i] != dcoef[k][j-1][i]))
						fym = 2.0*dcoef[k][j][i]*jump[k][jm][i] / (dcoef[k][j][i] + dcoef[k][j-1][i]) / (y_s[k][j][i]+y_s[k][j-1][i]);
					else fym = 0;
					if ((k != ngz-1) && (dcoef[k][j][i] != dcoef[k+1][j][i]))
						fzp = 2.0*dcoef[k][j][i]*jump[kp][j][i] / (dcoef[k][j][i] + dcoef[k+1][j][i]) / (z_t[k][j][i]+z_t[k+1][j][i]);
					else fzp = 0;
					if ((k != 0) && (dcoef[k][j][i] != dcoef[k-1][j][i]))
						fzm = 2.0*dcoef[k][j][i]*jump[km][j][i] / (dcoef[k][j][i] + dcoef[k-1][j][i]) / (z_t[k][j][i]+z_t[k-1][j][i]);
					else fzm = 0;


					bvec[k][j][i] = rhs[k][j][i] - fxp - fxm - fyp - fym - fzp - fzm;
				}
			}
		}

		for (i = 0; i < 6; i++) {
			ierr = DMDAVecRestoreArray(slv->dm, met->jac_v[i], &jac[i]);CHKERRQ(ierr);
		}
		ierr = stella_dmap_restore(slv->dmap, &rhs);CHKERRQ(ierr);
		ierr = stella_dmap_restore(slv->dmap, &jump);CHKERRQ(ierr);
		ierr = DMDAVecRestoreArray(slv->dm, b, &bvec);CHKERRQ(ierr);
		ierr = DMDAVecRestoreArray(slv->dm, slv->level.ldcoef, &dcoef);CHKERRQ(ierr);
	}

	ierr = VecPointwiseMult(b, slv->level.metric->jac, b);CHKERRQ(ierr);

	ierr = stella_boundary_apply_rhs(slv->boundary, slv->dm, slv->rhs);CHKERRQ(ierr);

	ierr = VecPointwiseMult(slv->level.add_cont, slv->level.pm, slv->level.add_cont);CHKERRQ(ierr);
	ierr = VecAXPY(b, 1, slv->level.add_cont);CHKERRQ(ierr);
	ierr = VecPointwiseMult(b, slv->level.nscale, b);CHKERRQ(ierr);

	if (stella_log(slv, STELLA_LOG_PROBLEM)) {
		PetscViewer bout;
		PetscViewerASCIIOpen(slv->comm, "b.txt", &bout);
		ierr = VecView(b, bout);CHKERRQ(ierr);
	}

	return 0;
}


/**
 * Creates and initializes data needed for setting up elliptic discretization.
 */
static PetscErrorCode stella_setup(stella *slv, int offset[], int stride[])
{
	PetscErrorCode ierr;
	PetscInt xs, ys, zs, xm, ym, zm, ngx, ngy, ngz;
	stella_pc *shell;
	PCType pc_type;

	if (stella_log(slv, STELLA_LOG_STATUS)) {
		ierr = stella_io_print(slv->comm, "Setting up electric field solver");CHKERRQ(ierr);
	}

	slv->ts = 0;

	if (stella_log(slv, STELLA_LOG_RESIDUAL)) {
		ierr = petsc_options_set_value("-ksp_monitor_short", NULL);CHKERRQ(ierr);
	}

	ierr = DMDASetFieldName(slv->dm, 0,"potential");CHKERRQ(ierr);

	if (slv->grid.nd == 2) {
		ierr = DMDAGetCorners(slv->dm, &xs, &ys, 0, &xm, &ym, 0);CHKERRQ(ierr);
		slv->dmap = stella_dmap_create_2d(xs - offset[0], ys - offset[1], xm, ym, stride);
	} else if (slv->grid.nd == 3) {
		ierr = DMDAGetCorners(slv->dm, &xs, &ys, &zs, &xm, &ym, &zm);CHKERRQ(ierr);
		slv->dmap = stella_dmap_create_3d(xs - offset[0], ys - offset[1], zs - offset[2], xm, ym, zm, stride);
	} else {
		SETERRQ1(PETSC_COMM_WORLD, PETSC_ERR_SUP, "Unsupported dimmension: %d", slv->grid.nd);
	}
	ierr = DMDAGetInfo(slv->dm, 0, &ngx, &ngy, &ngz,0,0,0,0,0,0,0,0,0);CHKERRQ(ierr);

	ierr = KSPCreate(slv->comm, &slv->ksp);CHKERRQ(ierr);
	if (stella_log(slv, STELLA_LOG_EIGS)) {
		ierr = KSPSetComputeEigenvalues(slv->ksp, PETSC_TRUE);CHKERRQ(ierr);
	}

	ierr = KSPSetFromOptions(slv->ksp);CHKERRQ(ierr);

	if (!slv->options.algebraic) {
		ierr = KSPGetPC(slv->ksp, &slv->pc);
		ierr = PCSetType(slv->pc, PCSHELL);CHKERRQ(ierr);
		ierr = stella_pc_create(&shell);CHKERRQ(ierr);
		ierr = PCShellSetApply(slv->pc, stella_pc_apply);CHKERRQ(ierr);
		ierr = PCShellSetDestroy(slv->pc, stella_pc_destroy);CHKERRQ(ierr);
		ierr = PCShellSetContext(slv->pc, shell);CHKERRQ(ierr);
		ierr = PCShellSetName(slv->pc, "boxmg");CHKERRQ(ierr);
		ierr = PCShellSetSetUp(slv->pc, stella_pc_setup);CHKERRQ(ierr);
	}

	ierr = stella_fd_create(&slv->fd, STELLA_FD_STANDARD_O2);CHKERRQ(ierr);
	ierr = stella_operator_create(&slv->op, &slv->level, slv->fd, slv->grid.nd);CHKERRQ(ierr);
	ierr = stella_boundary_create(&slv->boundary, &slv->level,
	                              slv->dmap, &slv->state, slv->fd);CHKERRQ(ierr);
	slv->op->axisymmetric = slv->options.axisymmetric;
	slv->boundary->axisymmetric = slv->options.axisymmetric;

	ierr = DMCreateGlobalVector(slv->dm, &slv->level.dcoef);CHKERRQ(ierr);
	ierr = DMCreateLocalVector(slv->dm, &slv->level.ldcoef);CHKERRQ(ierr);
	ierr = DMCreateGlobalVector(slv->dm, &slv->level.bcoef);CHKERRQ(ierr);
	ierr = DMCreateLocalVector(slv->dm, &slv->level.lbcoef);CHKERRQ(ierr);
	ierr = DMCreateGlobalVector(slv->dm, &slv->level.add_cont);CHKERRQ(ierr);
	ierr = DMCreateLocalVector(slv->dm, &slv->level.ladd_cont);CHKERRQ(ierr);
	ierr = DMCreateGlobalVector(slv->dm, &slv->level.nscale);CHKERRQ(ierr);
	ierr = DMCreateGlobalVector(slv->dm, &slv->level.pm);CHKERRQ(ierr);
	ierr = VecSet(slv->level.pm, 1.0);CHKERRQ(ierr);
	ierr = VecSet(slv->level.add_cont, 0.0);CHKERRQ(ierr);
	ierr = VecSet(slv->level.nscale, 1.0);CHKERRQ(ierr);
	ierr = DMCreateGlobalVector(slv->dm, &slv->rhs);CHKERRQ(ierr);
	ierr = DMCreateGlobalVector(slv->dm, &slv->x);CHKERRQ(ierr);

	if (!slv->options.algebraic) {
		#ifdef WITH_BOXMG
		ef_bmg_mat *mat_ctx = (ef_bmg_mat*) malloc(sizeof(ef_bmg_mat));
		if (slv->grid.nd == 2) {
			mat_ctx->op2 = bmg2_operator_create(slv->grid.topo2);
			mat_ctx->nd = 2;
			ierr = MatCreateShell(slv->comm, xm*ym, xm*ym, ngx*ngy, ngx*ngy, mat_ctx, &slv->A);CHKERRQ(ierr);

		} else {
			ef_bmg_mat *mat_ctx = (ef_bmg_mat*) malloc(sizeof(ef_bmg_mat));
			mat_ctx->op3 = bmg3_operator_create(slv->grid.topo3);
			mat_ctx->nd = 3;
			ierr = MatCreateShell(slv->comm, xm*ym*zm, xm*ym*zm, ngx*ngy*ngz, ngx*ngy*ngz, mat_ctx, &slv->A);CHKERRQ(ierr);
		}
		ierr = MatShellSetOperation(slv->A, MATOP_MULT, (void(*)(void)) ef_bmg_mult);CHKERRQ(ierr);
		#endif
	} else {
		ierr = DMCreateMatrix(slv->dm, &slv->A);CHKERRQ(ierr);
	}

	if (stella_log(slv, STELLA_LOG_STATUS)) {
		ierr = stella_io_print(slv->comm, "set-up is finished");CHKERRQ(ierr);
	}
	return 0;
}


PetscErrorCode stella_init(stella **solver_ctx, MPI_Comm comm,
                           int nGlobal[], int nProcs[], int nLocal[],
                           int offset[], int stride[],
                           int cartCoord[], int periodic[], int periodic_storage,
                           int nd, int axisymmetric)
{
	PetscErrorCode ierr;

	stella *slv = (stella*) malloc(sizeof(stella));
	slv->comm = comm;
	slv->options.axisymmetric = axisymmetric;
	{
		PetscBool flg;
		ierr = petsc_options_has_name("-boxmg", &flg);CHKERRQ(ierr);
		if (flg) {
		#ifdef WITH_BOXMG
			slv->options.algebraic = 0;
		#endif
		} else {
			slv->options.algebraic = 1;
		}
	}

	stella_set_log(slv, STELLA_LOG_STATUS | STELLA_LOG_RESIDUAL);
	ierr = stella_grid_setup(&slv->grid, &slv->dm, &slv->comm, nGlobal, nProcs,
	                         nLocal, cartCoord, periodic, periodic_storage, nd);CHKERRQ(ierr);
	slv->level.dm = slv->dm;
	ierr = stella_setup(slv, offset, stride);CHKERRQ(ierr);

	*solver_ctx = slv;

	return 0;
}


static PetscErrorCode stella_updatebc(stella *slv)
{
	PetscErrorCode ierr;
	ierr = VecSet(slv->level.pm, 1.0);CHKERRQ(ierr);
	ierr = VecSet(slv->level.add_cont, 0.0);CHKERRQ(ierr);
	ierr = VecSet(slv->level.nscale, 1.0);CHKERRQ(ierr);

	ierr = stella_boundary_apply(slv->boundary, slv->A, slv->dm);CHKERRQ(ierr);

	if (slv->options.algebraic) {
		ierr = MatAssemblyBegin(slv->A, MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
		ierr = MatAssemblyEnd(slv->A, MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
		if (stella_log(slv, STELLA_LOG_PROBLEM)) {
			PetscViewer vout;
			ierr = PetscViewerBinaryOpen(slv->comm, "a.dat", FILE_MODE_WRITE, &vout);CHKERRQ(ierr);
			ierr = MatView(slv->A, vout);CHKERRQ(ierr);
			ierr = PetscViewerDestroy(&vout);CHKERRQ(ierr);
		}
	} else {
		#ifdef WITH_BOXMG
		if (stella_log(slv, STELLA_LOG_PROBLEM)) {
			stella_bmg_mat *ctx;
			ierr = MatShellGetContext(slv->A, (void**) &ctx);CHKERRQ(ierr);
			if (ctx->nd == 2)
				bmg2_operator_dump(ctx->op2);
			else
				bmg3_operator_dump(ctx->op3);
		}
		#endif
	}

	ierr = KSPSetOperators(slv->ksp, slv->A, slv->A);CHKERRQ(ierr);

	ierr = stella_setup_rhs(slv);CHKERRQ(ierr);

	return 0;
}


PetscErrorCode stella_solve(stella *slv)
{
	PetscErrorCode ierr;
	PetscInt i, j, k, xs, ys, zs, xm, ym, zm;
	PetscInt ngx,ngy,ngz;
	DMBoundaryType bx, by, bz;
	PetscBool boxmg_direct = 0;
	slv->ts++;

	ierr = petsc_options_has_name("-boxmg_direct", &boxmg_direct);CHKERRQ(ierr);

	ierr = stella_updatebc(slv);CHKERRQ(ierr);
	ierr = VecSet(slv->x, 0);CHKERRQ(ierr);
	if (!slv->options.algebraic && boxmg_direct) {
		PC bmg_pc;
		ierr = KSPGetPC(slv->ksp, &bmg_pc);CHKERRQ(ierr);
		ierr = PCApply(bmg_pc, slv->rhs, slv->x);CHKERRQ(ierr);
	} else {
		ierr = KSPSolve(slv->ksp, slv->rhs, slv->x);CHKERRQ(ierr);
	}

	if (slv->state.sol) {
	  ierr = stella_io_vtkwrite(slv->dm, slv->x, "phi", 0, slv->ts);CHKERRQ(ierr);
	  ierr = stella_io_vtkwrite(slv->dm, slv->sol, "sol", 0, 0);CHKERRQ(ierr);
	  Vec err;
	  ierr = VecDuplicate(slv->sol, &err);CHKERRQ(ierr);
	  ierr = VecWAXPY(err, -1, slv->x, slv->sol);CHKERRQ(ierr);
	  ierr = VecAbs(err);CHKERRQ(ierr);
	  ierr = stella_io_vtkwrite(slv->dm, err, "err", 0, 0);CHKERRQ(ierr);
	  ierr = VecDestroy(&err);CHKERRQ(ierr);
	}


	ierr = DMDAGetInfo(slv->dm, 0, &ngx, &ngy, &ngz,0,0,0,0,0, &bx, &by, &bz,0);CHKERRQ(ierr);

	ierr = DMDAGetCorners(slv->dm, &xs, &ys, &zs, &xm, &ym, &zm);CHKERRQ(ierr);

	if (slv->grid.nd == 2) {
		Vec loc_x;
		PetscScalar **xvec;
		double **phi;

		ierr = stella_dmap_get(slv->dmap, slv->state.phi, &phi);CHKERRQ(ierr);
		ierr = DMGetLocalVector(slv->dm, &loc_x);CHKERRQ(ierr);
		ierr = DMGlobalToLocalBegin(slv->dm, slv->x, INSERT_VALUES, loc_x);CHKERRQ(ierr);
		ierr = DMGlobalToLocalEnd(slv->dm, slv->x, INSERT_VALUES, loc_x);CHKERRQ(ierr);
		ierr = DMDAVecGetArray(slv->dm, loc_x, &xvec);CHKERRQ(ierr);

		for (j = ys; j < ys + ym; j++) {
			for (i = xs; i < xs + xm; i++) {
				phi[j][i] = xvec[j][i];
			}
		}

		if ((by == DM_BOUNDARY_PERIODIC) && (ys + ym == ngy) && slv->grid.overlap_periodic) {
			for (i = xs; i < xs + xm; i++) {
				phi[ngy][i] = xvec[ngy][i];
			}
		}

		if ((bx == DM_BOUNDARY_PERIODIC) && (xs + xm == ngx) && slv->grid.overlap_periodic) {
			for (j = ys; j < ys + ym; j++) {
				phi[j][ngx] = xvec[j][ngx];
			}
		}

		ierr = stella_dmap_restore(slv->dmap, &phi);CHKERRQ(ierr);
		ierr = DMDAVecRestoreArray(slv->dm, loc_x, &xvec);CHKERRQ(ierr);
		if (stella_log(slv, STELLA_LOG_VTK)) {
			ierr = DMLocalToGlobalBegin(slv->dm, loc_x, INSERT_VALUES, slv->x);CHKERRQ(ierr);
			ierr = DMLocalToGlobalEnd(slv->dm, loc_x, INSERT_VALUES, slv->x);CHKERRQ(ierr);
			ierr = stella_io_vtkwrite(slv->dm, slv->x, "phi", slv->grid.id, slv->ts);CHKERRQ(ierr);
		}
		ierr = DMRestoreLocalVector(slv->dm, &loc_x);CHKERRQ(ierr);
	} else if (slv->grid.nd == 3) {
		PetscScalar ***xvec;
		double ***phi;
		Vec loc_x;

		ierr = stella_dmap_get(slv->dmap, slv->state.phi, &phi);CHKERRQ(ierr);
		ierr = DMGetLocalVector(slv->dm, &loc_x);CHKERRQ(ierr);
		ierr = DMGlobalToLocalBegin(slv->dm, slv->x, INSERT_VALUES, loc_x);CHKERRQ(ierr);
		ierr = DMGlobalToLocalEnd(slv->dm, slv->x, INSERT_VALUES, loc_x);CHKERRQ(ierr);
		ierr = DMDAVecGetArray(slv->dm, loc_x, &xvec);CHKERRQ(ierr);

		for (k = zs; k < zs + zm; k++) {
			for (j = ys; j < ys + ym; j++) {
				for (i = xs; i < xs + xm; i++) {
					phi[k][j][i] = xvec[k][j][i];
				}
			}
		}

		if ((bz == DM_BOUNDARY_PERIODIC) && (zs + zm == ngz) && slv->grid.overlap_periodic) {
			for (j = ys; j < ys + ym; j++) {
				for (i = xs; i < xs + xm; i++) {
					phi[ngz][j][i] = xvec[ngz][j][i];
				}
			}
		}
		if ((by == DM_BOUNDARY_PERIODIC) && (ys + ym == ngy) && slv->grid.overlap_periodic) {
			for (k = zs; k < zs + zm; k++) {
				for (i = xs; i < xs + xm; i++) {
					phi[k][ngy][i] = xvec[k][ngy][i];
				}
			}
		}
		if ((bx == DM_BOUNDARY_PERIODIC) && (xs + xm == ngx) && slv->grid.overlap_periodic) {
			for (k = zs; k < zs + zm; k++) {
				for (j = ys; j < ys + ym; j++) {
					phi[k][j][ngx] = xvec[k][j][ngx];
				}
			}
		}

		ierr = stella_dmap_restore(slv->dmap, &phi);CHKERRQ(ierr);
		ierr = DMDAVecRestoreArray(slv->dm, loc_x, &xvec);CHKERRQ(ierr);
		if (stella_log(slv, STELLA_LOG_VTK)) {
			ierr = DMLocalToGlobalBegin(slv->dm, loc_x, INSERT_VALUES, slv->x);CHKERRQ(ierr);
			ierr = DMLocalToGlobalEnd(slv->dm, loc_x, INSERT_VALUES, slv->x);CHKERRQ(ierr);
			ierr = stella_io_vtkwrite(slv->dm, slv->x, "phi", slv->grid.id, slv->ts);CHKERRQ(ierr);
		}
		ierr = DMRestoreLocalVector(slv->dm, &loc_x);CHKERRQ(ierr);
	}

	if (stella_log(slv, STELLA_LOG_PROBLEM)) {
		PetscViewer vout;
		ierr = PetscViewerBinaryOpen(slv->comm, "sol.dat", FILE_MODE_WRITE, &vout);CHKERRQ(ierr);
		ierr = VecView(slv->x, vout);CHKERRQ(ierr);
		ierr = PetscViewerDestroy(&vout);CHKERRQ(ierr);
	}

	return 0;
}


PetscErrorCode stella_set_state(stella *slv, double phi[], double dcoef[],
                                double bcoef[], double jump[])
{
	PetscErrorCode ierr;

	slv->state.phi = phi;
	slv->state.dcoef = dcoef;
	slv->state.bcoef = bcoef;
	slv->state.jump = jump;
	slv->state.sol = NULL;

	if (slv->grid.nd == 2) {
		PetscScalar **dcoefvec;
		double **dcostella_state;
		PetscInt i, j;
		PetscInt xs, ys, xm, ym;
		ierr = DMDAGetCorners(slv->dm, &xs, &ys, 0, &xm, &ym, 0);CHKERRQ(ierr);
		ierr = stella_dmap_get(slv->dmap, slv->state.dcoef, &dcostella_state);CHKERRQ(ierr);
		ierr = DMDAVecGetArray(slv->dm, slv->level.dcoef, &dcoefvec);CHKERRQ(ierr);
		for (j = ys; j < ys + ym; j++) {
			for (i = xs; i < xs + xm; i++) {
				dcoefvec[j][i] = dcostella_state[j][i];
			}
		}
		ierr = stella_dmap_restore(slv->dmap, &dcostella_state);CHKERRQ(ierr);
		ierr = DMDAVecRestoreArray(slv->dm, slv->level.dcoef, &dcoefvec);CHKERRQ(ierr);
	} else if (slv->grid.nd == 3) {
		PetscScalar ***dcoefvec;
		double ***dcostella_state;
		PetscInt i, j, k;
		PetscInt xs, ys, zs, xm, ym, zm;
		ierr = DMDAGetCorners(slv->dm, &xs, &ys, &zs, &xm, &ym, &zm);CHKERRQ(ierr);
		ierr = stella_dmap_get(slv->dmap, slv->state.dcoef, &dcostella_state);CHKERRQ(ierr);
		ierr = DMDAVecGetArray(slv->dm, slv->level.dcoef, &dcoefvec);CHKERRQ(ierr);
		for (k = zs; k < zs + zm; k++) {
			for (j = ys; j < ys + ym; j++) {
				for (i = xs; i < xs + xm; i++) {
					dcoefvec[k][j][i] = dcostella_state[k][j][i];
				}
			}
		}
		ierr = stella_dmap_restore(slv->dmap, &dcostella_state);CHKERRQ(ierr);
		ierr = DMDAVecRestoreArray(slv->dm, slv->level.dcoef, &dcoefvec);CHKERRQ(ierr);

	}

	ierr = DMGlobalToLocalBegin(slv->dm, slv->level.dcoef, INSERT_VALUES, slv->level.ldcoef);CHKERRQ(ierr);
	ierr = DMGlobalToLocalEnd(slv->dm, slv->level.dcoef, INSERT_VALUES, slv->level.ldcoef);CHKERRQ(ierr);

  if (slv->grid.nd == 2) {
    PetscScalar **bcoefvec;
    double **bcostella_state;
    PetscInt i, j;
    PetscInt xs, ys, xm, ym;
    ierr = DMDAGetCorners(slv->dm, &xs, &ys, 0, &xm, &ym, 0);CHKERRQ(ierr);
    ierr = stella_dmap_get(slv->dmap, slv->state.bcoef, &bcostella_state);CHKERRQ(ierr);
    ierr = DMDAVecGetArray(slv->dm, slv->level.bcoef, &bcoefvec);CHKERRQ(ierr);
    for (j = ys; j < ys + ym; j++) {
        for (i = xs; i < xs + xm; i++) {
            bcoefvec[j][i] = bcostella_state[j][i];
        }
    }
    ierr = stella_dmap_restore(slv->dmap, &bcostella_state);CHKERRQ(ierr);
    ierr = DMDAVecRestoreArray(slv->dm, slv->level.bcoef, &bcoefvec);CHKERRQ(ierr);
  } else if (slv->grid.nd == 3) {
		PetscScalar ***bcoefvec;
		double ***bcostella_state;
		PetscInt i, j, k;
		PetscInt xs, ys, zs, xm, ym, zm;
		ierr = DMDAGetCorners(slv->dm, &xs, &ys, &zs, &xm, &ym, &zm);CHKERRQ(ierr);
		ierr = stella_dmap_get(slv->dmap, slv->state.bcoef, &bcostella_state);CHKERRQ(ierr);
		ierr = DMDAVecGetArray(slv->dm, slv->level.bcoef, &bcoefvec);CHKERRQ(ierr);
		for (k = zs; k < zs + zm; k++) {
			for (j = ys; j < ys + ym; j++) {
				for (i = xs; i < xs + xm; i++) {
					bcoefvec[k][j][i] = bcostella_state[k][j][i];
				}
			}
		}
		ierr = stella_dmap_restore(slv->dmap, &bcostella_state);CHKERRQ(ierr);
		ierr = DMDAVecRestoreArray(slv->dm, slv->level.bcoef, &bcoefvec);CHKERRQ(ierr);
  }

  ierr = DMGlobalToLocalBegin(slv->dm, slv->level.bcoef, INSERT_VALUES, slv->level.lbcoef);CHKERRQ(ierr);
  ierr = DMGlobalToLocalEnd(slv->dm, slv->level.bcoef, INSERT_VALUES, slv->level.lbcoef);CHKERRQ(ierr);

	return 0;
}


PetscErrorCode stella_set_sol(stella *slv, double sol[])
{
	PetscErrorCode ierr;

	ierr = DMCreateGlobalVector(slv->dm, &slv->sol);
	slv->state.sol = sol;

	if (slv->grid.nd == 2) {
		PetscScalar **svec;
		double **sol_state;
		PetscInt i, j;
		PetscInt xs, ys, xm, ym;
		ierr = DMDAGetCorners(slv->dm, &xs, &ys, 0, &xm, &ym, 0);CHKERRQ(ierr);
		ierr = stella_dmap_get(slv->dmap, slv->state.sol, &sol_state);CHKERRQ(ierr);
		ierr = DMDAVecGetArray(slv->dm, slv->sol, &svec);CHKERRQ(ierr);

		for (j = ys; j < ys + ym; j++) {
			for (i = xs; i < xs + xm; i++) {
				svec[j][i] = sol_state[j][i];
			}
		}
		ierr = stella_dmap_restore(slv->dmap, &sol_state);CHKERRQ(ierr);
		ierr = DMDAVecRestoreArray(slv->dm, slv->sol, &svec);CHKERRQ(ierr);
	} else if (slv->grid.nd == 3) {
		PetscScalar ***svec;
		double ***sol_state;
		PetscInt i, j, k;
		PetscInt xs, ys, zs, xm, ym, zm;
		ierr = DMDAGetCorners(slv->dm, &xs, &ys, &zs, &xm, &ym, &zm);CHKERRQ(ierr);
		ierr = stella_dmap_get(slv->dmap, slv->state.sol, &sol_state);CHKERRQ(ierr);
		ierr = DMDAVecGetArray(slv->dm, slv->sol, &svec);CHKERRQ(ierr);
		for (k = zs; k < zs + zm; k++) {
			for (j = ys; j < ys + ym; j++) {
				for (i = xs; i < xs + xm; i++) {
					svec[k][j][i] = sol_state[k][j][i];
				}
			}
		}
		ierr = stella_dmap_restore(slv->dmap, &sol_state);CHKERRQ(ierr);
		ierr = DMDAVecRestoreArray(slv->dm, slv->sol, &svec);CHKERRQ(ierr);
	}

	return 0;
}


PetscErrorCode stella_set_rhs(stella *slv, double rhs[])
{
	PetscErrorCode ierr;

	slv->state.rhs = rhs;
	ierr = stella_setup_rhs(slv);CHKERRQ(ierr);

	return ierr;
}


PetscErrorCode stella_set_grid(stella *slv, int is[], int ie[], int num_cells, double xyz[])
{
	PetscErrorCode ierr;
	PetscInt i, j, k;
	PetscInt xs, ys, zs, xm, ym, zm;
	DM cda;
	Vec gc;

	for (i = 0; i < slv->grid.nd; i++) {
		slv->grid.is[i] = is[i];
		slv->grid.ie[i] = ie[i];
	}
	slv->grid.num_cells = num_cells;
	slv->grid.xyz = xyz;

	ierr = DMGetCoordinateDM(slv->dm, &cda);CHKERRQ(ierr);
	ierr = DMGetCoordinates(slv->dm, &gc);CHKERRQ(ierr);
	ierr = DMDAGetCorners(cda, &xs, &ys, &zs, &xm, &ym, &zm);CHKERRQ(ierr);
	if (slv->grid.nd == 2) {
		DMDACoor2d **coors;
		double **x, **y;
		ierr = stella_dmap_get(slv->dmap, slv->grid.xyz, &x);CHKERRQ(ierr);
		ierr = stella_dmap_get(slv->dmap, slv->grid.xyz + num_cells, &y);CHKERRQ(ierr);
		ierr = DMDAVecGetArray(cda, gc, &coors);CHKERRQ(ierr);

		for (j = ys; j < ys + ym; j++) {
			for (i = xs; i < xs + xm; i++) {
				coors[j][i].x = x[j][i];
				coors[j][i].y = y[j][i];
			}
		}

		ierr = stella_dmap_restore(slv->dmap, &x);CHKERRQ(ierr);
		ierr = stella_dmap_restore(slv->dmap, &y);CHKERRQ(ierr);
		ierr = DMDAVecRestoreArray(cda, gc, &coors);CHKERRQ(ierr);
	} else if (slv->grid.nd == 3) {
		DMDACoor3d ***coors;
		double ***x, ***y, ***z;
		ierr = DMDAVecGetArray(cda, gc, &coors);CHKERRQ(ierr);
		ierr = stella_dmap_get(slv->dmap, slv->grid.xyz, &x);CHKERRQ(ierr);
		ierr = stella_dmap_get(slv->dmap, slv->grid.xyz + num_cells, &y);CHKERRQ(ierr);
		ierr = stella_dmap_get(slv->dmap, slv->grid.xyz + 2*num_cells, &z);CHKERRQ(ierr);

		for (k = zs; k < zs + zm; k++) {
			for (j = ys; j < ys + ym; j++) {
				for (i = xs; i < xs + xm; i++) {
					coors[k][j][i].x = x[k][j][i];
					coors[k][j][i].y = y[k][j][i];
					coors[k][j][i].z = z[k][j][i];
				}
			}
		}

		ierr = stella_dmap_restore(slv->dmap, &x);CHKERRQ(ierr);
		ierr = stella_dmap_restore(slv->dmap, &y);CHKERRQ(ierr);
		ierr = stella_dmap_restore(slv->dmap, &z);CHKERRQ(ierr);
		ierr = DMDAVecRestoreArray(cda, gc, &coors);CHKERRQ(ierr);

	}

	ierr = stella_metric_create(&slv->level.metric, slv->dm, slv->fd);CHKERRQ(ierr);

	if (stella_log(slv, STELLA_LOG_PROBLEM)) {
		PetscViewer vout;
		ierr = PetscViewerBinaryOpen(slv->comm, "grid.dat", FILE_MODE_WRITE, &vout);CHKERRQ(ierr);
		ierr = VecView(gc, vout);CHKERRQ(ierr);
		ierr = PetscViewerDestroy(&vout);CHKERRQ(ierr);
	}

	return 0;
}


PetscErrorCode stella_setup_op(stella *slv)
{
	PetscErrorCode ierr;
    Vec Vec_tmp;
    PetscScalar    **array;
    PetscInt         m,n;

	ierr = stella_operator_assemble(slv->op, slv->A, slv->dm);CHKERRQ(ierr);
	ierr = stella_boundary_apply(slv->boundary, slv->A, slv->dm);CHKERRQ(ierr);

	if (stella_log(slv, STELLA_LOG_STATUS)) {
		ierr = stella_io_print(slv->comm, "Matrix assembled");CHKERRQ(ierr);
	}

	if (slv->options.algebraic) {
		ierr = MatAssemblyBegin(slv->A, MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
		ierr = MatAssemblyEnd(slv->A, MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
		if (stella_log(slv, STELLA_LOG_PROBLEM)) {
			PetscViewer vout;
			ierr = PetscViewerBinaryOpen(slv->comm, "a.dat", FILE_MODE_WRITE, &vout);CHKERRQ(ierr);
			ierr = MatView(slv->A, vout);CHKERRQ(ierr);
			ierr = PetscViewerDestroy(&vout);CHKERRQ(ierr);
		}
	} else {
		#ifdef WITH_BOXMG
		if (stella_log(slv, STELLA_LOG_PROBLEM)) {
			stella_bmg_mat *ctx;
			ierr = MatShellGetContext(slv->A, (void**) &ctx);CHKERRQ(ierr);
			if (ctx->nd == 2)
				bmg2_operator_dump(ctx->op2);
			else
				bmg3_operator_dump(ctx->op3);
		}
		#endif
	}

	ierr = KSPSetOperators(slv->ksp, slv->A, slv->A);CHKERRQ(ierr);

	return 0;
}


int stella_log(stella *slv, stella_log_level level)
{
	return slv->options.log_level & level;
}


void stella_set_log(stella *slv, stella_log_level level)
{
	slv->options.log_level = level;
}


PetscErrorCode stella_cleanup(stella *slv)
{
	PetscErrorCode ierr;
	int i;

	ierr = stella_boundary_destroy(slv->boundary);CHKERRQ(ierr);
	ierr = stella_operator_destroy(slv->op);CHKERRQ(ierr); free(slv->op);
	ierr = stella_fd_destroy(slv->fd);CHKERRQ(ierr); free(slv->fd);
	stella_dmap_destroy(slv->dmap); free(slv->dmap);
	ierr = stella_metric_destroy(slv->level.metric);CHKERRQ(ierr);
	free(slv->level.metric);

	ierr = DMDestroy(&slv->dm);CHKERRQ(ierr);
	ierr = KSPDestroy(&slv->ksp);CHKERRQ(ierr);

	return 0;
}
