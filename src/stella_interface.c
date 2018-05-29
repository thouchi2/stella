#include <stdio.h>
#include <stdlib.h>

#include "stella_solver.h"
#include "stella_io.h"


void stella_set_patch(void **solver, int is[], int ie[], int *bc_type,
                   int *norm_dir, double **dirichlet)

{
	efs *slv = (efs*) *solver;

	stella_boundary_add(slv->boundary, (stella_bctype) *bc_type, *norm_dir,
	                is, ie, *dirichlet);
}


void stella_init(void **solverCtx, MPI_Fint *fcomm,
             int nGlobal[], int nProcs[], int nLocal[],
             int offset[], int stride[],
             int cartCoord[], int periodic[], int periodic_storage,
             int *nd, int *ng, int axisymmetric, PetscErrorCode *ierr)
{
	efs *slv;

	MPI_Comm *comm = malloc(sizeof(MPI_Comm));
	#ifdef C_CLIENT
	*comm = *fcomm;
	#else
	*comm = MPI_Comm_f2c(*fcomm);
	#endif
	*ierr = efs_create(&slv, *comm);
	slv->options.axisymmetric = axisymmetric;
	//efs_set_log(slv, EFS_LOG_STATUS | EFS_LOG_RESIDUAL | EFS_LOG_PROBLEM);
	efs_set_log(slv, EFS_LOG_STATUS);
	*ierr = stella_grid_setup(&slv->grid, &slv->dm, &slv->comm, nGlobal, nProcs,
	                      nLocal, cartCoord, periodic, periodic_storage, *nd);
	slv->level.dm = slv->dm;
	slv->grid.id = *ng;
	*ierr = efs_setup(slv, offset, stride);

	*solverCtx = slv;
}


void stella_set_state(void **solver, double **phi, double **dcoef,
                  double **bcoef, double **jump, PetscErrorCode *ierr)
{
	efs *slv = (efs*) *solver;

	*ierr = efs_set_state(slv, *phi, *dcoef, *bcoef, *jump);
}


void stella_set_rhs(void **solver, double **rhs, PetscErrorCode *ierr)
{
	efs *slv = (efs*) *solver;

	*ierr = efs_set_rhs(slv, *rhs);
}


void stella_set_sol(void **solver, double sol[], PetscErrorCode *ierr)
{
	efs *slv = (efs*) *solver;

	*ierr = efs_set_sol(slv, sol);
}


void stella_set_grid(void **solver, int is[], int ie[], int *nCells, double **xyz,
                  PetscErrorCode *ierr)
{
	efs *slv = (efs*) *solver;

	*ierr = efs_set_grid(slv, is, ie, *nCells, *xyz);
}


void stella_setup_op(void **solver, PetscErrorCode *ierr)
{
	efs *slv = (efs*) *solver;

	*ierr = efs_setup_op(slv);
}


void stella_solve(void **solver, PetscErrorCode *ierr)
{
	efs *slv = (efs*) *solver;

	*ierr = efs_solve(slv);
}


void stella_cleanup(void **solver, PetscErrorCode *ierr)
{
	efs *slv = (efs*) *solver;

	*ierr = efs_destroy(slv);
	free(slv);
}


void stella_plot_sol(void **solver, int *idx, PetscErrorCode *ierr)
{
	efs *slv = (efs*) *solver;

	*ierr = stella_io_vtkwrite(slv->dm, slv->sol, "solution", slv->grid.id, *idx);
}
