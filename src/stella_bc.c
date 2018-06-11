#include "stella_bc.h"
#include "stella_neumann.h"
#include "stella_dirichlet.h"


PetscErrorCode stella_bc_create(stella_bc **efbc, stella_bctype btype, char *norm_dir, double *values,
                                stella_level *level, stella_dmap *dmap,
                                stella_state *state, stella_fd *fd)
{
	PetscErrorCode ierr;
	int i;

	stella_bc *bc = (stella_bc*) malloc(sizeof(stella_bc));

	bc->level = level;
	bc->slv_dmap = dmap;
	bc->state = state;
	bc->fd = fd;
	bc->norm_dirs = norm_dir;
	bc->values = values;
	bc->axisymmetric = 0;
	bc->destroy = NULL;
	bc->symmetric = NULL;
	bc->apply = NULL;
	bc->apply_rhs = NULL;
	bc->btype = btype;

	if (btype == STELLA_DIRICHLET || btype == STELLA_SCHWARZ) {
		ierr = stella_dirichlet_create(bc);CHKERRQ(ierr);
	} else if (btype == STELLA_NEUMANN) {
		ierr = stella_neumann_create(bc);CHKERRQ(ierr);
	}

	*efbc = bc;

	return 0;
}


PetscErrorCode stella_bc_symmetric(stella_bc *bc, Mat A, DM da)
{
	PetscErrorCode ierr;

	if (bc->symmetric) {
		ierr = bc->symmetric(bc, A, da);CHKERRQ(ierr);
	}

	return 0;
}


PetscErrorCode stella_bc_apply(stella_bc *bc, Mat A, DM da)
{
	PetscErrorCode ierr;

	if (bc->apply) {
		ierr = bc->apply(bc, A, da);CHKERRQ(ierr);
	}

	return 0;
}


PetscErrorCode stella_bc_apply_rhs(stella_bc *bc, DM da, Vec rhs)
{
	PetscErrorCode ierr;

	if (bc->apply_rhs) {
		ierr = bc->apply_rhs(bc, da, rhs);CHKERRQ(ierr);
	}

	return 0;
}


PetscErrorCode stella_bc_destroy(stella_bc *bc)
{
	PetscErrorCode ierr;

	if (bc->destroy) {
		ierr = bc->destroy(bc);CHKERRQ(ierr);
	}

	return 0;
}
