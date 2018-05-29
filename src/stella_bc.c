#include "stella_bc.h"
#include "stella_neumann.h"
#include "stella_dirichlet.h"


PetscErrorCode stella_bc_create(stella_bc **efbc, stella_bctype btype, int norm_dir, int is[], int ie[],
                                double *dirichlet, stella_level *level, stella_dmap *dmap,
                                stella_state *state, stella_fd *fd)
{
	PetscErrorCode ierr;
	int i;

	stella_bc *bc = (stella_bc*) malloc(sizeof(stella_bc));
	bc->patch = (stella_patch*) malloc(sizeof(stella_patch));
	for (i = 0; i < 3; i++) {
		bc->patch->corners.is[i] = is[i] - 1;
		bc->patch->corners.ie[i] = ie[i] - 1;
		PetscInt ngx, ngy, ngz;
		ierr = DMDAGetInfo(level->dm, 0, &ngx, &ngy, &ngz,0,0,0,0,0,0,0,0,0);CHKERRQ(ierr);
		if (i == 0 && bc->patch->corners.ie[i] >= ngx) {
			bc->patch->corners.ie[i]--;
		}
		if (i == 1 && bc->patch->corners.ie[i] >= ngy) {
			bc->patch->corners.ie[i]--;
		}
		if (i == 2 && bc->patch->corners.ie[i] >= ngz) {
			bc->patch->corners.ie[i]--;
		}
	}
	bc->patch->bc_type = btype;
	bc->patch->norm_dir = norm_dir;
	bc->patch->dirichlet = dirichlet;
  for (i = 0; i < 3; i++) {
    bc->patch->offset[i] = 0;
  }
  for (i = 0; i < 3; i++) {
    bc->patch->stride[i] = ie[i] - is[i] + 1;
  }
	bc->level = level;
	bc->slv_dmap = dmap;
	bc->state = state;
	bc->fd = fd;
	bc->axisymmetric = 0;
	bc->destroy = NULL;
	bc->symmetric = NULL;
	bc->apply = NULL;
	bc->apply_rhs = NULL;
	bc->btype = btype;

	if (btype == STELLA_DIRICHLET_G || btype == STELLA_DIRICHLET_P || btype == STELLA_SCHWARZ) {
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
	free(bc->patch);

	return 0;
}
