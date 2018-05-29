#include <stdlib.h>

#include "stella_boundary.h"

#define min(a,b) ((a) <= (b) ? (a) : (b))
#define max(a,b) ((a) >= (b) ? (a) : (b))

static int overlaps(stella_patch *bc1, stella_patch *bc2) {

  int both_neumann;
  int same_norm_dir;
  int result;

  both_neumann = bc1->bc_type == STELLA_NEUMANN && bc2->bc_type == STELLA_NEUMANN;
  same_norm_dir = bc1->norm_dir == bc2->norm_dir;

  result = (!both_neumann || same_norm_dir) &&
    bc1->corners.ie[0] >= bc2->corners.is[0] &&
    bc1->corners.is[0] <= bc2->corners.ie[0] &&
    bc1->corners.ie[1] >= bc2->corners.is[1] &&
    bc1->corners.is[1] <= bc2->corners.ie[1] &&
    bc1->corners.ie[2] >= bc2->corners.is[2] &&
    bc1->corners.is[2] <= bc2->corners.ie[2];

  return result;

}


static void intersect(stella_patch *bc1, stella_patch *bc2, stella_corner *intersection) {

  intersection->is[0] = max(bc1->corners.is[0], bc2->corners.is[0]);
  intersection->is[1] = max(bc1->corners.is[1], bc2->corners.is[1]);
  intersection->is[2] = max(bc1->corners.is[2], bc2->corners.is[2]);
  intersection->ie[0] = min(bc1->corners.ie[0], bc2->corners.ie[0]);
  intersection->ie[1] = min(bc1->corners.ie[1], bc2->corners.ie[1]);
  intersection->ie[2] = min(bc1->corners.ie[2], bc2->corners.ie[2]);

}


static int points_removed(int nd, stella_patch *patch, stella_corner *intersection, int dir) {

  int result;
  int otherdir1, otherdir2;
  int i;
  int N[3];
  int cut_lower, cut_upper;
  int removed_in_dir;

  if (nd == 2) {
    otherdir1 = (dir + 1) % 2;
  } else {
    otherdir1 = (dir + 1) % 3;
    otherdir2 = (dir + 2) % 3;
  }

  for (i = 0; i < nd; ++i) {
    N[i] = patch->corners.ie[i] - patch->corners.is[i] + 1;
  }

  cut_lower = !(patch->corners.is[dir] < intersection->is[dir] &&
    patch->corners.ie[dir] == intersection->ie[dir]);
  cut_upper = !(patch->corners.is[dir] == intersection->is[dir] &&
    patch->corners.ie[dir] > intersection->ie[dir]);

  removed_in_dir = 0;
  if (cut_lower) {
    removed_in_dir += intersection->ie[dir] - patch->corners.is[dir] + 1;
  }
  if (cut_upper) {
    removed_in_dir += patch->corners.ie[dir] - intersection->is[dir] + 1;
  }
  removed_in_dir = min(removed_in_dir, N[dir]);

  if (nd == 2) {
    result = removed_in_dir * N[otherdir1];
  } else {
    result = removed_in_dir * N[otherdir1] * N[otherdir2];
  }

  return result;

}


int bc_priority(stella_patch *bc) {

  int result;

  switch (bc->bc_type) {
  case (STELLA_NEUMANN):
    result = 0;
    break;
  case (STELLA_DIRICHLET_P):
    result = 1;
    break;
  case (STELLA_DIRICHLET_G):
    result = 1;
    break;
  case (STELLA_SCHWARZ):
    result = 2;
    break;
  default:
    result = INT_MAX;
    break;
  }

  return result;

}


static void fix_corner(int nd, stella_patch *bc1, stella_patch *bc2) {

  stella_corner intersection;
  int i;
  int cost;
  int min_cost;
  stella_patch *patch_to_cut;
  int cut_dir;
  int cut_lower, cut_upper;

  intersect(bc1, bc2, &intersection);

  min_cost = INT_MAX;
  patch_to_cut = NULL;
  cut_dir = 0;

  if (bc_priority(bc1) <= bc_priority(bc2)) {
    for (i = 0; i < nd; ++i) {
      cost = points_removed(nd, bc1, &intersection, i);
      if (cost < min_cost) {
        min_cost = cost;
        patch_to_cut = bc1;
        cut_dir = i;
      }
    }
  }

  if (bc_priority(bc1) >= bc_priority(bc2)) {
    for (i = 0; i < nd; ++i) {
      cost = points_removed(nd, bc2, &intersection, i);
      if (cost < min_cost) {
        min_cost = cost;
        patch_to_cut = bc2;
        cut_dir = i;
      }
    }
  }

  i = cut_dir;

  cut_lower = !(patch_to_cut->corners.is[i] < intersection.is[i] &&
    patch_to_cut->corners.ie[i] == intersection.ie[i]);
  cut_upper = !(patch_to_cut->corners.is[i] == intersection.is[i] &&
    patch_to_cut->corners.ie[i] > intersection.ie[i]);

  if (cut_lower) {
    patch_to_cut->offset[i] = intersection.ie[i]+1 - patch_to_cut->corners.is[i];
    patch_to_cut->corners.is[i] = intersection.ie[i]+1;
  }
  if (cut_upper) {
    patch_to_cut->corners.ie[i] = intersection.is[i]-1;
  }

}


static void handle_corner_overlap(stella_boundary *bnd, stella_bc *bc)
{
	int i;
	for (i = 0; i < bnd->num_bc; i++) {
		if (overlaps(bc->patch, bnd->bc[i]->patch)) {
      fix_corner(bnd->slv_dmap->nd, bc->patch, bnd->bc[i]->patch);
		}
	}
}


PetscErrorCode stella_boundary_create(stella_boundary **stella_bnd, stella_level *level,
                                  stella_dmap *dmap, stella_state *state, stella_fd *fd)
{
	stella_boundary *bnd = (stella_boundary*) malloc(sizeof(stella_boundary));

	bnd->bc = NULL;
	bnd->level = level;
	bnd->slv_dmap = dmap;
	bnd->state = state;
	bnd->fd = fd;
	bnd->num_bc = 0;
	bnd->axisymmetric = 0;

	*stella_bnd = bnd;
	return 0;
}


PetscErrorCode stella_boundary_apply(stella_boundary *bnd, Mat A, DM da)
{
	PetscErrorCode ierr;
	MatType mtype;
	PetscBool is_boxmg;
	int i;

	ierr = MatGetType(A, &mtype);CHKERRQ(ierr);
	ierr = PetscStrcmp(mtype, MATSHELL, &is_boxmg);CHKERRQ(ierr);

	ierr = DMGlobalToLocalBegin(da, bnd->level->add_cont, INSERT_VALUES, bnd->level->ladd_cont);CHKERRQ(ierr);
	ierr = DMGlobalToLocalEnd(da, bnd->level->add_cont, INSERT_VALUES, bnd->level->ladd_cont);CHKERRQ(ierr);

	for (i = 0; i < bnd->num_bc; i++) {
		if (bnd->bc[i]->patch->bc_type == STELLA_NEUMANN) {
			ierr = stella_bc_apply(bnd->bc[i], A, da);CHKERRQ(ierr);
		}
	}
	for (i = 0; i < bnd->num_bc; i++) {
		if (bnd->bc[i]->patch->bc_type == STELLA_DIRICHLET_G ||
		    bnd->bc[i]->patch->bc_type == STELLA_DIRICHLET_P ||
		    bnd->bc[i]->patch->bc_type == STELLA_SCHWARZ) {
			ierr = stella_bc_apply(bnd->bc[i], A, da);CHKERRQ(ierr);
		}
	}

	if (!is_boxmg) {
		for (i = 0; i < bnd->num_bc; i++) {
			ierr = stella_bc_symmetric(bnd->bc[i], A, da);CHKERRQ(ierr);
		}
	}

	ierr = DMLocalToGlobalBegin(da, bnd->level->ladd_cont, ADD_VALUES, bnd->level->add_cont);CHKERRQ(ierr);
	ierr = DMLocalToGlobalEnd(da, bnd->level->ladd_cont, ADD_VALUES, bnd->level->add_cont);CHKERRQ(ierr);

	return 0;
}


PetscErrorCode stella_boundary_apply_rhs(stella_boundary *bnd, DM da, Vec rhs)
{
	PetscErrorCode ierr;
	int i;

	for (i = 0; i < bnd->num_bc; i++) {
		ierr = stella_bc_apply_rhs(bnd->bc[i], da, rhs);CHKERRQ(ierr);
	}

	return 0;
}


PetscErrorCode stella_boundary_add(stella_boundary *bnd, stella_bctype btype, int norm_dir,
                               int is[], int ie[], double *dirichlet)
{
	PetscErrorCode ierr;
	int ind = bnd->num_bc;

	if ((bnd->num_bc % 10) == 0) {
		bnd->bc = realloc(bnd->bc, (bnd->num_bc + 10) * sizeof(stella_bc*));
	}

	ierr = stella_bc_create(&bnd->bc[ind], btype, norm_dir, is, ie, dirichlet,
	                    bnd->level, bnd->slv_dmap, bnd->state, bnd->fd);CHKERRQ(ierr);
	bnd->bc[ind]->axisymmetric = bnd->axisymmetric;
  handle_corner_overlap(bnd, bnd->bc[ind]);
	bnd->num_bc++;

	return 0;
}


PetscErrorCode stella_boundary_destroy(stella_boundary *bnd)
{
	PetscErrorCode ierr;
	int i;

	for (i = 0; i < bnd->num_bc; i++) {
		ierr = stella_bc_destroy(bnd->bc[i]);CHKERRQ(ierr);
	}

	return 0;
}
