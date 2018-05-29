#include <stdlib.h>
#include <stdio.h>

#include "boundary.h"

#ifndef max
#define max(a,b) ((a) > (b) ? (a) : (b))
#endif

#ifndef min
#define min(a,b) (((a) < (b)) ? (a) : (b))
#endif

int contains(int *proc_is, int *proc_ie, int *hole_is, int *hole_ie, int nd)
{
	int flag = 1;
	int i;
	for (i = 0; i < nd; i++) {
		flag = flag && (hole_is[i] < proc_ie[i]);
		flag = flag && (hole_ie[i] >= proc_is[i]);
	}

	return flag;
}


blist *boundary_create(grid *grd, problem *pb)
{
	int i, j;

	blist *lst = (blist*) malloc(sizeof(blist));
	lst->len = 4;
	if (pb->id == ELECTRODE) lst->len += pb->nholes;
	if (grd->nd == 3) lst->len += 2;
	lst->v = (boundary*) malloc(lst->len*sizeof(boundary));

	lst->len = 0;
	if (grd->nd == 3) {
		for (i = 0; i < 6; i++) {
			if (i == NORTH
			    && grd->cart_coord[1] == grd->num_procs[1]-1
			    && (!grd->periodic[1])) {
				lst->v[lst->len].dirichlet = (double*) malloc(grd->num_pts*sizeof(double));
				int ii,kk;
				for (kk = 0; kk < grd->nz; kk++) {
					for (ii = 0; ii < grd->nx; ii++) {
						int ind = kk*grd->nx*grd->ny+(grd->ny-1)*grd->nx+ii;
						lst->v[lst->len].dirichlet[kk*grd->nx+ii] =
							pb->sol(grd->x[ind], grd->y[ind], grd->z[ind]);
					}
				}

				lst->v[lst->len].type = pb->boundary[i];

				lst->v[lst->len].is[0] = grd->is[0];
				lst->v[lst->len].ie[0] = grd->ie[0];
				lst->v[lst->len].is[2] = grd->is[2];
				lst->v[lst->len].ie[2] = grd->ie[2];
				lst->v[lst->len].is[1] = grd->ie[1];
				lst->v[lst->len].ie[1] = grd->ie[1];
				lst->v[lst->len].norm_dir = -2;

				lst->len++;
			} else if (i == SOUTH
			           && grd->cart_coord[1] == 0
			           && (!grd->periodic[1])) {
				lst->v[lst->len].dirichlet = (double*) malloc(grd->num_pts*sizeof(double));
				int ii,kk;
				for (kk = 0; kk < grd->nz; kk++) {
					for (ii = 0; ii < grd->nx; ii++) {
						int ind = kk*grd->nx*grd->ny + ii;
						lst->v[lst->len].dirichlet[kk*grd->nx+ii] =
							pb->sol(grd->x[ind], grd->y[ind], grd->z[ind]);
					}
				}

				lst->v[lst->len].type = pb->boundary[i];

				lst->v[lst->len].is[0] = grd->is[0];
				lst->v[lst->len].ie[0] = grd->ie[0];
				lst->v[lst->len].is[2] = grd->is[2];
				lst->v[lst->len].ie[2] = grd->ie[2];
				lst->v[lst->len].is[1] = grd->is[1];
				lst->v[lst->len].ie[1] = grd->is[1];
				lst->v[lst->len].norm_dir = 2;

				lst->len++;
			} else if (i == EAST
			           && grd->cart_coord[0] == grd->num_procs[0]-1
			           && (!grd->periodic[0])) {
				lst->v[lst->len].dirichlet = (double*) malloc(grd->num_pts*sizeof(double));

				int jj,kk;
				for (kk = 0; kk < grd->nz; kk++) {
					for (jj = 0; jj < grd->ny; jj++) {
						int ind = kk*grd->nx*grd->ny + jj*grd->nx + grd->nx-1;
						lst->v[lst->len].dirichlet[kk*grd->ny+jj] =
							pb->sol(grd->x[ind], grd->y[ind], grd->z[ind]);
					}
				}

				lst->v[lst->len].type = pb->boundary[i];

				lst->v[lst->len].is[0] = grd->ie[0];
				lst->v[lst->len].ie[0] = grd->ie[0];
				lst->v[lst->len].is[1] = grd->is[1];
				lst->v[lst->len].ie[1] = grd->ie[1];
				lst->v[lst->len].is[2] = grd->is[2];
				lst->v[lst->len].ie[2] = grd->ie[2];
				lst->v[lst->len].norm_dir = -1;

				lst->len++;
			} else if (i == WEST
			           && grd->cart_coord[0] == 0
			           && (!grd->periodic[0])) {
				lst->v[lst->len].dirichlet = (double*) malloc(grd->num_pts*sizeof(double));

				int jj, kk;
				for (kk = 0; kk < grd->nz; kk++) {
					for (jj = 0; jj < grd->ny; jj++) {
						int ind = kk*grd->nx*grd->ny + jj*grd->nx;
						lst->v[lst->len].dirichlet[kk*grd->ny+jj] =
							pb->sol(grd->x[ind], grd->y[ind], grd->z[ind]);
					}
				}

				lst->v[lst->len].type = pb->boundary[i];

				lst->v[lst->len].is[0] = grd->is[0];
				lst->v[lst->len].ie[0] = grd->is[0];
				lst->v[lst->len].is[1] = grd->is[1];
				lst->v[lst->len].ie[1] = grd->ie[1];
				lst->v[lst->len].is[2] = grd->is[2];
				lst->v[lst->len].ie[2] = grd->ie[2];
				lst->v[lst->len].norm_dir = 1;

				lst->len++;
			} else if (i == FRONT
			           && grd->cart_coord[2] == grd->num_procs[2] - 1
			           && (!grd->periodic[2])) {
				lst->v[lst->len].dirichlet = (double*) malloc(grd->num_pts*sizeof(double));
				int ii, jj;
				for (jj = 0; jj < grd->ny; jj++) {
					for (ii = 0; ii < grd->nx; ii++) {
						int ind = (grd->nz-1)*grd->nx*grd->ny + jj*grd->nx + ii;
						lst->v[lst->len].dirichlet[jj*grd->nx + ii] =
							pb->sol(grd->x[ind], grd->y[ind], grd->z[ind]);
					}
				}

				lst->v[lst->len].type = pb->boundary[i];

				lst->v[lst->len].is[0] = grd->is[0];
				lst->v[lst->len].ie[0] = grd->ie[0];
				lst->v[lst->len].is[1] = grd->is[1];
				lst->v[lst->len].ie[1] = grd->ie[1];
				lst->v[lst->len].is[2] = grd->ie[2];
				lst->v[lst->len].ie[2] = grd->ie[2];
				lst->v[lst->len].norm_dir = -3;

				lst->len++;
			} else if (i == BACK
			           && grd->cart_coord[2] == 0
			           && (!grd->periodic[2])) {
				lst->v[lst->len].dirichlet = (double*) malloc(grd->num_pts*sizeof(double));
				int ii, jj;
				for (jj = 0; jj < grd->ny; jj++) {
					for (ii = 0; ii < grd->nx; ii++) {
						int ind = jj*grd->nx + ii;
						lst->v[lst->len].dirichlet[jj*grd->nx + ii] =
							pb->sol(grd->x[ind], grd->y[ind], grd->z[ind]);
					}
				}

				lst->v[lst->len].type = pb->boundary[i];

				lst->v[lst->len].is[0] = grd->is[0];
				lst->v[lst->len].ie[0] = grd->ie[0];
				lst->v[lst->len].is[1] = grd->is[1];
				lst->v[lst->len].ie[1] = grd->ie[1];
				lst->v[lst->len].is[2] = grd->is[2];
				lst->v[lst->len].ie[2] = grd->is[2];
				lst->v[lst->len].norm_dir = 3;

				lst->len++;
			}
		}
	} else {
		for (i = 0; i < 4; i++) {
			if (i == NORTH
			    && grd->cart_coord[1] == grd->num_procs[1]-1
			    && (!grd->periodic[1])) {
				lst->v[lst->len].dirichlet = (double*) malloc(grd->num_pts
				                                              *sizeof(double));
				for (j = 0; j < grd->nx; j++) {
					lst->v[lst->len].dirichlet[j] = pb->sol(
						grd->x[(grd->ny-1)*grd->nx + j],
						grd->y[(grd->ny-1)*grd->nx + j], 0);
				}

				lst->v[lst->len].type = pb->boundary[i];

				lst->v[lst->len].is[0] = grd->is[0];
				lst->v[lst->len].ie[0] = grd->ie[0];
				lst->v[lst->len].is[1] = grd->ie[1];
				lst->v[lst->len].ie[1] = grd->ie[1];
				lst->v[lst->len].is[2] = 0;
				lst->v[lst->len].ie[2] = 0;
				lst->v[lst->len].norm_dir = -2;

				lst->len++;
			} else if (i == SOUTH
			           && grd->cart_coord[1] == 0
			           && (!grd->periodic[1])) {
				lst->v[lst->len].dirichlet = (double*) malloc(grd->num_pts*sizeof(double));

				for (j = 0; j < grd->nx; j++)
					lst->v[lst->len].dirichlet[j] = pb->sol(grd->x[j], grd->y[j], 0);

				lst->v[lst->len].type = pb->boundary[i];

				lst->v[lst->len].is[0] = grd->is[0];
				lst->v[lst->len].ie[0] = grd->ie[0];
				lst->v[lst->len].is[1] = grd->is[1];
				lst->v[lst->len].ie[1] = grd->is[1];
				lst->v[lst->len].is[2] = 0;
				lst->v[lst->len].ie[2] = 0;
				lst->v[lst->len].norm_dir = 2;

				lst->len++;
			} else if (i == EAST
			           && grd->cart_coord[0] == grd->num_procs[0]-1
			           && (!grd->periodic[0])) {
				lst->v[lst->len].dirichlet = (double*) malloc(grd->num_pts*sizeof(double));

				for (j = 0; j < grd->ny; j++)
					lst->v[lst->len].dirichlet[j] = pb->sol(grd->x[j*grd->nx + grd->nx - 1], grd->y[j*grd->nx + grd->nx - 1], 0);

				lst->v[lst->len].type = pb->boundary[i];

				lst->v[lst->len].is[0] = grd->ie[0];
				lst->v[lst->len].ie[0] = grd->ie[0];
				lst->v[lst->len].is[1] = grd->is[1];
				lst->v[lst->len].ie[1] = grd->ie[1];
				lst->v[lst->len].is[2] = 0;
				lst->v[lst->len].ie[2] = 0;
				lst->v[lst->len].norm_dir = -1;

				lst->len++;
			} else if (i == WEST
			           && grd->cart_coord[0] == 0
			           && (!grd->periodic[0])) {
				lst->v[lst->len].dirichlet = (double*) malloc(grd->num_pts*sizeof(double));
				for (j = 0; j < grd->ny; j++)
					lst->v[lst->len].dirichlet[j] = pb->sol(grd->x[j*grd->nx], grd->y[j*grd->nx], 0);

				lst->v[lst->len].type = pb->boundary[i];

				lst->v[lst->len].is[0] = grd->is[0];
				lst->v[lst->len].ie[0] = grd->is[0];
				lst->v[lst->len].is[1] = grd->is[1];
				lst->v[lst->len].ie[1] = grd->ie[1];
				lst->v[lst->len].is[2] = 0;
				lst->v[lst->len].ie[2] = 0;
				lst->v[lst->len].norm_dir = 1;

				lst->len++;
			}
		}
	}

	if (pb->id == ELECTRODE) {
		for (i = 0; i < pb->nholes; ++i) {
			boundary *bnd = lst->v + lst->len;
			int el_is[3];
			int el_ie[3];

			bnd->type = DIRICHLET;

			int k;
			for (k = 0; k < grd->nd; k++) {
				double esize = pb->holes[i].rel_size[k] * grd->num_global[k];
				double offset = pb->holes[i].rel_offset[k] * grd->num_global[k];
				el_is[k] = offset - esize * 0.5;
				el_ie[k] = offset + esize * 0.5;
			}

			if (contains(grd->is, grd->ie,
			             el_is, el_ie, grd->nd)) {
				int ne[3];
				ne[2] = 1;
				for (k = 0; k < grd->nd; k++) {
					bnd->is[k] = max(el_is[k], grd->is[k]);
					bnd->ie[k] = min(el_ie[k], grd->ie[k]);
					ne[k] = bnd->ie[k] - bnd->is[k] + 1;
				}

				printf("Electrode: (%d %d %d) -> (%d %d %d)\n",
				       bnd->is[0], bnd->is[1], bnd->is[2],
				       bnd->ie[0], bnd->ie[1], bnd->ie[2]);

				bnd->dirichlet = (double*) malloc(ne[0]*ne[1]*ne[2]*sizeof(double));

				if (grd->nd == 2) {
					int ii,jj;
					for (jj = 0; jj < ne[1]; jj++) {
						for (ii = 0; ii < ne[0]; ii++) {
							int ind = (bnd->is[1] - grd->is[1] + jj)*grd->nx + (bnd->is[0] - grd->is[0] + ii);
							bnd->dirichlet[jj*ne[0] + ii] = pb->sol(grd->x[ind], grd->y[ind], 0);
							/* if (i-4) */
							/* 	bnd->dirichlet[jj*nbx + ii] = 1.0; */
							/* else */
							/* 	bnd->dirichlet[jj*nbx + ii] = -1.0; */
						}
					}
				} else {
					int ii,jj,kk;
					for (kk = 0; kk < ne[2]; kk++) {
						for (jj = 0; jj < ne[1]; jj++) {
							for (ii = 0; ii < ne[0]; ii++) {
								int ind = (((bnd->is[2] - grd->is[2] + kk)*grd->nx*grd->ny) +
								           ((bnd->is[1] - grd->is[1] + jj)*grd->nx) +
								           ((bnd->is[0] - grd->is[0] + ii)));
								bnd->dirichlet[kk*ne[0]*ne[1] + jj*ne[0] + ii] =
									pb->sol(grd->x[ind], grd->y[ind], grd->z[ind]);
							}
						}
					}
				}

				lst->v[lst->len].norm_dir = 1;

				lst->len++;
			}
		}
	}


	lst->v = realloc(lst->v, lst->len*sizeof(boundary));

	return lst;
}


void boundary_destroy(blist *lst)
{
	int i;

	for (i = 0; i < lst->len; i++) {
		free(lst->v[i].dirichlet);
	}

	free(lst->v);
}
