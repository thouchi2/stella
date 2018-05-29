#include <stdlib.h>

#include "stella_fd.h"


static PetscErrorCode set_second_order(stella_fd *fd)
{
	PetscErrorCode ierr;

	ierr = stella_stencil_create(&fd->first[-1], 3, 0);CHKERRQ(ierr);
	ierr = stella_stencil_create(&fd->second[-1], 4, 0);CHKERRQ(ierr);

	ierr = stella_stencil_create(&fd->first[0], 3, 1);CHKERRQ(ierr);
	ierr = stella_stencil_create(&fd->second[0], 3, 1);CHKERRQ(ierr);

	ierr = stella_stencil_create(&fd->first[1], 3, 2);CHKERRQ(ierr);
	ierr = stella_stencil_create(&fd->second[1], 4, 3);CHKERRQ(ierr);

	fd->first[0]->v[-1] = -0.5;
	fd->first[0]->v[0] = 0;
	fd->first[0]->v[1] = 0.5;

	fd->first[-1]->v[0] = -3*.5;
	fd->first[-1]->v[1] = 4*.5;
	fd->first[-1]->v[2] = -.5;

	fd->first[1]->v[0] = 3*.5;
	fd->first[1]->v[-1] = -4*.5;
	fd->first[1]->v[-2] = .5;

	fd->second[0]->v[-1] = 1;
	fd->second[0]->v[0] = -2;
	fd->second[0]->v[1] = 1;

	fd->second[-1]->v[0] = 2;
	fd->second[-1]->v[1] = -5;
	fd->second[-1]->v[2] = 4;
	fd->second[-1]->v[3] = -1;

	fd->second[1]->v[0] = 2;
	fd->second[1]->v[-1] = -5;
	fd->second[1]->v[-2] = 4;
	fd->second[1]->v[-3] = -1;

	return 0;
}


PetscErrorCode stella_fd_create(stella_fd **effd, stella_fd_type type)
{
	PetscErrorCode ierr;

	stella_fd *fd = (stella_fd*) malloc(sizeof(stella_fd));

	fd->first = (stella_stencil**) malloc(3*sizeof(stella_stencil*));
	fd->second = (stella_stencil**) malloc(3*sizeof(stella_stencil*));
	fd->first++;
	fd->second++;

	if (strcmp(type, STELLA_FD_STANDARD_O2) == 0) {
		ierr = set_second_order(fd);CHKERRQ(ierr);
	}

	*effd = fd;
	return 0;
}


PetscErrorCode stella_fd_destroy(stella_fd *fd)
{
	PetscErrorCode ierr;
	int i;

	fd->first--;
	fd->second--;

	for(i = 0; i < 3; i++) {
		ierr = stella_stencil_destroy(fd->first[i]);CHKERRQ(ierr);
		free(fd->first[i]);
		ierr = stella_stencil_destroy(fd->second[i]);CHKERRQ(ierr);
		free(fd->second[i]);
	}

	free(fd->first);
	free(fd->second);

	return 0;
}
