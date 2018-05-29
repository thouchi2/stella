#ifndef STELLA_FD_H
#define STELLA_FD_H

#include <petscsys.h>

#include "stella_stencil.h"

/**
 * file: stella_fd.h
 *
 * Class that is used to hold finite difference
 * coefficients on a grid with unit spacing
 */

typedef const char* stella_fd_type;

#define STELLA_FD_STANDARD_O2 "standard second order"

/**
 * Contains finite difference stencil
 * coefficients.
 */
typedef struct {
	stella_stencil **first; /**< first derivative finite difference stencil */
	stella_stencil **second;/**< second derivative finite difference stencil */
	stella_fd_type type;
} stella_fd;


PetscErrorCode stella_fd_create(stella_fd **fd, stella_fd_type type);


PetscErrorCode stella_fd_destroy(stella_fd *fd);


#endif
