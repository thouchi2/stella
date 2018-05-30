#include "stella_mat.h"


PetscErrorCode stella_bmg_SetValuesStencil(Mat mat, PetscInt m, const MatStencil idxm[], PetscInt n,
                                           const MatStencil idxn[], const PetscScalar v[],
                                           InsertMode addv)
{
	#ifdef WITH_BOXMG
	PetscInt ierr;
	stella_bmg_mat *ctx;
	ierr = MatShellGetContext(mat, (void**) &ctx);

	if (ctx->nd == 2) {
		unsigned int i, j;

		grid_coord_2d *coords = (grid_coord_2d*) malloc(n*m*sizeof(grid_coord_2d));

		for (i = 0; i < m; i++) {
			for (j = 0; j < n; j++) {
				// This logic may be incorrect
				if (idxm[i].i-1 == idxn[j].i &&
				    idxm[i].j == idxn[j].j) coords[i*n+j].dir = BMG2_W;
				if (idxm[i].i+1 == idxn[j].i &&
				    idxm[i].j == idxn[j].j) coords[i*n+j].dir = BMG2_E;
				if (idxm[i].i == idxn[j].i &&
				    idxm[i].j-1 == idxn[j].j) coords[i*n+j].dir = BMG2_S;
				if (idxm[i].i == idxn[j].i &&
				    idxm[i].j+1 == idxn[j].j) coords[i*n+j].dir = BMG2_N;
				if (idxm[i].i == idxn[j].i &&
				    idxm[i].j == idxn[j].j) coords[i*n+j].dir = BMG2_C;
				if (idxm[i].i-1 == idxn[j].i &&
				    idxm[i].j-1 == idxn[j].j) coords[i*n+j].dir = BMG2_SW;
				if (idxm[i].i+1 == idxn[j].i &&
				    idxm[i].j-1 == idxn[j].j) coords[i*n+j].dir = BMG2_SE;
				if (idxm[i].i-1 == idxn[j].i &&
				    idxm[i].j+1 == idxn[j].j) coords[i*n+j].dir = BMG2_NW;
				if (idxm[i].i+1 == idxn[j].i &&
				    idxm[i].j+1 == idxn[j].j) coords[i*n+j].dir = BMG2_NE;
				coords[i*n+j].i = idxm[i].i;
				coords[i*n+j].j = idxm[i].j;
			}
		}

		bmg2_operator_set(ctx->op2, n*m, coords, (double*)v);

		free(coords);
	} else {
		unsigned int i, j;

		grid_coord_3d *coords = (grid_coord_3d*) malloc(n*m*sizeof(grid_coord_3d));

		for (i = 0; i < m; i++) {
			for (j = 0; j < n; j++) {
				coords[i*n+j].i = idxm[i].i;
				coords[i*n+j].j = idxm[i].j;
				coords[i*n+j].k = idxm[i].k;

				if (idxm[i].i == idxn[j].i &&
				    idxm[i].j == idxn[j].j &&
				    idxm[i].k == idxn[j].k) coords[i*n+j].dir = BMG3_P;
				else if (idxm[i].i-1 == idxn[j].i &&
				    idxm[i].j == idxn[j].j &&
				    idxm[i].k == idxn[j].k) coords[i*n+j].dir = BMG3_PW;
				else if (idxm[i].i+1 == idxn[j].i &&
				    idxm[i].j == idxn[j].j &&
				    idxm[i].k == idxn[j].k) { // PE
					coords[i*n+j].dir = BMG3_PW;
					coords[i*n+j].i++;
				}
				else if (idxm[i].i == idxn[j].i &&
				    idxm[i].j-1 == idxn[j].j &&
				    idxm[i].k == idxn[j].k) coords[i*n+j].dir = BMG3_PS;
				else if (idxm[i].i == idxn[j].i &&
				    idxm[i].j+1 == idxn[j].j &&
				    idxm[i].k == idxn[j].k) { // PN
					coords[i*n+j].dir = BMG3_PS;
					coords[i*n+j].j++;
				}
				else if (idxm[i].i-1 == idxn[j].i &&
				    idxm[i].j-1 == idxn[j].j &&
				    idxm[i].k == idxn[j].k) coords[i*n+j].dir = BMG3_PSW;
				else if (idxm[i].i+1 == idxn[j].i &&
				    idxm[i].j+1 == idxn[j].j &&
				    idxm[i].k == idxn[j].k) { // PNE
					coords[i*n+j].dir = BMG3_PSW;
					coords[i*n+j].i++;
					coords[i*n+j].j++;
				}
				else if (idxm[i].i+1 == idxn[j].i &&
				    idxm[i].j-1 == idxn[j].j &&
				    idxm[i].k == idxn[j].k) { // PSE
					coords[i*n+j].dir = BMG3_PNW;
					coords[i*n+j].i++;
				}
				else if (idxm[i].i-1 == idxn[j].i &&
				    idxm[i].j+1 == idxn[j].j &&
				    idxm[i].k == idxn[j].k) { // PNW
					coords[i*n+j].dir = BMG3_PNW;
					coords[i*n+j].j++;
				}
				else if (idxm[i].i == idxn[j].i &&
				    idxm[i].j == idxn[j].j &&
				    idxm[i].k-1 == idxn[j].k) coords[i*n+j].dir = BMG3_B;
				else if (idxm[i].i == idxn[j].i &&
				    idxm[i].j == idxn[j].j &&
				    idxm[i].k+1 == idxn[j].k) { // T
					coords[i*n+j].dir = BMG3_B;
					coords[i*n+j].k++;
				}
				else if (idxm[i].i-1 == idxn[j].i &&
				    idxm[i].j == idxn[j].j &&
				    idxm[i].k-1 == idxn[j].k) coords[i*n+j].dir = BMG3_BW;
				else if (idxm[i].i+1 == idxn[j].i &&
				    idxm[i].j == idxn[j].j &&
				    idxm[i].k+1 == idxn[j].k) { // TNE
					coords[i*n+j].dir = BMG3_BW;
					coords[i*n+j].i++;
					coords[i*n+j].k++;
				}
				else if (idxm[i].i == idxn[j].i &&
				    idxm[i].j-1 == idxn[j].j &&
				    idxm[i].k-1 == idxn[j].k) coords[i*n+j].dir = BMG3_BS;
				else if (idxm[i].i == idxn[j].i &&
				    idxm[i].j+1 == idxn[j].j &&
				    idxm[i].k+1 == idxn[j].k){ // TN
					coords[i*n+j].dir = BMG3_BS;
					coords[i*n+j].j++;
					coords[i*n+j].k++;
				}
				else if (idxm[i].i-1 == idxn[j].i &&
				    idxm[i].j-1 == idxn[j].j &&
				    idxm[i].k-1 == idxn[j].k) coords[i*n+j].dir = BMG3_BSW;
				else if (idxm[i].i+1 == idxn[j].i &&
				    idxm[i].j+1 == idxn[j].j &&
				    idxm[i].k+1 == idxn[j].k){ // TNE
					coords[i*n+j].dir = BMG3_BSW;
					coords[i*n+j].i++;
					coords[i*n+j].j++;
					coords[i*n+j].k++;
				}
				else if (idxm[i].i == idxn[j].i &&
				    idxm[i].j+1 == idxn[j].j &&
				    idxm[i].k-1 == idxn[j].k){ // BN
					coords[i*n+j].dir = BMG3_BN;
					coords[i*n+j].j++;
				}
				else if (idxm[i].i == idxn[j].i &&
				    idxm[i].j-1 == idxn[j].j &&
				    idxm[i].k+1 == idxn[j].k){ // TS
					coords[i*n+j].dir = BMG3_BN;
					coords[i*n+j].k++;
				}
				else if (idxm[i].i+1 == idxn[j].i &&
				    idxm[i].j == idxn[j].j &&
				    idxm[i].k-1 == idxn[j].k){ // BE
					coords[i*n+j].dir = BMG3_BE;
					coords[i*n+j].i++;
				}
				else if (idxm[i].i-1 == idxn[j].i &&
				    idxm[i].j == idxn[j].j &&
				    idxm[i].k+1 == idxn[j].k){ // TW
					coords[i*n+j].dir = BMG3_BE;
					coords[i*n+j].k++;
				}
				else if (idxm[i].i+1 == idxn[j].i &&
				    idxm[i].j-1 == idxn[j].j &&
				    idxm[i].k-1 == idxn[j].k){ // BSE
					coords[i*n+j].dir = BMG3_BSE;
					coords[i*n+j].i++;
				}
				else if (idxm[i].i-1 == idxn[j].i &&
				    idxm[i].j+1 == idxn[j].j &&
				    idxm[i].k+1 == idxn[j].k){ // TNW
					coords[i*n+j].dir = BMG3_BSE;
					coords[i*n+j].j++;
					coords[i*n+j].k++;
				}
				else if (idxm[i].i-1 == idxn[j].i &&
				    idxm[i].j+1 == idxn[j].j &&
				    idxm[i].k-1 == idxn[j].k){ // BNW
					coords[i*n+j].dir = BMG3_BNW;
					coords[i*n+j].j++;
				}
				else if (idxm[i].i+1 == idxn[j].i &&
				    idxm[i].j-1 == idxn[j].j &&
				    idxm[i].k+1 == idxn[j].k){ // TSE
					coords[i*n+j].dir = BMG3_BNW;
					coords[i*n+j].i++;
					coords[i*n+j].k++;
				}
				else if (idxm[i].i+1 == idxn[j].i &&
				    idxm[i].j+1 == idxn[j].j &&
				    idxm[i].k-1 == idxn[j].k){ // BNE
					coords[i*n+j].dir = BMG3_BNE;
					coords[i*n+j].i++;
					coords[i*n+j].j++;
				}
				else if (idxm[i].i-1 == idxn[j].i &&
				    idxm[i].j-1 == idxn[j].j &&
				    idxm[i].k+1 == idxn[j].k){ // TSW
					coords[i*n+j].dir = BMG3_BNE;
					coords[i*n+j].k++;
				} else {
					SETERRQ(PETSC_COMM_WORLD, PETSC_ERR_LIB, "(Cedar) Invalid stencil direction");
				}
			}
		}
		bmg3_operator_set(ctx->op3, n*m, coords, (double*)v);

		free(coords);
	}
	#endif
	return 0;
}

PetscErrorCode stella_bmg_mult(Mat A, Vec x, Vec b)
{
	#ifdef WITH_BOXMG
	PetscErrorCode ierr;
	stella_bmg_mat *ctx;
	double *barr;
	const double *xarr;

	ierr = MatShellGetContext(A, (void**) &ctx);CHKERRQ(ierr);
	ierr = VecGetArray(b, &barr);CHKERRQ(ierr);
	ierr = VecGetArrayRead(x, &xarr);CHKERRQ(ierr);

	if (ctx->nd == 2)
		bmg2_operator_apply(ctx->op2, xarr, barr);
	else
		bmg3_operator_apply(ctx->op3, xarr, barr);

	ierr = VecRestoreArrayRead(x, &xarr);CHKERRQ(ierr);
	ierr = VecRestoreArray(b, &barr);CHKERRQ(ierr);
	#endif

	return 0;
}
