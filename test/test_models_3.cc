#include <limits.h>
#include <math.h>
#include <cmath>
#include <mpi.h>
#include <vector>
#include "gtest/gtest.h"

#include "problem.h"
#include "map.h"
#include "grid.h"
#include "solver.h"


namespace {


static double estimate_order(std::vector<double> & hs, std::vector<double> & norms)
{
	return std::log(norms[1] / norms[0]) / std::log(hs[1] / hs[0]);
}

static double norm(double v[], int len)
{
	int i;

	double curr = 0;
	for (i = 0; i < len; i++) {
		if (fabs(v[i]) > curr) curr = fabs(v[i]);
	}

	return curr;
}


static double mpi_norm(double v[], int len)
{
	int rank;
	double n;

	double mynorm = norm(v, len);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);

	MPI_Reduce(&mynorm, &n, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
	if (rank == 0) return n;
	else return 0;
}


TEST(Models3, UniformSine)
{
	std::vector<double> norms;
	std::vector<double> hs;
	int nx = 20;
	double tol = 1e-5;

	for (int k = 0; k < 2; k++) {
		int i, ierr;

		hs.push_back(1.0/(nx-1));

		grid *grd = grid_create(0, 1, nx,
		                        0, 1, nx,
		                        0, 1, nx);
		problem *pb = problem_create(SIN, 3, 0);
		solver *sol = solver_create(grd, pb);

		ierr = solver_init(sol, grd);
		ASSERT_EQ(ierr, 0);
		ierr = solver_run(sol);
		ASSERT_EQ(ierr, 0);

		double *u = new double[grd->num_pts];
		double *diff = new double[grd->num_pts];
		grid_eval(grd, pb->sol, u);

		for (i = 0; i < grd->num_pts; i++)
			diff[i] = sol->state->phi[i] - u[i];

		double nrm = mpi_norm(diff, grd->num_pts);
		norms.push_back(nrm);

		problem_destroy(pb); free(pb);
		ierr = solver_destroy(sol); free(sol);
		ASSERT_EQ(ierr, 0);
		grid_destroy(grd); free(grd);
		delete[] u;
		delete[] diff;

		nx = (nx-1)*2 + 1;
	}

	int rank;
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	if (rank == 0) {
		double order = estimate_order(hs, norms);
		ASSERT_GT(order, 2 - tol);
	}
}


TEST(Models3, UniformTSine)
{
	std::vector<double> norms;
	std::vector<double> hs;
	int nx = 53;
	double tol = 1e-3;

	for (int k = 0; k < 2; k++) {
		int i, ierr;

		hs.push_back(1.0/(nx-1));

		grid *grd = grid_create(0, 1, nx,
		                        0, 1, nx,
		                        0, 1, nx);
		problem *pb = problem_create(TSINE, 3, 0);
		solver *sol = solver_create(grd, pb);

		ierr = solver_init(sol, grd);
		ASSERT_EQ(ierr, 0);
		ierr = solver_run(sol);
		ASSERT_EQ(ierr, 0);

		double *u = new double[grd->num_pts];
		double *diff = new double[grd->num_pts];
		grid_eval(grd, pb->sol, u);

		for (i = 0; i < grd->num_pts; i++)
			diff[i] = sol->state->phi[i] - u[i];

		double nrm = mpi_norm(diff, grd->num_pts);
		norms.push_back(nrm);

		problem_destroy(pb); free(pb);
		ierr = solver_destroy(sol); free(sol);
		ASSERT_EQ(ierr, 0);
		grid_destroy(grd); free(grd);
		delete[] u;
		delete[] diff;

		nx = (nx-1)*2 + 1;
	}

	int rank;
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	if (rank == 0) {
		double order = estimate_order(hs, norms);
		ASSERT_GT(order, 2 - tol);
	}
}


TEST(Models3, UniformMixed)
{
	std::vector<double> norms;
	std::vector<double> hs;
	int nx = 40;
	double tol = 1e-2;

	for (int k = 0; k < 2; k++) {
		int i, ierr;

		hs.push_back(1.0/(nx-1));

		grid *grd = grid_create(0, 1, nx,
		                        0, 1, nx,
		                        0, 1, nx);
		problem *pb = problem_create(MIXED, 3, 0);
		solver *sol = solver_create(grd, pb);

		ierr = solver_init(sol, grd);
		ASSERT_EQ(ierr, 0);
		ierr = solver_run(sol);
		ASSERT_EQ(ierr, 0);

		double *u = new double[grd->num_pts];
		double *diff = new double[grd->num_pts];
		grid_eval(grd, pb->sol, u);

		for (i = 0; i < grd->num_pts; i++)
			diff[i] = sol->state->phi[i] - u[i];

		double nrm = mpi_norm(diff, grd->num_pts);
		norms.push_back(nrm);

		problem_destroy(pb); free(pb);
		ierr = solver_destroy(sol); free(sol);
		ASSERT_EQ(ierr, 0);
		grid_destroy(grd); free(grd);
		delete[] u;
		delete[] diff;

		nx = (nx-1)*2 + 1;
	}

	int rank;
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	if (rank == 0) {
		double order = estimate_order(hs, norms);
		ASSERT_GT(order, 2 - tol);
	}
}


TEST(Models3, Periodic)
{
	std::vector<double> norms;
	std::vector<double> hs;
	int nx = 40;
	double tol = 1e-2;

	for (int k = 0; k < 2; k++) {
		int i, ierr;

		hs.push_back(1.0/(nx-1));

		grid *grd = grid_create(0, 1, nx,
		                        0, 1, nx,
		                        0, 1, nx);

		grd->periodic[1] = 1;
		map *mp = map_create(MAP_POLAR, grd->nd);
		grid_apply_map(grd, mp);

		problem *pb = problem_create(PERIODIC, 3, mp->id);
		solver *sol = solver_create(grd, pb);

		ierr = solver_init(sol, grd);
		ASSERT_EQ(ierr, 0);
		ierr = solver_run(sol);
		ASSERT_EQ(ierr, 0);

		double *u = new double[grd->num_pts];
		double *diff = new double[grd->num_pts];
		grid_eval(grd, pb->sol, u);

		for (i = 0; i < grd->num_pts; i++)
			diff[i] = sol->state->phi[i] - u[i];

		double nrm = mpi_norm(diff, grd->num_pts);
		norms.push_back(nrm);

		map_destroy(mp); free(mp);
		problem_destroy(pb); free(pb);
		ierr = solver_destroy(sol); free(sol);
		ASSERT_EQ(ierr, 0);
		grid_destroy(grd); free(grd);
		delete[] u;
		delete[] diff;

		nx = (nx-1)*2 + 1;
	}

	int rank;
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	if (rank == 0) {
		double order = estimate_order(hs, norms);
		ASSERT_GT(order, 2 - tol);
	}
}

TEST(Models3, Jump)
{
	std::vector<double> norms;
	std::vector<double> hs;
	std::array<int, 2> nvals = {40, 80};
	double tol = 1e-1;

	for (int k = 0; k < 2; k++) {
		int i, ierr;

		auto nx = nvals[k];

		hs.push_back(2.0/(nx-1));

		grid *grd = grid_create(-1, 1, nx,
		                        -1, 1, nx,
		                        -1, 1, nx);
		problem *pb = problem_create(JUMP, 3, 0);
		solver *sol = solver_create(grd, pb);

		ierr = solver_init(sol, grd);
		ASSERT_EQ(ierr, 0);
		ierr = solver_run(sol);
		ASSERT_EQ(ierr, 0);

		double *u = new double[grd->num_pts];
		double *diff = new double[grd->num_pts];
		grid_eval(grd, pb->sol, u);

		for (i = 0; i < grd->num_pts; i++)
			diff[i] = sol->state->phi[i] - u[i];

		double nrm = mpi_norm(diff, grd->num_pts);
		norms.push_back(nrm);

		problem_destroy(pb); free(pb);
		ierr = solver_destroy(sol); free(sol);
		ASSERT_EQ(ierr, 0);
		grid_destroy(grd); free(grd);
		delete[] u;
		delete[] diff;
	}

	int rank;
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	if (rank == 0) {
		double order = estimate_order(hs, norms);
		ASSERT_GT(order, 2 - tol);
	}
}

}  // namespace
