Single Block Testing Framework
==============================

MPI framework for setting up model problems used in testing the Efield solver
and discretization from PlasComCM.  The solver is designed to solve
`- \div \eps \grad \phi = rhs` where `\eps` may be a grid-aligned discontinuity.
The framework includes:

- Creating distributed structured grid
- Defining analytical mappings that can be applied to a structured grid
- Boundary condition setup
  + Neumann
  + Dirichlet
- Interface condition setup

Building
--------
    mkdir build
	cmake ..
	make

Running
-------
	mpirun -np 4 ./tools/driver-sb

Command Line Options
--------------------

| Option             | Description                                          |
|--------------------|------------------------------------------------------|
| `-l p`             | List Model Problems                                  |
| `-l m`             | List Mappings                                        |
| `-x <nx>`          | Specify number of grid points in first dimension     |
| `-m <map key>`     | Specify the mapping                                  |
| `-p <problem key>` | Specify the model problem                            |
| `-d <n>`           | Specify number of grid points in each dimension (2D) |
| `-s <n>`           | Specify number of grid points in each dimension (3D) |
