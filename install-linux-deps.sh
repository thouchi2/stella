sudo apt-get update -qq
sudo apt-get install -y cmake gcc gfortran mpich doxygen git
export TMPDIR=/tmp

# Go get PETSc and build it.
git clone https://bitbucket.org/petsc/petsc petsc
pushd $PETSC_DIR
./configure --with-mpi=1 --with-debug=$DEBUG --with-shared-libraries=1 --download-hypre --downloadmpich
make
popd
