sudo apt-get update -qq
sudo apt-get install -y cmake gcc gfortran mpich doxygen git
export TMPDIR=/tmp

ln -s /usr/local/lib/gcc/5/libgfortran.dylib /usr/local/lib/libgfortran.dylib
ln -s /usr/local/lib/gcc/5/libgfortran.a /usr/local/lib/libgfortran.a

# Go get PETSc and build it.
git clone https://bitbucket.org/petsc/petsc petsc
pushd $PETSC_DIR
./configure --with-debug=yes --with-shared-libraries=1 --download-hypre --download-mpich --download-f2cblaslapack=1
make
make test
popd
