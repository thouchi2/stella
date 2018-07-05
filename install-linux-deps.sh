sudo apt-get update -qq
sudo apt-get install -y cmake gcc gfortran mpich doxygen git
export TMPDIR=/tmp

# Go get PETSc and build it.
git clone https://bitbucket.org/petsc/petsc petsc
pushd $PETSC_DIR
./configure --with-debug=yes --with-shared-libraries=1 --download-hypre --download-mpich --download-f2cblaslapack=1
make
popd

wget -q --no-check-certificate http://www.mpich.org/static/downloads/3.2/mpich-3.2.tar.gz
tar -xzf mpich-3.2.tar.gz
cd mpich-3.2
mkdir build && cd build
../configure CC=$PRK_CC CXX=$PRK_CXX --disable-fortran --prefix=$TRAVIS_ROOT/mpich
make -j4
