
sudo apt-get install -y cmake gcc mpich doxygen git
export TMPDIR=/tmp

# Go get PETSc and build it.
PN=petsc
PV=3.7.7
wget http://ftp.mcs.anl.gov/pub/petsc/release-snapshots/${PN}-${PV}.tar.gz
tar xfz ${PN}-${PV}.tar.gz
mv ${PN}-${PV} petsc
cd petsc
python2 ./configure --with-cc=mpicc --with-cxx=mpicxx --with-fc=mpif90 --with-debug=yes --with-shared-libraries=1 --download-hypre --download-f2cblaslapack=1
make
cd ..
