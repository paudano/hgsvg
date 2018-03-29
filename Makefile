all: mcutils/src/samToBed blasr/alignment/bin/blasr

mcutils/src/samToBed:
	cd mcutils/src && make

blasr/alignment/bin/blasr:
	cd blasr && make HDF5INCLUDEDIR=$(abspath ../hdf5/build/include) HDF5LIBDIR=$(abspath ../hdf5/build/lib) CPP="g++ -std=c++03 " -j 8
