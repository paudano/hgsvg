all: mcutils/src/samToBed blasr/alignment/bin/blasr

mcutils/src/samToBed:
	cd mcutils/src && make

blasr/alignment/bin/blasr:
	cd blasr && make HDF5INCLUDEDIR=$(PWD)/../include HDF5LIBDIR=$(PWD)/../lib CPP="g++ -static -std=c++03 " -j 8
