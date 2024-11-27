LIB_MKL_THREADED = -L$(MKLROOT)/lib/intel64 -lmkl_intel_lp64 -lmkl_gnu_thread -lmkl_core -lm -ldl
LIB_MKL_SEQ = -L$(MKLROOT)/lib/intel64 -lmkl_intel_lp64 -lmkl_sequential -lmkl_core -lm -ldl


#####  MUMPS
# MUMPSDIR = /fserver/mawussiz/Software/MUMPS_5.2.1
# LMETISDIR = /fserver/mawussiz/Software/metis_dir
MUMPSDIR = /usr
LMETISDIR = /usr
SCOTCHDIR = /usr
# MPIDIR = $(I_MPI_ROOT)/intel64
MPIDIR = /usr
LIB_SCALAP  = -L$(MKLROOT)/lib/intel64 -lmkl_scalapack_lp64 -lmkl_blacs_intelmpi_lp64 -lmkl_intel_lp64 -lmkl_gnu_thread -lmkl_core
LIB_MPI = -L$(MPIDIR)/lib  -lmpi
LIB_MUMPS = -L$(MUMPSDIR)/lib -ldmumps -lsmumps -lmumps_common  -lpord
LIB_METIS    = -L$(LMETISDIR)/lib -lmetis
INC_METIS    = -I$(LMETISDIR)/include
LIB_SCOTCH	 = -L$(SCOTCHDIR)/lib -lesmumps -lscotch -lscotcherr

## SuperLU
# LIB_SUPERLU_THREADED = -L/fserver/mawussiz/Software/superlu_mt_dir/lib -lsuperlu_mt_OPENMP
# LIB_SUPERLU_SEQ = -L/fserver/mawussiz/Software/superlu_dir/lib -lsuperlu_5.1

# LIB_SUPERLU_THREADED = -L/usr/lib -lsuperlu_mt_OPENMP
LIB_SUPERLU_THREADED = -L/usr/lib -lsuperlu_mt_PTHREAD
LIB_SUPERLU_SEQ = -L/usr/lib -lsuperlu

LIB_LOAD_MATRIX  = -L./ReadMatrixMarket -lloadmatrix

LIB_BENCHMARK = -L./benchmark/build/src -lbenchmark -lpthread


LIBRARY_SEQ =  $(LIB_SUPERLU_SEQ) $(LIB_LOAD_MATRIX) $(LIB_MKL_THREADED)
LIBRARY_PARDISO =  $(LIB_LOAD_MATRIX) $(LIB_MKL_THREADED)	
LIBRARY_THREADED =  $(LIB_SUPERLU_THREADED) $(LIB_LOAD_MATRIX) $(LIB_MKL_SEQ)


INC_LOAD_MATRIX = -I./ReadMatrixMarket
INC_BENCHMARK = -I./benchmark/include
INC_PARDISO      = -I$(MKLROOT)/include $(INC_LOAD_MATRIX)
# INC_SEQ          = -I/fserver/mawussiz/Software/superlu_dir/include $(INC_LOAD_MATRIX)
# INC_THREADED     = -I/fserver/mawussiz/Software/superlu_mt_dir/include  $(INC_LOAD_MATRIX)
INC_SEQ          = -I/usr/include/superlu $(INC_LOAD_MATRIX)
INC_THREADED     = -I/usr/include/superlu_mt  $(INC_LOAD_MATRIX)
INC_SUITESPARSE = -I/fserverlts/mawussiz/Software/SuiteSparse/include
INC_MPI = -I$(MPIDIR)/include
INC_MUMPS = -I$(MUMPSDIR)/include $(INC_MPI) -I$(MKLROOT)/include $(INC_LOAD_MATRIX) $(INC_METIS)

CC = g++
CFLAGS	= -fopenmp -fPIC -std=c++20 -Wno-deprecated -Wall -O0 -g

all: test

test: test_superLU_MT test_superILU test_pardiso  test_mumps test_spmv bench bench_spmv

test_mumps: test_mumps.o
	mpifort $(CFLAGS)  -o $@ $^ $(LIB_MUMPS) $(LIB_SCOTCH) $(LIB_METIS)   $(LIB_SCALAP)  $(LIB_MPI)  $(LIB_LOAD_MATRIX) -fopenmp -lgfortran   -lstdc++ -ldl -g

test_superILU: test_superILU.o 
	 $(CC) $(CFLAGS)  -o $@ $^ $(LIBRARY_SEQ)


test_pardiso: test_pardiso.o
	$(CC) -m64 $(CFLAGS)  -o $@ $^ $(LIBRARY_PARDISO)

test_spmv: test_spmv.o
	$(CC) -m64 $(CFLAGS)  -o $@ $^ $(LIBRARY_PARDISO)

test_superLU_MT: test_superLU_MT.o 
	$(CC) $(CFLAGS)  -o $@ $^ $(LIBRARY_THREADED)

bench: bench.o
	$(CC) $(CFLAGS) -o $@ $^ $(LIB_BENCHMARK)

bench_spmv: bench_spmv.o
	$(CC) -m64 $(CFLAGS) -o $@ $^ $(LIBRARY_PARDISO) $(LIB_BENCHMARK)

test_mumps.o: test_mumps.c
	mpicxx -DAdd_ -m64  $(CFLAGS) $(INC_MUMPS) -fPIC -c $^ -o $@

test_superILU.o: test_superILU.c
	$(CC) $(CFLAGS) $(INC_SEQ) -fPIC -c $^ -o $@

test_pardiso.o: test_pardiso.c
	$(CC) $(CFLAGS) $(INC_PARDISO) -fPIC -c $^ -o $@

test_spmv.o: test_spmv.c
	$(CC) $(CFLAGS) $(INC_PARDISO) -fPIC -c $^ -o $@

test_superLU_MT.o: test_superLU_MT.c
	$(CC) $(CFLAGS) $(INC_THREADED) -fPIC -c $^ -o $@

bench.o: bench.cc
	$(CC) $(CFLAGS) $(INC_BENCHMARK) -fPIC -c $^ -o $@

bench_spmv.o: bench_spmv.cc
	$(CC) $(CFLAGS) $(INC_PARDISO) $(INC_BENCHMARK) -fPIC -c $^ -o $@


clean:
	rm -f  *.o *~ test_superLU_MT  test_pardiso test_mumps test_spmv test_superILU bench bench_spmv
