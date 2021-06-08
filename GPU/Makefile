LIB_LOAD_MATRIX  = -L../ReadMatrixMarket -lloadmatrix
INC_LOAD_MATRIX = -I../ReadMatrixMarket

CC = nvcc
CFLAGS =   -arch=sm_70 -O3  #-ftz=true  # V100

all: test_spmv  test_ilu0


test_spmv: spmv.o
	 $(CC) $(CFLAGS)  -o $@ $^ $(LIB_LOAD_MATRIX) -lcublas -lcusparse -lcudart -L$(MKLROOT)/lib -lmkl_intel_lp64 -lmkl_sequential -lmkl_core -lm -ldl -lgomp


test_ilu0: ilu0.o
	 $(CC) $(CFLAGS)  -o $@ $^ $(LIB_LOAD_MATRIX) -lcublas -lcusparse -lcudart  -lm -ldl -lgomp

spmv.o : spmv.cu
	$(CC) $(CFLAGS) $(INC_LOAD_MATRIX) -c $^ -o $@

ilu0.o : ilu0.cu
	$(CC) $(CFLAGS) $(INC_LOAD_MATRIX) -c $^ -o $@

clean:
	rm -f test_spmv  test_ilu0  *~ *.o
