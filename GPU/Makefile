CSR_LIB = -L/home/mawussi/Software/readmatrixmarket2csr/ -lcsr 
CSR_INC  = -I/home/mawussi/Software/readmatrixmarket2csr

CC = nvcc
CFLAGS =   -arch=sm_70 -O3  #-ftz=true  # V100
#CFLAGS =   -arch=sm_70 -O3   # V100
#CFLAGS =   -arch=sm_60 -O3 -G # P100

all: spmv solve  ilu0



spmv: spmv.o
	 $(CC) $(CFLAGS)  -o $@ $^ $(CSR_LIB) -lcublas -lcusparse -lcudart -L$(MKLROOT)/lib -lmkl_intel_lp64 -lmkl_sequential -lmkl_core -lm -ldl -lgomp


ilu0: ilu0.o
	 $(CC) $(CFLAGS)  -o $@ $^ $(CSR_LIB) -lcublas -lcusparse -lcudart -L$(MKLROOT)/lib -lmkl_intel_lp64 -lmkl_sequential -lmkl_core -lm -ldl -lgomp

solve: solve.o
	 $(CC) $(CFLAGS)  -o $@ $^ $(CSR_LIB) -lcublas -lcusparse -lcudart -L$(MKLROOT)/lib -lmkl_intel_lp64 -lmkl_sequential -lmkl_core -lm -ldl -lgomp


spmv.o : spmv.cu
	$(CC) $(CFLAGS) $(CSR_INC) -c $^ -o $@

ilu0.o : ilu0.cu
	$(CC) $(CFLAGS) $(CSR_INC) -c $^ -o $@


solve.o : solve.cu
	$(CC) $(CFLAGS) $(CSR_INC) -c $^ -o $@
clean:
	rm -f spmv  ilu0  solve *.o *~
