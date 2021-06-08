* * *

## About ReSpaSol  

This repository is a collection of codes and scripts to evaluate the performance benefit
of using reduced precision in parallel sparse linear systems solvers. The associated
manuscript is available at [http://eprints.maths.manchester.ac.uk/2811/](http://eprints.maths.manchester.ac.uk/2811).


## Dependencies

### Parse direct solvers 
To use the routines provided here, the parallel sparse solvers evaluated must be downloaded and compiled.
These include:

1. MUMPS 5.2.1 or recent version available at [http://mumps.enseeiht.fr/](http://mumps.enseeiht.fr)
2. SuperLU_MT 3.1  available at [https://github.com/group-gu/superlu-mt](https://github.com/group-gu/superlu-mt)
3. SuperLU 5.2.1 available at [https://github.com/xiaoyeli/superlu](https://github.com/xiaoyeli/superlu)
4. Pardiso from MKL 2020 or a recent version available at [http://scc.ustc.edu.cn/zlsc/intel/2020/mkl/ps2020/get_started.htm](http://scc.ustc.edu.cn/zlsc/intel/2020/mkl/ps2020/get_started.htm)


### Sparse matrices
The matrices used in the experimentations are selected  from the [https://sparse.tamu.edu/](SuiteSparse Matrix Collection).
We use the matrix market format. The matrix can be downloaded directly from the SuiteSparse Matrix Collection.

| Matrices |Link to download Matrix Market formats |
| :---: | :---: | 
| 2cubes_sphere | https://suitesparse-collection-website.herokuapp.com/MM/Um/2cubes_sphere.tar.gz |
| :---: | :---: |
| ASIC_320ks | https://suitesparse-collection-website.herokuapp.com/MM/Sandia/ASIC_320ks.tar.gz   |