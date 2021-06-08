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
The matrices are divided in two groups. The first 21 matrices are from the medium size group with 700, 000 to 5, 000, 000 nonzero elements. It takes a few
seconds on average to factorize these matrices using 10 cores. The second group contains larger matrices with 7,000,000
to 64,000,000 nonzeros and it takes on average a few minutes to factorize most of the matrices in this
group using 10 cores.

#### Moderate size matrices 
| Matrices |Links to download Matrix Market formats |
| :--- | :--- | 
|2cubes_sphere | https://suitesparse-collection-website.herokuapp.com/MM/Um/2cubes_sphere.tar.gz |
|ASIC_320ks | https://suitesparse-collection-website.herokuapp.com/MM/Sandia/ASIC_320ks.tar.gz   |
|Baumann|https://suitesparse-collection-website.herokuapp.com/MM/Watson/Baumann.tar.gz |
|cfd2|https://suitesparse-collection-website.herokuapp.com/MM/Rothberg/cfd2.tar.gz |
|crashbasis |https://suitesparse-collection-website.herokuapp.com/MM/QLi/crashbasis.tar.gz|
|ct20stif |https://suitesparse-collection-website.herokuapp.com/MM/Boeing/ct20stif.tar.gz |
|dc1 |https://suitesparse-collection-website.herokuapp.com/MM/IBM_EDA/dc1.tar.gz |
|Dubcova3 |https://suitesparse-collection-website.herokuapp.com/MM/UTEP/Dubcova3.tar.gz|
|ecology2 | https://suitesparse-collection-website.herokuapp.com/MM/McRae/ecology2.tar.gz|
|FEM_3D_thermal2| https://suitesparse-collection-website.herokuapp.com/MM/Botonakis/FEM_3D_thermal2.tar.gz|
|G2_circuit |https://suitesparse-collection-website.herokuapp.com/MM/AMD/G2_circuit.tar.gz|
|Goodwin_095| https://suitesparse-collection-website.herokuapp.com/MM/Goodwin/Goodwin_095.tar.gz|
|matrix-new_3|https://suitesparse-collection-website.herokuapp.com/MM/Schenk_IBMSDS/matrix-new_3.tar.gz |
|offshore |https://suitesparse-collection-website.herokuapp.com/MM/Um/offshore.tar.gz|
|para-10 |https://suitesparse-collection-website.herokuapp.com/MM/Schenk_ISEI/para-10.tar.gz|
|parabolic_fem|https://suitesparse-collection-website.herokuapp.com/MM/Wissgott/parabolic_fem.tar.gz |
|ss1|https://suitesparse-collection-website.herokuapp.com/MM/VLSI/ss1.tar.gz |
|stomach|https://suitesparse-collection-website.herokuapp.com/MM/Norris/stomach.tar.gz|
|thermomech_TK|https://suitesparse-collection-website.herokuapp.com/MM/Botonakis/thermomech_TK.tar.gz|
|tmt_unsym|https://suitesparse-collection-website.herokuapp.com/MM/CEMW/tmt_unsym.tar.gz|
|xenon2|https://suitesparse-collection-website.herokuapp.com/MM/Ronis/xenon2.tar.gz|