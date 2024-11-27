#!/bin/bash

Matrices1="2cubes_sphere \
  ASIC_320ks \
  Baumann \
  cfd2 \
  crashbasis \
  ct20stif \
  dc1 \
  Dubcova3 \
  ecology2 \
  FEM_3D_thermal2 \
  G2_circuit \
  Goodwin_095 \
  matrix-new_3 \
  offshore \
  para-10 \
  parabolic_fem \
  ss1 \
  stomach \
  thermomech_TK \
  tmt_unsym \
  xenon2"
Matrices2=" \
  af_shell10 \
  af_shell2 \
  atmosmodd \
  atmosmodl \
  cage13 \
  CurlCurl_2 \
  dielFilterV2real \
  Geo_1438 \
  Hook_1498 \
  ML_Laplace \
  nlpkkt80 \
  Serena \
  Si87H76 \
  StocF-1465 \
  Transport"

for val in $Matrices1; do
	echo $val
	for n in {0..10}; do
		OMP_NUM_THREADS=4 taskset -c 0,1,2,3 ./test_pardiso ./matrices/moderate/$val.mtx ./bench_pardiso_v1.csv
	done
	echo ""
done

for val in $Matrices2; do
	echo $val
	for n in {0..10}; do
		OMP_NUM_THREADS=4 taskset -c 0,1,2,3 ./test_pardiso ./matrices/big/$val.mtx ./bench_pardiso_v1.csv
	done
	echo ""
done

# for n in {0..10}; do
# 	OMP_NUM_THREADS=4 taskset -c 0,1,2,3 ./test_pardiso ./matrices/moderate/ASIC_320ks.mtx ./bench_pardiso_v1.csv
# done
