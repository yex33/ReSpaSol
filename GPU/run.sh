#!/bin/bash
# Define a string variable with a value
Matrices1="2cubes_sphere  ASIC_320ks  Baumann  cfd2  crashbasis  ct20stif  dc1  Dubcova3  FEM_3D_thermal2  G2_circuit  Goodwin_095  parabolic_fem  ss1  stomach  thermomech_TK  xenon2"
Matrices2="af_shell2  atmosmodd  cage13  CurlCurl_2  ecology2  matrix-new_3  ML_Laplace  offshore  para-10  Si87H76  tmt_unsym"
Matrices3="af_shell10  atmosmodl  dielFilterV2real  Geo_1438  Hook_1498  nlpkkt80  Serena  StocF-1465  Transport"


# Iterate the string variable using for loop
for val in $Matrices1; do
    for VARIABLE in 1 2;  do
        echo $val
        ./$1 ../matrices/moderate/$val.mtx
    done
    echo ""
done


for val in $Matrices2; do
    for VARIABLE in 1 2;  do
        echo $val
        ./$1 ../matrices/big/$val.mtx
    done
    echo ""
done


for val in $Matrices3; do
    for VARIABLE in 1 2;  do
        echo $val
        ./$1 /scratch/shared/mawussi/big2/$val.mtx
    done
    echo ""
done
