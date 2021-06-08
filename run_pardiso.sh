#!/bin/bash
# Define a string variable with a value

Matrices1="2cubes_sphere  ASIC_320ks Baumann cfd2 crashbasis ct20stif dc1 Dubcova3 ecology2  FEM_3D_thermal2 G2_circuit Goodwin_095 matrix-new_3 offshore para-10  parabolic_fem  ss1  stomach thermomech_TK tmt_unsym  xenon2"

Matrices2="af_shell10 af_shell2 atmosmodd atmosmodl cage13 CurlCurl_2 dielFilterV2real  Geo_1438 Hook_1498  ML_Laplace nlpkkt80  Serena  Si87H76  StocF-1465  Transport"


# Iterate the string variable using for loop

for val in $Matrices1; do
        echo $val
        ./test_pardiso ../sparse_multiprecision_bench/matrices/moderate/$val.mtx
        echo ""
        echo "********************************************"
        echo ""
done


for val in $Matrices2; do
        echo $val
        ./test_pardiso  ../sparse_multiprecision_bench/matrices/big/$val.mtx
        echo ""
        echo "********************************************"
        echo ""
done

