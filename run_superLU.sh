#!/bin/bash
# Define a string variable with a value

Matrices1="2cubes_sphere  ASIC_320ks Baumann cfd2  dc1 crashbasis ct20stif  Dubcova3 ecology2 FEM_3D_thermal2 G2_circuit Goodwin_095 offshore para-10  parabolic_fem  stomach thermomech_TK tmt_unsym  xenon2"

Matrices2="af_shell10 af_shell2 CurlCurl_2"


# Iterate the string variable using for loop

for val in $Matrices1; do
        echo $val
        ./test_superLU_MT ../sparse_multiprecision_bench/matrices/moderate/$val.mtx
        echo ""
        echo "********************************************"
        echo ""
done


for val in $Matrices2; do
        echo $val
        ./test_superLU_MT  ../sparse_multiprecision_bench/matrices/big/$val.mtx
        echo ""
        echo "********************************************"
        echo ""
done

