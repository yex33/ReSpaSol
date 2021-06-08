#!/bin/bash

#Download matrices
wget https://suitesparse-collection-website.herokuapp.com/MM/Um/2cubes_sphere.tar.gz
wget https://suitesparse-collection-website.herokuapp.com/MM/Sandia/ASIC_320ks.tar.gz
wget https://suitesparse-collection-website.herokuapp.com/MM/Watson/Baumann.tar.gz
wget https://suitesparse-collection-website.herokuapp.com/MM/Rothberg/cfd2.tar.gz
wget https://suitesparse-collection-website.herokuapp.com/MM/QLi/crashbasis.tar.gz
wget https://suitesparse-collection-website.herokuapp.com/MM/Boeing/ct20stif.tar.gz
wget https://suitesparse-collection-website.herokuapp.com/MM/IBM_EDA/dc1.tar.gz
wget https://suitesparse-collection-website.herokuapp.com/MM/UTEP/Dubcova3.tar.gz
wget https://suitesparse-collection-website.herokuapp.com/MM/McRae/ecology2.tar.gz
wget https://suitesparse-collection-website.herokuapp.com/MM/Botonakis/FEM_3D_thermal2.tar.gz
wget https://suitesparse-collection-website.herokuapp.com/MM/AMD/G2_circuit.tar.gz
wget https://suitesparse-collection-website.herokuapp.com/MM/Goodwin/Goodwin_095.tar.gz
wget https://suitesparse-collection-website.herokuapp.com/MM/Schenk_IBMSDS/matrix-new_3.tar.gz
wget https://suitesparse-collection-website.herokuapp.com/MM/Um/offshore.tar.gz
wget https://suitesparse-collection-website.herokuapp.com/MM/Schenk_ISEI/para-10.tar.gz
wget https://suitesparse-collection-website.herokuapp.com/MM/Wissgott/parabolic_fem.tar.gz
wget https://suitesparse-collection-website.herokuapp.com/MM/VLSI/ss1.tar.gz
wget https://suitesparse-collection-website.herokuapp.com/MM/Norris/stomach.tar.gz
wget https://suitesparse-collection-website.herokuapp.com/MM/Botonakis/thermomech_TK.tar.gz
wget https://suitesparse-collection-website.herokuapp.com/MM/CEMW/tmt_unsym.tar.gz
wget https://suitesparse-collection-website.herokuapp.com/MM/Ronis/xenon2.tar.gz

#Unzip the files
tar -xvf 2cubes_sphere.tar.gz
tar -xvf ASIC_320ks.tar.gz
tar -xvf Baumann.tar.gz
tar -xvf cfd2.tar.gz
tar -xvf crashbasis.tar.gz
tar -xvf ct20stif.tar.gz
tar -xvf dc1.tar.gz
tar -xvf Dubcova3.tar.gz
tar -xvf ecology2.tar.gz
tar -xvf FEM_3D_thermal2.tar.gz
tar -xvf G2_circuit.tar.gz
tar -xvf Goodwin_095.tar.gz
tar -xvf matrix-new_3.tar.gz
tar -xvf offshore.tar.gz
tar -xvf para-10.tar.gz
tar -xvf parabolic_fem.tar.gz
tar -xvf ss1.tar.gz
tar -xvf stomach.tar.gz
tar -xvf thermomech_TK.tar.gz
tar -xvf tmt_unsym.tar.gz
tar -xvf xenon2.tar.gz
