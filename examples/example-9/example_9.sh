#!/usr/bin/env bash

nthreads=10
wtbexec="../../bin/wtb.x"

mkdir ./out/

cat > input_gaas_opt.dat << EOF
NTHREADS= 1
 
SYSDIM= "3D"
DFT= "V"

OUTPUT= "./out/"
CALC_DATA= "./out/"
PARAMS_FILE= "GaAs_tb.dat"
KPATH_FILE= "gaas-kpoints.dat"
KPATH_BSE= "gaas-kpoints-bse.dat" 


MESH_TYPE= "RK3D"
RK= 80

BANDS= T
BSE= T
SPEC= T
DTDIAG= T

COULOMB_POT= V3D
NBANDSC= 1
NBANDSV= 3

CSHIFT= 0.08
ENSPECI= 0.0
ENSPECF= 4.0

BSE_BND= F
BSE_WF= F
BERRY_EXC= F
EXC_WF_I= 1
EXC_WF_F= 1

      	
EOF



$wtbexec < input_gaas_opt.dat





