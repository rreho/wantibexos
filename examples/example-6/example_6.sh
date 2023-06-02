#!/usr/bin/env bash

nthreads=10
wtbexec="../../bin/wtb.x"

mkdir ./out/

cat > input_mos2_opt.dat << EOF
NTHREADS= $nthreads
SYSDIM= "2D"
DFT= "V"

OUTPUT= "./out/"
CALC_DATA= "./out/"
PARAMS_FILE= "tb_mos2.dat"                                                             
PARAMS_MMN= "hBN.mmn"                                                    
 

MESH_TYPE= "RK2D"
RK= 36

BSE= T
SPEC= T
DTDIAG= T

BSE_BND= T
BSE_WF= T
BERRY_EXC= T
EXC_WF_I= 1
EXC_WF_F= 1
KPATH_FILE= "tmd-kpoints.dat"  
KPATH_BSE= "tmd-kpoints-bse.dat" 

COULOMB_POT= V2DT2
NBANDSC= 1
NBANDSV= 1
LC= 8.00

CSHIFT= 0.08
ENSPECI= 0.0
ENSPECF= 4.0
	
EOF



$wtbexec < input_mos2_opt.dat





