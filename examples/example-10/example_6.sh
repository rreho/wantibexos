#!/usr/bin/env bash

nthreads=1
wtbexec="../../bin/wtb.x"

mkdir ./out/

cat > input_bi2se3_opt.dat << EOF
NTHREADS= $nthreads
SYSDIM= "3D"
DFT= "O"

OUTPUT= "./out/"
CALC_DATA= "./out/"
PARAMS_FILE= "bi2se3_hr.dat"                                                             
PARAMS_MMN="bi2se3.mmn"                                                    
 
NGX= 3
NGY= 3
NGZ= 3

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

COULOMB_POT= V3DL
NBANDSC= 1
NBANDSV= 1
LC= 8.00

CSHIFT= 0.08
ENSPECI= 0.0
ENSPECF= 4.0
	
EOF

$wtbexec < input_bi2se3_opt.dat
