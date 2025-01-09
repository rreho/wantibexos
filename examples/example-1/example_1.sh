#!/usr/bin/env bash

orbw="c1 c1-pz c1-sp2 c2 c-pz"
orbw=($orbw)

nthreads=4
wtbexec="../../bin/wtb.x"

#c1 - carbon 1 all orbitals contribution
#c1-pz - carbon 1 pz orbitals contribution
#c1-sp2 - carbon 1 sp2 orbitals contribution
#c2 - carbon 2 all orbitals contribution (pz only)
#c-pz carbon 1 + carbon 2 pz orbitals contribution


cat > c1.dat << EOF
1
1
1
1
0
EOF

cat > c1-pz.dat << EOF
0
0
0
1
0
EOF

cat > c1-sp2.dat << EOF
1
1
1
0
0
EOF

cat > c2.dat << EOF
0
0
0
0
1
EOF

cat > cpz.dat << EOF
0
0
0
1
1
EOF


cat > tmd-kpoints.dat << EOF
6
50
0.0 0.0 0.0               !Gamma
0.6667  -0.3333  0.0      !K
0.6667  -0.3333  0.0      !K
0.3333   0.3333 0.0       !K'
0.3333   0.3333 0.0       !K'
0.0 0.0 0.0               !Gamma
EOF

mkdir ./out

for (( i = 0; i<${#orbw[@]}; i++ )); do	

orbwl="${orbw[$i]}.dat"

mkdir ./out/$orbwl

cat > input-$orbwl.dat << EOF
NTHREADS= $nthreads
SYSDIM= "2D"
DFT= "V"

OUTPUT= "/out/"
!OUTPUT= "./out/$orbwl/"
CALC_DATA= "./out/"
PARAMS_FILE= "tb_graphene.dat"                                                             
KPATH_FILE= "tmd-kpoints.dat"                                                                                                                 
KPATH_BSE= "tmd-kpoints-bse.dat"

!ORB_W= "$orbwl" 

BANDS= T

BSE= T
SPEC= T
DTDIAG= T
NGX=6
NGY=6
NGZ=1
BSE_BND= T
BSE_WF= T
BERRY_EXC= F
EXC_WF_I= 1
EXC_WF_F= 1

COULOMB_POT= V2DT2
NBANDSC= 1
NBANDSV= 1
LC= 8.00

CSHIFT= 0.08
ENSPECI= 0.0
ENSPECF= 4.0
	
EOF


$wtbexec < input-$orbwl.dat



echo $i/${#orbw[@]}

done

