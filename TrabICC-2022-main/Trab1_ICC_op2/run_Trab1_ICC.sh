#!/bin/bash

TAM="32 33 64 65 128 129 256 257 512 1000"
METRICA="L3 L2CACHE FLOPS_DP"

mkdir Output_Trab1_ICC
for k in $TAM
do
    echo ${k}
    for j in $METRICA
    do
        likwid-perfctr -C 3 -g ${j} ./invmat -r ${k} -i 10 > Output_Trab1_ICC/OP2_${k}_${j}_SemOtimiz.log
    done
done