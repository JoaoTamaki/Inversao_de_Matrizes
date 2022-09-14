#!/bin/bash

TAM="10"
METRICA="L3 L2CACHE FLOPS_DP"

mkdir Output_Trab1_ICC
for k in $TAM
do
    for j in $METRICA
    do
        likwid-perfctr -C 3 -g ${j} ./invmat -r ${k} -i 10 > Output_Trab1_ICC/OP1_${k}_${j}_SemOtimiz.log
    done
done