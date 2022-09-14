#!/bin/bash

TAM="10"
METRICA="L3 L2CACHE FLOPS_DP"

for k in $TAM
do
    for j in $METRICA
    do
        likwid-perfctr -C 3 -g ${j} -m ${i} ./invmat -r ${k} -i 10 > OP1_${k}_${j}_Otimiz.log
    done
done