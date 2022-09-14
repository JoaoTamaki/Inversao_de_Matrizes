#!/bin/bash

TAM="32 33 64 65 128 129 256 257 512 1000"
METRICA="L3 L2CACHE FLOPS_DP"

for k in $TAM
do
    echo ${k}
    for j in $METRICA
    do
        likwid-perfctr -C 3 -g ${j} -m ${i} ./invmat -r ${k} -i 10 > OP1_${k}_${j}_Otimiz.log
    done
done