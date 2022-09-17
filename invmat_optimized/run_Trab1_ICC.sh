#!/bin/bash

TRAB="1 2"
OP="1 2"
TAM="32 33 64 65 128 129 256 257 512 1000 2000 4000 6000 10000"
METRICA="L3 L2CACHE FLOPS_DP"
for t in $TRAB
do
        for i in $OP
        do
                for k in $TAM
                do
                        echo ${k}
                        for j in $METRICA
                        do
                                likwid-perfctr -f -C 3 -g ${j} ./Trab${t}_ICC_op${i}/invmat -r ${k} -i 10 > Output_Trab${t}_op${i}_ICC/OP${i}_${k}_${j}_T${t}.log
                        done
                done
        done
done