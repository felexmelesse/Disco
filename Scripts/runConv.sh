#!/bin/bash

make

NR=(0016 0032 0064 0128 0256 0512)

for nr in "${NR[@]}"
do
    sed -i "s/^Num_R\s.*$/Num_R ${nr}/" in.par
    sed -i "s/^Num_Checkpoints\s.*$/Num_Checkpoints 0/" in.par
    if [ "$nr" -lt 100 ]; then
        ./disco
    else
        mpiexec -np 3 ./disco
    fi
    mv output.h5 output.$nr.h5
done

python Python/alfvenwaveAnalysis.py output.*.h5
