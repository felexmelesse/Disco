#!/bin/bash

make

# NR=(0016 0032 0064 0128 0256)
# NR=(0032 0048 0064 0096 0128)
# NR=(0096 0128 0192 0384 0256 0512 0768 1024 2048)
NR=(4096 2048 1024 512)

for nr in "${NR[@]}"
do
    if [ "$(uname)" == "Darwin" ]; then
        sed -e "s/^Num_R[[:blank:]].*$/Num_R ${nr}/" -i '' in.par
        sed -e "s/^Num_Checkpoints[[:blank:]].*$/Num_Checkpoints 0/" -i '' in.par
    else
        sed -i "s/^Num_R\s.*$/Num_R ${nr}/" in.par
        sed -i "s/^Num_Checkpoints\s.*$/Num_Checkpoints 0/" in.par
    fi

    mpiexec -np 12 disco
    mv output.h5 output.$nr.h5
done
