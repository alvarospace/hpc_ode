#!/bin/bash

program=./main

input=("input/res_gri_2000.csv" "input/res_est_2000.csv" "input/res_h2_2000.csv")
mechanisms=("gri30.yaml" "mechanism/ESTiMatE-Mech_04042022.xml" "mechanism/h2_li_2004.xml")
output=out.csv

append=0

arrThreads=(1 2 4 8 14 28)


for j in {0..2}
do
    for threads in "${arrThreads[@]}"
    do
        for i in 1
        do
            # From 50 to 1000 size package
            pack=$(( i * 65 ))

            echo "Mechanism = ${input[$j]}; Threads = ${threads};  Package Size = ${pack}"

            # Run program (first time append=0)
            OMP_NUM_THREADS=${threads} ${program} ${input[${j}]} ${output} ${mechanisms[${j}]} ${pack} ${append} >> log.out

            # Remove output
            rm ${output}

            append=1

        done
    done
done
