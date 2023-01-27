#!/bin/bash

# Script to test program performance (parameters)

program=./main

inputFile=res.csv
outputFile=out.csv

append=0

arrThreads=(1 2 4 8 16 32)

for threads in "${arrThreads[@]}"
do
	for i in {1..20}
	do
		# From 5 to 100 size package
		pack=$(( i * 5 ))
		
		echo "Threads = ${threads}/32;  Package Size = ${pack}/100"

		# Run program (first time append=0)
		OMP_NUM_THREADS=${threads} ${program} ${inputFile} ${outputFile} ${pack} ${append} >> log.out

		# Remove output
		rm ${outputFile}

		append=1

	done
done
