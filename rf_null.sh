#!/bin/bash

#$ -q long
#$ -pe smp 1

#run permutation tests for betaCV validation in random forest network

module load python/3.5.2
cd ~/network_validation/null_dist

for run in {1..1000}
do
	python Network_Valid_null_rf10000.py
done

