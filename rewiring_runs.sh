#!/bin/csh

#$ -q long
#$ -pe smp 4

#module load R/3.4.0-intel-17.1
#setenv R_LIBS ~/rlibs

#loop for rewiring networks permutations

set i = 1
while ($i <= 20)
	cd ~/piggyBac_rewire_perm
	Rscript rewiring_perms_DTWMIC_CRC.R
	@ i++
end
