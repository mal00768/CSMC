#!/bin/bash

#MOAB -N test
#MOAB -l nodes=4:ppn=4
#MOAB -j oe
#MOAB -l walltime=4:00:00
#MOAB -q backfill

module load gnu-openmpi
cd $PBS_O_WORKDIR
mpirun ./a.out 2levels.inpt

