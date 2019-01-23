#!/bin/bash -l 
#SBATCH -A C3SE2018-1-20
#SBATCH -N 1 
#SBATCH -n 16 
#SBATCH -J meneco_yli 
#SBATCH -t 01:00:00 
 
export USR_DIR=/c3se/NOBACKUP/groups/c3-c3se605-15-2 
module load intel Python 
 
cd $USR_DIR/rhto/
 
meneco.py -d rhto.xml -r yeastGEM_820.xml -s menecoSeeds.sbml -t menecoTargets.sbml > meneco.out