#!/bin/bash -l 
#SBATCH -A C3SE2018-1-20
#SBATCH -N 1 
#SBATCH -n 16 
#SBATCH -J meneco_yli 
#SBATCH -t 01:00:00 
 
export USR_DIR=/c3se/NOBACKUP/groups/c3-c3se605-15-2 
module load intel Python 
 
cd $USR_DIR/rhto/
 
python /c3se/users/eduardk/Hebbe/pyasp-1.4.2/meneco-1.5.0/meneco.py -d ../../scrap/rhtoDraft.xml -r ../../ComplementaryData/yeastGEM_820.xml -s ../../ComplementaryData/reconstruction/menecoSeeds.sbml -t ../../ComplementaryData/reconstruction/menecoTargets.sbml > ../../ComplementaryData/reconstruction/meneco.out