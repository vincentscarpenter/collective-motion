#!/bin/bash
# Standard output and error:
#SBATCH -o ./80-p_r-1.0e-3_64x64x260_pr.out.%j
#SBATCH -e ./80-p_r-1.0e-3_64x64x260_pr.err.%j
#SBATCH -D ./
#SBATCH -J 80-p_r-1.0e-3_64x64x260_pr
#SBATCH --nodes=4
#SBATCH --tasks-per-node=40
#SBATCH --cpus-per-task=1
#SBATCH --mail-type=all
#SBATCH --mail-user=carpenter@mpia.de
# Wall clock limit
#SBATCH --time=06:00:00

pc_run -f $PENCIL_HOME/config/hosts/isaac/isaac1.bc.rzg.mpg.de.conf
