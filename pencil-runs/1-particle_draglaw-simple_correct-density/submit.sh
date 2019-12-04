#!/bin/bash
# Standard output and error:
#SBATCH -o ./1-p_dls_cd.out.%j
#SBATCH -e ./1-p_dls_cd.err.%j
#SBATCH -D ./
#SBATCH -J 1-p_dls_cd
#SBATCH --nodes=2
#SBATCH --tasks-per-node=40
#SBATCH --cpus-per-task=1
#SBATCH --mail-type=end
#SBATCH --mail-user=carpenter@mpia.de
# Wall clock limit
#SBATCH --time=04:00:00

pc_run -f $PENCIL_HOME/config/hosts/isaac/isaac1.bc.rzg.mpg.de.conf
