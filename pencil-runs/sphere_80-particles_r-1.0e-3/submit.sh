#!/bin/bash
# Standard output and error:
#SBATCH -o ./s_80-p_r-1.0e-3.out.%j
#SBATCH -e ./s_80-p_r-1.0e-3.err.%j
#SBATCH -D ./
#SBATCH -J s_80-p_r-1.0e-3
#SBATCH --nodes=2
#SBATCH --tasks-per-node=40
#SBATCH --cpus-per-task=1
#SBATCH --mail-type=end
#SBATCH --mail-user=carpenter@mpia.de
# Wall clock limit
#SBATCH --time=08:00:00

pc_run -f $PENCIL_HOME/config/hosts/isaac/isaac1.bc.rzg.mpg.de.conf
