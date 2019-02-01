#!/bin/bash
# Standard output and error:
#SBATCH -o ./1x1_br.out.%j
#SBATCH -e ./1x1_br.err.%j
#SBATCH -D ./
#SBATCH -J 1x1_br
#SBATCH --nodes=4
#SBATCH --tasks-per-node=32
#SBATCH --cpus-per-task=1
#SBATCH --mail-type=end
#SBATCH --mail-user=carpenter@mpia.de
# Wall clock limit
#SBATCH --time=12:00:00

pc_run -f $PENCIL_HOME/config/hosts/isaac/isaac1.bc.rzg.mpg.de.conf
