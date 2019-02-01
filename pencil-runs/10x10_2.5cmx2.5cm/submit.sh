#!/bin/bash
# Standard output and error:
#SBATCH -o ./10x10_2.5x2.5_br.out.%j
#SBATCH -e ./10x10_2.5x2.5_br.err.%j
#SBATCH -D ./
#SBATCH -J 10x10_2.5x2.5_br
#SBATCH --nodes=8
#SBATCH --tasks-per-node=32
#SBATCH --cpus-per-task=1
#SBATCH --mail-type=end
#SBATCH --mail-user=carpenter@mpia.de
# Wall clock limit
#SBATCH --time=24:00:00

pc_run -f $PENCIL_HOME/config/hosts/isaac/isaac1.bc.rzg.mpg.de.conf
