#!/bin/bash
#SBATCH --time=168:00:00
#SBATCH --mem=10000M
#SBATCH --account=scjp1
#SBATCH --output=%u-%x-%j.log
#SBATCH --error=%u-%x-%j.err
#SBATCH --mail-user=ccrober@umich.edu
#SBATCH --mail-type=END,FAIL
#SBATCH --signal=B:TERM@60

module load openjdk/18.0.1.1
module load singularity/4.1.3



