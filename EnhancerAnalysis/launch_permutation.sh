#!/bin/bash
#SBATCH --partition=F16
#SBATCH --cpus-per-task=1
#SBATCH --nodes=1
#SBATCH --mem=2G
#SBATCH --array=0-4

date
seed_start=(32000 34000 36000 38000 40000)
seed_end=(33999 35999 37999 39999 41999)
Rscript /hot/code/jlivingstone/GitHub/uclahs-cds/project-CancerBiology-OutlierAnalysis/EnhancerAnalysis/permutation_test.R \
        --seed.start ${seed_start[$SLURM_ARRAY_TASK_ID]} \
        --seed.end ${seed_end[$SLURM_ARRAY_TASK_ID]} \
        --elite 0
date
