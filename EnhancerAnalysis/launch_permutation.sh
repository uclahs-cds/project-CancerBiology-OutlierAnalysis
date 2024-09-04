#!/bin/bash
#SBATCH --partition=F16
#SBATCH --cpus-per-task=1
#SBATCH --nodes=1
#SBATCH --mem=2G
#SBATCH --array=1-9
#SBATCH --nodelist=F16-5

date
seed_start=(90000 91000 92000 93000 94000 95000 96000 97000 98000 99000)
seed_end=(90999 91999 92999 93999 94999 95999 96999 97999 98999 99999)
Rscript /hot/code/jlivingstone/GitHub/uclahs-cds/project-CancerBiology-OutlierAnalysis/EnhancerAnalysis/permutation_test.R \
        --seed.start ${seed_start[$SLURM_ARRAY_TASK_ID]} \
        --seed.end ${seed_end[$SLURM_ARRAY_TASK_ID]} \
        --elite 0
date
