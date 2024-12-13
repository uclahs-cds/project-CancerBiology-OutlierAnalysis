#!/bin/bash
#SBATCH --job-name=outseekr_hatzis
#SBATCH --partition=F72
#SBATCH --exclusive

date
Rscript /hot/code/jlivingstone/GitHub/uclahs-cds/project-CancerBiology-OutlierAnalysis/test_OutSeekR/run_outseekr_datasets.R \
    --dataset Hatzis.1 \
    --workers 30
date
