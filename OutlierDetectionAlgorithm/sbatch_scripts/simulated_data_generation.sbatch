#!/bin/bash

#SBATCH -J simulate_data%A_%a
#SBATCH -o simulate_data%A_%a.out
#SBATCH -e simulate_data%A_%a.err
#SBATCH -p F32
#SBATCH -c 10
#SBATCH --mem=20G
#SBATCH --array=1-10

source /etc/profile.d/module.sh
module load R

dataset='BRCA_EU'
working_directory='/hot/users/jlivingstone/outlier/run_method'
distrbution_identification='/hot/users/jlivingstone/outlier/run_method/2023-11-20_Distribution_Identification_short_BRCA_EU.rda'

Rscript /hot/users/jlivingstone/bl-github/project-CancerBiology-OutlierAnalysis/OutlierDetectionAlgorithm/4.Simulated_Data_generation_2.R \
--dataset.name $dataset \
--working.directory $working_directory \
--distribution.identification.file $distrbution_identification \
--simulated.data.file /hot/users/jlivingstone/outlier/run_method/2023-11-21_Simulated_data_generation_1_BRCA_EU.${SLURM_ARRAY_TASK_ID}.rda
