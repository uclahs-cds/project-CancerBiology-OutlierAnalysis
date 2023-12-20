#!/bin/bash

#SBATCH -J calculate_pvalue%A_%a
#SBATCH -o calculate_pvalue%A_%a.out
#SBATCH -e calculate_pvalue%A_%a.err
#SBATCH -p F32
#SBATCH -c 20
#SBATCH --mem=32G
#SBATCH --array=1-18

source /etc/profile.d/module.sh
module load R

dataset='BRCA_EU'
working_directory='/hot/users/jlivingstone/outlier/run_method/remove_four_patients'
rank_file='/hot/user/jlivingstone/outlier/run_method/remove_four_patients/2023-12-19_Final_outlier_rank_bic_BRCA_EU.4.short.rda'
combined_file='/hot/user/jlivingstone/outlier/run_method/remove_four_patients/2023-12-20_Simulated_Data_5method_combine_BRCA_EU.4.short.rda'
patients_to_remove=4

Rscript /hot/code/jlivingstone/GitHub/uclahs-cds/project-CancerBiology-OutlierAnalysis/OutlierDetectionAlgorithm/7.Significant_Outlier_Detection.R \
--dataset.name $dataset \
--working.directory $working_directory \
--rank.file $rank_file \
--combined.file $combined_file \
--row.chunks 1000 \
--matrix.chunk ${SLURM_ARRAY_TASK_ID} \
--patients.to.remove $patients_to_remove
