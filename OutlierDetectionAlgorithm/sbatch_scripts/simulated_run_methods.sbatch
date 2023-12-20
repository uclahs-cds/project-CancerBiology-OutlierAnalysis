#!/bin/bash

#SBATCH -J run_simulate%A_%a
#SBATCH -o run_simulate%A_%a.out
#SBATCH -e run_simulate%A_%a.err
#SBATCH -p F16
#SBATCH -c 8
#SBATCH --mem=16G
#SBATCH --array=3-10

source /etc/profile.d/module.sh
module load R

dataset='BRCA_EU'
working_directory='/hot/users/jlivingstone/outlier/run_method/remove_four_patients'
distribution_file='/hot/users/jlivingstone/outlier/run_method/2023-11-20_Distribution_Identification_short_BRCA_EU.rda'
simulated_data_file="/hot/users/jlivingstone/outlier/run_method/2023-11-22_Simulated_data_generation_2_BRCA_EU.${SLURM_ARRAY_TASK_ID}.rda"
patients_to_remove=4

Rscript /hot/code/jlivingstone/GitHub/uclahs-cds/project-CancerBiology-OutlierAnalysis/OutlierDetectionAlgorithm/5.Simulated_Data_5method.R \
--dataset.name $dataset \
--working.directory $working_directory \
--distribution.identification.file $distribution_file \
--simulated.data.file $simulated_data_file \
--patients.to.remove $patients_to_remove
