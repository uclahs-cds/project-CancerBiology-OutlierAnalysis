#!/bin/bash

#SBATCH -J run_outlier%A_%a
#SBATCH -o run_outlier%A_%a.out
#SBATCH -e run_outlier%A_%a.err
#SBATCH -p F16
#SBATCH -c 8
#SBATCH --mem=16G

source /etc/profile.d/module.sh
module load R

dataset='BRCA_EU'
working_directory='/hot/users/jlivingstone/outlier/run_method/remove_four_patients'
matrix_file='/hot/users/jlivingstone/outlier/NikZainal_2016/original/SupplementaryTable7Transcriptomic342.txt'
patients_to_remove=4

# initial run with all patients - variables are from export.file
Rscript /hot/code/jlivingstone/GitHub/uclahs-cds/project-CancerBiology-OutlierAnalysis/OutlierDetectionAlgorithm/1.Outlier_Detection.R \
--dataset.name $dataset \
--working.directory $working_directory \
--data.matrix.file $matrix_file \
--patients.to.remove $patients_to_remove
