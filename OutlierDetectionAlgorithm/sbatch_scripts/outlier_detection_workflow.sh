#!/bin/bash

# fill in variable names in script FIRST
sbatch /hot/code/jlivingstone/GitHub/uclahs-cds/project-CancerBiology-OutlierAnalysis/OutlierDetectionAlgorithm/sbatch_scripts/run_outlier_detection.sbatch 

# Scripts 2., 3., & 4. only need to be run once to generate simulated data
Rscript /hot/code/jlivingstone/GitHub/uclahs-cds/project-CancerBiology-OutlierAnalysis/OutlierDetectionAlgorithm/2.Distribution_Identification.R \
--dataset.name BRCA_EU \
--working.directory /hot/user/jlivingstone/outlier/run_method \
--outlier.rank.file /hot/users/jlivingstone/outlier/run_method/2023-11-20_BRCA-EU_final_outlier_rank_bic.long.rda

Rscript /hot/code/jlivingstone/GitHub/uclahs-cds/project-CancerBiology-OutlierAnalysis/OutlierDetectionAlgorithm/3.Simulated_Data_generation_1.R \
--dataset.name BRCA_EU \
--working.directory /hot/users/jlivingstone/outlier/run_method \
--outlier.rank.file /hot/users/jlivingstone/outlier/run_method/2023-11-20_BRCA-EU_final_outlier_rank_bic.long.rda \
--ntimes 10

# fill in variable names in script FIRST
sbatch /hot/code/jlivingstone/GitHub/uclahs-cds/project-CancerBiology-OutlierAnalysis/OutlierDetectionAlgorithm/sbatch_scripts/simulated_data_generation.sbatch

# fill in variable names in script FIRST
sbatch /hot/code/jlivingstone/GitHub/uclahs-cds/project-CancerBiology-OutlierAnalysis/OutlierDetectionAlgorithm/sbatch_scripts/simulated_run_methods.sbatch

Rscript /hot/code/jlivingstone/GitHub/uclahs-cds/project-CancerBiology-OutlierAnalysis/OutlierDetectionAlgorithm/6.Simulated_Data_5method_combine.R \
--dataset.name BRCA_EU \
--working.directory /hot/user/jlivingstone/outlier/run_method/remove_four_patients \
--file.date 2023-12-19 \
--patients.to.remove 4

# fill in variable names in script FIRST
# array values will depend on chunks chosen (nrow(expression.matrix) / matrix.chunks)
sbatch /hot/code/jlivingstone/GitHub/uclahs-cds/project-CancerBiology-OutlierAnalysis/OutlierDetectionAlgorithm/sbatch_scripts/calculate_pvalue.sbatch

Rscript /hot/code/jlivingstone/GitHub/uclahs-cds/project-CancerBiology-OutlierAnalysis/OutlierDetectionAlgorithm/8.Significant_Outlier_Pvalue_Calculation.R \
--dataset.name BRCA_EU \
--working.directory /hot/user/jlivingstone/outlier/run_method/remove_four_patients \
--patients.to.remove 4

Rscript /hot/code/jlivingstone/GitHub/uclahs-cds/project-CancerBiology-OutlierAnalysis/OutlierDetectionAlgorithm/9.Patient_Outlier_Identification.R \
--dataset.name BRCA_EU \
--working.directory /hot/users/jlivingstone/outlier/run_method \
--outlier.rank.file /hot/users/jlivingstone/outlier/run_method/2023-11-20_BRCA-EU_final_outlier_rank_bic.short.rda \
--qvalue.cutoff 0.01
