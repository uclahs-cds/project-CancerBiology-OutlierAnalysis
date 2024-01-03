### check_sample_id.R ##############################################################
# Description
# Need to check to see that samples from Nik-Zanial 2016 match those in ICGC-EU/UK

### HISTORY ########################################################################
# Version	Date		Developer	Comments
# 0.01		2023-10-31	jlivingstone	initial code

### PREAMBLE ####################################################################
library(BoutrosLab.plotting.general);

setwd('/hot/users/jlivingstone/outlier')

# read in Supplementary Table 1 from Nik-Zanial 2016
paper <- read.delim(
	file = file.path('/hot/users/jlivingstone/outlier/NikZainal_2016/original/', 'SupplementaryTable1.CLINICAL.PATHOLOGY.DATA.FREEZE.ANALYSIS.v4.032015.csv'),
	as.is = TRUE,
	skip = 1,
	sep = ','
	)
# n = 560
paper.samples <- paper$sample_name

# submitted_sample_id + submitted_sample_id
eu <- read.delim(
	file = file.path('/hot/ref/cohort/ICGC/BRCA/EU/original/', 'sample.BRCA-EU.tsv.gz'),
	as.is = TRUE
	)
length(intersect(paper.samples, unique(substr(eu$submitted_sample_id, 1, 6))))
# 264!
eu.samples <- intersect(paper.samples, unique(substr(eu$submitted_sample_id, 1, 6)))

uk <- read.delim(
	file = file.path('/hot/ref/cohort/ICGC/BRCA/UK/original/', 'sample.BRCA-UK.tsv.gz'),
	as.is = TRUE
	)
length(intersect(paper.samples, unique(substr(uk$submitted_sample_id, 1, 6))))
# 45
uk.samples <- intersect(paper.samples, unique(substr(uk$submitted_sample_id, 1, 6)))

# intersect(eu.samples, uk.samples)
# total overlap - I wonder if the mutation calls are the same ...
# use EU samples to begin with

# specimen_type + submitted_donor_id
eu.specimen <- read.delim(
	file = file.path('/hot/ref/cohort/ICGC/BRCA/EU/original/', 'specimen.BRCA-EU.tsv.gz'),
	as.is = TRUE
	)

# merge eu & eu.specimen and reduce to primary samples
eu.specimen.order <- eu.specimen[match(eu$submitted_specimen_id, eu.specimen$submitted_specimen_id),]
eu.merged <- cbind(eu, eu.specimen.order)

eu.primary <- eu.merged[-grep('Normal', eu.merged$specimen_type),]

# clinical + submitted_donor_id
eu.donor <- read.delim(
	file = file.path('/hot/ref/cohort/ICGC/BRCA/EU/original/', 'donor.BRCA-EU.tsv.gz'),
	as.is = TRUE
	
# now merge eu.primary and eu.donor to have both clinical and PD id
eu.donor.order <- eu.donor[match(eu.primary$submitted_donor_id, eu.donor$submitted_donor_id), ]
eu.all <- cbind(eu.primary, eu.donor.order)
eu.all$sample_name <- substr(eu.all$submitted_sample_id, 1, 6)

# now match paper clinical and icgc clinical, n = 264
samples <- intersect(eu.all$sample_name, paper$sample_name)
paper.clinical <- paper[match(samples, paper$sample_name), ]
eu.clinical <- eu.all[match(samples, eu.all$sample_name), ]

# check age
age.check <- data.frame(
	paper.id = paper.clinical$sample_name,
	paper.age = paper.clinical$donor_age_at_diagnosis,
	icgc.id = eu.clinical$submitted_sample_id,
	icgc.age = eu.clinical$donor_age_at_diagnosis,
	stringsAsFactors = FALSE
	)

table(age.check$paper.age == age.check$icgc.age, useNA = 'ifany')
#FALSE  TRUE  <NA> 
#   12   244     8 
# FALSE is due to 'over_80' classification, two patients are off by 1 year

# check stage
stage.check <- data.frame(
	paper.id = paper.clinical$sample_name,
	paper.stage = paper.clinical$T_stage,
	paper.met = paper.clinical$M_stage,
	paper.node = paper.clinical$N_stage,
	icgc.id = eu.clinical$submitted_sample_id,
	icgc.stage = substr(eu.clinical$tumour_stage, 1, 2),
	icgc.node = substr(eu.clinical$tumour_stage, 3, 4),
	icgc.met = substr(eu.clinical$tumour_stage, 5, 6),
	stringsAsFactors = FALSE
	)

table(stage.check$paper.stage == stage.check$icgc.stage, useNA = 'ifany')
# FALSE  TRUE 
#   80   184 
table(stage.check$paper.node == stage.check$icgc.node, useNA = 'ifany')
# due to missing data in icgc data 'Tx|Mx|Nx' vs ''
table(stage.check$paper.met == stage.check$icgc.met, useNA = 'ifany')
# one 'extra' difference due to MX vs Mx

grade.check <- data.frame(
	paper.id = paper.clinical$sample_name,
	paper.grade = c('G1', 'G2', 'G3', '')[match(paper.clinical$tumour_grade, c('I', 'II', 'III', 'no_data_supplied'))],
	icgc.id = eu.clinical$submitted_sample_id,
	icgc.grade = eu.clinical$tumour_grade,
	stringsAsFactors = FALSE
	)
table(grade.check$paper.grade == grade.check$icgc.grade)
# TRUE 
# 264 

## Conclusion, they are the same samples ##
