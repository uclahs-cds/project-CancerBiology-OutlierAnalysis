## bedr_search_plots.R #############################################################
# Description
# overlap SV breakpoints with GeneHancer-v5.18 regions

### HISTORY #########################################################################
# Version	Date		Developer	Comments
# 0.01		2023-11-08	jlivingstone	initial code

#####################################################################################
library(BoutrosLab.plotting.general)
library(BoutrosLab.utilities)

if (dataset == 'icgc') {
	overlap.matrix <- read.delim(
		file = '/hot/user/jlivingstone/outlier/2023-11-07_outlier_GeneHancer_ICGC_segment_overlap.txt',
		as.is = TRUE
		)
	filename <- generate.filename('outlier', 'ICGC_overlap_recurrcent_histogram', 'png')

	patient.ylimit <- c(0, 35)
	patient.yat <- seq(0, 35, 10)
	patient.xlimit <- c(0, 120000)
	patient.xat <- seq(0, 120000, 10000)

	gh.ylimit <- c(0, 20000)
	gh.yat <- seq(0, 20000, 2000)
	gh.xlimit <- c(0, 250)
	gh.xat <- seq(0, 250, 50)
	gh.breaks <- 20
	}
if (dataset == 'cpcgene') {
	overlap.matrix <- read.delim(
		file = '/hot/user/jlivingstone/outlier/2023-11-03_outlier_GeneHancer_CPCG_OS_overlap.txt',
		as.is = TRUE
		)
	filename <- generate.filename('outlier', 'CPCG_overlap_recurrcent_histogram', 'png')

	gh.ylimit <- c(0, 100000)
	gh.yat <- seq(0, 100000, 10000)
	gh.xlimit <- c(0, 150)
	gh.xat <- seq(0, 150, 50)
	gh.breaks <- 150
	}

patient <- data.frame(
	number = colSums(overlap.matrix)
	)
ghid <- data.frame(
	number = rowSums(overlap.matrix, na.rm = TRUE)
	)

# number of enhancers that overlap with deletions
patient.plot <- create.histogram(
	x = patient$number,
	filename = generate.filename('outlier', 'overlap_patient_recurrence', 'png'),
	xaxis.cex = 1,
	yaxis.cex = 1,
	xlab.cex = 1,
	ylab.cex = 1,
	xlab.label = 'Enhancers that overlap deletions',
	ylab.label = 'Count',
	ylimit = patient.ylimit,
	yat = patient.yat,
	xlimit = patient.xlimit,
	xat = patient.xat,
	xaxis.rot = 45,
	breaks = 20,
	type = 'count',
	resolution = 200
	)

ghid.plot <- create.histogram(
	x = ghid$number,
	filename = generate.filename('outlier', 'overlap_ghid_recurrence', 'png'),
	xaxis.cex = 1,
	yaxis.cex = 1,
	xlab.cex = 1,
	ylab.cex = 1,
	xlab.label = 'Per enhancer number of patients with deletion',
	ylab.label = 'Count',
	ylimit = gh.ylimit,
	yat = gh.yat,
	xlimit = gh.xlimit,
	xat = gh.xat,
	xaxis.rot = 45,
	breaks = gh.breaks,
	type = 'count',
	resolution = 200
	)
	
create.multipanelplot(
	plot.objects = list(patient.plot, ghid.plot),
	filename = filename,
	resolution = 200,
	layout.width = 2,
	layout.height = 1,
	x.spacing = 0.5,
	ylab.axis.padding = 6.5,
	width = 12,
	height = 8
	)
