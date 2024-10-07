### HISTORY ######################################################################
# This script analyzes RNA abundance of genes on chromosome 6 in an outlier patient
# with ecDNA. It calculates robust z-scores for the genes and visualizes the results
# using a bar plot.
# Date: 2024-08-13

### PREAMBLE ####################################################################
# Load necessary libraries
library(BoutrosLab.plotting.general);
library(BoutrosLab.utilities);

# Source the helper library
source(here::here('common_functions.R'));

# Load the datafile
load(file.path(get.outlier.data.dir(), '2024-10-03_Figure1_2_3_4_min_input.rda'));

load.multiple.computed.variables(c(
    'brca.outlier.symbol',
    'brca.cnv.chr.new.gis.fpkm.order.match'
    ));

brca.cnv.chr.new.gis.fpkm.order.match.chr6 <- brca.cnv.chr.new.gis.fpkm.order.match[brca.cnv.chr.new.gis.fpkm.order.match.chr$chromosome == 'chr6', ];
brca.cnv.chr.new.gis.fpkm.order.match.chr6.all.gene.location <- brca.cnv.chr.new.gis.fpkm.order.match.chr[brca.cnv.chr.new.gis.fpkm.order.match.chr$chromosome %in% 'chr6', ];

# Filter FPKM data for genes on chromosome 6
fpkm.chr6 <- fpkm.tumor.symbol.filter.brca[
    match(brca.cnv.chr.new.gis.fpkm.order.match.chr6$Hugo_Symbol, fpkm.tumor.symbol.filter.brca$Symbol),
    patient.part.brca
    ];

# Calculate robust z-scores for each gene
fpkm.chr6.median.z <- NULL;
for (i in 1:nrow(fpkm.chr6)) {
    gene.median <- median(as.numeric(fpkm.chr6[i, ]));
    gene.mad <- mad(as.numeric(fpkm.chr6[i, ]));
    robust.z <- (as.numeric(fpkm.chr6[i, ]) - gene.median) / gene.mad;
    fpkm.chr6.median.z <- rbind(fpkm.chr6.median.z, robust.z);
    }

# Set row and column names for the z-score matrix
rownames(fpkm.chr6.median.z) <- rownames(fpkm.chr6);
colnames(fpkm.chr6.median.z) <- colnames(fpkm.chr6);


ecdna.patient <- 'TCGA.A2.A3XX.01A';

# Prepare data for plotting
fpkm.chr6.median.z.bar <- data.frame(
    cbind(fpkm.chr6.median.z[, ecdna.patient, drop = FALSE], sample = c(1:nrow(fpkm.chr6.median.z)))
    );




# Get ecDNA data
outlier.patient.tag.01.brca.patient.sum <- apply(outlier.patient.tag.01.brca, 2, sum);
# Read in the amplicon data
brca.amplicon$sample_barcode <- gsub('-', '.', brca.amplicon$sample_barcode);

# Match amplicon data with outlier patient data
brca.amplicon.match <- brca.amplicon[brca.amplicon$sample_barcode %in% substr(names(outlier.patient.tag.01.brca.patient.sum), 1, 15),];
brca.amplicon.match.ec <- brca.amplicon.match[brca.amplicon.match$amplicon_classification %in% "Circular",];

# Process each amplicon interval
brca.amplicon.match.ec.each <- NULL
for (i in 1:nrow(brca.amplicon.match.ec)) {
    original.df <- brca.amplicon.match.ec[i,];
    split_intervals <- strsplit(original.df$amplicon_intervals, ",");
    num_rows <- sapply(split_intervals, length);
    split_df <- original.df[rep(seq_len(nrow(original.df)), num_rows), ];
    split_df$amplicon_intervals <- unlist(split_intervals);
    row.names(split_df) <- NULL;
    brca.amplicon.match.ec.each <- rbind(brca.amplicon.match.ec.each, split_df);
    }

# Split intervals and create a data frame with chromosome, start, and end
split.intervals <- strsplit(as.character(brca.amplicon.match.ec.each$amplicon_intervals), "[:-]")
brca.amplicon.match.ec.each.chr <- cbind(
    brca.amplicon.match.ec.each[,1:3],
    chr = sapply(split.intervals, "[[", 1),
    start = sapply(split.intervals, "[[", 2),
    end = sapply(split.intervals, "[[", 3)
    );

# Convert hg19 to hg38
library(rtracklayer)
library(liftOver)
library(GenomicRanges)

hg19_location <- brca.amplicon.match.ec.each.chr[,4:6];
hg19_location[,2] <- as.numeric(hg19_location[,2]);
hg19_location[,3] <- as.numeric(hg19_location[,3]);
hg19_location[,1] <- paste("chr", hg19_location[,1], sep = '');
hg19_location.gr <- makeGRangesFromDataFrame(hg19_location);


unpacked.chain <- "hg19ToHg38.over.chain";
if (!file.exists(unpacked.chain)) {
    url2 <- "http://hgdownload.cse.ucsc.edu/goldenpath/hg19/liftOver/hg19ToHg38.over.chain.gz";
    packed.chain <- "hg19ToHg38.over.chain.gz";
    download.file(url2, packed.chain);
    system(paste0("gzip -d ", packed.chain));
}
ch <- import.chain(unpacked.chain);

hg38_location <- liftOver(hg19_location.gr, ch);

hg38_location.reduce <- NULL;
for (i in 1:length(hg38_location)) {
    combined_location <- reduce(hg38_location[[i]]);
    chr <- tryCatch(seqnames(combined_location)[1], error = function(e) NULL);
    start <- tryCatch(start(combined_location)[1], error = function(e) NULL);
    end <- tryCatch(end(combined_location)[length(end(combined_location))], error = function(e) NULL);
    if (!is.null(chr) && !is.null(start) && !is.null(end)) {
        combined_location_str <- paste(chr, start, end, sep = " ");
        } 
    else {
        combined_location_str <- c('0 0 0');
        }
    hg38_location.reduce <- rbind(hg38_location.reduce, combined_location_str);
    }

hg38_location.reduce.split <- strsplit(hg38_location.reduce[,1], " ");
hg38_location.reduce.split.t <- t(data.frame(hg38_location.reduce.split));
rownames(hg38_location.reduce.split.t) <- NULL;
colnames(hg38_location.reduce.split.t) <- c('chr', 'start', 'end');

brca.amplicon.match.ec.each.chr.hg38 <- cbind(brca.amplicon.match.ec.each.chr[,1:3], hg38_location.reduce.split.t);

# Patients only have ecDNA data
outlier.patient.tag.01.brca.ecDNA <- outlier.patient.tag.01.brca[,substr(colnames(outlier.patient.tag.01.brca), 1, 15) %in% unique(brca.amplicon.match.ec.each.chr.hg38$sample_barcode)];
outlier.patient.tag.01.brca.ecDNA.match <- outlier.patient.tag.01.brca.ecDNA[apply(outlier.patient.tag.01.brca.ecDNA, 1, sum) > 0,];


brca.cnv.chr.new.gis.fpkm.order.match.chr.outlier.location <- brca.cnv.chr.new.gis.fpkm.order.match.chr[brca.cnv.chr.new.gis.fpkm.order.match.chr$gene_name %in% brca.outlier.symbol,1:5];
brca.cnv.chr.new.gis.fpkm.order.match.chr.outlier.location.ecdna.match <- brca.cnv.chr.new.gis.fpkm.order.match.chr.outlier.location[brca.cnv.chr.new.gis.fpkm.order.match.chr.outlier.location$gene_name %in% fpkm.tumor.symbol.filter.brca[rownames(outlier.patient.tag.01.brca.ecDNA.match),'Symbol'],];

# Find if the outlier genes are from ecDNA
ecDNA.outlier.patient.location <- NULL;
for (i in 1:nrow(brca.amplicon.match.ec.each.chr.hg38)) {
    target.patient <- brca.amplicon.match.ec.each.chr.hg38$sample_barcode[i];
    target.chr <- brca.amplicon.match.ec.each.chr.hg38$chr[i];
    target.start <- brca.amplicon.match.ec.each.chr.hg38$start[i];
    target.end <- brca.amplicon.match.ec.each.chr.hg38$end[i];
    target.outlier <- rownames(outlier.patient.tag.01.brca.ecDNA.match)[outlier.patient.tag.01.brca.ecDNA.match[,substr(colnames(outlier.patient.tag.01.brca.ecDNA.match), 1, 15) %in% target.patient] == 1];
    target.symbol <- fpkm.tumor.symbol.filter.brca[rownames(fpkm.tumor.symbol.filter.brca) %in% target.outlier, 'Symbol'];
    target.symbol.location <- brca.cnv.chr.new.gis.fpkm.order.match.chr.outlier.location.ecdna.match[brca.cnv.chr.new.gis.fpkm.order.match.chr.outlier.location.ecdna.match$gene_name %in% target.symbol,];
    
    if (nrow(target.symbol.location) > 0) {
        for (j in 1:nrow(target.symbol.location)) {
            if (target.symbol.location$chromosome[j] == target.chr & target.symbol.location$start[j] >= target.start & target.symbol.location$end[j] <= target.end) {
                ecDNA.outlier <- cbind(brca.amplicon.match.ec.each.chr.hg38[i,], target.symbol.location[j,]);
                } 
            else {
                ecDNA.outlier <- NULL
                }
            ecDNA.outlier.patient.location <- rbind(ecDNA.outlier.patient.location, ecDNA.outlier);
            }
        }
    }



# Set bar colors, highlighting outlier genes in red
bar.colours <- rep('black', nrow(fpkm.chr6.median.z.bar));
outlier.gene.indices <- which(
    brca.cnv.chr.new.gis.fpkm.order.match.chr6$Hugo_Symbol %in%
        ecDNA.outlier.patient.location$gene_name[
            ecDNA.outlier.patient.location$sample_barcode %in% substr(ecdna.patient, 1, 15)
            ]
    );
bar.colours[outlier.gene.indices] <- 'red';

# Define ecDNA region start and end genes
chr6.ecdna.start <- as.numeric(ecDNA.outlier.patient.location[
    ecDNA.outlier.patient.location$sample_barcode %in% substr(ecdna.patient, 1, 15),
    ]$start[1]);
chr6.ecdna.end <- as.numeric(ecDNA.outlier.patient.location[
    ecDNA.outlier.patient.location$sample_barcode %in% substr(ecdna.patient, 1, 15),
    ]$end[1]);
chr6.ecdna.start.gene <- which(
    as.numeric(brca.cnv.chr.new.gis.fpkm.order.match.chr6.all.gene.location$start) > chr6.ecdna.start
    )[1];
chr6.ecdna.end.gene <- which(
    as.numeric(brca.cnv.chr.new.gis.fpkm.order.match.chr6.all.gene.location$end) < chr6.ecdna.end
    )[length(which(as.numeric(brca.cnv.chr.new.gis.fpkm.order.match.chr6.all.gene.location$end) < chr6.ecdna.end))];


### PLOTTING ####################################################################

ecdna.bar.chr6 <- create.barplot(
    formula = TCGA.A2.A3XX.01A ~ sample,
    data = fpkm.chr6.median.z.bar,
    main = expression('RNA abundance of genes on chr6 - outlier patient with ecDNA'),
    col = bar.colours,
    main.cex = 1.4,
    border.col = bar.colours,
    border.lwd = 1,
    xlab.label = NULL,
    xlab.cex = 0,
    ylab.label = expression('Robust z-score'),
    ylab.cex = 1.3,
    yaxis.cex = 1.1,
    xaxis.cex = 0,
    xaxis.fontface = 1,
    yaxis.fontface = 1,
    yaxis.tck = c(0.2, 0),
    xaxis.tck = c(0, 0),
    add.rectangle = TRUE,
    xleft.rectangle = c(chr6.ecdna.start.gene - 1),
    xright.rectangle = c(chr6.ecdna.end.gene + 1),
    ybottom.rectangle = -20,
    ytop.rectangle = 170,
    col.rectangle = 'grey',
    alpha.rectangle = 0.5,
    legend = list(
        inside = list(
            fun = draw.key,
            args = list(
                key = list(
                    points = list(
                        col = 'red',
                        pch = 22,
                        cex = 1.5,
                        fill = c('red')
                        ),
                    text = list(
                        lab = c('Outlier genes on ecDNA')
                        ),
                    padding.text = 1,
                    cex = 1
                    )
                ),
            x = 0.65,
            y = 0.95
            )
        )
    );


save.outlier.figure(
    ecdna.bar.chr6,
    c('Figure2d', 'ecdna', 'chr6', 'barplot'),
    width = 8,
    height = 3.5
    );

save.session.profile(file.path('output', 'Figure2d.txt'));
