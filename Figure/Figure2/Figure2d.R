### ECDNA ANALYSIS #########################################################
# This script analyzes RNA abundance of genes on chromosome 6 in an outlier patient 
# with ecDNA. It calculates robust z-scores for the genes and visualizes the results 
# using a bar plot. 
# Date: 2024-08-13

### Load Required Libraries ######################################################
library(BoutrosLab.utilities);
library(BoutrosLab.plotting.general);
# library(rtracklayer)
# library(liftOver)
# library(GenomicRanges)
# 
# 
# 
# 
# 
# ### Preprocess Amplicon Data #####################################################
# brca.amplicon$sample_barcode <- gsub('-', '.', brca.amplicon$sample_barcode);
# outlier.patient.tag.01.brca.sum <- apply(outlier.patient.tag.01.brca, 2, sum);
# 
# brca.amplicon.match <- brca.amplicon[
#     brca.amplicon$sample_barcode %in% substr(
#         names(outlier.patient.tag.01.brca.sum), 1, 15
#         ),
#     ];
# 
# brca.amplicon.match.ec <- brca.amplicon.match[
#     brca.amplicon.match$amplicon_classification == "Circular",
#     ];
# 
# 
# split.amplicon.intervals <- function(df) {
#     intervals <- strsplit(df$amplicon_intervals, ",")
#     num.intervals <- sapply(intervals, length)
#     
#     split.df <- df[rep(seq_len(nrow(df)), num.intervals), ]
#     split.df$amplicon_intervals <- unlist(intervals)
#     row.names(split.df) <- NULL
#     
#     return(split.df)
#     };
# 
# brca.amplicon.match.ec.each <- do.call(
#     rbind,
#     lapply(
#         seq_len(nrow(brca.amplicon.match.ec)),
#         function(i) split.amplicon.intervals(brca.amplicon.match.ec[i,])
#         )
#     );
# 
# # Extract Chromosome, Start, and End from Intervals
# split.intervals <- strsplit(
#     as.character(brca.amplicon.match.ec.each$amplicon_intervals),
#     "[:-]"
#     );
# 
# brca.amplicon.match.ec.each.chr <- cbind(
#     brca.amplicon.match.ec.each[,1:3],
#     chr = sapply(split.intervals, "[[", 1),
#     start = sapply(split.intervals, "[[", 2),
#     end = sapply(split.intervals, "[[", 3)
#     );
# 
# 
# hg19.location <- brca.amplicon.match.ec.each.chr[,4:6];
# hg19.location[,2] <- as.numeric(hg19.location[,2]);
# hg19.location[,3] <- as.numeric(hg19.location[,3]);
# hg19.location[,1] <- paste0("chr", hg19.location[,1]);
# hg19.location.gr <- makeGRangesFromDataFrame(hg19.location);
# 
# 
# chain.url <- "http://hgdownload.cse.ucsc.edu/goldenpath/hg19/liftOver/hg19ToHg38.over.chain.gz";
# chain.file <- "hg19ToHg38.over.chain.gz";
# download.file(chain.url, chain.file);
# system(paste0("gzip -d ", chain.file));
# ch <- import.chain(sub("\\.gz$", "", chain.file));
# 
# 
# hg38.location <- liftOver(hg19.location.gr, ch);
# 
# 
# process.lifted.location <- function(location) {
#     if (length(location) == 0) return(c("0", "0", "0"))
#     combined.location <- reduce(location)
#     c(
#         as.character(seqnames(combined.location)[1]),
#         as.character(start(combined.location)[1]),
#         as.character(end(combined.location)[length(end(combined.location))])
#         )
#     };
# 
# hg38.location.reduce <- t(sapply(hg38.location, process.lifted.location));
# colnames(hg38.location.reduce) <- c('chr', 'start', 'end');
# 
# 
# brca.amplicon.match.ec.each.chr.hg38 <- cbind(
#     brca.amplicon.match.ec.each.chr[,1:3],
#     hg38.location.reduce
#     );
# # Summing the outlier patient tags for each BRCA sample
# outlier.patient.tag.01.brca.patient.sum <- apply(outlier.patient.tag.01.brca, 2, sum);
# 
# # Extracting ecDNA related outlier patient tags for BRCA samples
# outlier.patient.tag.01.brca.ecDNA <- outlier.patient.tag.01.brca[
#     , substr(colnames(outlier.patient.tag.01.brca), 1, 15) %in% unique(brca.amplicon.match.ec.each.chr.hg38$sample_barcode)
#     ];
# 
# # Filtering BRCA samples with ecDNA outlier tags
# outlier.patient.tag.01.brca.ecDNA.match <- outlier.patient.tag.01.brca[
#     apply(outlier.patient.tag.01.brca.ecDNA, 1, sum) > 0,
#     ];
# 
# # Matching outlier locations in BRCA CNV data for specific gene names
# brca.cnv.chr.new.gis.fpkm.order.match.chr.outlier.location <- brca.cnv.chr.new.gis.fpkm.order.match.chr[
#     brca.cnv.chr.new.gis.fpkm.order.match.chr$gene_name %in% brca.outlier.symbol,
#     ];
# 
# # Matching outlier locations with ecDNA in BRCA CNV data
# brca.cnv.chr.new.gis.fpkm.order.match.chr.outlier.location.ecdna.match <- brca.cnv.chr.new.gis.fpkm.order.match.chr.outlier.location[
#     brca.cnv.chr.new.gis.fpkm.order.match.chr.outlier.location$gene_name %in% fpkm.tumor.symbol.filter.brca[
#         rownames(outlier.patient.tag.01.brca.ecDNA.match), 'Symbol'
#         ],
#     ];
# 
# 
# # Pre-compute patient barcodes
# patient_barcodes <- substr(colnames(outlier.patient.tag.01.brca.ecDNA.match), 1, 15);
# 
# 
# ecDNA.outlier.patient.location <- list();
# for (i in 1:nrow(brca.amplicon.match.ec.each.chr.hg38)) {
#     
#     # Extract target patient details from the current row
#     target.patient <- brca.amplicon.match.ec.each.chr.hg38$sample_barcode[i];
#     target.chr <- brca.amplicon.match.ec.each.chr.hg38$chr[i];
#     target.start <- brca.amplicon.match.ec.each.chr.hg38$start[i];
#     target.end <- brca.amplicon.match.ec.each.chr.hg38$end[i];
#     
#     target.outlier <- rownames(outlier.patient.tag.01.brca.ecDNA.match)[
#         outlier.patient.tag.01.brca.ecDNA.match[, patient_barcodes %in% target.patient] == 1
#         ];
#     
#     # Get the corresponding gene symbols for the outliers
#     target.symbol <- fpkm.tumor.symbol.filter.brca[target.outlier, 'Symbol', drop = FALSE];
#     
#     # Locate the symbols in the ecDNA matched location dataset
#     target.symbol.location <- brca.cnv.chr.new.gis.fpkm.order.match.chr.outlier.location.ecdna.match[
#         brca.cnv.chr.new.gis.fpkm.order.match.chr.outlier.location.ecdna.match$gene_name %in% target.symbol$Symbol, 
#         ];
#     
#     # Check if symbols are within the target range
#     matches <- target.symbol.location$chromosome == target.chr &
#         target.symbol.location$start >= target.start &
#         target.symbol.location$end <= target.end;
#     
#     if (any(matches)) {
#         
#         ecDNA.outlier <- cbind(
#             brca.amplicon.match.ec.each.chr.hg38[i, ],
#             target.symbol.location[matches, ]
#             );
#         
#         
#         ecDNA.outlier.patient.location[[length(ecDNA.outlier.patient.location) + 1]] <- ecDNA.outlier;
#         }
#     }
# 
# # Combine all results into a single dataframe
# ecDNA.outlier.patient.location <- do.call(rbind, ecDNA.outlier.patient.location)
# 
# 
# brca.cnv.chr.new.gis.fpkm.order.match.chr6 <- brca.cnv.chr.new.gis.fpkm.order.match[brca.cnv.chr.new.gis.fpkm.order.match.chr$chromosome == 'chr6',];
# brca.cnv.chr.new.gis.fpkm.order.match.chr6.all.gene.location <- brca.cnv.chr.new.gis.fpkm.order.match.chr[brca.cnv.chr.new.gis.fpkm.order.match.chr$chromosome %in% 'chr6', ];
# 
# # Filter FPKM data for genes on chromosome 6
# fpkm.chr6 <- fpkm.tumor.symbol.filter.brca[
#     match(brca.cnv.chr.new.gis.fpkm.order.match.chr6$Hugo_Symbol, fpkm.tumor.symbol.filter.brca$Symbol), 
#     patient.part.brca
#     ];
# 
# # Calculate robust z-scores for each gene
# fpkm.chr6.median.z <- NULL;
# for (i in 1:nrow(fpkm.chr6)) {
#     gene.median <- median(as.numeric(fpkm.chr6[i,]));
#     gene.mad <- mad(as.numeric(fpkm.chr6[i,]));
#     robust.z <- (as.numeric(fpkm.chr6[i,]) - gene.median) / gene.mad;
#     fpkm.chr6.median.z <- rbind(fpkm.chr6.median.z, robust.z);
#     }
# 
# # Set row and column names for the z-score matrix
# rownames(fpkm.chr6.median.z) <- rownames(fpkm.chr6);
# colnames(fpkm.chr6.median.z) <- colnames(fpkm.chr6);
# 
# 
# ecdna.patient <- 'TCGA.A2.A3XX.01A';
# 
# # Prepare data for plotting
# fpkm.chr6.median.z.bar <- data.frame(
#     cbind(fpkm.chr6.median.z[, ecdna.patient, drop = FALSE], sample = c(1:nrow(fpkm.chr6.median.z)))
#     );
# 
# 
# 
# # Set bar colors, highlighting outlier genes in red
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
    ecDNA.outlier.patient.location$sample_barcode %in% substr(ecdna.patient, 1, 15),]$start[1]
    );
chr6.ecdna.end <- as.numeric(ecDNA.outlier.patient.location[
    ecDNA.outlier.patient.location$sample_barcode %in% substr(ecdna.patient, 1, 15),]$end[1]
    );
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
    col.rectangle = "grey",
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

output.directory <- get0('output.directory', ifnotfound = 'figures');

# Save the plot as a PNG
BoutrosLab.plotting.general::write.plot(
    trellis.object = ecdna.bar.chr6, 
    filename = file.path(output.directory, 'Figure_2_d.png'),
    width = 8,
    height = 3.5,
    size.units = 'in',
    resolution = 1200
);
