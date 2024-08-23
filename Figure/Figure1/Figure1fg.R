### HISTORY #####################################################################
# This script performs meta-analysis for multiple datasets across chromosomes
# Date: 2024-08-13


library(metafor);

### DATA PREPARATION ############################################################

### 1. Chromosomal enrichment ###################################################

### 1. TCGA
# Get chromosomal location information for outlier genes
gene.list <- rownames(outlier.gene.fdr.01.brca);
gene.list.sub <- substr(gene.list, 1, 15);
ensembl <- biomaRt:::useEnsembl(biomart = "ensembl", 
                     dataset = "hsapiens_gene_ensembl", 
                     mirror = "useast");
ensembl <- biomaRt:::useDataset(dataset = "hsapiens_gene_ensembl", mart = ensembl);
gene.position <- biomaRt:::getBM(attributes = c('ensembl_gene_id', 'hgnc_symbol', 'chromosome_name',
                                                'start_position', 'end_position', 'band', "gene_biotype", "entrezgene_id"),
                                 filters = 'ensembl_gene_id', 
                                 values = gene.list.sub, 
                                 mart = ensembl);
gene.position.brca <- gene.position;

# Get chromosomal location information for all genes
fpkm.tumor.symbol.filter.max.brca <- apply(fpkm.tumor.symbol.filter.brca[,patient.part.brca], 1, max);
gene.list <- rownames(fpkm.tumor.symbol.filter.brca)[fpkm.tumor.symbol.filter.max.brca > 5];
gene.list.sub <- substr(gene.list, 1, 15);
ensembl <- biomaRt:::useEnsembl(biomart = "ensembl", 
                                dataset = "hsapiens_gene_ensembl", 
                                mirror = "useast");
ensembl <- biomaRt:::useDataset(dataset = "hsapiens_gene_ensembl", mart = ensembl);
gene.position <- biomaRt:::getBM(attributes = c('ensembl_gene_id', 'hgnc_symbol', 'chromosome_name',
                                                'start_position', 'end_position', 'band', "gene_biotype", "entrezgene_id"),
                                 filters = 'ensembl_gene_id', 
                                 values = gene.list.sub, 
                                 mart = ensembl);

gene.position.brca.all <- gene.position;


chr.position.brca <- data.frame(as.matrix(table(gene.position.brca$chromosome_name)))
chr.name <- c('1', '2', '3', '4', '5', '6', '7', '8', '9', '10', 
              '11', '12', '13', '14', '15', '16', '17', '18', '19', '20', 
              '21', '22', 'MT', 'X', 'Y');
chr.position.order.brca <- chr.position.brca[chr.name,,drop = FALSE ];
rownames(chr.position.order.brca) <- chr.name;
chr.position.order.brca[is.na(chr.position.order.brca$as.matrix.table.gene.position.brca.chromosome_name..),] <- 0;
chr.position.outlier.brca <- data.frame(cbind(chr = c(1:25), count = as.numeric(chr.position.order.brca[,1])));


# all transcripts from BRCA data
chr.position.brca.all <- data.frame(as.matrix(table(gene.position.brca.all$chromosome_name)))
chr.position.order.brca.all <- chr.position.brca.all[chr.name,,drop = FALSE ];
rownames(chr.position.order.brca.all) <- chr.name;
chr.position.order.brca.all[is.na(chr.position.order.brca.all$as.matrix.table.gene.position.brca.all.chromosome_name..),] <- 0;
chr.position.outlier.brca.all <- data.frame(cbind(chr = c(1:25), count = as.numeric(chr.position.order.brca.all[,1])));

# segment plot
p.value.chr.brca.fisher.sub <- NULL;
p.value.chr.brca.odd.sub <- NULL;
for (i in 1:25) {
   total.gene <- nrow(fpkm.tumor.symbol.filter.brca);
    chr.gene <- chr.position.outlier.brca.all$count[i] # number of genes on the chromosome of interest in the population
    total.outlier <- nrow(outlier.gene.fdr.01.brca) # sample size
    chr.outlier <- chr.position.outlier.brca$count[i] # number of genes on the chromosome of interest in the sample
     if (is.na(chr.outlier)) {
        chr.outlier <- 0;
        }

     p_value <- fisher.test(matrix(c(chr.outlier, total.outlier - chr.outlier, chr.gene - chr.outlier, total.gene - total.outlier - chr.gene + chr.outlier), nrow=2), alternative="two.sided")$p.value;
    p.value.chr.brca.fisher.sub <- c(p.value.chr.brca.fisher.sub, p_value);
     
     odd.ratio <- fisher.test(matrix(c(chr.outlier, total.outlier - chr.outlier, chr.gene - chr.outlier, total.gene - total.outlier - chr.gene + chr.outlier), nrow=2), alternative="two.sided")
     odd.ratio.ci <- c(odd.ratio$estimate ,odd.ratio$conf.int);
    p.value.chr.brca.odd.sub <- rbind(p.value.chr.brca.odd.sub, odd.ratio.ci);
    }



### 2. METABIRC
# Get chromosomal location information for outlier genes
gene.list <- substr(rownames(outlier.gene.fdr.01.meta), 1, nchar(rownames(outlier.gene.fdr.01.meta))-3)
ensembl <- biomaRt:::useEnsembl(biomart = "ensembl", 
                     dataset = "hsapiens_gene_ensembl", 
                     mirror = "uswest");
ensembl <- biomaRt:::useDataset(dataset = "hsapiens_gene_ensembl", mart = ensembl);
gene.position <- biomaRt:::getBM(attributes = c('ensembl_gene_id', 'hgnc_symbol', 'chromosome_name',
                                                'start_position', 'end_position', 'band', "gene_biotype", "entrezgene_id"),
                                 filters = "entrezgene_id", 
                                 values = gene.list, 
                                 mart = ensembl);
gene.position.meta <- gene.position;

# Get chromosomal location information for all genes
gene.list <- substr(rownames(fpkm.tumor.symbol.filter.meta.symbol), 1, nchar(rownames(fpkm.tumor.symbol.filter.meta.symbol))-3)
ensembl <- biomaRt:::useEnsembl(biomart = "ensembl", 
                     dataset = "hsapiens_gene_ensembl", 
                     mirror = "uswest");
ensembl <- biomaRt:::useDataset(dataset = "hsapiens_gene_ensembl", mart = ensembl);
gene.position <- biomaRt:::getBM(attributes = c('ensembl_gene_id', 'hgnc_symbol', 'chromosome_name',
                                                'start_position', 'end_position', 'band', "gene_biotype", "entrezgene_id"),
                                 filters = "entrezgene_id", 
                                 values = gene.list, 
                                 mart = ensembl);
gene.position.meta.all <- gene.position;


chr.position.meta.all <- data.frame(as.matrix(table(gene.position.meta.all$chromosome_name)))
chr.position.order.meta.all <- chr.position.meta.all[chr.name,,drop = FALSE ];
rownames(chr.position.order.meta.all) <- chr.name;
chr.position.order.meta.all[is.na(chr.position.order.meta.all$as.matrix.table.gene.position.meta.all.chromosome_name..),] <- 0;
chr.position.outlier.meta.all <- data.frame(cbind(chr = c(1:25), count = as.numeric(chr.position.order.meta.all[,1])));

# segment plot
p.value.chr.meta.fisher.sub <- NULL;
p.value.chr.meta.odd.sub <- NULL;
for (i in 1:25) {
   total.gene <- nrow(fpkm.tumor.symbol.filter.meta.symbol);
    chr.gene <- chr.position.outlier.meta.all$count[i] # number of genes on the chromosome of interest in the population
    total.outlier <- nrow(outlier.gene.fdr.01.meta) # sample size
    chr.outlier <- chr.position.outlier.meta$count[i] # number of genes on the chromosome of interest in the sample
     if (is.na(chr.outlier)) {
        chr.outlier <- 0;
        }

     p_value <- fisher.test(matrix(c(chr.outlier, total.outlier - chr.outlier, chr.gene - chr.outlier, total.gene - total.outlier - chr.gene + chr.outlier), nrow=2), alternative="two.sided")$p.value;
    p.value.chr.meta.fisher.sub <- c(p.value.chr.meta.fisher.sub, p_value);
     
     odd.ratio <- fisher.test(matrix(c(chr.outlier, total.outlier - chr.outlier, chr.gene - chr.outlier, total.gene - total.outlier - chr.gene + chr.outlier), nrow=2), alternative="two.sided")
     odd.ratio.ci <- c(odd.ratio$estimate ,odd.ratio$conf.int);
    p.value.chr.meta.odd.sub <- rbind(p.value.chr.meta.odd.sub, odd.ratio.ci);
    }



### 3. ISPY
# Get chromosomal location information for outlier genes
gene.list <- rownames(outlier.gene.fdr.01.ispy);
gene.position <- biomaRt:::getBM(attributes = c('ensembl_gene_id', 'hgnc_symbol', 'chromosome_name',
                                                'start_position', 'end_position', 'band', "gene_biotype", "entrezgene_id"),
                                 filters = "hgnc_symbol", 
                                 values = gene.list, 
                                 mart = ensembl);
gene.position.ispy <- gene.position;

# Get chromosomal location information for all genes
gene.list <- rownames(fpkm.tumor.symbol.filter.ispy);
gene.position <- biomaRt:::getBM(attributes = c('ensembl_gene_id', 'hgnc_symbol', 'chromosome_name',
                                                'start_position', 'end_position', 'band', "gene_biotype", "entrezgene_id"),
                                 filters = "hgnc_symbol", 
                                 values = gene.list, 
                                 mart = ensembl);
gene.position.ispy.all <- gene.position;

# segment plot
p.value.chr.ispy.fisher.sub <- NULL;
p.value.chr.ispy.odd.sub <- NULL;
for (i in 1:25) {
   total.gene <- nrow(fpkm.tumor.symbol.filter.ispy);
    chr.gene <- chr.position.outlier.ispy.all$count[i] # number of genes on the chromosome of interest in the population
    total.outlier <- nrow(outlier.gene.fdr.01.ispy) # sample size
    chr.outlier <- chr.position.outlier.ispy$count[i] # number of genes on the chromosome of interest in the sample
     if (is.na(chr.outlier)) {
        chr.outlier <- 0;
        }

     p_value <- fisher.test(matrix(c(chr.outlier, total.outlier - chr.outlier, chr.gene - chr.outlier, total.gene - total.outlier - chr.gene + chr.outlier), nrow=2), alternative="two.sided")$p.value;
    p.value.chr.ispy.fisher.sub <- c(p.value.chr.ispy.fisher.sub, p_value);
     
     odd.ratio <- fisher.test(matrix(c(chr.outlier, total.outlier - chr.outlier, chr.gene - chr.outlier, total.gene - total.outlier - chr.gene + chr.outlier), nrow=2), alternative="two.sided")
     odd.ratio.ci <- c(odd.ratio$estimate ,odd.ratio$conf.int);
    p.value.chr.ispy.odd.sub <- rbind(p.value.chr.ispy.odd.sub, odd.ratio.ci);
    }



### 4. MATADOR
# Get chromosomal location information for outlier genes
gene.list <- rownames(outlier.gene.fdr.01.matador);
gene.list.sub <- substr(gene.list, 1, 15);
ensembl <- biomaRt:::useEnsembl(biomart = "ensembl", 
                     dataset = "hsapiens_gene_ensembl", 
                     mirror = "useast");
ensembl <- biomaRt:::useDataset(dataset = "hsapiens_gene_ensembl", mart = ensembl);
gene.position <- biomaRt:::getBM(attributes = c('ensembl_gene_id', 'hgnc_symbol', 'chromosome_name',
                                                'start_position', 'end_position', 'band', "gene_biotype", "entrezgene_id"),
                                 filters = 'ensembl_gene_id', 
                                 values = gene.list.sub, 
                                 mart = ensembl);

# Get chromosomal location information for all genes
fpkm.tumor.symbol.filter.metador.max <- apply(fpkm.tumor.symbol.filter.metador, 1, max);
fpkm.tumor.symbol.filter.metador.max.filter <- fpkm.tumor.symbol.filter.metador[fpkm.tumor.symbol.filter.metador.max > 5,];
gene.list <- rownames(fpkm.tumor.symbol.filter.metador.max.filter);
gene.list.sub <- substr(gene.list, 1, 15);
gene.position <- biomaRt:::getBM(attributes = c('ensembl_gene_id', 'hgnc_symbol', 'chromosome_name',
                                                'start_position', 'end_position', 'band', "gene_biotype", "entrezgene_id"),
                                 filters = 'ensembl_gene_id', 
                                 values = gene.list.sub, 
                                 mart = ensembl);
gene.position.metador.all <- gene.position;

# segment plot
p.value.chr.metador.fisher.sub <- NULL;
p.value.chr.metador.odd.sub <- NULL;
for (i in 1:25) {
   total.gene <- nrow(fpkm.tumor.symbol.filter.metador);
    chr.gene <- chr.position.outlier.metador.all$count[i] # number of genes on the chromosome of interest in the population
    total.outlier <- nrow(outlier.gene.fdr.01.matador) # sample size
    chr.outlier <- chr.position.outlier.metador$count[i] # number of genes on the chromosome of interest in the sample
     if (is.na(chr.outlier)) {
        chr.outlier <- 0;
        }

     p_value <- fisher.test(matrix(c(chr.outlier, total.outlier - chr.outlier, chr.gene - chr.outlier, total.gene - total.outlier - chr.gene + chr.outlier), nrow=2), alternative="two.sided")$p.value;
    p.value.chr.metador.fisher.sub <- c(p.value.chr.metador.fisher.sub, p_value);
     
     odd.ratio <- fisher.test(matrix(c(chr.outlier, total.outlier - chr.outlier, chr.gene - chr.outlier, total.gene - total.outlier - chr.gene + chr.outlier), nrow=2), alternative="two.sided")
     odd.ratio.ci <- c(odd.ratio$estimate ,odd.ratio$conf.int);
    p.value.chr.metador.odd.sub <- rbind(p.value.chr.metador.odd.sub, odd.ratio.ci);
    }




### 5. ICGC BRCA-EU
gene.position.icgc.all.chr <-  gsub(":.*", "", fpkm.data.icgc$loc[as.numeric(rownames(fpkm.tumor.symbol.filter.icgc))]);
gene.position.icgc.all.chr.table <- data.frame(as.matrix(table(gene.position.icgc.all.chr)));
chr.position.outlier.icgc.all <- gene.position.icgc.all.chr.table[chr.name,,drop = FALSE ];
chr.position.outlier.icgc.all <- data.frame(cbind(chr = c(1:25), count = as.numeric(chr.position.outlier.icgc.all[,1])));

gene.position.icgc.chr <- gsub(":.*", "", fpkm.data.icgc$loc[as.numeric(outlier.gene.fdr.01.icgc$gene)]);
chr.position.icgc <- data.frame(as.matrix(table(gene.position.icgc.chr)))
chr.name <- c('1', '2', '3', '4', '5', '6', '7', '8', '9', '10', 
                '11', '12', '13', '14', '15', '16', '17', '18', '19', '20', 
                '21', '22', 'MT', 'X', 'Y');
chr.position.order.icgc <- chr.position.icgc[chr.name,,drop = FALSE ];
rownames(chr.position.order.icgc) <- chr.name;
chr.position.order.icgc[is.na(chr.position.order.icgc$as.matrix.table.fpkm.data.icgc.chr.outlier..),] <- 0;
chr.position.outlier.icgc <- data.frame(cbind(chr = c(1:25), count = as.numeric(chr.position.order.icgc[,1])));

# segment plot
p.value.chr.icgc.fisher.sub <- NULL;
p.value.chr.icgc.odd.sub <- NULL;
for (i in 1:25) {
   total.gene <- nrow(fpkm.tumor.symbol.filter.icgc);
    chr.gene <- chr.position.outlier.icgc.all$count[i] # number of genes on the chromosome of interest in the population
    total.outlier <- nrow(outlier.gene.fdr.01.icgc) # sample size
    chr.outlier <- chr.position.outlier.icgc$count[i] # number of genes on the chromosome of interest in the sample
     if (is.na(chr.outlier)) {
        chr.outlier <- 0;
        }

     p_value <- fisher.test(matrix(c(chr.outlier, total.outlier - chr.outlier, chr.gene - chr.outlier, total.gene - total.outlier - chr.gene + chr.outlier), nrow=2), alternative="two.sided")$p.value;
    p.value.chr.icgc.fisher.sub <- c(p.value.chr.icgc.fisher.sub, p_value);
     
     odd.ratio <- fisher.test(matrix(c(chr.outlier, total.outlier - chr.outlier, chr.gene - chr.outlier, total.gene - total.outlier - chr.gene + chr.outlier), nrow=2), alternative="two.sided")
     odd.ratio.ci <- c(odd.ratio$estimate ,odd.ratio$conf.int);
    p.value.chr.icgc.odd.sub <- rbind(p.value.chr.icgc.odd.sub, odd.ratio.ci);
    }




### Meta-analysis
p.value.chr.brca.odd.sub.df <- data.frame(p.value.chr.brca.odd.sub);
p.value.chr.meta.odd.sub.df <- data.frame(p.value.chr.meta.odd.sub);
p.value.chr.ispy.odd.sub.df <- data.frame(p.value.chr.ispy.odd.sub);
p.value.chr.metador.odd.sub.df <- data.frame(p.value.chr.metador.odd.sub);
p.value.chr.icgc.odd.sub.df <- data.frame(p.value.chr.icgc.odd.sub);


#   - use natural log
ln.odd.brca <- log(p.value.chr.brca.odd.sub.df$odds.ratio);
se.odd.brca <- (log(p.value.chr.brca.odd.sub.df$V3)-log(p.value.chr.brca.odd.sub.df$V2)) / 3.92;
ln.odd.meta <- log(p.value.chr.meta.odd.sub.df$odds.ratio);
se.odd.meta <- (log(p.value.chr.meta.odd.sub.df$V3)-log(p.value.chr.meta.odd.sub.df$V2)) / 3.92;
ln.odd.ispy <- log(p.value.chr.ispy.odd.sub.df$odds.ratio);
se.odd.ispy <- (log(p.value.chr.ispy.odd.sub.df$V3)-log(p.value.chr.ispy.odd.sub.df$V2)) / 3.92;
ln.odd.metador <- log(p.value.chr.metador.odd.sub.df$odds.ratio);
se.odd.metador <- (log(p.value.chr.metador.odd.sub.df$V3)-log(p.value.chr.metador.odd.sub.df$V2)) / 3.92;
ln.odd.icgc <- log(p.value.chr.icgc.odd.sub.df$odds.ratio);
se.odd.icgc <- (log(p.value.chr.icgc.odd.sub.df$V3)-log(p.value.chr.icgc.odd.sub.df$V2)) / 3.92;



chr.odd.se.5 <- list();

# Loop through each row of the chromosome odds ratio dataframe
for (i in 1:nrow(p.value.chr.brca.odd.sub.df)) {
    
    # Extract the odds ratios and standard errors for each dataset
    chr.odd <- c(
        ln.odd.brca[i],
        ln.odd.meta[i],
        ln.odd.ispy[i],
        ln.odd.metador[i],
        ln.odd.icgc[i]
        );
    chr.se <- c(
        se.odd.brca[i],
        se.odd.meta[i],
        se.odd.ispy[i],
        se.odd.metador[i],
        se.odd.icgc[i]
        );
    
    # Combine odds ratios and standard errors into a single data frame
    chr.all <- data.frame(cbind(chr.odd, chr.se));
    
    # Store the data frame in the list
    chr.odd.se.5[[i]] <- chr.all;
    
    }

metafor.chr.odd.ci.p.5 <- NULL;

# Loop through each chromosome's data for meta-analysis
for (i in 1:nrow(p.value.chr.brca.odd.sub.df)) {
    
    # Extract the sample data for the current chromosome
    chr.odd.se.sample <- chr.odd.se.5[[i]];
    
    # Remove infinite values from the sample data
    chr.odd.se.sample.inf <- chr.odd.se.sample[
        !is.infinite(chr.odd.se.sample$chr.odd) & !is.infinite(chr.odd.se.sample$chr.se)
        ];
    
    # Perform meta-analysis if there are more than one valid sample
    if (nrow(chr.odd.se.sample.inf) > 1) {
        
        metafor.chr <- rma.uni(
            yi = chr.odd, 
            sei = chr.se, 
            data = chr.odd.se.sample.inf, 
            method = 'DL'
            );
        
        metafor.chr.odd <- exp(metafor.chr$beta);
        metafor.chr.lower <- exp(metafor.chr$ci.lb);
        metafor.chr.upper <- exp(metafor.chr$ci.ub);
        metafor.chr.p <- metafor.chr$pval;
        
        metafor.all <- c(
            metafor.chr.odd, 
            metafor.chr.lower, 
            metafor.chr.upper, 
            metafor.chr.p
            );
        
        } 
    else {
        
        metafor.all <- c(NA, NA, NA, NA);
        
        }
    
    # Append the results to the matrix
    metafor.chr.odd.ci.p.5 <- rbind(metafor.chr.odd.ci.p.5, metafor.all);
    
    }

# Create a data frame to store meta-analysis results with labels
metafor.chr.odd.ci.p.data.5 <- data.frame(
    p.value = metafor.chr.odd.ci.p.5[1:24, 4],
    odd = metafor.chr.odd.ci.p.5[1:24, 1],
    ci.min = metafor.chr.odd.ci.p.5[1:24, 2], 
    ci.max = metafor.chr.odd.ci.p.5[1:24, 3]
    );

metafor.chr.odd.ci.p.data.label.5 <- data.frame(
    metafor.chr.odd.ci.p.data.5,
    labels = as.factor(paste('Chr', chr.name[1:24], sep = ''))
    );





### 2. Exon number/Gene length/GC content/RNA abundance #########################

### - 1. GC content
# 1. TCGA-BRCA
gene.list <- rownames(outlier.gene.fdr.01.brca);
gene.list.sub <- substr(gene.list, 1, 15);
gene.gc.brca <- gene.gc[match(gene.list.sub, gene.gc$Gene.stable.ID), ];
gene.gc.brca.content <- as.numeric(gene.gc.brca$Gene...GC.content);
outlier.gene.fdr.01.brca.gc <- cbind(
    outlier.gene.fdr.01.brca,                                             
    GC.content = gene.gc.brca.content
    );

fpkm.tumor.symbol.filter.max.brca <- apply(fpkm.tumor.symbol.filter.brca[, patient.part.brca], 1, max);
fpkm.tumor.symbol.filter.brca.max.5 <- fpkm.tumor.symbol.filter.brca[fpkm.tumor.symbol.filter.max.brca > 5, ];
gene.list <- rownames(fpkm.tumor.symbol.filter.brca.max.5);
gene.list.sub <- substr(gene.list, 1, 15);
gene.gc.brca.all <- gene.gc[match(gene.list.sub, gene.gc$Gene.stable.ID), ];
gene.gc.brca.content.all <- as.numeric(gene.gc.brca.all$Gene...GC.content);
outlier.gene.fdr.all.brca.gc <- cbind(
    fpkm.tumor.symbol.filter.brca.max.5,                                
    GC.content = gene.gc.brca.content.all
    );

gene.position.all.gc.non.brca <- outlier.gene.fdr.all.brca.gc[!(rownames(outlier.gene.fdr.all.brca.gc) %in% rownames(outlier.gene.fdr.01.brca)), ];
gene.position.all.gc.outlier.brca <- outlier.gene.fdr.01.brca.gc;
gc.num.non.brca <- gene.position.all.gc.non.brca$GC.content;
gc.num.outlier.brca <- gene.position.all.gc.outlier.brca$GC.content;

gc.box.brca <- data.frame(cbind(
    gc.content = c(gc.num.non.brca, gc.num.outlier.brca),
    status = c(rep('non', length(gc.num.non.brca)),
               rep('out', length(gc.num.outlier.brca)))
    ));
gc.box.brca <- na.omit(gc.box.brca);    
gc.box.brca$gc.content <- as.numeric(gc.box.brca$gc.content);

# 2. METABRIC
gene.list <- rownames(outlier.gene.fdr.01.meta);
gene.list.sub <- substr(gene.list, 1, nchar(gene.list) - 3);
gene.position.meta.order <- gene.position.meta[match(gene.list.sub, gene.position.meta$entrezgene_id), ];
outlier.gene.fdr.01.meta.ensg <- cbind(
    outlier.gene.fdr.01.meta,
    ensembl = gene.position.meta.order$ensembl_gene_id
    );

gene.gc.meta <- gene.gc[match(gene.position.meta.order[gene.position.meta.order$entrezgene_id %in% gene.list.sub, ]$ensembl_gene_id, gene.gc$Gene.stable.ID), ];
gene.gc.meta.content <- as.numeric(gene.gc.meta$Gene...GC.content);
outlier.gene.fdr.01.meta.gc <- cbind(
    na.omit(outlier.gene.fdr.01.meta.ensg), 
    GC.content = gene.gc.meta.content
    );

gene.list <- rownames(fpkm.tumor.symbol.filter.meta.symbol);
gene.list.sub <- substr(gene.list, 1, nchar(gene.list) - 3);
gene.position.meta.order.all <- gene.position.meta.all[match(gene.list.sub, gene.position.meta.all$entrezgene_id), ];
outlier.gene.fdr.all.meta.ensg <- cbind(
    fpkm.tumor.symbol.filter.meta.symbol,
    ensembl = gene.position.meta.order.all$ensembl_gene_id
    );
gene.gc.meta.all <- gene.gc[match(gene.position.meta.order.all$ensembl_gene_id, gene.gc$Gene.stable.ID), ];
gene.gc.meta.content.all <- as.numeric(gene.gc.meta.all$Gene...GC.content);
outlier.gene.fdr.all.meta.gc <- cbind(
    na.omit(outlier.gene.fdr.all.meta.ensg),
    GC.content = gene.gc.meta.content.all
    );

gene.position.all.gc.non.meta <- outlier.gene.fdr.all.meta.gc[!(rownames(outlier.gene.fdr.all.meta.gc) %in% rownames(outlier.gene.fdr.01.meta.gc)), ];
gene.position.all.gc.outlier.meta <- outlier.gene.fdr.01.meta.gc;
gc.num.non.meta <- gene.position.all.gc.non.meta$GC.content;
gc.num.outlier.meta <- gene.position.all.gc.outlier.meta$GC.content;

gc.box.meta <- data.frame(cbind(
    gc.content = c(gc.num.non.meta, gc.num.outlier.meta),
    status = c(rep('non', length(gc.num.non.meta)),
               rep('out', length(gc.num.outlier.meta)))
    ));
gc.box.meta <- na.omit(gc.box.meta);    
gc.box.meta$gc.content <- as.numeric(gc.box.meta$gc.content);

# 3. I-SPY-2
gene.list <- rownames(outlier.gene.fdr.01.ispy);
gene.position.ispy.order <- gene.position.ispy[match(gene.list, gene.position.ispy$hgnc_symbol), ];
outlier.gene.fdr.01.ispy.ensg <- cbind(
    outlier.gene.fdr.01.ispy,
    ensembl = gene.position.ispy.order$ensembl_gene_id
    );
gene.gc.ispy <- gene.gc[match(gene.position.ispy.order$ensembl_gene_id, gene.gc$Gene.stable.ID), ];
gene.gc.ispy.content <- as.numeric(gene.gc.ispy$Gene...GC.content);
outlier.gene.fdr.01.ispy.gc <- cbind(
    outlier.gene.fdr.01.ispy.ensg,
    GC.content = gene.gc.ispy.content
    );

gene.list <- rownames(fpkm.tumor.symbol.filter.ispy);
gene.position.ispy.order.all <- gene.position.ispy.all[match(gene.list, gene.position.ispy.all$hgnc_symbol), ];
outlier.gene.fdr.all.ispy.ensg <- cbind(
    fpkm.tumor.symbol.filter.ispy,
    ensembl = gene.position.ispy.order.all$ensembl_gene_id
    );
gene.gc.ispy.all <- gene.gc[match(gene.position.ispy.order.all$ensembl_gene_id, gene.gc$Gene.stable.ID), ];
gene.gc.ispy.content.all <- as.numeric(gene.gc.ispy.all$Gene...GC.content);
outlier.gene.fdr.all.ispy.gc <- cbind(
    outlier.gene.fdr.all.ispy.ensg,
    GC.content = gene.gc.ispy.content.all
    );

gene.position.all.gc.non.ispy <- outlier.gene.fdr.all.ispy.gc[!(rownames(outlier.gene.fdr.all.ispy.gc) %in% rownames(outlier.gene.fdr.01.ispy.gc)), ];
gene.position.all.gc.outlier.ispy <- outlier.gene.fdr.01.ispy.gc;
gc.num.non.ispy <- gene.position.all.gc.non.ispy$GC.content;
gc.num.outlier.ispy <- gene.position.all.gc.outlier.ispy$GC.content;

gc.box.ispy <- data.frame(cbind(
    gc.content = c(gc.num.non.ispy, gc.num.outlier.ispy),
    status = c(rep('non', length(gc.num.non.ispy)),
               rep('out', length(gc.num.outlier.ispy)))
    ));

# 4. MATADOR
gene.list <- rownames(outlier.gene.fdr.01.matador);
gene.list.metadordor.sub <- substr(gene.list, 1, 15);
gene.gc.metador <- gene.gc[match(gene.list.metadordor.sub, gene.gc$Gene.stable.ID), ];
gene.gc.metador.content <- as.numeric(gene.gc.metador$Gene...GC.content);
outlier.gene.fdr.01.matador.gc <- cbind(
    outlier.gene.fdr.01.matador,
    GC.content = gene.gc.metador.content
    );

gene.list <- rownames(fpkm.tumor.symbol.filter.metador.symbol);
gene.list.metadordor.sub.all <- substr(gene.list, 1, 15);
gene.gc.metador.all <- gene.gc[match(gene.list.metadordor.sub.all, gene.gc$Gene.stable.ID), ];
gene.gc.metador.content.all <- as.numeric(gene.gc.metador.all$Gene...GC.content);
outlier.gene.fdr.all.matador.gc <- cbind(
    fpkm.tumor.symbol.filter.metador.symbol,
    GC.content = gene.gc.metador.content.all
    );

gene.position.all.gc.non.metador <- outlier.gene.fdr.all.matador.gc[!(rownames(outlier.gene.fdr.all.matador.gc) %in% rownames(outlier.gene.fdr.01.matador.gc)), ];
gene.position.all.gc.outlier.metador <- outlier.gene.fdr.01.matador.gc;
gc.num.non.metador <- gene.position.all.gc.non.metador$GC.content;
gc.num.outlier.metador <- gene.position.all.gc.outlier.metador$GC.content;

gc.box.metador <- data.frame(cbind(
    gc.content = c(gc.num.non.metador, gc.num.outlier.metador),
    status = c(rep('non', length(gc.num.non.metador)),
               rep('out', length(gc.num.outlier.metador)))
    ));
gc.box.metador <- na.omit(gc.box.metador);    
gc.box.metador$gc.content <- as.numeric(gc.box.metador$gc.content);

# 5. ICGC
gene.list.icgc <- fpkm.data.icgc$Ensembl[as.numeric(outlier.gene.fdr.01.icgc$gene)];
gene.gc.icgc <- gene.gc[match(gene.list.icgc, gene.gc$Gene.stable.ID), ];
gene.gc.icgc.content <- as.numeric(gene.gc.icgc$Gene...GC.content);
outlier.gene.fdr.01.icgc.gc <- cbind(
    outlier.gene.fdr.01.icgc,
    GC.content = gene.gc.icgc.content
    );

gene.list.icgc.sub.all <- fpkm.data.icgc$Ensembl[as.numeric(rownames(fpkm.tumor.symbol.filter.icgc))];
gene.gc.icgc.all <- gene.gc[match(gene.list.icgc.sub.all, gene.gc$Gene.stable.ID), ];
gene.gc.icgc.content.all <- as.numeric(gene.gc.icgc.all$Gene...GC.content);
outlier.gene.fdr.all.icgc.gc <- cbind(
    fpkm.tumor.symbol.filter.icgc,
    GC.content = gene.gc.icgc.content.all
    );

gene.position.all.gc.non.icgc <- outlier.gene.fdr.all.icgc.gc[!(rownames(outlier.gene.fdr.all.icgc.gc) %in% rownames(outlier.gene.fdr.01.icgc.gc)), ];
gene.position.all.gc.outlier.icgc <- outlier.gene.fdr.01.icgc.gc;
gc.num.non.icgc <- gene.position.all.gc.non.icgc$GC.content;
gc.num.outlier.icgc <- gene.position.all.gc.outlier.icgc$GC.content;

gc.box.icgc <- data.frame(cbind(
    gc.content = c(gc.num.non.icgc, gc.num.outlier.icgc),
    status = c(rep('non', length(gc.num.non.icgc)),
               rep('out', length(gc.num.outlier.icgc)))
    ));
gc.box.icgc <- na.omit(gc.box.icgc);    
gc.box.icgc$gc.content <- as.numeric(gc.box.icgc$gc.content);



# Meta-analysis

brca.out.mean <- mean(gc.box.brca$gc.content[gc.box.brca$status == 'out']);
brca.out.sd <- sd(gc.box.brca$gc.content[gc.box.brca$status == 'out']);
brca.out.n <- length(gc.box.brca$gc.content[gc.box.brca$status == 'out']);
brca.non.mean <- mean(gc.box.brca$gc.content[gc.box.brca$status == 'non']);
brca.non.sd <- sd(gc.box.brca$gc.content[gc.box.brca$status == 'non']);
brca.non.n <- length(gc.box.brca$gc.content[gc.box.brca$status == 'non']);

meta.out.mean <- mean(gc.box.meta$gc.content[gc.box.meta$status == 'out']);
meta.out.sd <- sd(gc.box.meta$gc.content[gc.box.meta$status == 'out']);
meta.non.mean <- mean(gc.box.meta$gc.content[gc.box.meta$status == 'non']);
meta.non.sd <- sd(gc.box.meta$gc.content[gc.box.meta$status == 'non']);
meta.out.n <- length(gc.box.meta$gc.content[gc.box.meta$status == 'out']);
meta.non.n <- length(gc.box.meta$gc.content[gc.box.meta$status == 'non']);

ispy.out.mean <- mean(gc.box.ispy$gc.content[gc.box.ispy$status == 'out']);
ispy.out.sd <- sd(gc.box.ispy$gc.content[gc.box.ispy$status == 'out']);
ispy.non.mean <- mean(gc.box.ispy$gc.content[gc.box.ispy$status == 'non']);
ispy.non.sd <- sd(gc.box.ispy$gc.content[gc.box.ispy$status == 'non']);
ispy.out.n <- length(gc.box.ispy$gc.content[gc.box.ispy$status == 'out']);
ispy.non.n <- length(gc.box.ispy$gc.content[gc.box.ispy$status == 'non']);

metador.out.mean <- mean(gc.box.metador$gc.content[gc.box.metador$status == 'out']);
metador.out.sd <- sd(gc.box.metador$gc.content[gc.box.metador$status == 'out']);
metador.non.mean <- mean(gc.box.metador$gc.content[gc.box.metador$status == 'non']);
metador.non.sd <- sd(gc.box.metador$gc.content[gc.box.metador$status == 'non']);
metador.out.n <- length(gc.box.metador$gc.content[gc.box.metador$status == 'out']);
metador.non.n <- length(gc.box.metador$gc.content[gc.box.metador$status == 'non']);

icgc.out.mean <- mean(gc.box.icgc$gc.content[gc.box.icgc$status == 'out']);
icgc.out.sd <- sd(gc.box.icgc$gc.content[gc.box.icgc$status == 'out']);
icgc.non.mean <- mean(gc.box.icgc$gc.content[gc.box.icgc$status == 'non']);
icgc.non.sd <- sd(gc.box.icgc$gc.content[gc.box.icgc$status == 'non']);
icgc.out.n <- length(gc.box.icgc$gc.content[gc.box.icgc$status == 'out']);
icgc.non.n <- length(gc.box.icgc$gc.content[gc.box.icgc$status == 'non']);


smd.gc.matrix.5 <- data.frame(cbind(
    study = c(1, 2, 3, 4, 5), # 1: TCGA-BRCA, 2: METABRIC, 3: I-SPY2, 4: MATADOR, 5: ICGC
    
    # Mean of the gc number of outlier genes
    mean.out = c(brca.out.mean, meta.out.mean, ispy.out.mean, metador.out.mean, icgc.out.mean), 
    # SD of the gc number of outlier genes
    sd.out = c(brca.out.sd, meta.out.sd, ispy.out.sd, metador.out.sd, icgc.out.sd),
    # Number of outlier genes
    n.out = c(brca.out.n, meta.out.n, ispy.out.n, metador.out.n, icgc.out.n),
    
    # Mean of the gc number of non-outlier genes
    mean.non = c(brca.non.mean, meta.non.mean, ispy.non.mean, metador.non.mean, icgc.non.mean),
    # SD of the gc number of non-outlier genes
    sd.non = c(brca.non.sd, meta.non.sd, ispy.non.sd, metador.non.sd, icgc.non.sd),
    # Number of non-outlier genes
    n.non = c(brca.non.n, meta.non.n, ispy.non.n, metador.non.n, icgc.non.n)
    ));



escalc.gc.matrix.5 <- escalc(measure = "SMD", 
       m1i = mean.out, 
       m2i = mean.non, 
       sd1i = sd.out,  
       sd2i = sd.non, 
       n1i = n.out, 
       n2i = n.non, 
       data = smd.gc.matrix.5, 
       append = TRUE);


smd.metafor.gc.5 <- rma.uni(yi = yi, vi = vi, data=escalc.gc.matrix.5, 
                          # test="knha", 
                          method = 'DL');






### - 2. Gene length

### 1. TCGA-BRCA

# Calculate gene lengths for the sample and population gene sets
outlier.length.brca <- gene.position.brca$end_position - gene.position.brca$start_position + 1;
gene.length.brca <- gene.position.brca.all$end_position - gene.position.brca.all$start_position + 1;

# Create a subset of non-outlier genes
gene.position.all.non.brca <- gene.position.brca.all[!(gene.position.brca.all$ensembl_gene_id %in% gene.position.brca.1$ensembl_gene_id), ];
gene.non.length.brca <- gene.position.all.non.brca$end_position - gene.position.all.non.brca$start_position + 1;

# Create a data frame for the lengths and status of genes
length.box.brca <- data.frame(cbind(
    length.content = c(gene.non.length.brca, outlier.length.brca),
    status = c(rep('non', length(gene.non.length.brca)),
               rep('out', length(outlier.length.brca)))
    ));
length.box.brca <- na.omit(length.box.brca);    
length.box.brca$length.content <- as.numeric(length.box.brca$length.content);


### 2. METABRIC

# Calculate gene lengths for METABRIC sample and population gene sets
gene.position.meta.length <- cbind(gene.position.meta,
                                   length = gene.position.meta$end_position - gene.position.meta$start_position + 1);
outlier.length.meta <- gene.position.meta$end_position - gene.position.meta$start_position + 1;
gene.length.meta <- gene.position.meta.all$end_position - gene.position.meta.all$start_position + 1;

# Create a subset of non-outlier genes
gene.position.meta.all.non <- gene.position.meta.all[!(gene.position.meta.all$ensembl_gene_id %in% gene.position.meta.length$ensembl_gene_id), ];
gene.non.length.meta <- gene.position.meta.all.non$end_position - gene.position.meta.all.non$start_position + 1;

# Create a data frame for the lengths and status of genes
length.box.meta <- data.frame(cbind(
    length.content = c(gene.non.length.meta, outlier.length.meta),
    status = c(rep('non', length(gene.non.length.meta)),
               rep('out', length(outlier.length.meta)))
    ));
length.box.meta <- na.omit(length.box.meta);    
length.box.meta$length.content <- as.numeric(length.box.meta$length.content);


### 3. I-SPY-2

# Calculate gene lengths for I-SPY-2 sample and population gene sets
gene.position.ispy.length <- cbind(gene.position.ispy,
                                   length = gene.position.ispy$end_position - gene.position.ispy$start_position + 1);

# Create vectors for the lengths of the sample and population gene sets
outlier.length.ispy <- gene.position.ispy.length$length;
gene.length.ispy <- gene.position.ispy.all$end_position - gene.position.ispy.all$start_position + 1;

# Create a subset of non-outlier genes
gene.position.ispy.all.non <- gene.position.ispy.all[!(gene.position.ispy.all$ensembl_gene_id %in% gene.position.ispy.length$ensembl_gene_id), ];
gene.non.length.ispy <- gene.position.ispy.all.non$end_position - gene.position.ispy.all.non$start_position + 1;

# Create a data frame for the lengths and status of genes
length.box.ispy <- data.frame(cbind(
    length.content = c(gene.non.length.ispy, outlier.length.ispy),
    status = c(rep('non', length(gene.non.length.ispy)),
               rep('out', length(outlier.length.ispy)))
    ));
length.box.ispy <- na.omit(length.box.ispy);    
length.box.ispy$length.content <- as.numeric(length.box.ispy$length.content);


### 4. MATADOR

# Calculate gene lengths for MATADOR sample and population gene sets
gene.position.metador.length <- cbind(gene.position.metador,
                                      length = gene.position.metador$end_position - gene.position.metador$start_position + 1);

outlier.length.metador <- gene.position.metador.length$length;
gene.length.metador <- gene.position.metador.all$end_position - gene.position.metador.all$start_position + 1;

# Create a subset of non-outlier genes
gene.position.metador.all.non <- gene.position.metador.all[!(gene.position.metador.all$ensembl_gene_id %in% gene.position.metador.length$ensembl_gene_id), ];
gene.non.length.metador <- gene.position.metador.all.non$end_position - gene.position.metador.all.non$start_position + 1;

# Create a data frame for the lengths and status of genes
length.box.metador <- data.frame(cbind(
    length.content = c(gene.non.length.metador, outlier.length.metador),
    status = c(rep('non', length(gene.non.length.metador)),
               rep('out', length(outlier.length.metador)))
    ));
length.box.metador <- na.omit(length.box.metador);    
length.box.metador$length.content <- as.numeric(length.box.metador$length.content);


### 5. ICGC

# Extract the start and end positions for ICGC non-outlier genes
first_number.icgc.non <- gsub("^.*:(\\d+)-.*$", "\\1", fpkm.data.icgc$loc[as.numeric(outlier.gene.fdr.all.icgc$gene)[-(as.numeric(outlier.gene.fdr.01.icgc$gene))]]);
second_number.icgc.non <- gsub("^.*-(\\d+)$", "\\1", fpkm.data.icgc$loc[as.numeric(outlier.gene.fdr.all.icgc$gene)[-(as.numeric(outlier.gene.fdr.01.icgc$gene))]]);

# Extract the start and end positions for ICGC outlier genes
first_number.icgc.out <- gsub("^.*:(\\d+)-.*$", "\\1", fpkm.data.icgc$loc[as.numeric(outlier.gene.fdr.01.icgc$gene)]);
second_number.icgc.out <- gsub("^.*-(\\d+)$", "\\1", fpkm.data.icgc$loc[as.numeric(outlier.gene.fdr.01.icgc$gene)]);

# Calculate the lengths of non-outlier and outlier genes
outlier.length.icgc <- as.numeric(second_number.icgc.out) - as.numeric(first_number.icgc.out) + 1;
gene.non.length.icgc <- as.numeric(second_number.icgc.non) - as.numeric(first_number.icgc.non) + 1;

# Create a data frame for the lengths and status of genes
length.box.icgc <- data.frame(cbind(
    length.content = c(gene.non.length.icgc, outlier.length.icgc),
    status = c(rep('non', length(gene.non.length.icgc)),
               rep('out', length(outlier.length.icgc)))
    ));
length.box.icgc <- na.omit(length.box.icgc);    
length.box.icgc$length.content <- as.numeric(length.box.icgc$length.content);



brca.out.mean <- mean(length.box.brca$length.content[length.box.brca$status == 'out']);
brca.out.sd <- sd(length.box.brca$length.content[length.box.brca$status == 'out']);
brca.out.n <- length(length.box.brca$length.content[length.box.brca$status == 'out']);
brca.non.mean <- mean(length.box.brca$length.content[length.box.brca$status == 'non']);
brca.non.sd <- sd(length.box.brca$length.content[length.box.brca$status == 'non']);
brca.non.n <- length(length.box.brca$length.content[length.box.brca$status == 'non']);

meta.out.mean <- mean(length.box.meta$length.content[length.box.meta$status == 'out']);
meta.out.sd <- sd(length.box.meta$length.content[length.box.meta$status == 'out']);
meta.non.mean <- mean(length.box.meta$length.content[length.box.meta$status == 'non']);
meta.non.sd <- sd(length.box.meta$length.content[length.box.meta$status == 'non']);
meta.out.n <- length(length.box.meta$length.content[length.box.meta$status == 'out']);
meta.non.n <- length(length.box.meta$length.content[length.box.meta$status == 'non']);

ispy.out.mean <- mean(length.box.ispy$length.content[length.box.ispy$status == 'out']);
ispy.out.sd <- sd(length.box.ispy$length.content[length.box.ispy$status == 'out']);
ispy.non.mean <- mean(length.box.ispy$length.content[length.box.ispy$status == 'non']);
ispy.non.sd <- sd(length.box.ispy$length.content[length.box.ispy$status == 'non']);
ispy.out.n <- length(length.box.ispy$length.content[length.box.ispy$status == 'out']);
ispy.non.n <- length(length.box.ispy$length.content[length.box.ispy$status == 'non']);

metador.out.mean <- mean(length.box.metador$length.content[length.box.metador$status == 'out']);
metador.out.sd <- sd(length.box.metador$length.content[length.box.metador$status == 'out']);
metador.non.mean <- mean(length.box.metador$length.content[length.box.metador$status == 'non']);
metador.non.sd <- sd(length.box.metador$length.content[length.box.metador$status == 'non']);
metador.out.n <- length(length.box.metador$length.content[length.box.metador$status == 'out']);
metador.non.n <- length(length.box.metador$length.content[length.box.metador$status == 'non']);

icgc.out.mean <- mean(length.box.icgc$length.content[length.box.icgc$status == 'out']);
icgc.out.sd <- sd(length.box.icgc$length.content[length.box.icgc$status == 'out']);
icgc.non.mean <- mean(length.box.icgc$length.content[length.box.icgc$status == 'non']);
icgc.non.sd <- sd(length.box.icgc$length.content[length.box.icgc$status == 'non']);
icgc.out.n <- length(length.box.icgc$length.content[length.box.icgc$status == 'out']);
icgc.non.n <- length(length.box.icgc$length.content[length.box.icgc$status == 'non']);


smd.length.matrix.5 <- data.frame(cbind(
    study = c(1, 2, 3, 4, 5), # 1: TCGA-BRCA, 2: METABRIC, 3: I-SPY2, 4: MATADOR, 5: ICGC
    
    # Mean of the length number of outlier genes
    mean.out = c(brca.out.mean, meta.out.mean, ispy.out.mean, metador.out.mean, icgc.out.mean), 
    # SD of the length number of outlier genes
    sd.out = c(brca.out.sd, meta.out.sd, ispy.out.sd, metador.out.sd, icgc.out.sd),
    # Number of outlier genes
    n.out = c(brca.out.n, meta.out.n, ispy.out.n, metador.out.n, icgc.out.n),
    
    # Mean of the length number of non-outlier genes
    mean.non = c(brca.non.mean, meta.non.mean, ispy.non.mean, metador.non.mean, icgc.non.mean),
    # SD of the length number of non-outlier genes
    sd.non = c(brca.non.sd, meta.non.sd, ispy.non.sd, metador.non.sd, icgc.non.sd),
    # Number of non-outlier genes
    n.non = c(brca.non.n, meta.non.n, ispy.non.n, metador.non.n, icgc.non.n)
    ));


escalc.length.matrix.5 <- escalc(
    measure = "SMD", 
    m1i = mean.out, 
    m2i = mean.non, 
    sd1i = sd.out,  
    sd2i = sd.non, 
    n1i = n.out, 
    n2i = n.non, 
    data = smd.length.matrix.5, 
    append = TRUE);


smd.metafor.length.5 <- rma.uni(yi = yi, vi = vi, data=escalc.length.matrix.5, method = 'DL');




### - 3. Exon number
library(TxDb.Hsapiens.UCSC.hg38.knownGene);

## setup transcriptDb
txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene;
## get exon locations for each gene
exons <- exonsBy(txdb,'gene');

exons.reduce <- reduce(exons);

## print number of exons for each gene 
exon.num <- sapply(exons.reduce, length);



# 1. TCGA-BRCA
#   - get entrez id
exon.num.order.brca <- exon.num[match(gene.position.brca.all.entrez$entrezgene_id, names(exon.num))];

gene.position.all.exon.brca <- data.frame(cbind(
    gene.position.brca.all.entrez,
    exon.num = exon.num.order.brca
    ));

gene.position.all.exon.non.brca <- gene.position.all.exon.brca[!(gene.position.all.exon.brca$ensembl_gene_id %in% gene.position.brca.entrez$ensembl_gene_id),];
gene.position.all.exon.outlier.brca <- gene.position.all.exon.brca[gene.position.all.exon.brca$ensembl_gene_id %in% gene.position.brca.entrez$ensembl_gene_id,];
exon.num.non.brca <- gene.position.all.exon.non.brca$exon.num;
exon.num.outlier.brca <- gene.position.all.exon.outlier.brca$exon.num;

exon.box.brca <- data.frame(cbind(
    exon.content = c(na.omit(exon.num.non.brca), na.omit(exon.num.outlier.brca)),
    status = c(rep('non', length(na.omit(exon.num.non.brca))),
               rep('out', length(na.omit(exon.num.outlier.brca))))
    ));
exon.box.brca <- na.omit(exon.box.brca);    
exon.box.brca$exon.content <- as.numeric(exon.box.brca$exon.content);


# 2. METABRIC
exon.num.order.meta <- exon.num[match(gene.position.meta.all$entrezgene_id, names(exon.num))];
gene.position.all.exon.meta <- data.frame(cbind(
    gene.position.meta.all,
    exon.num = exon.num.order.meta
    ));

gene.position.all.exon.non.meta <- gene.position.all.exon.meta[!(gene.position.all.exon.meta$ensembl_gene_id %in% gene.position.meta.order$ensembl_gene_id),];
gene.position.all.exon.outlier.meta <- gene.position.all.exon.meta[gene.position.all.exon.meta$ensembl_gene_id %in% gene.position.meta.order$ensembl_gene_id,];
exon.num.non.meta <- gene.position.all.exon.non.meta$exon.num;
exon.num.outlier.meta <- gene.position.all.exon.outlier.meta$exon.num;

exon.box.meta <- data.frame(cbind(
    exon.content = c(na.omit(exon.num.non.meta), na.omit(exon.num.outlier.meta)),
    status = c(rep('non', length(na.omit(exon.num.non.meta))),
               rep('out', length(na.omit(exon.num.outlier.meta))))
    ));
exon.box.meta <- na.omit(exon.box.meta);    
exon.box.meta$exon.content <- as.numeric(exon.box.meta$exon.content);




# 3. ISPY
exon.num.order.ispy <- exon.num[match(gene.position.ispy.all$entrezgene_id, names(exon.num))];
gene.position.all.exon.ispy <- data.frame(cbind(
    gene.position.ispy.all,
    exon.num = exon.num.order.ispy
    ));

gene.position.all.exon.non.ispy <- gene.position.all.exon.ispy[!(gene.position.all.exon.ispy$ensembl_gene_id %in% gene.position.ispy.order$ensembl_gene_id),];
gene.position.all.exon.outlier.ispy <- gene.position.all.exon.ispy[gene.position.all.exon.ispy$ensembl_gene_id %in% gene.position.ispy.order$ensembl_gene_id,];
exon.num.non.ispy <- gene.position.all.exon.non.ispy$exon.num;
exon.num.outlier.ispy <- gene.position.all.exon.outlier.ispy$exon.num;

exon.box.ispy <- data.frame(cbind(
    exon.content = c(na.omit(exon.num.non.ispy), na.omit(exon.num.outlier.ispy)),
    status = c(rep('non', length(na.omit(exon.num.non.ispy))),
               rep('out', length(na.omit(exon.num.outlier.ispy))))
    ));
exon.box.ispy <- na.omit(exon.box.ispy);    
exon.box.ispy$exon.content <- as.numeric(exon.box.ispy$exon.content);


# 4. MATADOR
exon.num.order.metador <- exon.num[match(gene.position.metador.all.entrez$entrezgene_id, names(exon.num))];

gene.position.all.exon.metador <- data.frame(cbind(
    gene.position.metador.all.entrez,
    exon.num = exon.num.order.metador
    ));

gene.position.all.exon.non.metador <- gene.position.all.exon.metador[!(gene.position.all.exon.metador$ensembl_gene_id %in% gene.position.metador.entrez$ensembl_gene_id),];
gene.position.all.exon.outlier.metador <- gene.position.all.exon.metador[gene.position.all.exon.metador$ensembl_gene_id %in% gene.position.metador.entrez$ensembl_gene_id,];
exon.num.non.metador <- gene.position.all.exon.non.metador$exon.num;
exon.num.outlier.metador <- gene.position.all.exon.outlier.metador$exon.num;

exon.box.metador <- data.frame(cbind(
    exon.content = c(na.omit(exon.num.non.metador), na.omit(exon.num.outlier.metador)),
    status = c(rep('non', length(na.omit(exon.num.non.metador))),
               rep('out', length(na.omit(exon.num.outlier.metador))))
    ));
exon.box.metador <- na.omit(exon.box.metador);    
exon.box.metador$exon.content <- as.numeric(exon.box.metador$exon.content);




# 5. ICGC
gene.list <- fpkm.data.icgc$Ensembl[as.numeric(outlier.gene.fdr.all.icgc$gene)];
ensembl <- biomaRt:::useEnsembl(biomart = "ensembl", 
                     dataset = "hsapiens_gene_ensembl", 
                     mirror = "useast");
ensembl <- biomaRt:::useDataset(dataset = "hsapiens_gene_ensembl", mart = ensembl);
gene.position.entrez <- biomaRt:::getBM(attributes = c('ensembl_gene_id', 'hgnc_symbol', 'chromosome_name',
                                                'start_position', 'end_position', 'band', "gene_biotype", "entrezgene_id"),
                                 filters = 'ensembl_gene_id', 
                                 values = gene.list, 
                                 mart = ensembl);
gene.position.icgc.all.entrez <- gene.position.entrez;

gene.list <- fpkm.data$Ensembl[as.numeric(outlier.gene.fdr.01.icgc$gene)]
gene.position.entrez <- biomaRt:::getBM(attributes = c('ensembl_gene_id', 'hgnc_symbol', 'chromosome_name',
                                                'start_position', 'end_position', 'band', "gene_biotype", "entrezgene_id"),
                                 filters = 'ensembl_gene_id', 
                                 values = gene.list, 
                                 mart = ensembl);
gene.position.icgc.entrez <- gene.position.entrez;

exon.num.order.icgc <- exon.num[match(gene.position.icgc.all.entrez$entrezgene_id, names(exon.num))];
gene.position.all.exon.icgc <- data.frame(cbind(
    gene.position.icgc.all.entrez,
    exon.num = exon.num.order.icgc
    ));

gene.position.all.exon.non.icgc <- gene.position.all.exon.icgc[!(gene.position.all.exon.icgc$ensembl_gene_id %in% gene.position.icgc.entrez$ensembl_gene_id),];
gene.position.all.exon.outlier.icgc <- gene.position.all.exon.icgc[gene.position.all.exon.icgc$ensembl_gene_id %in% gene.position.icgc.entrez$ensembl_gene_id,];
exon.num.non.icgc <- gene.position.all.exon.non.icgc$exon.num;
exon.num.outlier.icgc <- gene.position.all.exon.outlier.icgc$exon.num;

exon.box.icgc <- data.frame(cbind(
    exon.content = c(na.omit(exon.num.non.icgc), na.omit(exon.num.outlier.icgc)),
    status = c(rep('non', length(na.omit(exon.num.non.icgc))),
               rep('out', length(na.omit(exon.num.outlier.icgc))))
    ));
exon.box.icgc <- na.omit(exon.box.icgc);    
exon.box.icgc$exon.content <- as.numeric(exon.box.icgc$exon.content);


# Meta-analysis

brca.out.mean <- mean(exon.box.brca$exon.content[exon.box.brca$status == 'out']);
brca.out.sd <- sd(exon.box.brca$exon.content[exon.box.brca$status == 'out']);
brca.out.n <- length(exon.box.brca$exon.content[exon.box.brca$status == 'out']);
brca.non.mean <- mean(exon.box.brca$exon.content[exon.box.brca$status == 'non']);
brca.non.sd <- sd(exon.box.brca$exon.content[exon.box.brca$status == 'non']);
brca.non.n <- length(exon.box.brca$exon.content[exon.box.brca$status == 'non']);

meta.out.mean <- mean(exon.box.meta$exon.content[exon.box.meta$status == 'out']);
meta.out.sd <- sd(exon.box.meta$exon.content[exon.box.meta$status == 'out']);
meta.non.mean <- mean(exon.box.meta$exon.content[exon.box.meta$status == 'non']);
meta.non.sd <- sd(exon.box.meta$exon.content[exon.box.meta$status == 'non']);
meta.out.n <- length(exon.box.meta$exon.content[exon.box.meta$status == 'out']);
meta.non.n <- length(exon.box.meta$exon.content[exon.box.meta$status == 'non']);

ispy.out.mean <- mean(exon.box.ispy$exon.content[exon.box.ispy$status == 'out']);
ispy.out.sd <- sd(exon.box.ispy$exon.content[exon.box.ispy$status == 'out']);
ispy.non.mean <- mean(exon.box.ispy$exon.content[exon.box.ispy$status == 'non']);
ispy.non.sd <- sd(exon.box.ispy$exon.content[exon.box.ispy$status == 'non']);
ispy.out.n <- length(exon.box.ispy$exon.content[exon.box.ispy$status == 'out']);
ispy.non.n <- length(exon.box.ispy$exon.content[exon.box.ispy$status == 'non']);

metador.out.mean <- mean(exon.box.metador$exon.content[exon.box.metador$status == 'out']);
metador.out.sd <- sd(exon.box.metador$exon.content[exon.box.metador$status == 'out']);
metador.non.mean <- mean(exon.box.metador$exon.content[exon.box.metador$status == 'non']);
metador.non.sd <- sd(exon.box.metador$exon.content[exon.box.metador$status == 'non']);
metador.out.n <- length(exon.box.metador$exon.content[exon.box.metador$status == 'out']);
metador.non.n <- length(exon.box.metador$exon.content[exon.box.metador$status == 'non']);

icgc.out.mean <- mean(exon.box.icgc$exon.content[exon.box.icgc$status == 'out']);
icgc.out.sd <- sd(exon.box.icgc$exon.content[exon.box.icgc$status == 'out']);
icgc.non.mean <- mean(exon.box.icgc$exon.content[exon.box.icgc$status == 'non']);
icgc.non.sd <- sd(exon.box.icgc$exon.content[exon.box.icgc$status == 'non']);
icgc.out.n <- length(exon.box.icgc$exon.content[exon.box.icgc$status == 'out']);
icgc.non.n <- length(exon.box.icgc$exon.content[exon.box.icgc$status == 'non']);



smd.exon.matrix.5 <- data.frame(cbind(
    study = c(1, 2, 3, 4, 5), # 1: TCGA-BRCA, 2: METABRIC, 3: I-SPY2, 4: MATADOR, 5: ICGC
    
    # Mean of the exon number of outlier genes
    mean.out = c(brca.out.mean, meta.out.mean, ispy.out.mean, metador.out.mean, icgc.out.mean), 
    # SD of the exon number of outlier genes
    sd.out = c(brca.out.sd, meta.out.sd, ispy.out.sd, metador.out.sd, icgc.out.sd),
    # Number of outlier genes
    n.out = c(brca.out.n, meta.out.n, ispy.out.n, metador.out.n, icgc.out.n),
    
    # Mean of the exon number of non-outlier genes
    mean.non = c(brca.non.mean, meta.non.mean, ispy.non.mean, metador.non.mean, icgc.non.mean),
    # SD of the exon number of non-outlier genes
    sd.non = c(brca.non.sd, meta.non.sd, ispy.non.sd, metador.non.sd, icgc.non.sd),
    # Number of non-outlier genes
    n.non = c(brca.non.n, meta.non.n, ispy.non.n, metador.non.n, icgc.non.n)
    ));


# 2. metafor
escalc.exon.matrix.5 <- escalc(
    measure = "SMD", 
    m1i = mean.out, 
    m2i = mean.non, 
    sd1i = sd.out,  
    sd2i = sd.non, 
    n1i = n.out, 
    n2i = n.non, 
    data = smd.exon.matrix.5, 
    append = TRUE);


smd.metafor.exon.5 <- rma.uni(yi = yi, vi = vi, data=escalc.exon.matrix.5,  method = 'DL');





### - 4) RNA abundance

# 1. TCGA-BRCA
outlier.rna.brca <- apply(fpkm.tumor.symbol.filter.brca[rownames(outlier.gene.fdr.01.brca),patient.part.brca], 1, median);

gene.non.rna.brca <- apply(fpkm.tumor.symbol.filter.brca[!(rownames(fpkm.tumor.symbol.filter.brca) %in% rownames(outlier.gene.fdr.01.brca)), patient.part.brca], 1, median)

rna.box.brca <- data.frame(cbind(
    rna.content = c(gene.non.rna.brca, outlier.rna.brca),
    status = c(rep('non', length(gene.non.rna.brca)),
               rep('out', length(outlier.rna.brca)))
    ));
rna.box.brca <- na.omit(rna.box.brca);    
rna.box.brca$rna.content <- as.numeric(rna.box.brca$rna.content);


# 2. METABRIC
outlier.rna.meta <- apply(fpkm.tumor.symbol.filter.meta[rownames(outlier.gene.fdr.01.meta),patient.part.meta], 1, median);

gene.non.rna.meta <- apply(fpkm.tumor.symbol.filter.meta[!(rownames(fpkm.tumor.symbol.filter.meta) %in% rownames(outlier.gene.fdr.01.meta)), patient.part.meta], 1, median)

rna.box.meta <- data.frame(cbind(
    rna.content = c(gene.non.rna.meta, outlier.rna.meta),
    status = c(rep('non', length(gene.non.rna.meta)),
               rep('out', length(outlier.rna.meta)))
    ));
rna.box.meta <- na.omit(rna.box.meta);    
rna.box.meta$rna.content <- as.numeric(rna.box.meta$rna.content);


# 3. I-SPY-2
outlier.rna.ispy <- apply(fpkm.tumor.symbol.filter.ispy[rownames(gene.rank.order.cosine.observed.p.value.fdr.01.ispy),patient.part.ispy], 1, median);

gene.non.rna.ispy <- apply(fpkm.tumor.symbol.filter.ispy[!(rownames(fpkm.tumor.symbol.filter.ispy) %in% rownames(gene.rank.order.cosine.observed.p.value.fdr.01.ispy)), patient.part.ispy], 1, median)

rna.box.ispy <- data.frame(cbind(
    rna.content = c(gene.non.rna.ispy, outlier.rna.ispy),
    status = c(rep('non', length(gene.non.rna.ispy)),
               rep('out', length(outlier.rna.ispy)))
    ));
rna.box.ispy <- na.omit(rna.box.ispy);    
rna.box.ispy$rna.content <- as.numeric(rna.box.ispy$rna.content);


# 4. MATADOR
outlier.rna.metador <- apply(fpkm.tumor.symbol.filter.metador[rownames(gene.rank.order.cosine.observed.p.value.fdr.01.metador),patient.part.metador], 1, median);

gene.non.rna.metador <- apply(fpkm.tumor.symbol.filter.metador[!(rownames(fpkm.tumor.symbol.filter.metador) %in% rownames(gene.rank.order.cosine.observed.p.value.fdr.01.metador)), patient.part.metador], 1, median)

rna.box.metador <- data.frame(cbind(
    rna.content = c(gene.non.rna.metador, outlier.rna.metador),
    status = c(rep('non', length(gene.non.rna.metador)),
               rep('out', length(outlier.rna.metador)))
    ));
rna.box.metador <- na.omit(rna.box.metador);    
rna.box.metador$rna.content <- as.numeric(rna.box.metador$rna.content);


# 5. ICGC
outlier.rna.icgc <- apply(fpkm.tumor.symbol.filter.icgc[gene.rank.order.cosine.observed.p.value.max.filter.fdr.01.icgc$gene,patient.part.icgc], 1, median);

gene.non.rna.icgc <- apply(fpkm.tumor.symbol.filter.icgc[!(rownames(fpkm.tumor.symbol.filter.icgc) %in% gene.rank.order.cosine.observed.p.value.max.filter.fdr.01.icgc$gene), patient.part.icgc], 1, median)

rna.box.icgc <- data.frame(cbind(
    rna.content = c(gene.non.rna.icgc, outlier.rna.icgc),
    status = c(rep('non', length(gene.non.rna.icgc)),
               rep('out', length(outlier.rna.icgc)))
    ));
rna.box.icgc <- na.omit(rna.box.icgc);    
rna.box.icgc$rna.content <- as.numeric(rna.box.icgc$rna.content);




brca.out.mean <- mean(rna.box.brca$rna.content[rna.box.brca$status == 'out']);
brca.out.sd <- sd(rna.box.brca$rna.content[rna.box.brca$status == 'out']);
brca.out.n <- length(rna.box.brca$rna.content[rna.box.brca$status == 'out']);
brca.non.mean <- mean(rna.box.brca$rna.content[rna.box.brca$status == 'non']);
brca.non.sd <- sd(rna.box.brca$rna.content[rna.box.brca$status == 'non']);
brca.non.n <- length(rna.box.brca$rna.content[rna.box.brca$status == 'non']);

meta.out.mean <- mean(rna.box.meta$rna.content[rna.box.meta$status == 'out']);
meta.out.sd <- sd(rna.box.meta$rna.content[rna.box.meta$status == 'out']);
meta.non.mean <- mean(rna.box.meta$rna.content[rna.box.meta$status == 'non']);
meta.non.sd <- sd(rna.box.meta$rna.content[rna.box.meta$status == 'non']);
meta.out.n <- length(rna.box.meta$rna.content[rna.box.meta$status == 'out']);
meta.non.n <- length(rna.box.meta$rna.content[rna.box.meta$status == 'non']);

ispy.out.mean <- mean(rna.box.ispy$rna.content[rna.box.ispy$status == 'out']);
ispy.out.sd <- sd(rna.box.ispy$rna.content[rna.box.ispy$status == 'out']);
ispy.non.mean <- mean(rna.box.ispy$rna.content[rna.box.ispy$status == 'non']);
ispy.non.sd <- sd(rna.box.ispy$rna.content[rna.box.ispy$status == 'non']);
ispy.out.n <- length(rna.box.ispy$rna.content[rna.box.ispy$status == 'out']);
ispy.non.n <- length(rna.box.ispy$rna.content[rna.box.ispy$status == 'non']);

metador.out.mean <- mean(rna.box.metador$rna.content[rna.box.metador$status == 'out']);
metador.out.sd <- sd(rna.box.metador$rna.content[rna.box.metador$status == 'out']);
metador.non.mean <- mean(rna.box.metador$rna.content[rna.box.metador$status == 'non']);
metador.non.sd <- sd(rna.box.metador$rna.content[rna.box.metador$status == 'non']);
metador.out.n <- length(rna.box.metador$rna.content[rna.box.metador$status == 'out']);
metador.non.n <- length(rna.box.metador$rna.content[rna.box.metador$status == 'non']);

icgc.out.mean <- mean(rna.box.icgc$rna.content[rna.box.icgc$status == 'out']);
icgc.out.sd <- sd(rna.box.icgc$rna.content[rna.box.icgc$status == 'out']);
icgc.non.mean <- mean(rna.box.icgc$rna.content[rna.box.icgc$status == 'non']);
icgc.non.sd <- sd(rna.box.icgc$rna.content[rna.box.icgc$status == 'non']);
icgc.out.n <- length(rna.box.icgc$rna.content[rna.box.icgc$status == 'out']);
icgc.non.n <- length(rna.box.icgc$rna.content[rna.box.icgc$status == 'non']);



smd.rna.matrix.5 <- data.frame(cbind(
    study = c(1, 2, 3, 4, 5), # 1: TCGA-BRCA, 2: METABRIC, 3: I-SPY2, 4: MATADOR, 5: ICGC
    
    # Mean of the rna number of outlier genes
    mean.out = c(brca.out.mean, meta.out.mean, ispy.out.mean, metador.out.mean, icgc.out.mean), 
    # SD of the rna number of outlier genes
    sd.out = c(brca.out.sd, meta.out.sd, ispy.out.sd, metador.out.sd, icgc.out.sd),
    # Number of outlier genes
    n.out = c(brca.out.n, meta.out.n, ispy.out.n, metador.out.n, icgc.out.n),
    
    # Mean of the rna number of non-outlier genes
    mean.non = c(brca.non.mean, meta.non.mean, ispy.non.mean, metador.non.mean, icgc.non.mean),
    # SD of the rna number of non-outlier genes
    sd.non = c(brca.non.sd, meta.non.sd, ispy.non.sd, metador.non.sd, icgc.non.sd),
    # Number of non-outlier genes
    n.non = c(brca.non.n, meta.non.n, ispy.non.n, metador.non.n, icgc.non.n)
    ));


escalc.rna.matrix.5 <- escalc(
    measure = "SMD", 
    m1i = mean.out, 
    m2i = mean.non, 
    sd1i = sd.out,  
    sd2i = sd.non, 
    n1i = n.out, 
    n2i = n.non, 
    data = smd.rna.matrix.5, 
    append = TRUE);


smd.metafor.rna.5 <- rma.uni(yi = yi, vi = vi, data=escalc.rna.matrix.5,  method = 'DL');





### All results ###


metafor.smd.all.matrix.5 <- data.frame(cbind(
    p.value = c(smd.metafor.exon.5$pval,
                smd.metafor.length.5$pval,
                smd.metafor.gc.5$pval,
                smd.metafor.rna.5$pval),
    odd = c(exp(smd.metafor.exon.5$beta),
            exp(smd.metafor.length.5$beta),
            exp(smd.metafor.gc.5$beta),
            exp(smd.metafor.rna.5$beta)),
    ci.min = c(exp(smd.metafor.exon.5$ci.lb),
               exp(smd.metafor.length.5$ci.lb),
               exp(smd.metafor.gc.5$ci.lb),
               exp(smd.metafor.rna.5$ci.lb)),
    ci.max = c(exp(smd.metafor.exon.5$ci.ub),
               exp(smd.metafor.length.5$ci.ub),
               exp(smd.metafor.gc.5$ci.ub),
               exp(smd.metafor.rna.5$ci.ub))
    ));

metafor.smd.all.matrix.label.5 <- data.frame(cbind(
    metafor.smd.all.matrix.5,
    labels = as.factor(c('Exon number', 'Gene length', 'GC content', 'RNA abundance'))
    ));

# FDR correction for the SMD matrix
metafor.smd.all.matrix.label.fdr.5 <- p.adjust(metafor.smd.all.matrix.label.5$p.value, method = 'BH');



### PLOTTING RESULTS ############################################################

# Perform FDR correction on the p-values
fdr.metafor.schr <- p.adjust(metafor.chr.odd.ci.p.data.label.rev.5$p.value, method = 'BH');
dot.colours <- vector(length=nrow(metafor.chr.odd.ci.p.data.label.rev.5));
dot.colours <- rep('grey70',nrow(metafor.chr.odd.ci.p.data.label.rev.5));
dot.colours[fdr.metafor.schr.5 < 0.01 & metafor.chr.odd.ci.p.data.label.rev.5$odd < 1] <- 'dodgerblue2';
dot.colours[fdr.metafor.schr.5 < 0.01 & metafor.chr.odd.ci.p.data.label.rev.5$odd > 1] <- 'red';

metafor.all.segplot.multi.5 <- BoutrosLab.plotting.general::create.segplot(
    formula = labels ~ log2(ci.min) + log2(ci.max),
    data = metafor.chr.odd.ci.p.data.label.rev.5,
    centers = log2(metafor.chr.odd.ci.p.data.label.rev.5$odd),
    main.cex = 0,
    ylab.label = expression('Chromosome'),
    yaxis.lab = metafor.chr.odd.ci.p.data.label.rev.5$labels,
    yaxis.fontface = 1,
    xlab.cex = 0,
    xlab.label = NULL,
    ylab.cex = 1.3,
    yaxis.cex = 1,
    xaxis.cex = 0,
    xlimits = c(-1.1, 1.3),
    xaxis.lab = c('0.5', '1', '2'),
    xat = seq(-1, 1, 1),
    xaxis.fontface = 1,
    yaxis.tck = c(0.2, 0),
    xaxis.tck = c(0.2, 0),
    segments.col = dot.colours,
    abline.v = 0,
    abline.lty = 3,
    add.rectangle = TRUE,
    xleft.rectangle = -1.7,
    xright.rectangle = 13,
    ybottom.rectangle = seq(1.5, 23.5, 2),
    ytop.rectangle = seq(2.5, 24.5, 2),
    # set rectangle colour
    col.rectangle = "grey",
    alpha.rectangle = 0.25,
    disable.factor.sorting = TRUE
    )
metafor.all.segplot.multi.5;


dot.colours <- vector(length=nrow(metafor.smd.all.matrix.label.5));
dot.colours <- rep('grey70',nrow(metafor.smd.all.matrix.label.5));
dot.colours[metafor.smd.all.matrix.label.fdr.5 < 0.01 & metafor.smd.all.matrix.label.5$odd < 1] <- 'dodgerblue2';
dot.colours[metafor.smd.all.matrix.label.fdr.5 < 0.01 & metafor.smd.all.matrix.label.5$odd > 1] <- 'red';
metafor.smd.all.matrix.label.segplot.multi.5 <- BoutrosLab.plotting.general::create.segplot(
    formula = labels ~ log2(ci.min) + log2(ci.max),
    data = metafor.smd.all.matrix.label.5,
    # add middle dots
    centers = log2(metafor.smd.all.matrix.label.5$odd),
    main.cex = 1.4,
    ylab.label = NULL,
    yaxis.lab = metafor.smd.all.matrix.label.5$labels,
    yaxis.fontface = 1,
    xlab.cex = 1.3,
    xlab.label = expression('Odds Ratio'),
    ylab.cex = 1.3,
    yaxis.cex = 1.1,
    xaxis.cex = 1,
    xlimits = c(-1.1, 1.3),
    xaxis.lab = c('0.5', '1', '2'),
    xat = seq(-1, 1, 1),
    xaxis.fontface = 1,
    yaxis.tck = c(0.2, 0),
    xaxis.tck = c(0.2, 0),
    segments.col = dot.colours,
    abline.v = 0,
    abline.lty = 3,
    add.rectangle = TRUE,
    xleft.rectangle = -1.7,
    xright.rectangle = 13,
    ybottom.rectangle = seq(1.5, 23.5, 2),
    ytop.rectangle = seq(2.5, 24.5, 2),
    # set rectangle colour
    col.rectangle = "grey",
    # set rectangle alpha (transparency)
    alpha.rectangle = 0.25,
    disable.factor.sorting = TRUE
)
metafor.smd.all.matrix.label.segplot.multi.5;



metafor.multi.chr.smd.5 <- create.multipanelplot(
    list(metafor.all.segplot.multi.5, metafor.smd.all.matrix.label.segplot.multi.5),
    resolution = 300,
    layout.height = 2,
    layout.width = 1,
    ylab.cex = 0,
    xlab.cex = 0,
    layout.skip = c(FALSE, FALSE),
    plot.objects.heights = c(5, 2),
    ylab.axis.padding = -7,
    xlab.axis.padding = 1,
    y.spacing = -3.5,
    bottom.padding = -0.5,
    top.padding = -0.5
    );

metafor.multi.chr.smd.5;


# Save the multi plot as a PDF
pdf(
    file = generate.filename(
        'metafor.multi.chr.smd.5', 
        'multipanel', 
        'pdf'
        ), 
    width = 5, 
    height = 6.5
    );
metafor.multi.chr.smd.5;
dev.off();

# Save the multi plot as a PNG
png(
    file = generate.filename(
        'metafor.multi.chr.smd.5', 
        'multipanel', 
        'png'
        ), 
    width = 5, 
    height = 6.5,
    unit = 'in', 
    res = 1200
    );
metafor.multi.chr.smd.5;
dev.off();
