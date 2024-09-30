### HISTORY #####################################################################
# This script processes p-value data from multiple datasets, combines them
# using Fisher's method, and calculates the false discovery rate (FDR) for
# further analysis. The results are then used in downstream analysis, including
# gene set enrichment analysis (GSEA).
# The plot is processed through Cytoscape.
# Date: 2024-08-23

### PREAMBLE ####################################################################
# Load necessary libraries
library(BoutrosLab.utilities);
library(dplyr);
library(tidyr);
library(poolr);

# Source the helper library
source(here::here('common_functions.R'));

# Load the datafile
load(file.path(get.outlier.data.dir(), '2024-09-10_Figure1.rda'));

outlier.gene.fdr.all.icgc$Symbol <- fpkm.data.icgc$Name[as.numeric(outlier.gene.fdr.all.icgc$gene)];
outlier.gene.fdr.all.ispy$Symbol <- rownames(outlier.gene.fdr.all.ispy);
outlier.gene.fdr.all.meta$Symbol <- fpkm.tumor.symbol.filter.meta.symbol[rownames(outlier.gene.fdr.all.meta), ]$Symbol;
outlier.gene.fdr.all.matador$Symbol <- substr(rownames(outlier.gene.fdr.all.matador), 17, nchar(rownames(outlier.gene.fdr.all.matador)));
outlier.gene.fdr.all.brca$Symbol <- fpkm.tumor.symbol.filter.brca[rownames(outlier.gene.fdr.all.brca), ]$Symbol;


icgc.all.fdr <- outlier.gene.fdr.all.icgc[, c('fdr', 'Symbol')];
ispy.all.fdr <- outlier.gene.fdr.all.ispy[, c('fdr', 'Symbol')];
meta.all.fdr <- outlier.gene.fdr.all.meta[, c('fdr', 'Symbol')];
metador.all.fdr <- outlier.gene.fdr.all.matador[, c('fdr', 'Symbol')];
brca.all.fdr <- outlier.gene.fdr.all.brca[, c('fdr', 'Symbol')];


colnames(icgc.all.fdr) <- c('fdr.icgc', 'Symbol');
colnames(ispy.all.fdr) <- c('fdr.ispy', 'Symbol');
colnames(meta.all.fdr) <- c('fdr.meta', 'Symbol');
colnames(metador.all.fdr) <- c('fdr.metador', 'Symbol');
colnames(brca.all.fdr) <- c('fdr.brca', 'Symbol');


brca.all.fdr$Symbol <- as.character(brca.all.fdr$Symbol);
metador.all.fdr$Symbol <- as.character(metador.all.fdr$Symbol);
meta.all.fdr$Symbol <- as.character(meta.all.fdr$Symbol);
ispy.all.fdr$Symbol <- as.character(ispy.all.fdr$Symbol);
icgc.all.fdr$Symbol <- as.character(icgc.all.fdr$Symbol);


combined.df <- icgc.all.fdr %>%
    rename(fdr.icgc = fdr.icgc) %>%
    full_join(ispy.all.fdr %>% rename(fdr.ispy = fdr.ispy), by = 'Symbol') %>%
    full_join(meta.all.fdr %>% rename(fdr.meta = fdr.meta), by = 'Symbol') %>%
    full_join(metador.all.fdr %>% rename(fdr.metador = fdr.metador), by = 'Symbol') %>%
    full_join(brca.all.fdr %>% rename(fdr.brca = fdr.brca), by = 'Symbol');

icgc.all.pvalue <- outlier.gene.fdr.all.icgc[, c('obs.p.value', 'Symbol')];
ispy.all.pvalue <- outlier.gene.fdr.all.ispy[, c('new.p.value', 'Symbol')];
meta.all.pvalue <- outlier.gene.fdr.all.meta[, c('new.p.value', 'Symbol')];
metador.all.pvalue <- outlier.gene.fdr.all.matador[, c('new.p.value', 'Symbol')];
brca.all.pvalue <- outlier.gene.fdr.all.brca[, c('new.p.value', 'Symbol')];


colnames(icgc.all.pvalue) <- c('pvalue.icgc', 'Symbol');
colnames(ispy.all.pvalue) <- c('pvalue.ispy', 'Symbol');
colnames(meta.all.pvalue) <- c('pvalue.meta', 'Symbol');
colnames(metador.all.pvalue) <- c('pvalue.metador', 'Symbol');
colnames(brca.all.pvalue) <- c('pvalue.brca', 'Symbol');


combined.df <- icgc.all.pvalue %>%
    rename(pvalue.icgc = pvalue.icgc) %>%
    full_join(ispy.all.pvalue %>% rename(pvalue.ispy = pvalue.ispy), by = 'Symbol') %>%
    full_join(meta.all.pvalue %>% rename(pvalue.meta = pvalue.meta), by = 'Symbol') %>%
    full_join(metador.all.pvalue %>% rename(pvalue.metador = pvalue.metador), by = 'Symbol') %>%
    full_join(brca.all.pvalue %>% rename(pvalue.brca = pvalue.brca), by = 'Symbol');


all.pvalue.df <- combined.df[, c(1, 3, 4, 5, 6, 2)];


combine.fisher.pvalue.all <- apply(all.pvalue.df[, 1:5], 1, function(x) {
    fisher(na.omit(x))$p
    });
combine.fisher.pvalue.all.fdr <- p.adjust(combine.fisher.pvalue.all, method = 'BH');


kegg.Pathways <- msigdbr(species = 'Homo sapiens', category = 'C2', subcategory = 'CP:KEGG');
kegg.List.name <- split(kegg.Pathways$gene_symbol, kegg.Pathways$gs_name);


names(combine.fisher.pvalue.all.fdr) <- all.pvalue.df$Symbol;
combine.fisher.pvalue.all.fdr.sort <- sort(combine.fisher.pvalue.all.fdr);
combine.fisher.pvalue.all.fdr.sort.log <- -log10(combine.fisher.pvalue.all.fdr.sort);

# GSEA
set.seed(42);
fgsea.Results.kegg.name.p.combine <- fgsea(
    pathways = kegg.List.name,
    stats = combine.fisher.pvalue.all.fdr.sort.log,
    minSize = 3,
    maxSize = 1000,
    scoreType = 'pos'
    );


top.Pathways.kegg.name.p.combine <- fgsea.Results.kegg.name.p.combine[order(fgsea.Results.kegg.name.p.combine$padj, decreasing = FALSE), ];
top.Pathways.kegg.name.p.combine.output <- data.frame(
    'GOID.PathwayID' = paste('KEGG:', substr(top.Pathways.kegg.name.p.combine$pathway, 6, nchar(top.Pathways.kegg.name.p.combine$pathway)), sep = ''),
    'P.value' = top.Pathways.kegg.name.p.combine$padj
    );
top.Pathways.kegg.name.p.combine.output.05 <- top.Pathways.kegg.name.p.combine.output[top.Pathways.kegg.name.p.combine.output$P.value < 0.05, ];


write.table(
    top.Pathways.kegg.name.p.combine.output.05,
    file = 'combine.fisher.pvalue.all.fdr.output.05.txt',
    sep = '\t',
    row.names = FALSE,
    col.names = TRUE
    );

save.session.profile(file.path('output', 'Figure1h.txt'));
