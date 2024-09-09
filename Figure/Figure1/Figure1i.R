### HISTORY #####################################################################
# This script generates a Manhattan plot of all genes from five datasets, 
# highlighting the distribution of outlier genes across chromosomes.
# The analysis is connected to Figure 1f,g.
# Date: 2024-08-16

library(BoutrosLab.plotting.general);




library(dplyr)
library(tidyr)
library(poolr)

outlier.gene.fdr.all.icgc.symbol <- outlier.gene.fdr.all.icgc;
outlier.gene.fdr.all.icgc.symbol$Symbol <- fpkm.data.icgc$Name[as.numeric(outlier.gene.fdr.all.icgc.symbol$gene)];

outlier.gene.fdr.all.ispy.symbol <- outlier.gene.fdr.all.ispy;
outlier.gene.fdr.all.ispy.symbol$Symbol <- rownames(outlier.gene.fdr.all.ispy.symbol);

outlier.gene.fdr.all.meta.symbol <- outlier.gene.fdr.all.meta;
outlier.gene.fdr.all.meta.symbol$Symbol <- fpkm.tumor.symbol.filter.meta.symbol[rownames(outlier.gene.fdr.all.meta.symbol),]$Symbol;

outlier.gene.fdr.all.matador.symbol <- outlier.gene.fdr.all.matador;
outlier.gene.fdr.all.matador.symbol$Symbol <- substr(rownames(outlier.gene.fdr.all.matador.symbol), 17, nchar(rownames(outlier.gene.fdr.all.matador.symbol)))

outlier.gene.fdr.all.brca.symbol <- outlier.gene.fdr.all.brca;
outlier.gene.fdr.all.brca.symbol$Symbol <- fpkm.tumor.symbol.filter.brca[rownames(outlier.gene.fdr.all.brca.symbol),]$Symbol



# p-value combine and then multiple testing correction
icgc.all.pvalue <- outlier.gene.fdr.all.icgc.symbol[, c('obs.p.value', 'Symbol')];
ispy.all.pvalue <- outlier.gene.fdr.all.ispy.symbol[, c('new.p.value', 'Symbol')];
meta.all.pvalue <- outlier.gene.fdr.all.meta.symbol[, c('new.p.value', 'Symbol')];
metador.all.pvalue <- outlier.gene.fdr.all.matador.symbol[, c('new.p.value', 'Symbol')];
brca.all.pvalue <- outlier.gene.fdr.all.brca.symbol[, c('new.p.value', 'Symbol')];

colnames(icgc.all.pvalue) <- c("pvalue_icgc", "Symbol");
colnames(ispy.all.pvalue) <- c("pvalue_ispy", "Symbol");
colnames(meta.all.pvalue) <- c("pvalue_meta", "Symbol");
colnames(metador.all.pvalue) <- c("pvalue_metador", "Symbol");
colnames(brca.all.pvalue) <- c("pvalue_brca", "Symbol");


combined_df <- icgc.all.pvalue %>%
    rename(pvalue_icgc = pvalue_icgc) %>%
    full_join(ispy.all.pvalue %>% rename(pvalue_ispy = pvalue_ispy), by = "Symbol") %>%
    full_join(meta.all.pvalue %>% rename(pvalue_meta = pvalue_meta), by = "Symbol") %>%
    full_join(metador.all.pvalue %>% rename(pvalue_metador = pvalue_metador), by = "Symbol") %>%
    full_join(brca.all.pvalue %>% rename(pvalue_brca = pvalue_brca), by = "Symbol")


all.pvalue.df <- combined_df[,c(1, 3, 4, 5, 6, 2)];
combine.fisher.pvalue.all <- apply(all.pvalue.df[,1:5], 1, function(x) { fisher(na.omit(x))$p});
combine.fisher.pvalue.all.fdr <- p.adjust(combine.fisher.pvalue.all, method = 'BH');


names(combine.fisher.pvalue.all.fdr) <- all.pvalue.df$Symbol;
combine.fisher.pvalue.all.fdr.sort <- sort(combine.fisher.pvalue.all.fdr);
combine.fisher.pvalue.all.fdr.sort.log <- -log10(combine.fisher.pvalue.all.fdr.sort);



gene.position.ispy.all.location <- gene.position.ispy.all[,2:5];
gene.position.meta.all.location <- gene.position.meta.all[,2:5];
gene.position.brca.all.location <- gene.position.brca.all[,2:5];
gene.position.metador.all.location <- gene.position.metador.all[,2:5];
gene.position.icgc.all.location <- gene.position.icgc.all[,2:5];

all.gene.location <- rbind(
    gene.position.meta.all.location,
    gene.position.brca.all.location,
    gene.position.ispy.all.location,
    gene.position.metador.all.location,
    gene.position.icgc.all.location
    );

all.gene.location.filter <- all.gene.location[!grepl("^CHR", all.gene.location$chromosome_name), ]
all.gene.location.filter.nodup <- all.gene.location.filter[!duplicated(all.gene.location.filter$hgnc_symbol), ]


chr.name <- c('1', '2', '3', '4', '5', '6', '7', '8', '9', '10', '11', '12', '13', '14', '15', '16', '17', '18', '19', '20', '21', '22', 'MT', 'X', 'Y');

all.gene.location.filter.nodup.order <- all.gene.location.filter.nodup[match(names(combine.fisher.pvalue.all.fdr.sort.log), all.gene.location.filter.nodup$hgnc_symbol),];

all.gene.location.filter.nodup.order$chromosome_name[all.gene.location.filter.nodup.order$chromosome_name == 'X'] <- 23;
all.gene.location.filter.nodup.order$chromosome_name[all.gene.location.filter.nodup.order$chromosome_name == 'Y'] <- NA;
all.gene.location.filter.nodup.order$chromosome_name[all.gene.location.filter.nodup.order$chromosome_name == 'MT'] <- NA;
all.gene.location.filter.nodup.order.fdr <- data.frame(cbind(
    all.gene.location.filter.nodup.order,
    fdr = combine.fisher.pvalue.all.fdr.sort.log)
    );
all.gene.location.filter.nodup.order.fdr.na <- na.omit(all.gene.location.filter.nodup.order.fdr);
    


all.gene.location.filter.nodup.order.fdr.na <- all.gene.location.filter.nodup.order.fdr.na[order(as.numeric(all.gene.location.filter.nodup.order.fdr.na$start_position)),];
all.gene.location.filter.nodup.order.fdr.na.chr <- all.gene.location.filter.nodup.order.fdr.na[order(as.numeric(all.gene.location.filter.nodup.order.fdr.na$chromosome_name)),]
all.gene.location.filter.nodup.order.fdr.na.chr <- all.gene.location.filter.nodup.order.fdr.na.chr[all.gene.location.filter.nodup.order.fdr.na.chr$chromosome_name %in% c(1:23),];
all.gene.location.filter.nodup.order.fdr.na.chr <- all.gene.location.filter.nodup.order.fdr.na.chr[!duplicated(all.gene.location.filter.nodup.order.fdr.na.chr$hgnc_symbol), ]

 # set up chromosome covariate colours to use for chr covariate, below
chr.colours <- force.colour.scheme(all.gene.location.filter.nodup.order.fdr.na.chr$chromosome_name, scheme = 'chromosome');


chr.n.genes <- numeric(23);
chr.tck <- numeric(23);
chr.pos.genes <- numeric(23);
chr.break <- c(0, numeric(23));

for (i in 1:23) {
    n <- sum(all.gene.location.filter.nodup.order.fdr.na.chr$chromosome_name == i);
    chr.n.genes[i] <- n;
    chr.break[i + 1] <- n + chr.break[i];
    chr.pos.genes[i] <- floor(chr.n.genes[i] / 2);
    chr.tck[i] <- chr.pos.genes[i] + which(all.gene.location.filter.nodup.order.fdr.na.chr$chromosome_name == i)[1];
    };

all.gene.location.filter.nodup.order.fdr.na.chr$ind <- seq_len(nrow(all.gene.location.filter.nodup.order.fdr.na.chr));

# Generate Manhattan plot
outlier.manhattan <- create.manhattanplot(
    formula = fdr ~ ind,
    data = all.gene.location.filter.nodup.order.fdr.na.chr,
    main = expression('All genes'),
    main.cex = 1.5,
    xlab.label = expression('Chromosome'),
    ylab.label = expression(FDR),
    xat = chr.tck,
    yaxis.tck = c(0.2, 0),
    xaxis.lab = c(1:22, 'X'),
    xaxis.tck = 0,
    xaxis.cex = 0.9,
    yaxis.cex = 1.1,
    xlab.cex = 1.3,
    ylab.cex = 1.3,
    yat = seq(0, 15, 5),
    yaxis.lab = expression(10^0, 10^-5, 10^-10, 10^-15),
    ylimits = c(-0.1, 18),
    xlimits = c(-max(chr.break) / 50, max(chr.break) * 1.02),
    col = chr.colours,
    pch = 20,
    cex = 0.7,
    add.rectangle = TRUE,
    xleft.rectangle = chr.break[seq(1, length(chr.break) - 1, 2)],
    ybottom.rectangle = -1,
    xright.rectangle = chr.break[c(seq(2, length(chr.break) - 1, 2), 24)],
    ytop.rectangle = 19,
    col.rectangle = 'grey',
    alpha.rectangle = 0.25,
    description = 'Manhattan plot created using BoutrosLab.plotting.general',
    resolution = 200
    );


# Save the plot as a PNG
output.directory <- get0('output.directory', ifnotfound = 'figures');

write.plot(
    trellis.object = outlier.manhattan,
    filename = file.path(output.directory, 'Figure_1_i.png'),
    width = 11,
    height = 4.5,
    size.units = 'in',
    resolution = 1200
);
