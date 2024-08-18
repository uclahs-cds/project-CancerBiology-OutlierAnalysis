### HISTORY #####################################################################
# This script generates a Manhattan plot of all genes from five datasets, 
# highlighting the distribution of outlier genes across chromosomes.
# Date: 2024-08-16


# Chromosome name vector
chr.name <- c('1', '2', '3', '4', '5', '6', '7', '8', '9', '10', '11', '12', '13', '14', '15', '16', '17', '18', '19', '20', '21', '22', 'MT', 'X', 'Y');

# Function to process gene position data
process.gene.position <- function(gene.position, gene.rank, chromosome.col = "chromosome_name", start.col = "start_position") {
    gene.position.order <- gene.position[match(gene.rank$gene_id, gene.position$gene_id), ];
    gene.position.order[[chromosome.col]][gene.position.order[[chromosome.col]] == 'X'] <- 24;
    gene.position.order[[chromosome.col]][gene.position.order[[chromosome.col]] == 'Y'] <- NA;
    gene.position.order[[chromosome.col]][gene.position.order[[chromosome.col]] == 'MT'] <- 23;

    gene.position.order.fdr <- data.frame(
        gene.position.order,
        fdr = gene.rank$fdr
        );

    gene.position.order.fdr.na <- na.omit(gene.position.order.fdr);
    gene.position.order.fdr.na <- gene.position.order.fdr.na[order(as.numeric(gene.position.order.fdr.na[[start.col]])), ];
    gene.position.order.fdr.na.chr <- gene.position.order.fdr.na[order(as.numeric(gene.position.order.fdr.na[[chromosome.col]])), ];
    gene.position.order.fdr.na.chr <- gene.position.order.fdr.na.chr[gene.position.order.fdr.na.chr[[chromosome.col]] %in% c(1:24), ];

    return(gene.position.order.fdr.na.chr);
    };

# Process data for each dataset
gene.position.brca.all.order.fdr.na.chr <- process.gene.position(gene.position.brca.all, gene.rank.order.cosine.observed.p.value.max.filter.1.brca);
gene.position.meta.all.order.fdr.na.chr <- process.gene.position(gene.position.meta.all, gene.rank.order.cosine.observed.p.value.meta);
gene.position.ispy.all.order.fdr.na.chr <- process.gene.position(gene.position.ispy.all, gene.rank.order.cosine.observed.p.value.ispy);
gene.position.metador.all.order.fdr.na.chr <- process.gene.position(gene.position.metador.all, gene.rank.order.cosine.observed.p.value.metador);

# Process ICGC data
locations <- fpkm.data.icgc$loc[as.numeric(rownames(gene.rank.order.cosine.observed.p.value.icgc))];

split.location <- function(location) {
    parts <- strsplit(location, ":|\\-")[[1]];
    return(c(parts[1], parts[2], parts[3]));
    };

split.locations <- t(sapply(locations, split.location));
location.df <- data.frame(split.locations, stringsAsFactors = FALSE);
names(location.df) <- c("Chromosome", "Start", "End");

gene.position.icgc.all.order.fdr <- data.frame(
    symbol = fpkm.data.icgc$Name[as.numeric(rownames(gene.rank.order.cosine.observed.p.value.icgc))],
    location.df,
    fdr = gene.rank.order.cosine.observed.p.value.icgc$fdr
    );

gene.position.icgc.all.order.fdr.na.chr <- process.gene.position(gene.position.icgc.all.order.fdr, gene.rank.order.cosine.observed.p.value.icgc, "Chromosome", "Start");

# Combine all datasets
all.gene.location <- rbind(
    gene.position.meta.all[, 2:5],
    gene.position.brca.all[, 2:5],
    gene.position.ispy.all[, 2:5],
    gene.position.metador.all[, 2:5],
    gene.position.icgc.all[, 2:5]
    );

all.gene.location.filter <- all.gene.location[!grepl("^CHR", all.gene.location$chromosome_name), ];
all.gene.location.filter.nodup <- all.gene.location.filter[!duplicated(all.gene.location.filter$hgnc_symbol), ];

all.gene.location.filter.nodup.order <- all.gene.location.filter.nodup[match(names(combine.fisher.pvalue.all.fdr.sort.log), all.gene.location.filter.nodup$hgnc_symbol), ];
all.gene.location.filter.nodup.order$chromosome_name[all.gene.location.filter.nodup.order$chromosome_name == 'X'] <- 23;
all.gene.location.filter.nodup.order$chromosome_name[all.gene.location.filter.nodup.order$chromosome_name %in% c('Y', 'MT')] <- NA;

all.gene.location.filter.nodup.order.fdr <- data.frame(
    all.gene.location.filter.nodup.order, 
    fdr = combine.fisher.pvalue.all.fdr.sort.log
    );

all.gene.location.filter.nodup.order.fdr.na.chr <- process.gene.position(all.gene.location.filter.nodup.order.fdr, data.frame(fdr = combine.fisher.pvalue.all.fdr.sort.log));
all.gene.location.filter.nodup.order.fdr.na.chr <- all.gene.location.filter.nodup.order.fdr.na.chr[!duplicated(all.gene.location.filter.nodup.order.fdr.na.chr$hgnc_symbol), ];

# Set up chromosome colors and positions
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


# Save the manhattan plot as a PDF
pdf(
    file = generate.filename(
        'combine_outlier', 
        'manhattan', 
        'pdf'
        ), 
    width = 11, 
    height = 4.5
    );
outlier.manhattan;
dev.off();

# Save the manhattan plot as a PNG
png(
    file = generate.filename(
        'combine_outlier', 
        'manhattan', 
        'png'
        ), 
    width = 11, 
    height = 4.5,
    unit = 'in', 
    res = 1200
    );
outlier.manhattan;
dev.off();



