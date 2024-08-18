### HISTORY #####################################################################
# This script analyzes RNAi gene dependency scores for outlier genes identified 
# in CCLE. 
# Date: 2024-08-16


# Read RNAi gene dependency file
gene.rnai <- read.delim2(
    file ='/CCLE/RNAi/D2_combined_gene_dep_scores.csv', 
    header =TRUE, 
    row.names =1, 
    sep =','
    );

# Match breast sample information
sample.info.breast.match <- sample.info.breast[
    match(colnames(fpkm.tumor.symbol.filter), rownames(sample.info.breast)),];

# Filter and transpose the RNAi data for breast samples
gene.rnai.breast.t <- gene.rnai[, na.omit(match(sample.info.breast.match$CCLE_Name, colnames(gene.rnai)))];
gene.rnai.breast.num <- apply(gene.rnai.breast.t,1,function(x){as.numeric(x)});
gene.rnai.breast.t.num <- data.frame(t(gene.rnai.breast.num));
rownames(gene.rnai.breast.t.num)<- gsub("\\s*\\(.*\\)$","", rownames(gene.rnai.breast.t));
colnames(gene.rnai.breast.t.num)<- rownames(sample.info.breast.match)[
    match(colnames(gene.rnai.breast.t), sample.info.breast.match$CCLE_Name)];

# Handle duplicate gene names
row_names <- row.names(gene.rnai.breast.t.num);
gene.rnai.breast.t.num.dup.gene <- data.frame();

for(i in 1:nrow(gene.rnai.breast.t.num)){
    row_name_parts <- unlist(strsplit(row_names[i],"&"));
    
    for(j in 1:length(row_name_parts)){
        new_row <- gene.rnai.breast.t.num[i,];
        row.names(new_row)<- row_name_parts[j];
        gene.rnai.breast.t.num.dup.gene <- rbind(gene.rnai.breast.t.num.dup.gene, new_row);
    }}

# Filter for FDR < 0.05 and match gene names
sample.outlier.05.overlap.row <- sample.outlier.05[
    rownames(gene.rank.p.value.one.gene.p0.fdr.05),];
gene.rnai.breast.t.num.match.05 <- gene.rnai.breast.t.num.dup.gene[
    rownames(gene.rnai.breast.t.num.dup.gene)%in% gsub("\\..*$","", rownames(gene.rank.p.value.one.gene.p0.fdr.05)),];
gene.rnai.breast.t.num.match.05.na <- gene.rnai.breast.t.num.match.05;
sample.outlier.05.overlap.na <- sample.outlier.05.overlap[
    match(rownames(gene.rnai.breast.t.num.match.05.na), gsub("\\..*$","", rownames(sample.outlier.05.overlap))),];
sample.outlier.05.overlap.na <- sample.outlier.05.overlap.na[, colnames(gene.rnai.breast.t.num.match.05.na)];
sample.outlier.05.overlap.na.rnai <- sample.outlier.05.overlap.na[, colnames(gene.rnai.breast.t.num.match.05.na)];
rownames(gene.rnai.breast.t.num.match.05.na)<- rownames(sample.outlier.05.overlap.na);

# Sum and filter out genes with no outliers
sample.outlier.05.overlap.na.sum <- apply(sample.outlier.05.overlap.na,1,sum);
sample.outlier.05.overlap.na <- sample.outlier.05.overlap.na[
    sample.outlier.05.overlap.na.sum > 0,];
gene.rnai.breast.t.num.match.05.na <- gene.rnai.breast.t.num.match.05.na[
    sample.outlier.05.overlap.na.sum > 0,];



outlier.gene.rnai.matrix.05 <-NULL;
for(i in 1:nrow(sample.outlier.05.overlap.na)){
    each.list <- rnai.quantile.05[[i]];
    each.list.data <- data.frame(
        gene = rep(rownames(each.list), ncol(each.list)),
        quantile = as.numeric(each.list),
        rnai = as.numeric(outlier.gene.rnai.score.05[[i]]),
        sample = colnames(each.list));
    
    outlier.gene.rnai.matrix.05 <- rbind(outlier.gene.rnai.matrix.05, each.list.data);
    }


outlier.gene.rnai.matrix.05.ess <- outlier.gene.rnai.matrix.05[
    outlier.gene.rnai.matrix.05$quantile == 0& outlier.gene.rnai.matrix.05$rnai <-0.5,];


outlier.gene.rnai.matrix.05.overlap <- outlier.gene.rnai.matrix.05[
    sub("\\..*","", outlier.gene.rnai.matrix.05$gene)%in% five.data.outlier.symbol,];

# Calculate mean RNAi scores and differences
outlier.gene.rnai.score.05.mean <- sapply(outlier.gene.rnai.score.05,function(x){
    mean(na.omit(unlist(x)));  
    });
outlier.gene.rnai.score.05.mean <- data.frame(rnai = outlier.gene.rnai.score.05.mean);
rownames(outlier.gene.rnai.score.05.mean)<- rownames(sample.outlier.05.overlap.na);

nonoutlier.gene.rnai.score.05.mean <- sapply(nonoutlier.gene.rnai.score.05,function(x){
    mean(na.omit(unlist(x)));  
    });
nonoutlier.gene.rnai.score.05.mean <- data.frame(rnai = nonoutlier.gene.rnai.score.05.mean);
rownames(nonoutlier.gene.rnai.score.05.mean)<- rownames(sample.outlier.05.overlap.na);

# Perform Wilcoxon test on RNAi scores
p.rnai.05 <- wilcox.test(
    outlier.gene.rnai.score.05.mean$rnai, 
    nonoutlier.gene.rnai.score.05.mean$rnai, 
    alternative ="two.sided", conf.int =TRUE);

# Create difference matrix for RNAi scores
gene.rnai.diff.matrix.05 <- data.frame(
    out = outlier.gene.rnai.score.05.mean$rnai,
    non = nonoutlier.gene.rnai.score.05.mean$rnai,
    diff = outlier.gene.rnai.score.05.mean$rnai - nonoutlier.gene.rnai.score.05.mean$rnai,
    symbol = sub("\\..*","", rownames(outlier.gene.rnai.score.05.mean)));
rownames(gene.rnai.diff.matrix.05)<- rownames(outlier.gene.rnai.score.05.mean);

# Set colors for the scatter plot
dot.colours <-rep('grey30', nrow(gene.rnai.diff.matrix.05));
dot.colours[gene.rnai.diff.matrix.05$diff <-0.5]<-'dodgerblue3';
dot.colours[gene.rnai.diff.matrix.05$diff >0.5]<-'red2';


interesting.points <- gene.rnai.diff.matrix.05.overlap$diff <-0.5;
text.x <- na.omit(gene.rnai.diff.matrix.05.overlap$non[interesting.points]);
text.y <- na.omit(gene.rnai.diff.matrix.05.overlap$out[interesting.points]);
text.labels <- na.omit(gene.rnai.diff.matrix.05.overlap$symbol[interesting.points]);

# Create the scatter plot
gene.scatter.05.minus.overlap.label <- create.scatterplot(
    formula = diff ~ non,
    data = gene.rnai.diff.matrix.05.overlap,
    col = dot.colours,
    alpha =.6,
    ylimits =c(-1.3,1.05),
    xlimits =c(-1.4,0.5),
    xaxis.fontface =1,
    yaxis.fontface =1,
    yaxis.tck =c(0.2,0),
    xaxis.tck =c(0.2,0),
    add.grid =FALSE,
    cex =1,
    xaxis.cex =1,
    yaxis.cex =1,
    xlab.cex =1.2,
    ylab.cex =1.2,
    main.cex =1.5,
    left.padding =0,
    main =expression('Gene effect score of the overlap outliers'),
    xlab.label =expression(paste('Mean of gene rnai score of non-outlier samples')),
    ylab.label =expression('Difference of gene effect score'),
    add.text =TRUE,
    text.x = text.x,
    text.y = text.y,
    text.labels = text.labels,
    text.fontface =1,
    abline.h =c(0,-0.5),
    abline.col =c('grey30','dodgerblue3'),
    abline.lwd =c(1.5,2),
    abline.lty =c(1,3)
    );



# Save the scatter plot as a PDF
pdf(
    file = generate.filename(
        'gene_dependency_diff_rnai', 
        'scatter', 
        'pdf'
        ), 
    width = 6, 
    height = 5
    );
gene.scatter.05.minus.overlap.label;
dev.off();

# Save the scatter plot as a PNG
png(
    file = generate.filename(
        'gene_dependency_diff_rnai', 
        'scatter', 
        'png'
        ), 
    width = 6, 
    height = 5,
    unit = 'in', 
    res = 1200
    );
gene.scatter.05.minus.overlap.label;
dev.off();




