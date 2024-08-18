### HISTORY #####################################################################
# This script processes and analyzes drug response data from Sanger.
# Date: 2024-08-16


# Load required libraries
library(tidyr);
library(reshape2);
library(dplyr);

# Read data files
sanger.drug <- read.delim2(
    file ='/CCLE/sanger/sanger-dose-response.csv', 
    header =TRUE, 
    sep =','
    );
sanger.drug.vi <- read.delim2(
    file ='/CCLE/sanger/sanger-viability.csv', 
    header =TRUE, 
    sep =','
    );
sanger.drug.info <- read.delim2(
    file ='/CCLE/sanger/screened_compounds_rel_8.5.csv', 
    header =TRUE, 
    sep =','
    );

# Data matching 
sanger.drug.match <- sanger.drug[sanger.drug$ARXSPAN_ID %in% colnames(sample.outlier.05.overlap.na),];
sanger.drug.match.name.unique <- unique(sanger.drug.match$DRUG_NAME);

# Separate and clean data
sanger.drug.match.dup <- separate_rows(sanger.drug.match, DRUG_NAME, sep =", ");
sanger.drug.match.dup <- data.frame(sanger.drug.match.dup);
sanger.drug.match.dup$IC50_PUBLISHED <- as.numeric(sanger.drug.match.dup$IC50_PUBLISHED);

sanger.drug.info.dup <- separate_rows(sanger.drug.info, TARGET, sep =", ");
sanger.drug.info.dup <- data.frame(sanger.drug.info.dup);
sanger.drug.info.dup <- separate_rows(sanger.drug.info.dup, TARGET, sep =",");
sanger.drug.info.dup <- data.frame(sanger.drug.info.dup);

# Filter CCLE and tissue-specific data
sanger.drug.info.dup.ccle <- sanger.drug.info.name.dup[sanger.drug.info.name.dup$merged.target.sanger.depmap %in% fdr.05.symbol,];
sanger.drug.info.dup.tissue <- sanger.drug.info.name.dup[sanger.drug.info.name.dup$merged.target.sanger.depmap %in% five.data.outlier.symbol,];

# Get more drug tartet data from depmap data
depmap.drug.info.match.sanger <- depmap.drug.info[depmap.drug.info$Drug.Name %in% unique(sanger.drug.match.dup$DRUG_NAME),];
depmap.drug.info.match.sanger.dup <- separate_rows(depmap.drug.info.match.sanger, repurposing_target, sep =", ");
depmap.drug.info.match.sanger.dup <- data.frame(depmap.drug.info.match.sanger.dup);

# Define functions to check capitalization
is_first_letter_capital <- function(word){return(grepl("^[A-Z]", word));
    };

contains_lowercase <- function(word){return(grepl("[a-z]", word));
    };

capital_status <- sapply(sanger.drug.info.name.dup$merged.target.sanger.depmap, contains_lowercase);

# Data matching and merging
depmap.drug.info$UpperDrugName <- toupper(depmap.drug.info$Drug.Name);
sanger_upper_drug_names <- toupper(sanger.drug.info.name$DRUG_NAME);
sanger_upper_drug_ids <- toupper(sanger.drug.info.name$drug.name.id);

# Match Sanger and DepMap drug targets
overlap.sanger.depmap.drug.target <-NULL;
for(i in 1:nrow(sanger.drug.info.name)){
    current_targets <-c(sanger_upper_drug_names[i], sanger_upper_drug_ids[i]);
    matched_indices <- which(depmap.drug.info$UpperDrugName %in% current_targets);
    overlap.sanger.depmap.drug.target <-c(
        overlap.sanger.depmap.drug.target,
        depmap.drug.info$repurposing_target[matched_indices][1]
        );
    };

sanger.drug.info.name$depmap.drug.target <- overlap.sanger.depmap.drug.target;
merged.target.sanger.depmap <- paste(sanger.drug.info.name$TARGET, sanger.drug.info.name$depmap.drug.target, sep =', ');
sanger.drug.info.name$merged.target.sanger.depmap <- merged.target.sanger.depmap;

# Separate and restructure data
sanger.drug.info.name.dup <- separate_rows(sanger.drug.info.name, merged.target.sanger.depmap, sep =", ");
sanger.drug.info.name.dup <- data.frame(sanger.drug.info.name.dup);

# Process IC50 and Z-score data
sanger.drug.match.dup.ic50 <- dcast(sanger.drug.match.dup, DRUG_NAME ~ ARXSPAN_ID, value.var ="IC50_PUBLISHED", fun.aggregate = mean);
rownames(sanger.drug.match.dup.ic50)<- sanger.drug.match.dup.ic50$DRUG_NAME;
sanger.drug.match.dup.ic50 <- sanger.drug.match.dup.ic50[,2:ncol(sanger.drug.match.dup.ic50)];

sanger.drug.match.dup$Z_SCORE_PUBLISHED <-as.numeric(sanger.drug.match.dup$Z_SCORE_PUBLISHED);
sanger.drug.match.dup.zscore <- dcast(sanger.drug.match.dup, DRUG_NAME ~ ARXSPAN_ID, value.var ="Z_SCORE_PUBLISHED", fun.aggregate = mean);
rownames(sanger.drug.match.dup.zscore)<- sanger.drug.match.dup.zscore$DRUG_NAME;
sanger.drug.match.dup.zscore <- sanger.drug.match.dup.zscore[,2:ncol(sanger.drug.match.dup.zscore)];

# Process sample outlier data
sample.outlier.05.overlap.na.samger.match <- sample.outlier.05.overlap.na[, colnames(sanger.drug.match.dup.ic50)];
rownames(sample.outlier.05.overlap.na.samger.match)<- unlist(sapply(rownames(sample.outlier.05.overlap.na.samger.match),function(x) sub("\\..*","", x)));
sample.outlier.05.overlap.na.samger.match.dup.sum <- apply(sample.outlier.05.overlap.na.samger.match,1,sum);
sample.outlier.05.overlap.na.samger.match.dup.filter <- sample.outlier.05.overlap.na.samger.match[sample.outlier.05.overlap.na.samger.match.dup.sum >0,];

# Process breast cancer related data
sanger.drug.breast.match.out <-list();
sanger.drug.breast.match.non <-list();
sanger.drug.breast.match.info <-list();

for(i in1:nrow(sample.outlier.05.overlap.na.samger.match.dup.filter)){
    patient.status <- sample.outlier.05.overlap.na.samger.match.dup.filter[i,];
    drug.target <- depmap.drug.info.match.sanger.dup[depmap.drug.info.match.sanger.dup$repurposing_target %in% rownames(patient.status),];
    out.value <- sanger.drug.match.dup.ic50[match(drug.target$Drug.Name, rownames(sanger.drug.match.dup.ic50)), colnames(patient.status)[patient.status ==1], drop =FALSE];
    non.value <- sanger.drug.match.dup.ic50[match(drug.target$Drug.Name, rownames(sanger.drug.match.dup.ic50)), colnames(patient.status)[patient.status ==0], drop =FALSE];
    
    sanger.drug.breast.match.out[[i]]<- out.value;
    sanger.drug.breast.match.non[[i]]<- non.value;
    sanger.drug.breast.match.info[[i]]<- drug.target;
    };

# Combine and calculate mean values for matched data
sanger.drug.breast.match.out.df <- bind_rows(sanger.drug.breast.match.out);
sanger.drug.breast.match.non.df <- bind_rows(sanger.drug.breast.match.non);
sanger.drug.breast.match.info.df <- bind_rows(sanger.drug.breast.match.info);

sanger.drug.breast.match.out.list.mean <- lapply(sanger.drug.breast.match.out,function(x){mean(na.omit(as.numeric(unlist(x))))});
sanger.drug.breast.match.non.list.mean <- lapply(sanger.drug.breast.match.non,function(x){
    apply(x, 2, function(y) {mean(na.omit(as.numeric(y)))});
    });

# Create final result dataframe
sanger.drug.breast.match.out.list.mean.unlist <- unlist(sanger.drug.breast.match.out.list.mean);
sanger.drug.breast.match.non.list.mean.df <- bind_rows(sanger.drug.breast.match.non.list.mean);
sanger.drug.breast.match.out.non.mean.df <- data.frame(out = sanger.drug.breast.match.out.list.mean.unlist,
                                                       sanger.drug.breast.match.non.list.mean.df);
rownames(sanger.drug.breast.match.out.non.mean.df) <- rownames(sample.outlier.05.overlap.na.samger.match.dup.filter);

# Final summary of mean values
sanger.drug.breast.match.out.mean <- apply(sanger.drug.breast.match.out.df, 1, function(x){mean(na.omit(as.numeric(x)))}); 
sanger.drug.breast.match.non.mean <- apply(sanger.drug.breast.match.non.df, 1, function(x){mean(na.omit(as.numeric(x)))}); 
sanger.drug.breast.match.mean.df <- data.frame(
    out = unlist(sanger.drug.breast.match.out.mean),
    non = unlist(sanger.drug.breast.match.non.mean),
    minus = unlist(sanger.drug.breast.match.out.mean)- unlist(sanger.drug.breast.match.non.mean)); 
sanger.drug.breast.match.mean.df <- cbind(sanger.drug.breast.match.mean.df, sanger.drug.breast.match.info.df); 

# Group and summarize by repurposing target
sanger.drug.breast.match.mean.df.merge <- sanger.drug.breast.match.mean.df %>%
    group_by(repurposing_target)%>%
    summarise(
        out_avg = mean(out),
        non_avg = mean(non),
        minus_avg = mean(minus),
        screen_first = first(screen),
        dose_first = first(dose),
        Drug_Name_first = first(Drug.Name),
        MOA_first = first(MOA),
        IDs_first = first(IDs),
        Synonyms_first = first(Synonyms),
        .groups ='drop');

# Z-score analysis for matched breast cancer data
sanger.zscore.drug.breast.match.out <-list(); 
sanger.zscore.drug.breast.match.non <-list(); 
sanger.zscore.drug.breast.match.info <-list(); 

for(i in1:nrow(sample.outlier.05.overlap.na.samger.match.dup.filter)){ 
    patient.status <- sample.outlier.05.overlap.na.samger.match.dup.filter[i,]; 
    drug.target <- depmap.drug.info.match.sanger.dup[depmap.drug.info.match.sanger.dup$repurposing_target %in% rownames(patient.status),]; 
    out.value <- sanger.drug.match.dup.zscore[match(drug.target$Drug.Name, rownames(sanger.drug.match.dup.zscore)), colnames(patient.status)[patient.status ==1], drop =FALSE]; 
    non.value <- sanger.drug.match.dup.zscore[match(drug.target$Drug.Name, rownames(sanger.drug.match.dup.zscore)), colnames(patient.status)[patient.status ==0], drop =FALSE]; 
    
    sanger.zscore.drug.breast.match.out[[i]]<- out.value; 
    sanger.zscore.drug.breast.match.non[[i]]<- non.value; 
    sanger.zscore.drug.breast.match.info[[i]]<- drug.target; 
    }; 

# Combine and summarize Z-score data
sanger.zscore.drug.breast.match.out.df <- bind_rows(sanger.zscore.drug.breast.match.out); 
sanger.zscore.drug.breast.match.non.df <- bind_rows(sanger.zscore.drug.breast.match.non); 
sanger.zscore.drug.breast.match.info.df <- bind_rows(sanger.zscore.drug.breast.match.info); 

sanger.zscore.drug.breast.match.out.list.mean <- lapply(sanger.zscore.drug.breast.match.out,function(x){mean(na.omit(as.numeric(unlist(x))))}); 
sanger.zscore.drug.breast.match.non.list.mean <- lapply(sanger.zscore.drug.breast.match.non,function(x){
    apply(x,2,function(y){mean(na.omit(as.numeric(y)))})});


# Unlist and bind Z-score data
sanger.zscore.drug.breast.match.out.list.mean.unlist <- unlist(sanger.zscore.drug.breast.match.out.list.mean);
sanger.zscore.drug.breast.match.non.list.mean.df <- bind_rows(sanger.zscore.drug.breast.match.non.list.mean);
sanger.zscore.drug.breast.match.out.non.mean.df <- data.frame(
    out = sanger.zscore.drug.breast.match.out.list.mean.unlist,
    sanger.zscore.drug.breast.match.non.list.mean.df
    );
rownames(sanger.zscore.drug.breast.match.out.non.mean.df)<- rownames(sample.outlier.05.overlap.na.samger.match.dup.filter);

# Calculate mean Z-scores
sanger.zscore.drug.breast.match.out.mean <- apply(sanger.zscore.drug.breast.match.out.df,1,function(x){mean(na.omit(as.numeric(x)))});
sanger.zscore.drug.breast.match.non.mean <- apply(sanger.zscore.drug.breast.match.non.df,1,function(x){mean(na.omit(as.numeric(x)))});

# Combine Z-score mean data
sanger.zscore.drug.breast.match.mean.df <- data.frame(
    out = unlist(sanger.zscore.drug.breast.match.out.mean),
    non = unlist(sanger.zscore.drug.breast.match.non.mean),
    minus = unlist(sanger.zscore.drug.breast.match.out.mean)- unlist(sanger.zscore.drug.breast.match.non.mean));
sanger.zscore.drug.breast.match.mean.df <- cbind(sanger.zscore.drug.breast.match.mean.df, sanger.zscore.drug.breast.match.info.df);

# Group by repurposing target and summarize
sanger.zscore.drug.breast.match.mean.df.merge <- sanger.zscore.drug.breast.match.mean.df %>%
    group_by(repurposing_target)%>%
    summarise(
        out_avg = mean(out),
        non_avg = mean(non),
        minus_avg = mean(minus),
        screen_first = first(screen),
        dose_first = first(dose),
        Drug_Name_first = first(Drug.Name),
        MOA_first = first(MOA),
        IDs_first = first(IDs),
        Synonyms_first = first(Synonyms),
        .groups ='drop');
sanger.zscore.drug.breast.match.mean.df.merge <- data.frame(sanger.zscore.drug.breast.match.mean.df.merge);

# Do not merge genes, keep each Z-score separately
sanger.zscore.drug.breast.match.out.non.each.df <- data.frame(
    out = apply(sanger.zscore.drug.breast.match.out.df, 1, function(x){na.omit(as.numeric(x))}),
    non = sanger.zscore.drug.breast.match.non.df
    );
rownames(sanger.zscore.drug.breast.match.out.non.each.df)<- rownames(sanger.zscore.drug.breast.match.out.df);

# Prepare data for boxplot
z.score.box.each <- data.frame(
    value = c(sanger.zscore.drug.breast.match.out.non.each.df$out, 
              unlist(sanger.zscore.drug.breast.match.out.non.each.df[,2:ncol(sanger.zscore.drug.breast.match.out.non.each.df)])),
    status = c(rep('a', nrow(sanger.zscore.drug.breast.match.out.non.each.df)),
               rep('b',((ncol(sanger.zscore.drug.breast.match.out.non.each.df)-1)*nrow(sanger.zscore.drug.breast.match.out.non.each.df))))
    );

# Perform Wilcoxon test on Z-scores
p.zscore <- wilcox.test(sanger.zscore.drug.breast.match.out.non.each.df$out, unlist(sanger.zscore.drug.breast.match.out.non.each.df[,2:ncol(sanger.zscore.drug.breast.match.out.non.each.df)]), alternative ="two.sided", conf.int =TRUE);

text.pvalue <- display.statistical.result(
    x = p.zscore$p.value,
    statistic.type ='p',
    symbol =' = ');

# Create legend key for the plot
key <-list(
    text =list(
        lab = text.pvalue, 
        cex =1), 
    x =0.25,
    y =0.95);

# Create Z-score boxplot
zscore.box.sanger <- BoutrosLab.plotting.general::create.boxplot(
    formula = value ~ status,
    data = z.score.box.each,
    main =expression('Z-score'),
    outlier =TRUE,
    add.stripplot =TRUE,
    add.rectangle =TRUE,
    xleft.rectangle =1.5,
    xright.rectangle =3,
    ybottom.rectangle =-3,
    ytop.rectangle =5,
    col.rectangle ="grey",
    alpha.rectangle =0.25,
    main.cex =1.5,
    xlab.label =NULL,
    xlab.cex =0,
    ylab.label =expression('z-score'),
    ylab.cex =1.3,
    yaxis.cex =1.1,
    xaxis.cex =1.1,
    xaxis.lab =c('Outlier','Non-outlier'),
    xaxis.fontface =1,
    yaxis.fontface =1,
    yaxis.tck =c(0.2,0),
    xaxis.tck =c(0.2,0),
    xaxis.rot =90,
    ylimits =c(-3,4.5),
    yat = seq(-2,4,1),
    sample.order ="none",
    legend =list(
        inside =list(
            fun = draw.key,
            args =list(
                key = key
                ),
            x =0.5,
            y =0.97,
            corner =c(0,1))),
    points.pch =16,
    points.cex =1,
    points.col ='grey70',
    lwd =1.2,
    col =c('red3','dodgerblue3'),
    alpha =0.25
    );




# Save the box plot as a PDF
pdf(
    file = generate.filename(
        'sanger_drug', 
        'box', 
        'pdf'
        ), 
    width = 4, 
    height = 6
    );
zscore.box.sanger;
dev.off();

# Save the box plot as a PNG
png(
    file = generate.filename(
        'sanger_drug', 
        'box', 
        'png'
        ), 
    width = 4, 
    height = 6,
    unit = 'in', 
    res = 1200
    );
zscore.box.sanger;
dev.off();

