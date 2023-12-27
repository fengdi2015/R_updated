if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("CNTools")


BiocManager::install("survminer")
BiocManager::install("RTCGAToolbox")
### getTCGA !!!
getTCGA <- function(disease="GBM", data.type="RNASeq2", type="", filter="Y", p=getOption("mc.cores", 2), clinical=FALSE, cvars="OS"){
  
  # Check if it's valid datatype
  data.good <-  c("RNASeq2", "RNASeq", "miRNASeq", "CNA_SNP", "CNV_SNP", "CNA_CGH", 
                  "Methylation", "Mutation", "mRNA_Array", "miRNA_Array")
  if(! (data.type %in% data.good)){
    message("Error: Not recognized datatype for Firehose\n")
    dat <- NULL
    return(dat)
  }
  
  # Parse links
  ldoc <- tryCatch({
    fileURL="https://gdac.broadinstitute.org/runs/stddata__latest/"
    ldoc <- htmlParse(rawToChar(GET(fileURL)$content))
    
    #ldoc <- XML::htmlTreeParse(sub("s", "", fileURL), useInternalNodes = F)
  } ,error = function(e){
    ldoc = NULL
  })
  
  if(is.null(ldoc)){
    message("Error: Problem connect to Firehose. Please ensure Internet connection is working.\n")
    dat <- NULL
    return(dat)
  }
  
  # parse to get data set and URL for each datasets
  datasets <- XML::xpathSApply(ldoc, "//a[contains(@href, 'Standardized+Data+Run+Release+Notes')]", XML::xmlValue)
  
  if(!(disease %in% datasets)){
    message("Error: Not recognized Disease Abbreviation for Firehose\n")
    dat <- NULL
    return(dat)
  }
  
  if(Sys.getenv("TAR") == "" | Sys.getenv("R_GZIPCMD") ==""){
    message("Error: TAR is not installed in the system. Data unzip failed.\n")
    dat <- NULL
    return(dat)
  }
  
  dataset <- disease
  llinks = unlist(XML::xpathApply(ldoc, "//a[@href]", XML::xmlGetAttr, "href"))
  dlinks = llinks[grepl(paste("/data/", dataset, "/", sep=""), llinks)]
  #ddoc = XML::htmlTreeParse(dlinks, useInternalNodes = T)
  ddoc= htmlParse(rawToChar(GET(dlinks)$content))
  # -----------------
  # Get data
  if(data.type=="RNASeq2"){
    dat <- RNASeqV2(ddoc=ddoc, dlinks=dlinks, dataset=dataset)
    
    if(is.null(dat)){
      return(dat)
    }
    if(!clinical){
      return(list(dat=dat, clinical=NULL, merged.dat=NULL))
    }
    
    if(clinical){
      gdats <- list()
      cli <- Clinical(ddoc=ddoc, dlinks=dlinks, dataset=dataset)
      mdat <- NULL
      if(!is.null(cli)){
        # combine clinical and data together
        mdat <- MatrixMerge(dat=dat, cli=cli, cvars=cvars)
      }
      return(list(dat=dat, clinical=cli, merged.dat=mdat))
    }
  }
  
  # ----------
  if(data.type=="RNASeq"){
    if(type==""){
      dat <- RNASeq(ddoc=ddoc, dlinks=dlinks, dataset=dataset)
    }
    else{
      dat <- RNASeq(ddoc=ddoc, dlinks=dlinks, dataset=dataset, type=type)
    }
    
    if(is.null(dat)){
      return(dat)
    }
    if(!clinical){
      return(list(dat=dat, clinical=NULL, merged.dat=NULL))
    }
    
    if(clinical){
      gdats <- list()
      cli <- Clinical(ddoc=ddoc, dlinks=dlinks, dataset=dataset)
      mdat <- NULL
      if(!is.null(cli)){
        # combine clinical and data together
        mdat <- MatrixMerge(dat=dat, cli=cli, cvars=cvars)
      }
      return(list(dat=dat, clinical=cli, merged.dat=mdat))
    }
  }
  # ----------
  if(data.type=="miRNASeq"){
    if(type==""){
      dat <- miRNASeq(ddoc=ddoc, dlinks=dlinks, dataset=dataset)
    }
    else{
      dat <- miRNASeq(ddoc=ddoc, dlinks=dlinks, dataset=dataset, type=type)
    }
    if(is.null(dat)){
      return(dat)
    }
    if(!clinical){
      return(list(dat=dat, clinical=NULL, merged.dat=NULL))
    }
    
    if(clinical){
      gdats <- list()
      cli <- Clinical(ddoc=ddoc, dlinks=dlinks, dataset=dataset)
      mdat <- NULL
      if(!is.null(cli)){
        # combine clinical and data together
        mdat <- MatrixMerge(dat=dat, cli=cli, cvars=cvars)
      }
      return(list(dat=dat, clinical=cli, merged.dat=mdat))
    }
  }
  # ----------
  if(data.type=="CNA_SNP"){
    dat <- CNA_SNP(ddoc=ddoc, dlinks=dlinks, dataset=dataset, filter=filter)
    
    if(is.null(dat)){
      return(dat)
    }
    if(!clinical){
      return(list(dat=dat, clinical=NULL, merged.dat=NULL))
    }
    
    if(clinical){
      gdats <- list()
      cli <- Clinical(ddoc=ddoc, dlinks=dlinks, dataset=dataset)
      mdat = NULL
      if(!is.null(cli)){
        # combine clinical and data together
        mdat <- MatrixMerge(dat=dat, cli=cli, cvars=cvars)
      }
      return(list(dat=dat, clinical=cli, merged.dat=mdat))
    }
  }
  # ----------
  if(data.type=="CNV_SNP"){
    dat <- CNV_SNP(ddoc=ddoc, dlinks=dlinks, dataset=dataset, filter=filter)
    
    if(is.null(dat)){
      return(dat)
    }
    if(!clinical){
      return(list(dat=dat, clinical=NULL, merged.dat=NULL))
    }
    
    if(clinical){
      gdats <- list()
      cli <- Clinical(ddoc=ddoc, dlinks=dlinks, dataset=dataset)
      mdat = NULL
      if(!is.null(cli)){
        # combine clinical and data together
        mdat <- MatrixMerge(dat=dat, cli=cli, cvars=cvars)
      }
      return(list(dat=dat, clinical=cli, merged.dat=mdat))
    }
  }
  # ----------
  if(data.type=="CNA_Seq"){
    dat <- CNA_Seq(ddoc=ddoc, dlinks=dlinks, dataset=dataset)
    
    if(is.null(dat)){
      return(dat)
    }
    if(!clinical){
      return(dat)
    }
    
    if(clinical){
      gdats <- list()
      cli <- Clinical(ddoc=ddoc, dlinks=dlinks, dataset=dataset)
      
      if(!is.null(cli)){
        # combine clinical and data together
        mdat <- MatrixMerge(dat=dat, cli=cli, cvars=cvars)
      }
      return(list(dat=dat, clinical=cli, merged.dat=mdat))
    }
  }
  # ----------
  if(data.type=="CNA_CGH"){
    if(type==""){
      dat <- CNA_CGH(ddoc=ddoc, dlinks=dlinks, dataset=dataset, filter=filter)
    }
    else{
      dat <- CNA_CGH(ddoc=ddoc, dlinks=dlinks, dataset=dataset, type=type, filter=filter)
    }
    
    if(is.null(dat)){
      return(dat)
    }
    if(!clinical){
      return(list(dat=dat, clinical=NULL, merged.dat=NULL))
    }
    
    if(clinical){
      gdats <- list()
      cli <- Clinical(ddoc=ddoc, dlinks=dlinks, dataset=dataset)
      mdat = NULL
      if(!is.null(cli)){
        # combine clinical and data together
        mdat <- MatrixMerge(dat=dat, cli=cli, cvars=cvars)
      }
      return(list(dat=dat, clinical=cli, merged.dat=mdat))
    }
  }
  # ----------
  if(data.type=="Methylation"){
    
    if(type==""){
      dat <- Methylation(ddoc=ddoc, dlinks=dlinks, dataset=dataset, p=p)
    }
    else{
      dat <- Methylation(ddoc=ddoc, dlinks=dlinks, dataset=dataset, p=p, type=type)
    }
    
    if(is.null(dat)){
      return(dat)
    }
    if(!clinical){
      return(list(dat=as.matrix(dat[, -c(1,2,3)]), cpgs=dat[,1:3], clinical=NULL, merged.dat=NULL))
    }
    
    if(clinical){
      gdats <- list()
      cli <- Clinical(ddoc=ddoc, dlinks=dlinks, dataset=dataset)
      mdat <- NULL
      if(!is.null(cli)){
        # combine clinical and data together
        mdat <- MatrixMerge(dat=dat, cli=cli, cvars=cvars)
      }
      return(list(dat=as.matrix(dat[, -c(1,2,3)]), cpgs=dat[,1:3], clinical=cli, merged.dat= mdat))
    }
  }
  # ----------
  if(data.type=="Mutation"){
    if(type==""){
      dat <- Mutation(ddoc=ddoc, dlinks=dlinks, dataset=dataset)
    }
    else{
      dat <- Mutation(ddoc=ddoc, dlinks=dlinks, dataset=dataset, type=type)
    }
    if(is.null(dat)){
      return(dat)
    }
    if(!clinical){
      return(list(dat=dat, clinical=NULL, merged.dat=NULL))
    }
    
    if(clinical){
      gdats <- list()
      cli <- Clinical(ddoc=ddoc, dlinks=dlinks, dataset=dataset)
      mdat <- NULL
      
      if(!is.null(cli)){
        # combine clinical and data together
        mdat <- MatrixMerge(dat=dat, cli=cli, cvars=cvars)
      }
      return(list(dat=dat, clinical=cli, merged.dat=mdat))
    }
  }
  # ----------
  if(data.type=="mRNA_Array"){
    if(type==""){
      dat <- mRNA_Array(ddoc=ddoc, dlinks=dlinks, dataset=dataset)
    }
    else{
      dat <- mRNA_Array(ddoc=ddoc, dlinks=dlinks, dataset=dataset, type=type)
    }
    
    if(is.null(dat)){
      return(dat)
    }
    if(!clinical){
      return(list(dat=dat, clinical=NULL, merged.dat=NULL))
    }
    
    if(clinical){
      gdats <- list()
      cli <- Clinical(ddoc=ddoc, dlinks=dlinks, dataset=dataset)
      mdat <- NULL
      if(!is.null(cli)){
        # combine clinical and data together
        mdat <- MatrixMerge(dat=dat, cli=cli, cvars=cvars)
      }
      return(list(dat=dat, clinical=cli, merged.dat=mdat))
    }
  }
  # ----------
  if(data.type=="miRNA_Array"){
    dat <- miRNA_Array(ddoc=ddoc, dlinks=dlinks, dataset=dataset)
    
    if(is.null(dat)){
      return(dat)
    }
    if(!clinical){
      return(list(dat=as.matrix(dat), clinical=NULL, merged.dat=NULL))
    }
    
    if(clinical){
      gdats <- list()
      cli <- Clinical(ddoc=ddoc, dlinks=dlinks, dataset=dataset)
      mdat <- NULL
      if(is.null(cli)){
        message("Error: clinical data not available.")
      }
      if(!is.null(cli)){
        # combine clinical and data together
        mdat <- MatrixMerge(dat=dat, cli=cli, cvars=cvars)
      }
      return(list(dat=as.matrix(dat), clinical=cli, merged.dat=as.matrix(mdat)))
    }
  }
}
#
#

# BiocManager::install("CNTools") # Needed for TCGA2STAT
# remotes::install_github("cran/TCGA2STAT")
library(TCGA2STAT)
library(ggplot2)
library(cowplot)
library(rmarkdown)
library(reshape2)
if (!require(survplot)) {
  install.packages("misc/survplot_0.0.7.tar.gz", repos = NULL, type="source")
}
source("Supplemental_R_script_1.R")
load_data <- function(disease = cancer, data.type = data.type, type = type, data_dir = data_dir, force_reload = FALSE) {
  FILE = paste0(data_dir, "/mtx_", disease, "_", data.type, "_", type, ".rda") # R object with data
  if (all(file.exists(FILE), !(force_reload))) {
    # If the data has been previously saved, load it
    load(file = FILE)
  } else {
    # If no saved data exists, get it from the remote source
    mtx <- getTCGA(disease = disease, data.type = data.type, type = type, clinical = TRUE)
    save(file = FILE, list = c("mtx")) # Save it
  }
  return(mtx)
}
summarize_data <- function(mtx = mtx) {
  print(paste0("Dimensions of expression matrix, genex X patients: ", paste(dim(mtx$dat), collapse = " ")))
  print(paste0("Dimensions of clinical matrix, patients X parameters: ", paste(dim(mtx$clinical), collapse = " ")))
  print(paste0("Dimensions of merged matrix, patients X parameters + genes: ", paste(dim(mtx$merged.dat), collapse = " ")))
  print("Head of the merged matrix")
  print(mtx$merged.dat[1:5, 1:10])
  print("Head of the clinical matrix")
  print(mtx$clinical[1:5, 1:7])
  print("List of clinical values, and frequency of each variable: ")
  clin_vars <- apply(mtx$clinical, 2, function(x) length(table(x[ !(is.na(x) & x != "" )]))) %>% as.data.frame()
  # Filter clinical variables to have at least 2, but no more than 10 categories,
  # And they are not dates
  clin_vars <- clin_vars[ as.numeric(clin_vars$.) > 1 & as.numeric(clin_vars$.) < 10 & !grepl("years|days|date|vital|OS|RFS|TIME|sample_type", rownames(clin_vars), perl = TRUE) , , drop = FALSE]
  print(kable(clin_vars))
  return(rownames(clin_vars))
}

# A function to create expression matrix
make_expression_matrix <- function(mtx = mtx, disease = cancer, data.type = data.type, type = type, results_dir = results_dir) {
  # transposed expression matrix (genes start at column 4) 
  mtx.expression <- mtx$merged.dat[, 4:ncol(mtx$merged.dat) ] %>% t 
  # Set column names as patient IDs
  colnames(mtx.expression) <- mtx$merged.dat$bcr 
  # Set row names as probe IDs
  rownames(mtx.expression) <- colnames(mtx$merged.dat)[ 4:ncol(mtx$merged.dat) ] 
  # Save gzipped matrix
  fileName.gz <- gzfile(paste0(results_dir, "/mtx_", disease, "_", data.type, "_", type, "_1expression.txt.gz"), "w")
  write.table(mtx.expression, fileName.gz, sep = ";", quote = FALSE)
  close(fileName.gz)
}
summarize_data <- function(mtx = mtx) {
  print(paste0("Dimensions of expression matrix, genex X patients: ", paste(dim(mtx$dat), collapse = " ")))
  print(paste0("Dimensions of clinical matrix, patients X parameters: ", paste(dim(mtx$clinical), collapse = " ")))
  print(paste0("Dimensions of merged matrix, patients X parameters + genes: ", paste(dim(mtx$merged.dat), collapse = " ")))
  print("Head of the merged matrix")
  print(mtx$merged.dat[1:5, 1:10])
  print("Head of the clinical matrix")
  print(mtx$clinical[1:5, 1:7])
  print("List of clinical values, and frequency of each variable: ")
  clin_vars <- apply(mtx$clinical, 2, function(x) length(table(x[ !(is.na(x) & x != "" )]))) %>% as.data.frame()
  # Filter clinical variables to have at least 2, but no more than 10 categories,
  # And they are not dates
  clin_vars <- clin_vars[ as.numeric(clin_vars$.) > 1 & as.numeric(clin_vars$.) < 10 & !grepl("years|days|date|vital|OS|RFS|TIME|sample_type", rownames(clin_vars), perl = TRUE) , , drop = FALSE]
  print(kable(clin_vars))
  return(rownames(clin_vars))
}

data_dir = "/home/di/Documents/Data/GenomeRunner/TCGAsurvival/data" # Mikhail
# data_dir = "/Users/stevenmeas/TCGAsurvival/data" # Steve
# Cancer types: http://www.liuzlab.org/TCGA2STAT/CancerDataChecklist.pdf
# Data types: http://www.liuzlab.org/TCGA2STAT/DataPlatforms.pdf
# Expression types: http://www.liuzlab.org/TCGA2STAT/DataPlatforms.pdf
# Clinical values: http://www.liuzlab.org/TCGA2STAT/ClinicalVariables.pdf
# General settings
data.type = "RNASeq2"
type = "" 
# data.type = "miRNASeq" # miRNASeq - "count" for raw read counts (default); "rpmmm" for normalized read counts
# type = "rpmmm"
# data.type = "Mutation"
# type = "somatic"
# type="all"
# data.type = "Methylation"
# type = "450K"
# clinical = TRUE
# Select cancer type
cancer = "PAAD" # Breast cancer
# cancer = "OV"   # Ovarian cancer
# cancer = "LIHC" # Liver hepatocellular carcinoma
# cancer = "HNSC" # Head and Neck squamous cell carcinoma
# cancer = "SARC" # Sarcoma
# cancer = "PAAD" # Pancreatic cancer
# cancer = "LUAD" # Lung adenocarcinoma
# cancer = "LUSC" # Lung squamous cell carcinoma
# cancer = "GBMLGG" # Glioma
# Select genes of interest
selected_genes = "LOXL2"
# Remove previous results
#system(paste0("rm -r ", selected_genes, ".", cancer, ".Analysis*"))
# Censor KM plots at certain duration
max_days <- 365 * 5 # Cutoff for maximum survival time, days x years 
# If true, Analysis 3 (survival differences in subcategories) will be run for each cancer
subcategories_in_all_cancers <- TRUE
######################
mtx <- load_data(disease = cancer, data.type = data.type, type = type, data_dir = data_dir, force_reload = FALSE)
clinical_annotations <- summarize_data(mtx = mtx)
# Prepare expression data
expr <- mtx$merged.dat[ , 4:ncol(mtx$merged.dat)] %>% as.matrix
# Filter out low expressed genes
# Should be more than 90% of non-zero values
ff <- genefilter::pOverA(p = 0.9, A = 0, na.rm = TRUE) 
expr <- expr[, apply(expr, 2, ff)] 
expr <- data.frame(AffyID = mtx$merged.dat$bcr, expr, stringsAsFactors = FALSE)

# Prepare clinical data
clin <- mtx$merged.dat[, 1:3]
colnames(clin)[1] <- "AffyID"
##
system("mkdir res") # Create results folder
# Censor KM plots at certain duration
index <- clin$OS <= max_days & !is.na(clin$OS)
clin <- clin[index, ]
expr <- expr[index, ]
# Run survival analysis for selected genes
kmplot(expr, clin, event_index=2, time_index=3,  affyid = selected_genes, auto_cutoff="true", transform_to_log2 = TRUE, cancer_type = cancer, fileType = "png", use_survminer = TRUE)
dev.off()
# Rename the results folder
system(paste0("mv res ", selected_genes, ".", cancer, ".Analysis1"))
###############
# All cancers with RNASeq2 data
all_cancers = c("ACC", "BLCA", "BRCA", "CESC", "CHOL", "COAD", "COADREAD", "DLBC", "ESCA", "GBM", "GBMLGG", "HNSC", "KICH", "KIPAN", "KIRC", "KIRP", "LGG", "LIHC", "LUAD", "LUSC", "MESO", "OV", "PAAD", "PCPG", "PRAD", "READ", "SARC", "SKCM", "STAD", "STES", "TGCT", "THCA", "THYM", "UCEC", "UCS", "UVM") # "LAML", 
data.type = "RNASeq2"; type = "" 
system("mkdir res") # Create results folder
for (cancer_type in all_cancers) {
  print(paste0("Processing cancer ", cancer_type))
  mtx <- load_data(disease = cancer_type, data.type = data.type, type = type, data_dir = data_dir, force_reload = FALSE)
  # Prepare expression data
  expr <- mtx$merged.dat[ , 4:ncol(mtx$merged.dat)] %>% as.matrix
  # Filter out low expressed genes
  # Should be more than 90% of non-zero values
  # ff <- genefilter::pOverA(p = 0.9, A = 0, na.rm = TRUE) 
  # expr <- expr[, apply(expr, 2, ff)] 
  expr <- data.frame(AffyID = mtx$merged.dat$bcr, expr, stringsAsFactors = FALSE)
  # Prepare clinical data
  clin <- mtx$merged.dat[, 1:3]
  colnames(clin)[1] <- "AffyID"
  # Censor KM plots at certain duration
  index <- clin$OS <= max_days & !is.na(clin$OS)
  clin <- clin[index, ]
  expr <- expr[index, ]
  # Run survival analysis for selected genes
  kmplot(expr, clin, event_index=2, time_index=3,  affyid = selected_genes, auto_cutoff="true", transform_to_log2 = TRUE, cancer_type = cancer_type, fileType = "png", use_survminer = TRUE)
}
### Plot the results of one gene across all cancers
# Read in analysis natrix
mtx <- read.table("res/global_stats.txt", sep = "\t", header = TRUE, stringsAsFactors = FALSE, fill = TRUE)
# Add -log10-transformed p-value
mtx <- mtx %>% mutate(log10.pval = -1 * log10(p.adjust(p.value, method = "BH")))
# Print into PNG
mtx_to_plot <- mtx %>% subset(Gene == selected_genes)
mtx_to_plot <- mtx_to_plot[order(mtx_to_plot$log10.pval), ]
mtx_to_plot$Cancer <- factor(mtx_to_plot$Cancer, levels = mtx_to_plot$Cancer)
ggplot(mtx_to_plot, aes(x = Cancer, y = log10.pval)) + 
  geom_bar(stat = "identity") + 
  theme(legend.position="none") +
  labs(x="Cancer", y="-log10(p-value)") +
  coord_flip() +
  theme_cowplot() 
ggsave(paste0("res/", selected_genes, "_all_TCGA_cancers.png"), width = 5, height = 6.5)
# Rename the results folder
system(paste0("mv res ", selected_genes, ".", cancer, ".Analysis2"))
########## analysis 3
System("mkdir res") # Create results folder
# If true, analyze survival differences in subcategories in all cancers
if (subcategories_in_all_cancers) {
  cancer_type <- c("ACC", "BLCA", "BRCA", "CESC", "CHOL", "COAD", "COADREAD", "DLBC", "ESCA", "GBM", "GBMLGG", "HNSC", "KICH", "KIPAN", "KIRC", "KIRP", "LGG", "LIHC", "LUAD", "LUSC", "MESO", "OV", "PAAD", "PCPG", "PRAD", "READ", "SARC", "SKCM", "STAD", "STES", "TGCT", "THCA", "THYM", "UCEC", "UCS", "UVM") # "LAML", 
} else {
  cancer_type <- cancer # Selected in settings
}
data.type = "RNASeq2"; type = "" 
for (cancer_type in cancer_type) {
  print(paste0("Processing cancer ", cancer_type))
  mtx <- load_data(disease = cancer_type, data.type = data.type, type = type, data_dir = data_dir, force_reload = FALSE)
  if (cancer_type == "BRCA") {
    # BRCA-specific - replace original annotations with XENA	
    mtx$clinical <- read.csv("data.TCGA/XENA_classification.csv", row.names = 1)
  }
  clinical_annotations <- summarize_data(mtx = mtx)
  # Prepare expression data
  expr <- mtx$merged.dat[ , 4:ncol(mtx$merged.dat)] %>% as.matrix
  # Filter out low expressed genes
  # Should be more than 90% of non-zero values
  # ff <- genefilter::pOverA(p = 0.9, A = 0, na.rm = TRUE) 
  # expr <- expr[, apply(expr, 2, ff)] 
  expr <- data.frame(AffyID = mtx$merged.dat$bcr, expr, stringsAsFactors = FALSE)
  # Prepare clinical data
  clin <- mtx$merged.dat[, 1:3]
  colnames(clin)[1] <- "AffyID"
  # Censor KM plots at certain duration
  index <- clin$OS <= max_days & !is.na(clin$OS)
  clin <- clin[index, ]
  expr <- expr[index, ]
  # Full clinical information
  
  # clin_full <- mtx$clinical
  # Match to the order of small clinical annitation	
  # clin_full <- clin_full[rownames(clin_full) %in% clin$AffyID, ]	
  # clin_full <- clin_full[match(clin$AffyID, rownames(clin_full)), ]
  
  clin_full <- mtx$clinical[rownames(mtx$clinical) %in% clin$AffyID, ]
  clin_full <- clin_full[match(expr$AffyID, rownames(clin_full)), ]
  all.equal(rownames(clin_full), expr$AffyID)
  
  if (cancer_type == "OV") {
    # OV - specific Prepare extra clinical annotation
    clin_extra <- read.table("data.TCGA/TCGA_489_UE.k4.txt", sep = "\t", header = TRUE)
    clin_extra$Sample.ID <- sapply(clin_extra$Sample.ID, function(x) strsplit(x, ".", fixed = TRUE)[[1]][1:3] %>% paste(., collapse = "-")) # Shorten sample IDs
    # Subset to common samples
    common_AffyID <- intersect(rownames(clin_full), clin_extra$Sample.ID)
    expr <- expr[ expr$AffyID %in% common_AffyID, ]
    clin <- clin[ clin$AffyID %in% common_AffyID, ]
    clin_full <- clin_full[ rownames(clin_full) %in% common_AffyID, ]
    clin_extra <- clin_extra[ clin_extra$Sample.ID %in% common_AffyID, ]
    # Make clinical annotation the same order
    clin <- clin[ match(expr$AffyID, clin$AffyID), ]
    clin_extra <- clin_extra[ match(expr$AffyID, clin_extra$Sample.ID), ]
    all.equal(expr$AffyID, clin$AffyID, rownames(clin_full), clin_extra$Sample.ID)
    # Add extra clinical annotation
    clin_full <- data.frame(clin_full, subtype = clin_extra$k4)
    clinical_annotations <- c( clinical_annotations, "subtype")
  }
  
  # For each clinical annotation
  for (annotation in clinical_annotations) { 
    # Get the number of patients per category in the current annotation
    annotations <- table(clin_full[, annotation], useNA = "no") 
    # How many categories in the current annotation
    num_of_annot_categories <- length(annotations) 
    # How many categories to select at one time
    for (num_of_selected_categories in 1:num_of_annot_categories) {
      # Select patients annotated with this categories
      patients <- rownames(clin_full)[ clin_full[, annotation] %in% names(annotations)[ num_of_selected_categories] ]
      # Get their index in the clin and expr matrixes
      index_patients <- which(clin$AffyID %in% patients)
      # If the number of patients annotated with the combination of categories is large enough, proceed
      if (length(index_patients) > 40) {
        print(paste("Processing annotation:", annotation, 
                    ", categories:", names(annotations)[ num_of_selected_categories],
                    ", number of patients:", length(index_patients)))
        # Get a subset of clinical information for these patients
        clin_selected <- clin[ index_patients, ]
        # Get a subset of expression information for these patients
        expr_selected <- expr[ index_patients, ]
        # For this subset of expression, filter out low expressed genes
        index_genes <- apply(expr_selected %>% dplyr::select(-AffyID), 2, ff) # index of expression values to keep
        expr_selected <- cbind(expr_selected$AffyID, dplyr::select(expr_selected, -AffyID)[, index_genes]) # patient IDs and expression values to keep
        # Perform actual survival analysis
        kmplot(expr_selected, clin_selected, event_index=2, time_index=3,  affyid = selected_genes, auto_cutoff="true", transform_to_log2 = TRUE, cancer_type = paste(c(cancer_type, annotation, names(annotations)[ num_of_selected_categories] ), collapse = "-"), fileType = "png", use_survminer = TRUE)
      }
    }
  }
}
# Rename the results folder
system(paste0("mv res ", selected_genes, ".", cancer, ".Analysis3"))

### Analysis 4: Selected genes, selected cancers, all combinations of clinical annotations
system("mkdir res") # Create results folder
clin_full <- mtx$clinical # Full clinical information
# For each clinical annotation
for (annotation in clinical_annotations) { 
  # Get the number of patients per category in the current annotation
  annotations <- table(clin_full[, annotation], useNA = "no") 
  # How many categories in the current annotation
  num_of_annot_categories <- length(annotations) 
  # How many categories to select at one time
  for (num_of_selected_categories in 1:num_of_annot_categories) {
    # All combinations of categories, m categories selected at one time
    combination_of_categories <- combn(x = names(annotations), m = num_of_selected_categories)
    # For each combination of categories (column)
    for (combination in 1:ncol(combination_of_categories)) {
      # Select patients annotated with this combination of categories
      patients <- rownames(clin_full)[ clin_full[, annotation] %in% combination_of_categories[, combination]]
      # Get their index in the clin and expr matrixes
      index_patients <- which(clin$AffyID %in% patients)
      # If the number of patients annotated with the combination of categories is large enough, proceed
      if (length(index_patients) > 40) {
        print(paste("Processing annotation:", annotation, 
                    ", categories:", paste(combination_of_categories[, combination], collapse = ","),
                    ", number of patients:", length(index_patients)))
        # Get a subset of clinical information for these patients
        clin_selected <- clin[ index_patients, ]
        # Get a subset of expression information for these patients
        expr_selected <- expr[ index_patients, ]
        # For this subset of expression, filter out low expressed genes
        ff <- genefilter::pOverA(p = 0.9, A = 0, na.rm = TRUE) 
        expr <- expr[, apply(expr, 2, ff)] 
        index_genes <- apply(expr_selected %>% dplyr::select(-AffyID), 2, ff) # index of expression values to keep
        expr_selected <- cbind(expr_selected$AffyID, select(expr_selected, -AffyID)[, index_genes]) # patient IDs and expression values to keep
        # Perform actual survival analysis
        kmplot(expr_selected, clin_selected, event_index=2, time_index=3,  affyid = selected_genes, auto_cutoff="true", transform_to_log2 = TRUE, cancer_type = paste(c(cancer, annotation, combination_of_categories[, combination]), collapse = "-"), fileType = "png", use_survminer = FALSE)
      }
    }
  }
}
# Rename the results folder
system("mv res res.genes.Analysis4")
```

```{r}
if (cancer %in% c("BRCA","OV")){
  continue_analysis <- TRUE
} else {
  continue_analysis <- FALSE
}
opts_chunk$set(eval = continue_analysis, include = continue_analysis)
```
### Analysis 5: Clinical-centric analysis. Selected cancer, selected clinical subcategory, survival difference between all pairs of subcategories

```{r}
system("mkdir -p res")
# Select cancer and load expression data

# cancer <- "BRCA"
mtx <- load_data(disease = cancer, data.type = data.type, type = type, data_dir = data_dir, force_reload = FALSE)
clinical_annotations <- summarize_data(mtx = mtx)
# Prepare expression data
expr <- mtx$merged.dat[ , 4:ncol(mtx$merged.dat)] %>% as.matrix
# Filter out low expressed genes
# Should be more than 90% of non-zero values
ff <- genefilter::pOverA(p = 0.9, A = 0, na.rm = TRUE) 
expr <- expr[, apply(expr, 2, ff)] 
expr <- data.frame(AffyID = mtx$merged.dat$bcr, expr, stringsAsFactors = FALSE)

clinical_annotations_selected <- "pathologyMstage" # Default selected clinical annotation

# Prepare clinical data
if (cancer == "BRCA") {
  # BRCA-specific - replace original annotations with XENA	
  mtx$clinical <- read.csv("data.TCGA/XENA_classification.csv", row.names = 1)
  clin <- mtx$merged.dat[, 1:3]
  colnames(clin)[1] <- "AffyID"
  # Full clinical information
  clin_full <- mtx$clinical
  # Match to the order of small clinical annitation
  clin_full <- clin_full[rownames(clin_full) %in% clin$AffyID, ]
  clin_full <- clin_full[match(clin$AffyID, rownames(clin_full)), ]
  all.equal(expr$AffyID, rownames(clin_full))
  # Select clinical category and subcategories
  clinical_annotations <- summarize_data(mtx = mtx)
  clinical_annotations_selected <- "PAM50Call_RNAseq" # For PAM50 classifier
  clin <- cbind(clin, clin_full[, clinical_annotations_selected, drop = FALSE])
} else {
  clin <- mtx$merged.dat[, 1:3]
  colnames(clin)[1] <- "AffyID"
  # Full clinical information
  clin_full <- mtx$clinical
  # Match to the order of small clinical annitation
  clin_full <- clin_full[rownames(clin_full) %in% clin$AffyID, ]
  clin_full <- clin_full[match(clin$AffyID, rownames(clin_full)), ]
  all.equal(expr$AffyID, rownames(clin_full))
}


if (cancer == "OV") {
  # OV-specific Prepare extra clinical annotation
  clin_extra <- read.table("data.TCGA/TCGA_489_UE.k4.txt", sep = "\t", header = TRUE)
  clin_extra$Sample.ID <- sapply(clin_extra$Sample.ID, function(x) strsplit(x, ".", fixed = TRUE)[[1]][1:3] %>% paste(., collapse = "-")) # Shorten sample IDs
  # Subset to common samples
  common_AffyID <- intersect(rownames(clin_full), clin_extra$Sample.ID)
  expr <- expr[ expr$AffyID %in% common_AffyID, ]
  clin <- clin[ clin$AffyID %in% common_AffyID, ]
  clin_full <- clin_full[ rownames(clin_full) %in% common_AffyID, ]
  clin_extra <- clin_extra[ clin_extra$Sample.ID %in% common_AffyID, ]
  # Make clinical annotation the same order
  clin <- clin[ match(expr$AffyID, clin$AffyID), ]
  clin_extra <- clin_extra[ match(expr$AffyID, clin_extra$Sample.ID), ]
  all.equal(expr$AffyID, clin$AffyID, rownames(clin_full), clin_extra$Sample.ID)
  # Add extra clinical annotation
  clin_full <- data.frame(clin_full, subtype = clin_extra$k4)
  
  # Select clinical category and subcategories
  clinical_annotations_selected <- "subtype" # For OV classifier
  clin <- cbind(clin, clin_full[, clinical_annotations_selected, drop = FALSE])
}

print(paste0("Number of patients in each subcategory, in the ", clinical_annotations_selected, " category"))
clin <- cbind(clin, clin_full[, clinical_annotations_selected, drop = FALSE])
table(clin[, clinical_annotations_selected]) 
# All pairs of subcategories
group_combinations <- combn(unique(clin[!is.na(clin[, clinical_annotations_selected]), clinical_annotations_selected]), 2) 
# Survival analysis for all group combinations
for (pair in 1:ncol(group_combinations)) {
  # KM analysis on subcategories
  kmplot.clin(clin = clin, event_index=2, time_index=3, clinical_annotations = clinical_annotations_selected, group1 = group_combinations[1, pair], group2 = group_combinations[2, pair], cancer_type = cancer, fileType = "png", use_survminer = TRUE)
}


# Plot expression boxplots for the selected gene in selected subcategories
# selected_genes = c("SDC1") # BRCA
# Prepare expression data
expr <- mtx$merged.dat[ , 4:ncol(mtx$merged.dat)] %>% as.matrix
# Filter out low expressed genes
# Should be more than 90% of non-zero values
# ff <- genefilter::pOverA(p = 0.9, A = 0, na.rm = TRUE) 
# expr <- expr[, apply(expr, 2, ff)] 
expr <- data.frame(AffyID = mtx$merged.dat$bcr, expr, stringsAsFactors = FALSE)
expr <- expr[ expr$AffyID %in% clin$AffyID, ]
expr <- expr[match(clin$AffyID, expr$AffyID), ]
# A matrix to plot
mtx_to_plot <- data.frame(Gene = log2(expr[, selected_genes] + 1), Clinical = clin[, clinical_annotations_selected])
mtx_to_plot <- mtx_to_plot[complete.cases(mtx_to_plot), ]
save(mtx_to_plot, file = "res/mtx_to_plot.rda")
# Plot and save
p <- ggplot(melt(mtx_to_plot, id.vars = "Clinical"), aes(x = Clinical, y = value, fill = Clinical)) +
  geom_boxplot() +
  ylab("log2 expression") + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
plot(p)
ggsave(filename = paste0("res/", cancer, "_", selected_genes, "_", clinical_annotations_selected, ".png"), p, width = 5, height = 4, device = "png")
system(paste0("mv res ", selected_genes, ".", cancer, ".Analysis5"))
```

```{r eval = FALSE}
## Analysis 6: Dimensionality reduction of a gene signature across all cancers using NMF, PCA, or FA
# For each cancer, extracts gene expression of a signature, reduces its dimensionality,
# plots a heatmap sorted by the first component, biplots, saves eigenvectors in files
# named after cancer, signature, method. They are used in `correlations.Rmd`
system("rm -r res")
system("mkdir res") # Create results folder
library(NMF)
library(ggfortify)
library(pheatmap)
signature      <- readLines("data.TCGA/EINAV_INTERFERON_SIGNATURE_IN_CANCER.txt") # Interferon signature
selected_genes <- "interferon_signature" 
signature=(strsplit(signature, ",")[[1]])
pdf(file = paste0("res/dimreduction_", selected_genes, ".pdf"))
# All cancers with RNASeq2 data
cancer_RNASeq2 = c("ACC", "BLCA", "BRCA", "CESC", "CHOL", "COAD", "COADREAD", "DLBC",
                   "ESCA", "GBM", "HNSC", "KICH", "KIPAN", 
                   "KIRC", "KIRP", "LIHC", "LUAD", "LUSC", 
                   "MESO", "OV", "PAAD", "PCPG", "PRAD", 
                   "READ", "SARC", "SKCM", "STAD", "TGCT", 
                   "THCA", "THYM", "UCEC", "UCS") # "GBMLGG", "LGG", 
cancer_RNASeq2 = c("ACC", "BLCA", "BRCA", "CESC", "CHOL", "COAD", "COADREAD", "DLBC",
                   "ESCA", "GBM", "HNSC", "KICH",  
                   "LIHC", "LUAD", "LUSC", 
                   "MESO", "OV", "PAAD", "PCPG", "PRAD", 
                   "READ", "SARC", "SKCM", "STAD", "TGCT", 
                   "THCA", "THYM", "UCEC", "UCS") # "GBMLGG", "LGG", 
data.type = "RNASeq2"; type = ""
#cancer_RNASeq2 = "PAAD"
# Process each cancer
for (cancer_type in cancer_RNASeq2) {
  print(paste0("Processing cancer ", cancer_type))
  mtx <- load_data(disease = cancer_type, data.type = data.type, type = type, data_dir = data_dir, force_reload = FALSE)
  # Prepare subset of expression data for a given signature
  signature=intersect(signature,colnames(mtx$merged.dat))
  mtx_selected <- t(log2(mtx$merged.dat[, signature]))
  colnames(mtx_selected) <- mtx$merged.dat$bcr
  # NMF 
  method <- "NMF"
  mtx_nmf <- nmf(t(mtx_selected), 3) # with three components
  pheatmap(mtx_selected[, order(mtx_nmf@fit@W[, 1])], cluster_cols = FALSE, treeheight_row = FALSE,   treeheight_col = FALSE, scale = "row", main = paste0(cancer_type, "_", method))
  mtx_reduced <- mtx_nmf@fit@W
  save(mtx_reduced, file = paste0("res/", cancer_type, "_", selected_genes, "_", method, ".Rda"))
  # PCA
  # https://cran.r-project.org/web/packages/ggfortify/vignettes/plot_pca.html
  method <- "PCA"
  mtx_pca <- prcomp(t(mtx_selected), scale. = TRUE, center = TRUE)
  pheatmap(mtx_selected[, order(mtx_pca$x[, 1])], cluster_cols = FALSE, treeheight_row = FALSE, treeheight_col = FALSE, scale = "row", main = paste0(cancer_type, "_", method))
  p <- autoplot(prcomp(t(mtx_selected), scale. = TRUE, center = TRUE), loadings = TRUE, loadings.label = TRUE) + ggtitle(paste0(cancer_type, "_", method))
  print(p)
  mtx_reduced <- mtx_pca$x
  save(mtx_reduced, file = paste0("res/", cancer_type, "_", selected_genes, "_", method, ".Rda"))
  # Factor analysis
  method <- "FA"
  mtx_fa <- factanal(t(mtx_selected), factors = 3, scores = "regression")
  pheatmap(mtx_selected[, order(mtx_fa$scores[, 1])], cluster_cols = FALSE, treeheight_row = FALSE, treeheight_col = FALSE, scale = "row", main = paste0(cancer_type, "_", method))
  p <- autoplot(mtx_fa,  loadings = TRUE, loadings.label = TRUE) + ggtitle(paste0(cancer_type, "_", method))
  print(p)
  mtx_reduced <- mtx_fa$scores
  save(mtx_reduced, file = paste0("res/", cancer_type, "_", selected_genes, "_", method, ".Rda"))
}
dev.off()
system("mv res res.genes.Analysis6666")
```

# References

