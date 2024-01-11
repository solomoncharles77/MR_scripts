# Run this script in R/4.3.1

# This script takes one argument
# - the path to the input hg19 gwas summary statistics file
# This script outputs a gwas summary statistics file with coordID and rsID

# Set global options and load libraries -----------------------------------
options(scipen = 999, "warnPartialMatchDollar"=TRUE)
# Load libraries
suppressPackageStartupMessages(library(optparse))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(TwoSampleMR))
suppressPackageStartupMessages(library(doMC))


# List options
option_list = list(
  make_option(c("-e", "--expogwas"), type="character", default=NULL, 
              help="path to exposure gwas summary statistics", metavar="character"),
  make_option(c("-f", "--outcgwas"), type="character", default=NULL, 
              help="path to outcome gwas summary statistics", metavar="character"),
  make_option(c("-t", "--expPhe"), type="character", default=NULL, 
              help="Trait identifier for expogwas e.g SBP", metavar="character"),
  make_option(c("-u", "--outPhe"), type="character", default=NULL, 
              help="Trait identifier for outcgwas e.g SCAD", metavar="character"),
  make_option(c("-c", "--cand"), type="character", default=NULL, 
              help="Path to noRHQ .txt file containing candidate gene list.", metavar="character"),
  make_option(c("-p", "--pcutexp"), type="numeric", default = 5e-03, 
              help="pvalue threshold for selecting exposure associations.", metavar="numeric"),
  make_option(c("-q", "--pcutout"), type="numeric", default = 1e-00, 
              help="pvalue threshold for selecting outcome associations.", metavar="numeric"),
  make_option(c("-s", "--sampsize"), type="numeric", default = NULL, 
              help="The sample size of the eQTL summary statistics.", metavar="numeric")
  
); 

opt_parser <-  OptionParser(option_list=option_list);
opt <-  parse_args(opt_parser);


# Data input error messages
if (is.null(opt$expogwas)){
  print_help(opt_parser)
  stop("Please provide path to exposure gwas summ stats.", call.=FALSE)
}

if (is.null(opt$outcgwas)){
  print_help(opt_parser)
  stop("Please provide path to outcome gwas summ stats.", call.=FALSE)
}

if (is.null(opt$expPhe)){
  print_help(opt_parser)
  stop("Please provide identifier for exposure gwas files.", call.=FALSE)
}

if (is.null(opt$cand)){
  print_help(opt_parser)
  stop("Please indicate ensemble ID for gene to be analysed.", call.=FALSE)
}

if (is.null(opt$sampsize)){
  print_help(opt_parser)
  stop("Please indicate the sample size of the eQTL.", call.=FALSE)
}

if (is.null(opt$outPhe)){
  print_help(opt_parser)
  stop("Please provide identifier for outcome gwas files.", call.=FALSE)
}


# # Assign input variables --------------------------------------------------
expoFileName <- paste0(opt$expogwas)
#expoFileName <- "/scratch/vasccell/cs806/pubAvail_QTL/GTEx/eQTL_catalogue/QTD000251_GTEx_heart_atrial_appendage_pval005.tsv.gz"
outcFileName <- paste0(opt$outcgwas)
#outcFileName <- "/scratch/vasccell/cs806/colocalization/cleanGWAS_Summary_Stats/GWAS_Sponteneous_Coronary_Artery_Dissection_Female_Adlam_2023_NatGenet_hg38.txt"
pcutoff1 <- as.numeric(opt$pcutexp)
pcutoff2 <- as.numeric(opt$pcutout)
expID <- paste0(opt$expPhe)
#expID <- "GTEX_Aorta"
outID <- paste0(opt$outPhe)
#outID <- "SCAD_Adlam2023"
cand <- paste0(opt$cand)
#cand <- "/alice-home/2/c/cs806/Documents/frailtyGWAS/frailtyMR_scripts/target1.txt"
ss <- as.numeric(opt$sampsize)
#ss <- 387


resDir1 <- paste0("mrResults/",  expID, "_", outID, "/")
if (!dir.exists(resDir1)) {
  dir.create(resDir1)
}

resDir2 <- paste0("mrFiles/",  expID, "_", outID, "/")
if (!dir.exists(resDir2)) {
  dir.create(resDir2)
}


gwas1 <- data.frame(fread(expoFileName))
gwas2 <- data.frame(fread(outcFileName))
cands <- read.table(cand)
invisible(gc())
cat("\n")
cat("\n")

# Subset SNPs that meet association threshold
gwas2_sig <- gwas2[which(gwas2$pvalue < pcutoff2), ]
gwas1_sig <- gwas1[which(gwas1$pvalue < pcutoff1), ]
gwas1_sig <- gwas1_sig[which(gwas1_sig$gene %in% cands$V1), ]
setnames(gwas1_sig, c("allele1", "allele2"),
         c("alt", "ref"), skip_absent = TRUE)

# How many genes have genome wide significant associations
topSigQTL <- gwas1_sig[!duplicated(gwas1_sig$gene), ]

mrRes <- foreach(rr = 1:nrow(topSigQTL)) %dopar% {
  cand <- topSigQTL[rr, "gene"]
  
  # Extract candidate gene
  eGene <- gwas1_sig[gwas1_sig$gene == cand, ]
  
  if (dim(eGene)[1] > 0) {
    eGene$n <- ss
  } else if (dim(eGene)[1] == 0) {
    cat("\n")
    cat("No exposure SNP made the pvalue threshold:", pcutoff1)
    cat("\n")
    cat("\n")
    cat("Consider adjusting the pvalue cutoffs")
    cat("\n")
  }
  
  if (dim(eGene)[1] == 0) {
    stop("Stopping the script.")
  }
  
  
  eGene$Phenotype <- paste0(cand, "_", expID)
  gwas1_match <- eGene[, c("snps", "rsid", "alt", "ref", "pvalue", "beta", "se", "n", "Phenotype")]
  colnames(gwas1_match) <- c("hg38_markername", "SNP", "effect_allele","other_allele", "pval", "beta", "se", "samplesize", "Phenotype")
  
  exposure_dat <- format_data(gwas1_match, type = "exposure")
  
  # Clumping data
  exp_clump <- clump_data(exposure_dat, clump_r2 = 0.001)
  
  
  # Harmonize
  gwas2_match <- gwas2_sig[gwas2_sig$permID %in% exp_clump$SNP,  ]
  gwas2_match <- gwas2_match[, c("hg38_markername", "permID", "ALT", "REF", "pvalue", "beta", "se", "n")]
  colnames(gwas2_match) <- c("hg38_markername", "SNP", "effect_allele","other_allele", "pval", "beta", "se", "samplesize")
  
  if (dim(gwas2_match)[1] >= 3) {
    gwas2_match$Phenotype <- outID
  } else if (dim(gwas2_match)[1] < 3) {
    cat("\n")
    cat("No overlapping IVs found after clumping, and at pvalue cutoff; exposure:", pcutoff1, " outcome:", pcutoff2)
    cat("\n")
    cat("\n")
    cat("Consider adjusting the pvalue cutoffs")
    cat("\n")
  }
  
  if (dim(gwas2_match)[1] < 3) {
    stop("Stopping the script.")
  }
  
  #Order
  exp_clump <- exp_clump[order(exp_clump$SNP), ]
  gwas2_match <- gwas2_match[order(gwas2_match$SNP), ]
  
  # Reformat
  outcome_dat <- format_data(gwas2_match, type = "outcome")
  
  # Export IVs
  # write.csv(exposure_dat, paste0("mrFiles/", expID, "_", cand, "_IVs_exposure.csv"), row.names = F)
  # write.csv(outcome_dat, paste0("mrFiles/", outID, "_", cand,  "_IVs_outcome.csv"), row.names = F)
  
  # exposure_dat <- read.csv(paste0("mrFiles/", expID, "_IVs_exposure.csv"))
  # outcome_dat <- read.csv(paste0("mrFiles/", outID, "_IVs_outcome.csv"))
  
  
  # Harmonisation
  dat <- harmonise_data(exposure_dat = exposure_dat, outcome_dat = outcome_dat)
  
  mr_report(dat,
            output_path = resDir1,
            author = "Charles Solomon",
            study = "Frailty_MR")
  
  write.csv(dat, paste0(resDir2, cand, "_" ,expID, "_against_", outID, "_clumped_IVs.csv"), row.names = F)
  
  # TODO
  # Output file to include ID of cand gene
  
  
}
