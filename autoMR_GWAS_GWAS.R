# Run this script in R/4.2.1

# This script takes one argument
# - the path to the input hg19 gwas summary statistics file
# This script outputs a gwas summary statistics file with coordID and rsID

# Set global options and load libraries -----------------------------------
options(scipen = 999, "warnPartialMatchDollar"=TRUE)
# Load libraries
suppressPackageStartupMessages(library(optparse))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(TwoSampleMR))


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
  make_option(c("-p", "--pcutexp"), type="numeric", default = 5e-08, 
              help="pvalue threshold for selecting exposure associations.", metavar="numeric"),
  make_option(c("-q", "--pcutout"), type="numeric", default = 1e-00, 
              help="pvalue threshold for selecting outcome associations.", metavar="numeric")
  
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

if (is.null(opt$outPhe)){
  print_help(opt_parser)
  stop("Please provide identifier for outcome gwas files.", call.=FALSE)
}


# # Assign input variables --------------------------------------------------
expoFileName <- paste0(opt$expogwas)
#expoFileName <- "/scratch/vasccell/cs806/colocalization/cleanGWAS_Summary_Stats/GWAS_Sex_Hormone-binding_Globulin_Levels_EurFemale_Ruth_2020_NatMed_hg38.txt"
outcFileName <- paste0(opt$outcgwas)
#outcFileName <- "/scratch/vasccell/cs806/colocalization/cleanGWAS_Summary_Stats/GWAS_Sponteneous_Coronary_Artery_Dissection_Female_Adlam_2023_NatGenet_hg38.txt"
pcutoff1 <- as.numeric(opt$pcutexp)
pcutoff2 <- as.numeric(opt$pcutout)
expID <- paste0(opt$expPhe)
#expID <- "SHBG_Ruth2020"
outID <- paste0(opt$outPhe)
#outID <- "SCAD_Adlam2023"

gwas1 <- data.frame(fread(expoFileName))
gwas2 <- data.frame(fread(outcFileName))
invisible(gc())
cat("\n")
cat("\n")

# Subset SNPs that meet association threshold
gwas1_sig <- gwas1[which(gwas1$pvalue < pcutoff1), ]
gwas2_sig <- gwas2[which(gwas2$pvalue < pcutoff2), ]


# Prep IVs from exposure GWAS
gwas1_sig <- gwas1_sig[, c("hg38_markername", "permID", "ALT", "REF", "eaf", "pvalue", "beta", "se", "n")]
colnames(gwas1_sig) <- c("hg38_markername", "SNP", "effect_allele","other_allele","eaf", "pval", "beta", "se", "samplesize")
gwas1_sig$Phenotype <- expID
exposure_dat <- format_data(gwas1_sig, type = "exposure")

# Clumping data
exp_clump <- clump_data(exposure_dat, clump_r2 = 0.001)


# Harmonize
gwas2_match <- gwas2_sig[gwas2_sig$permID %in% exp_clump$SNP,  ]
gwas2_match <- gwas2_match[, c("hg38_markername", "permID", "ALT", "REF", "eaf", "pvalue", "beta", "se", "n")]
colnames(gwas2_match) <- c("hg38_markername", "SNP", "effect_allele","other_allele","eaf", "pval", "beta", "se", "samplesize")

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
write.csv(exposure_dat, paste0("mrFiles/", expID, "_IVs_exposure.csv"), row.names = F)
write.csv(outcome_dat, paste0("mrFiles/", outID, "_IVs_outcome.csv"), row.names = F)

# exposure_dat <- read.csv(paste0("mrFiles/", expID, "_IVs_exposure.csv"))
# outcome_dat <- read.csv(paste0("mrFiles/", outID, "_IVs_outcome.csv"))


# Harmonisation
dat <- harmonise_data(exposure_dat = exposure_dat, outcome_dat = outcome_dat)

# mr_report(dat,
#           output_path = "mrResults",
#           author = "Charles Solomon",
#           study = "SCAD_MR")

mr_report(dat,
          output_path = "mrResults/reverse",
          author = "Charles Solomon",
          study = "SCAD_MR")

write.csv(dat, paste0("mrFiles/", expID, "_against_", outID, "_clumped_IVs.csv"), row.names = F)

