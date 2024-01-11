
# Run this script in R/4.2.1

library(data.table)
library(TwoSampleMR)
library(MendelianRandomization)
library(doMC)
library(ieugwasr)

# gwas1 <- data.frame(fread("/scratch/vasccell/cs806/colocalization/cleanGWAS_Summary_Stats/GWAS_Sponteneous_Coronary_Artery_Dissection_Female_Adlam_2023_NatGenet_hg38.txt"))
# gwas2 <- data.frame(fread("/scratch/vasccell/cs806/colocalization/cleanGWAS_Summary_Stats/GWAS_Sex_Hormone-binding_Globulin_Levels_EurFemale_Ruth_2020_NatMed_hg38.txt"))
# 
# gc()
# 
# pcutoff1 <- 1e-04
# pcutoff2 <- 5e-08
# 
# # Subset SNPs that meet association threshold
# gwas1_sig <- gwas1[which(gwas1$pvalue < pcutoff1), ]
# gwas2_sig <- gwas2[which(gwas2$pvalue < pcutoff2), ]
# 
# # Harmonize
# gwas1_match <- gwas1_sig[gwas1_sig$hg38_markername %in% gwas2_sig$hg38_markername, ]
# gwas2_match <- gwas2_sig[gwas2_sig$hg38_markername %in% gwas1_sig$hg38_markername, ]
# 
# gwas1_match <- gwas1_match[, c("hg38_markername", "permID", "ALT", "REF", "eaf", "pvalue", "beta", "se", "n")]
# gwas2_match <- gwas2_match[, c("hg38_markername", "permID", "ALT", "REF", "eaf", "pvalue", "beta", "se", "n")]
# 
# colnames(gwas1_match) <- c("hg38_markername", "SNP", "effect_allele","other_allele","eaf", "pval", "beta", "se", "samplesize")
# colnames(gwas2_match) <- c("hg38_markername", "SNP", "effect_allele","other_allele","eaf", "pval", "beta", "se", "samplesize")
# 
# gwas1_match$Phenotype <- "SCAD"
# gwas2_match$Phenotype <- "Horm"
# 
# #Order
# gwas1_match <- gwas1_match[order(gwas1_match$SNP), ]
# gwas2_match <- gwas2_match[order(gwas2_match$SNP), ]
# 
# write.csv(gwas1_match, "mrFiles/SCAD_Adlam2023_Female_IVs_outcome.csv", row.names = F)
# write.csv(gwas2_match, "mrFiles/SHBG_Ruth2020_Female_IVs_exposure.csv", row.names = F)

outc <- read.csv("mrFiles/SCAD_Adlam2023_Female_IVs_outcome.csv")
expo <- read.csv("mrFiles/SHBG_Ruth2020_Female_IVs_exposure.csv")

# Reformat
exposure_dat <- format_data(expo, type = "exposure")
outcome_dat <- format_data(outc, type = "outcome")

# Harmonisation
dat <- harmonise_data(exposure_dat = exposure_dat, outcome_dat = outcome_dat)

# Clumping data
dat_clump <- clump_data(dat[which(dat$mr_keep),], clump_r2 = 0.001)


dat_clump$R2 <- TwoSampleMR::get_r_from_bsen(b=dat_clump$beta.exposure,
                                             se=dat_clump$se.exposure,
                                             n=dat_clump$samplesize.exposure[1]) 

#dat_clump <- dat_clump[which(dat_clump$R2 > 0), ]

N <- dat_clump$samplesize.exposure[1]
k <- nrow(dat_clump)
dat_clump$F <- ((N-k-1)/(k)) * (dat_clump$R2/(1-dat_clump$R2))

#dat_clump <- dat_clump[which(dat_clump$F > 10), ]



# MR
mrResults <- mr(dat_clump)
mrScatterPlots <- mr_scatter_plot(mrResults, dat_clump)
mrPleio <- mr_pleiotropy_test(dat_clump)
mrHetero <- mr_heterogeneity(dat_clump)

res_single <- mr_singlesnp(dat_clump)
resForest <- forest_plot(mrResults)

# Direction
directiontest <- directionality_test(dat_clump)




### Perform SMR with MendelianRandomization package
mrIn <- dat_clump[, c("SNP", "beta.exposure", "se.exposure", "beta.outcome", "se.outcome")]
colnames(mrIn) <- c("SNP", "qtlBeta", "qtlSe", "gwasBeta", "gwasSe")
MRInputObject <- mr_input(snps = mrIn$SNP, bx = mrIn$qtlBeta, bxse = mrIn$qtlSe,
                          by = mrIn$gwasBeta, byse = mrIn$gwasSe)
MRAllObject_all <- mr_allmethods(MRInputObject, method = "all")
# mr_plot(MRAllObject_all)
mrForest <- mr_forest(MRInputObject)
mrFunnel <- mr_funnel(MRInputObject)
mr_loo(MRInputObject)


mr_presso_result <- run_mr_presso(dat_clump,          
                                  NbDistribution = 10000,          
                                  SignifThreshold = 0.05)

mr_presso_result

mr_report(dat_clump,
          output_path = "mrResults",
          author = "Charles Solomon",
          study = "")


