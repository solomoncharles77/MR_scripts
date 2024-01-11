library(MRInstruments)
library(TwoSampleMR)
library(data.table)
library(MendelianRandomization)


# Prep outcome data -------------------------------------------------------
gwas1 <- data.frame(fread("/scratch/vasccell/cs806/colocalization/cleanGWAS_Summary_Stats/GWAS_Sponteneous_Coronary_Artery_Dissection_Female_Adlam_2023_NatGenet_hg38.txt"))
pcutoff <- 5e-02

# Subset SNPs that meet association threshold
gwas1_sig <- gwas1[which(gwas1$pvalue < pcutoff), ]
#gwas1_sig <- gwas1
outcome_gwas <- gwas1_sig[, c("hg38_markername", "permID", "ALT", "REF", "eaf", "pvalue", "beta", "se", "n")]
colnames(outcome_gwas) <- c("hg38_markername", "SNP", "effect_allele","other_allele","eaf", "pval", "beta", "se", "samplesize")
outcome_dat <- format_data(outcome_gwas, type = "outcome")


# Get exposure data -------------------------------------------------------
data(gwas_catalog)
exposure_gwas <- subset(gwas_catalog, grepl("Menarche", Phenotype_simple))
unique(exposure_gwas$Author)
exposure_gwas <- subset(exposure_gwas, grepl("Day", Author))
exposure_gwas <- exposure_gwas[exposure_gwas$pval < 5*10^-8, ]
head(exposure_gwas[,c(7:12,18:21)])
exposure_dat <- format_data(exposure_gwas, type = "exposure")


# Harmonisation
dat <- harmonise_data(exposure_dat = exposure_dat, outcome_dat = outcome_dat)

# Clumping data
dat_clump <- clump_data(dat[which(dat$mr_keep),], clump_r2 = 0.0001)


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



ld_clump_local2 <- function (dat, clump_kb, clump_r2, clump_p, bfile, plink_bin) 
{
  shell <- ifelse(Sys.info()["sysname"] == "Windows", "cmd", 
                  "sh")
  fn <- tempfile()
  write.table(data.frame(SNP = dat[["rsid"]], P = dat[["pval"]]), 
              file = fn, row.names = F, col.names = T, quote = F)
  fun2 <- paste0(shQuote(plink_bin, type = shell), " --bfile ", 
                 shQuote(bfile, type = shell), " --clump ", shQuote(fn, 
                                                                    type = shell), " --clump-p1 ", clump_p, " --clump-r2 ", 
                 clump_r2, " --clump-kb ", clump_kb, " --out ", shQuote(fn, 
                                                                        type = shell))
  system(fun2)
  res <- read.table(paste(fn, "", sep = ""), header = T)
  unlink(paste(fn, "*", sep = ""))
  y <- subset(dat, !dat[["rsid"]] %in% res[["SNP"]])
  if (nrow(y) > 0) {
    message("Removing ", length(y[["rsid"]]), " of ", nrow(dat), 
            " variants due to LD with other variants or absence from LD reference panel")
  }
  return(subset(dat, dat[["rsid"]] %in% res[["SNP"]]))
}
