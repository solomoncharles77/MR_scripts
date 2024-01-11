library(openxlsx)

# arrange GTEx

# Organize and export gtexEQTL and CVD eCAVIAR ----------------------------
resDir <- "/scratch/vasccell/cs806/monocyteColoc/smrResults/resultTables/"
ll <- list.files(resDir, "csv")


# Export all SMR HEIDI Results
llall <- ll[grepl("_SMR_HEIDI.csv", ll)]
aa <- lapply(paste0(resDir, llall), read.csv)
llall <- sub("_SMR_HEIDI.csv", "", llall)
llall <- sub("20", "", llall)
llall <- sub("_CAD", "", llall)
llall <- sub("vanderHarst", "Harst", llall)
names(aa) <- llall
write.xlsx(aa, file = paste0(resDir, "all_SMR_Fairfax2014_Blueprint_Monocyte_Aragam2022CAD_vanderHarst2017CAD_Results.xlsx"))

# Export significant SMR results
smr <- ll[grepl("_SMR_HEIDI_sig_p_SMR.csv", ll)]
ss <- lapply(paste0(resDir, smr), read.csv)
names(ss) <- llall
write.xlsx(ss, file = paste0(resDir, "sig_SMR_Fairfax2014_Blueprint_Monocyte_Aragam2022CAD_vanderHarst2017CAD_Results.xlsx"))

# Export significant HEIDI results
heidi <- ll[grepl("_SMR_HEIDI_sig_p_HEIDI.csv", ll)]
hh <- lapply(paste0(resDir, heidi), read.csv)
names(hh) <- llall
write.xlsx(hh, file = paste0(resDir, "sig_HEIDI_Fairfax2014_Blueprint_Monocyte_Aragam2022CAD_vanderHarst2017CAD_Results.xlsx"))
