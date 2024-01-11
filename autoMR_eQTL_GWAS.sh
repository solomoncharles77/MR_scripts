#!/bin/bash

#PBS -N	MR_SCAD
#PBS -o MR_SCAD.output
#PBS -e MR_SCAD.error
#PBS -l walltime=95:00:00
#PBS -l vmem=80gb
#PBS -m bea
#PBS -M cs806@leicester.ac.uk
#PBS -l nodes=1:ppn=28


cd /lustre/alice3/scratch/vasccell/cs806/SCAD
module load R/4.3.1

Rscript svep1MR_scripts/autoMR_eQTL_GWAS.R -e ../pubAvail_QTL/GTEx/eQTL_catalogue/QTD000131_GTExartery_aorta_pval005.tsv.gz -f ../colocalization/cleanGWAS_Summary_Stats/GWAS_Systolic_Blood_Pressure_Evangelou_2018_Nature_hg38.txt -t GTEX_Aorta -u SBP_Evangelou2018 -c ENSG00000165124 -s 387
Rscript svep1MR_scripts/autoMR_eQTL_GWAS.R -e ../pubAvail_QTL/GTEx/eQTL_catalogue/QTD000136_GTEx_artery_coronary_pval005.tsv.gz -f ../colocalization/cleanGWAS_Summary_Stats/GWAS_Systolic_Blood_Pressure_Evangelou_2018_Nature_hg38.txt -t GTEX_CoronArt -u SBP_Evangelou2018 -c ENSG00000165124 -s 213

Rscript svep1MR_scripts/autoMR_eQTL_GWAS.R -e ../pubAvail_QTL/GTEx/eQTL_catalogue/QTD000131_GTExartery_aorta_pval005.tsv.gz -f ../colocalization/cleanGWAS_Summary_Stats/GWAS_Diastolic_Blood_Pressure_Evangelou_2018_Nature_hg38.txt -t GTEX_Aorta -u DBP_Evangelou2018 -c ENSG00000165124 -s 387
Rscript svep1MR_scripts/autoMR_eQTL_GWAS.R -e ../pubAvail_QTL/GTEx/eQTL_catalogue/QTD000136_GTEx_artery_coronary_pval005.tsv.gz -f ../colocalization/cleanGWAS_Summary_Stats/GWAS_Diastolic_Blood_Pressure_Evangelou_2018_Nature_hg38.txt -t GTEX_CoronArt -u DBP_Evangelou2018 -c ENSG00000165124 -s 213

Rscript svep1MR_scripts/autoMR_eQTL_GWAS.R -e ../pubAvail_QTL/GTEx/eQTL_catalogue/QTD000131_GTExartery_aorta_pval005.tsv.gz -f ../colocalization/cleanGWAS_Summary_Stats/GWAS_Pulse_Pressure_Evangelou_2018_Nature_hg38.txt -t GTEX_Aorta -u SBP_Evangelou2018 -c ENSG00000165124 -s 387
Rscript svep1MR_scripts/autoMR_eQTL_GWAS.R -e ../pubAvail_QTL/GTEx/eQTL_catalogue/QTD000136_GTEx_artery_coronary_pval005.tsv.gz -f ../colocalization/cleanGWAS_Summary_Stats/GWAS_Pulse_Pressure_Evangelou_2018_Nature_hg38.txt -t GTEX_CoronArt -u SBP_Evangelou2018 -c ENSG00000165124 -s 213


Rscript svep1MR_scripts/autoMR_eQTL_GWAS.R -e ../colocalization/eQTLData/imputeData_cisEQTL_2pc_10pf.txt -f ../colocalization/cleanGWAS_Summary_Stats/GWAS_Systolic_Blood_Pressure_Evangelou_2018_Nature_hg38.txt -t VSMC -u SBP_Evangelou2018 -c ENSG00000165124 -s 1499
