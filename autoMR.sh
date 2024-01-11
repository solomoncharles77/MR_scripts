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
module load R/4.2.1

Rscript svep1MR_scripts/autoMR.R -e ../colocalization/cleanGWAS_Summary_Stats/GWAS_SVEP1_OPGS000334_INTERVAL_Omicspred_hg38.txt_hg38.txt -f ../colocalization/cleanGWAS_Summary_Stats/GWAS_Systolic_Blood_Pressure_Evangelou_2018_Nature_hg38.txt -t SVEP1_OPGS000334 -u SBP_Evangelou2018 
Rscript svep1MR_scripts/autoMR.R -e ../colocalization/cleanGWAS_Summary_Stats/GWAS_SVEP1_OPGS000334_INTERVAL_Omicspred_hg38.txt_hg38.txt -f ../colocalization/cleanGWAS_Summary_Stats/GWAS_Diastolic_Blood_Pressure_Evangelou_2018_Nature_hg38.txt -t SVEP1_OPGS000334 -u DBP_Evangelou2018 
Rscript svep1MR_scripts/autoMR.R -e ../colocalization/cleanGWAS_Summary_Stats/GWAS_SVEP1_OPGS000334_INTERVAL_Omicspred_hg38.txt_hg38.txt -f ../colocalization/cleanGWAS_Summary_Stats/GWAS_Pulse_Pressure_Evangelou_2018_Nature_hg38.txt -t SVEP1_OPGS000334 -u PP_Evangelou2018 
