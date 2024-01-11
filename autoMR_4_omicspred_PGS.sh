#!/bin/bash

#PBS -N	MR_SVEP1
#PBS -o MR_SVEP1.output
#PBS -e MR_SVEP1.error
#PBS -l walltime=95:00:00
#PBS -l vmem=80gb
#PBS -m bea
#PBS -M cs806@leicester.ac.uk
#PBS -l nodes=1:ppn=28


cd /lustre/alice3/scratch/vasccell/cs806/svep1MR
module load R/4.2.1

Rscript svep1MR_scripts/autoMR_4_omicspred_PGS.R -e ../omicspred/mrFiles/SVEP1_OPGS000334_QTL_all.csv -f ../colocalization/cleanGWAS_Summary_Stats/GWAS_Systolic_Blood_Pressure_Evangelou_2018_Nature_hg38.txt -t SVEP1a -u SBP_Evangelou2018 
Rscript svep1MR_scripts/autoMR_4_omicspred_PGS.R -e ../omicspred/mrFiles/SVEP1_OPGS000426_QTL_all.csv -f ../colocalization/cleanGWAS_Summary_Stats/GWAS_Systolic_Blood_Pressure_Evangelou_2018_Nature_hg38.txt -t SVEP1a -u SBP_Evangelou2018 

Rscript svep1MR_scripts/autoMR_4_omicspred_PGS.R -e ../omicspred/mrFiles/SVEP1_OPGS000334_QTL_all.csv -f ../colocalization/cleanGWAS_Summary_Stats/GWAS_Diastolic_Blood_Pressure_Evangelou_2018_Nature_hg38.txt -t SVEP1a -u DBP_Evangelou2018 
Rscript svep1MR_scripts/autoMR_4_omicspred_PGS.R -e ../omicspred/mrFiles/SVEP1_OPGS000426_QTL_all.csv -f ../colocalization/cleanGWAS_Summary_Stats/GWAS_Diatolic_Blood_Pressure_Evangelou_2018_Nature_hg38.txt -t SVEP1a -u DBP_Evangelou2018 

Rscript svep1MR_scripts/autoMR_4_omicspred_PGS.R -e ../omicspred/mrFiles/SVEP1_OPGS000334_QTL_all.csv -f ../colocalization/cleanGWAS_Summary_Stats/GWAS_Pulse_Pressure_Evangelou_2018_Nature_hg38.txt -t SVEP1a -u PP_Evangelou2018 
Rscript svep1MR_scripts/autoMR_4_omicspred_PGS.R -e ../omicspred/mrFiles/SVEP1_OPGS000426_QTL_all.csv -f ../colocalization/cleanGWAS_Summary_Stats/GWAS_Pulse_Pressure_Evangelou_2018_Nature_hg38.txt -t SVEP1a -u PP_Evangelou2018 
