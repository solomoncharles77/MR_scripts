#!/bin/bash

#PBS -N	runSMR
#PBS -o runSMR.output
#PBS -e runSMR.error
#PBS -l walltime=03:00:00
#PBS -l vmem=80gb
#PBS -m bea
#PBS -M cs806@leicester.ac.uk
#PBS -l nodes=1:ppn=28

cd /lustre/alice3/scratch/vasccell/cs806/SCAD
module load R/4.1.0

# vsmc SCAD
/home/c/cs806/SMR/smr-1.3.1 --bfile /scratch/vasccell/cs806/colocalization/1000Genome/euroSamps1kGMerge --gwas-summary /scratch/vasccell/cs806/colocalization/zhu_SMR/GWAS_Sponteneous_Coronary_Artery_Dissection_BothSex_Adlam_2023_NatGenet_hg38_smrHEIDI.txt --beqtl-summary ../colocalization/eQTLData/imputeData_cisEQTL_2pc_10pf --out smrResults/vsmcEQTL_Adlam2023_SCAD_SMR_HEIDI
Rscript /scratch/vasccell/cs806/SCAD/SCAD_scripts/auto_process_SMR_HEIDI_Results.R smrResults/vsmcEQTL_Adlam2023_SCAD_SMR_HEIDI.smr

# vsmc CAD
/home/c/cs806/SMR/smr-1.3.1 --bfile /scratch/vasccell/cs806/colocalization/1000Genome/euroSamps1kGMerge --gwas-summary /scratch/vasccell/cs806/colocalization/zhu_SMR/GWAS_Coronary-Artery-Disease_Aragam_2022_NatGenet_hg38_smrHEIDI.txt --beqtl-summary ../colocalization/eQTLData/imputeData_cisEQTL_2pc_10pf --out smrResults/vsmcEQTL_Aragam2022_CAD_SMR_HEIDI
Rscript /scratch/vasccell/cs806/SCAD/SCAD_scripts/auto_process_SMR_HEIDI_Results.R smrResults/vsmcEQTL_Aragam2022_CAD_SMR_HEIDI.smr
