#!/bin/bash

# Call SNPs by chr in 60 samples in parallel
# Run this script from SNP_calling_pipeline_202205 directory using this command :
# parallel -a ./02_infos/chr_list.txt -k -j 10 srun -c 2 --mem=10G -p medium --time=7-00:00:00 -J 01_mpileup_{} -o log/01_mpileup_{}_%j.log 01_scripts/01_mpileup.sh {} &


# VARIABLES
REF="03_genome/genome.fasta"     # reference genome
BAM_LIST="02_infos/bam_list.txt" # list of bam files
CHR=$1                           # region or chromosome for this job 
OUT_DIR="05_calls"               # results directory
CPU=2                            # number of threads to use per job

# LOAD REQUIRED MODULES
module load bcftools/1.12

# 1. Call SNPs only : required -I option
bcftools mpileup -Ou -f $REF -r $CHR -b $BAM_LIST -I -a AD,DP,SP,ADF,ADR -q 5 --threads $CPU | bcftools call -a GP,GQ -mv -Oz --threads $CPU > $OUT_DIR/"$CHR".vcf.gz 

# Call SNPs AND indels : remove -I 
#bcftools mpileup -Ou -f $REF -r $CHR -b $BAM_LIST -a AD,DP,SP,ADF,ADR -q 5 --threads $CPU | bcftools call -a GP,GQ -mv -Oz --threads $CPU > $OUT_DIR/"$CHR"_SNPs_indels.vcf.gz 