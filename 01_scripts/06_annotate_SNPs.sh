#!/bin/bash

# Find intersection between SVs and genes
# /bin/sh 01_scripts/12_annotate_python_quantiles.sh &

# VARIABLES
GROUP=pop #the subgroup on whcih we are making the fst comparison -> it should be a file like GROUP.txt in the folder 02_info
POP_FILE1="02_infos/pop.txt" 

FST_DIR="09_fst"

GENOME="03_genome/genome.fasta"
ANGSD_DIR="08_angsd"

#FILT_VCF_FILE="$ANGSD_DIR/ALLDP1_MISS50_2all_maf0.05.vcf.gz"
FILT_VCF_FILE="$ANGSD_DIR/merged_NS30_AD1_2all.vcf.gz.maf0.05.vcf"
FILT_POP_DIR="$FST_DIR/bypop"

ID_SEX_POP="02_infos/ID_Sex_Pop_updated.txt"

NB_CPU=1

#VCF_DIR="11_filtered"


VCF_FST="$FST_DIR/$GROUP/"$(basename -s .vcf $FILT_VCF_FILE)".fst.vcf" # path to annotated VCF (the one not formatted for angsd) produced by 11_fst_by_group.sh 

DIST=1000 # the size of the window in which intersection will be done

ANNOT_DIR="10_annot"

FIXED_FST=0.5 # fixed Fst value to be considered

GENOME_ANNOT="03_genome/genome_annotation_table_simplified.tsv"


# LOAD REQUIRED MODULES
module load bedtools/2.30.0

mkdir $ANNOT_DIR/ALL_SNPs
if [[ ! -d $ANNOT_DIR/ALL_SNPs ]]
then
  mkdir $ANNOT_DIR/ALL_SNPs
fi
mkdir $ANNOT_DIR/outlier_SNPs_fixed"$FIXED_FST"
if [[ ! -d $ANNOT_DIR/outlier_SNPs_fixed"$FIXED_FST" ]]
then
  mkdir $ANNOT_DIR/outlier_SNPs_fixed"$FIXED_FST"
fi

# 1. Extract ALL genotyped and filtered SNPs
bcftools query -f '%CHROM\t%POS\t%END\t%ID\t%FST_RO_PU\n' $VCF_FST > $FST_DIR/$GROUP/"$(basename -s .vcf $VCF_FST)".bed
## Remove entries where no Fst have been calculated
less $FST_DIR/$GROUP/"$(basename -s .vcf $VCF_FST)".bed | awk '{ if ($5 != ".") { print } }' > $ANNOT_DIR/ALL_SNPs/all_fst_SNPs_ALL.bed

# Outlier SnP genes : will be used for "Eric" enrichment, i.e. outlier Fst SNP genes vs all known genes in the annotated genome
# 2. Extract SVs with Fst values over fixed threshold
less $ANNOT_DIR/ALL_SNPs/all_fst_SNPs_ALL.bed | awk -v val="$FIXED_FST" 'BEGIN{FS="\t"} $5 > val {print}' > $ANNOT_DIR/outlier_SNPs_fixed"$FIXED_FST"/outlier_SNPs_fixed"$FIXED_FST".bed 

# 3. Find overlap between these outlier SV genes and ALL known genes using Eric's script => for "Eric" enrichment
python 01_scripts/find_genes_overlapping_snps.py $GENOME_ANNOT $ANNOT_DIR/outlier_SNPs_fixed"$FIXED_FST"/outlier_SNPs_fixed"$FIXED_FST".bed $DIST $ANNOT_DIR/outlier_SNPs_fixed"$FIXED_FST"/outlier_SNPs_fixed"$FIXED_FST"_"$DIST".overlap 

cut -f1 $ANNOT_DIR/outlier_SNPs_fixed"$FIXED_FST"/outlier_SNPs_fixed"$FIXED_FST"_"$DIST".overlap | uniq > $ANNOT_DIR/outlier_SNPs_fixed"$FIXED_FST"/outlier_SNPs_fixed"$FIXED_FST"_"$DIST".overlap.IDs

# ALL SV genes, regardless of Fst : will be used for "Claire" enrichment, i.e. outlier Fst SV genes vs all SV genes, regardless of Fst
# 1. Find overlap between all SV genes and known genes in the genome
python 01_scripts/find_genes_overlapping_snps.py $GENOME_ANNOT $ANNOT_DIR/ALL_SNPs/all_fst_SNPs_ALL.bed $DIST $ANNOT_DIR/ALL_SNPs/all_fst_SNPs_ALL_"$DIST".overlap

cut -f1 $ANNOT_DIR/ALL_SVs/all_fst_SNPs_ALL_"$DIST".overlap | uniq > $ANNOT_DIR/ALL_SNPs/all_fst_SNPs_ALL_"$DIST".overlap.IDs

# All known genes IDs, SV or not : for "Eric" enrichment   ##SOME IDS ARE MISSING, USE all_genes.ids INSTEAD##
#less $GENOME_ANNOT | tail -n +2 | cut -f5 > 03_genome/"$(basename -s .tsv $GENOME_ANNOT)".all.IDs