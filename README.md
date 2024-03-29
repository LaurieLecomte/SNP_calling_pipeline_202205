# SNP calling pipeline (upgraded 2022 version)

This repository features the latest version of my SNP calling pipeline used in my M.Sc on Atlantic Salmon in Louis Bernatchez’s lab.

The pipeline calls SNPs and indels in 60 Atlantic Salmon genomes from 2 rivers (La Romaine and Puyjalon), using `bcftools mpileup`.

Scripts were adapted from previous work done by Clément Rougeux, Claire Mérot and Florent Sylvestre.


## Pipeline Overview

This pipeline features 4 main steps, SNP (and indel) calling, merging and filtering. 

The SNP calling step uses `bcftools mpileup` on 60 bam files, one for each individual (see [Raw Data](https://github.com/LaurieLecomte/SNP_calling_pipeline_202106#raw-data) below). For the sake of efficiency, this step has been parallelized by chromosome, e.g. SNPs were called across **all** samples at once on **1 chromosome at the time**, which yields 1 VCF file per chromosome.

All chromosomes' VCFs were then merged together in a single VCF file, which was then filtered.


## Prerequisites

### Files 

* A **reference genome** (`.fasta`) and its **index** (`.fai`) in `03_genome`
* **Bam files** for all samples and their index. These can be soft-linked in the `04_bam` folder for easier handling : if `$BAM_PATH` is the remote path to bam files use `for file in $(ls -1 $BAM_PATH/*); do ln -s $file ./04_bam; done`
* A **bam files list** in `02_infos`. This list can be generated with the following command, where `$BAM_DIR` is the path of the directory where bam files are located : `ls -1 $BAM_DIR/*.bam > 02_infos/bam_list.txt`
* A **chromosomes list** (or contigs, or sites) in `02_infos`. This list is used for parallelizing the SNP calling step. It can be produced from the indexed genome file ("$GENOME".fai) : `less "$GENOME".fai | cut -f1 > 02_infos/chr_list.txt`
* Optional : a list of samples IDs and their population (and/or sex) for popgen analysis, such as PCA or FST calculation, in `02_infos`. 

### Software Requirements
* `bcftools` (Danecek et al. 2021, https://pubmed.ncbi.nlm.nih.gov/33590861/) and dependencies : used version 1.12
* `GNU parallel` 
* Optional : `tmux` 


## Raw Data 

Adipose fin tissue has been sampled from 60 wild-born Atlantic Salmon in 2020 by the staff of LARSA (Laboratoire de recherche en science aquatique, Université Laval). DNA extractions were performed in september 2020 using Qiagen’s DNeasy Blood and Tissue kit. 

DNA samples were sent to Génome Québec in January 2021 for sequencing using Illumina’s NovaSeq6000 platform at a coverage of about 16X. 
Sequencing data was then processed using Éric Normandeau’s [wgs_sample_preparation pipeline](https://github.com/enormandeau/wgs_sample_preparation) in May 2021. 
