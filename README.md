# Project Description

A reproduction of the Wang et al. 2014 paper: “A comprehensive study design reveals treatment- and transcript abundance–dependent concordance between RNA-seq and microarray data”. 

# Contributors

Data Curator: Rebecca Spirgel

Programmer: Giulia Digiovanni

Analyst: Marissa Chiaradio

Biologist: Ledia Gebremedhin


# Repository Contents

Data Curator:

data_download.sh - A bash script to create a list of sample IDs from toxgroup 3

run_fastqc.qsub - A qsub job containing a bash script to run fastqc on the IDs in the folder

star_align.qsub - A qsub job containing a bash script to run STAR on paired end fastq.gz files with a shared ID (provided as input)

multiqc.qsub - A qsub job containing a bash script to run multiqc on the IDs in the folder


Analyst: 

analyst_script1.R - A limma analysis of microarray data.

analyst_rscript2.R - Computing concordance values for limma/deseq results. 

  Analyst output files: Contains several scatter plots, histograms, and barplots produced by R scripts. As well as csv files from limma results. 

