#########################
### SBS EXPLORER IN R ###
#########################

## RQ1 ###
# What is the strand directionality and CpG context of de-novo mutations? 
# I am using data from https://escholarship.org/uc/item/2b83c8gw
# very rare SNVs from TOPMed freeze 5, with allele frequency below 10^-4, as a proxy for germline mutations. 
# They then validated the inferred components using de novo mutations from trio studies.

### Notes to methods ###
# there are two possible starting points, either from vcf files of the matrices
# count each mutation directly into SBS96 or SBS192 bins (convert to frequencies)
# assign each mutation to leading or lagging based on BED overlap
# run signature extraction separately or jointly by strand class

#Minimum useful outputs
#A QC-annotated window or mutation table
#coordinates
#leading/lagging assignment
#overlap fraction
#ambiguous flag
#number of callable sites or context counts
#Leading vs lagging mutation spectrum
#SBS96 or SBS192 counts/rates
#enrichment or log2 fold-change
#plots
#Statistical comparison
#per-context enrichment tests
#overall cosine similarity / correlation between spectra
#Optional higher-level outputs
#PCA+ICA component intensities per window
#Correlation of component asymmetry with replication direction
#Signature extraction within leading and lagging subsets

#6. Variant filtering

#If starting from VCFs, decide:

#  SNVs only?
#  exclude indels?
#  exclude multiallelics or split them?
#  PASS only?
#  allele frequency threshold?
#  singleton only?
#  rare variants only?
#  
#  If you have VCFs
#Filter to high-quality rare SNVs
#Annotate each SNV with trinucleotide context
#Overlap each SNV with leading/lagging BED
#Count SBS96 or SBS192 separately for leading and lagging
#Normalize by callable context opportunities
#Compare spectra
#Optionally run signature extraction
#If you have TOPMed / gnomAD matrices
#Read matrix
#Read BED
#Assign each window a leading/lagging label or direction score
#Remove ambiguous / poorly covered windows
#Compare context-dependent rates between classes
#Run PCA+ICA on windows
#Correlate component intensities or paired asymmetry with direction score

#Then, if you want a vector “from the minus-strand point of view,” you use the reverse-complement-normalized representation consistently.

#But for most SBS analyses, that is unnecessary if you already normalized everything into pyrimidine-centered notation and stratified by replication class.

#The exact checks you need

#For this CpG-focused analysis, I would check:

#  Genome build matches between VCF, BED, and reference genome.
#Only SNVs are included.
#Single ALT only or properly split multiallelics.
#Reference allele matches reference genome at that position.
#No ambiguous BED overlaps unless you define a rule.
#Fork direction to plus leading/lagging mapping is explicitly documented.
#Counts and opportunities are both computed for C[C>G]G.
#If using precomputed matrices instead of VCFs, verify whether they are reference-oriented 192-context rates, because then the bookkeeping differs slightly from SBS96 counting.

#Mode 1: VCF input

#read SNVs
#extract trinucleotide context from reference genome
#build SBS192 and SBS96 labels
#assign each variant to a genomic window
#overlap with BED to get replication annotation

#Mode 2: matrix input

#read precomputed per-window SBS96 or SBS192 matrix
#overlap windows with BED
#use matrix directly

#So your code should always standardize inputs to a common structure like:

#  meta: chr, start, end, feature_label, direction_score
#mat192: windows x 192
#mat96: windows x 96


rm(list = ls(all = TRUE))

setwd("//wsl.localhost/Ubuntu/home/lucinka_rose/Marketa_group/science_paper")

###############################################################################
# PACKAGES
###############################################################################

pkgs <- c("BiocManager", "tidyverse", "dplyr", "tidyr",
          "ggplot2", "data.table", "zoo")
for (p in pkgs) {
  if (!require(p, character.only = TRUE, quietly = TRUE))
    install.packages(p)
  library(p, character.only = TRUE)
}

bioc_pkgs <- c("GenomicRanges", "IRanges", 
               "BSgenome.Hsapiens.UCSC.hg38", "VariantAnnotation")
for (p in bioc_pkgs) {
  if (!require(p, character.only = TRUE, quietly = TRUE))
    BiocManager::install(p)
  library(p, character.only = TRUE)
}

###############################################################################
# SETTINGS  -- only edit this block
###############################################################################

matrix_file      <- "gnomAD_100kb.txt"
bed_file         <- "final.hg38.bed"
output_dir       <- "strand_specific_output_gnomad"

# A TOPMed window is assigned to leading/lagging only if at least this
# fraction of its length is covered by a single-direction BED interval.
# Windows below this threshold are labelled "ambiguous" and excluded.
assign_threshold <- 0.5

# Pseudocount added before log2 to avoid log(0).
# 1e-6 is appropriate for mutation rates (small fractions ~0.001-0.01).
# Do NOT use 1e-10 or 1e-12 -- those amplify noise in near-zero channels.
PSEUDOCOUNT <- 1e-6

dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

#Even if you start from VCF, you eventually aggregate into that same structure.