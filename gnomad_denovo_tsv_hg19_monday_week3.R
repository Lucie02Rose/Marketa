library(data.table)
library(GenomicRanges)
library(VariantAnnotation)
library(Rsamtools)
library(Biostrings)
library(readr)
library(dplyr)
library(stringr)
library(purrr)
library(GenomeInfoDb)


###############################################################################
# Inputs
###############################################################################
setwd("//wsl.localhost/Ubuntu/home/lucinka_rose/Marketa_group/science_paper")

tsv_file <- "gnomad.exomes.v4.1.de_novo.high_quality_coding.tsv"
vcf_file <- "gnomad_denovo.vcf"
fasta_file <- "hg19.fa"
rep_bed_file <- "per_base_territories_20kb_line_numbers.bed"

##############################################################################

df <- read_tsv(tsv_file)

vcf_df <- df %>%
  mutate(
    # split locus into CHROM and POS
    CHROM = str_split(locus, ":", simplify = TRUE)[,1],
    POS   = str_split(locus, ":", simplify = TRUE)[,2],
    POS = as.integer(POS),
    # clean and split alleles
    alleles_clean = str_remove_all(alleles, "\\[|\\]|\""),
    alleles_list  = str_split(alleles_clean, ","),
    
    REF = map_chr(alleles_list, 1),
    ALT = map_chr(alleles_list, ~ paste(.x[-1], collapse = ","))
  ) %>%
  dplyr::select(CHROM, POS, REF, ALT)

# Step 1: prepare your data
vcf_df <- vcf_df %>%
  mutate(
    ID = ".",
    QUAL = ".",
    FILTER = "PASS",
    INFO = "."
  ) %>%
  dplyr::select(CHROM, POS, ID, REF, ALT, QUAL, FILTER, INFO)

# Step 2: write VCF header
writeLines(
  c(
    "##fileformat=VCFv4.2",
    "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO"
  ),
  "gnomad_denovo.vcf"
)

# Step 3: append data (NO column names!)
write.table(
  vcf_df,
  file = "gnomad_denovo.vcf",
  sep = "\t",
  row.names = FALSE,
  quote = FALSE,
  col.names = FALSE,   # <-- THIS is the key fix
  append = TRUE
)

###############################################################################
# Load replication-direction BED
# Expected columns: chr, start, end, isLeft, isRight, timing
# BED is 0-based half-open; GRanges is 1-based closed.
###############################################################################

bed_dt <- fread(rep_bed_file)

bed_dt[, win_id := .I]

win_gr <- GRanges(
  seqnames = bed_dt$chr,
  ranges   = IRanges(start = bed_dt$start + 1L, end = bed_dt$end)
)

mcols(win_gr)$win_id  <- bed_dt$win_id
mcols(win_gr)$isLeft  <- bed_dt$isLeft
mcols(win_gr)$isRight <- bed_dt$isRight
mcols(win_gr)$timing  <- bed_dt$timing

###############################################################################
# Load VCF and expand to one row per ALT allele
###############################################################################
vcf <- readVcf(vcf_file, genome = "hg19")

rr <- rowRanges(vcf)

mut_dt <- tibble(
  chr = as.character(seqnames(rr)),
  pos = start(rr),
  ref = as.character(ref(vcf)),
  alt = as.character(unlist(alt(vcf)))
)

# Keep only simple SNVs
mut_dt <- mut_dt %>%
  filter(
    nchar(ref) == 1,
    nchar(alt) == 1,
    ref %in% c("A", "C", "G", "T"),
    alt %in% c("A", "C", "G", "T"),
    ref != alt
  )

###############################################################################
# Extract trinucleotide context from hg19 FASTA
###############################################################################

fa <- FaFile(fasta_file)
open(fa)

# Need one base on each side of the SNV
tri_gr <- GRanges(
  seqnames = mut_dt$chr,
  ranges   = IRanges(start = mut_dt$pos - 1L, end = mut_dt$pos + 1L)
)

seqlevelsStyle(tri_gr) <- "RefSeq"
seqlevels(tri_gr)[1:5]
seqnames(fa)

tri_context <- toupper(as.character(getSeq(fa, tri_gr)))

mut_dt[, trinuc := tri_context]
mut_dt <- mut_dt[
  nchar(trinuc) == 3 &
    !grepl("N", trinuc)
]

# Optional sanity check: center base should equal REF
mut_dt <- mut_dt[substr(trinuc, 2, 2) == ref]

# SBS192 channel in + reference coordinates, e.g. TCG_T, CGA_A
mut_dt[, mut_type := paste0(trinuc, "_", alt)]

###############################################################################
# Assign mutations to replication windows
###############################################################################

mut_gr <- GRanges(
  seqnames = mut_dt$chr,
  ranges   = IRanges(start = mut_dt$pos, end = mut_dt$pos)
)

hits <- findOverlaps(mut_gr, win_gr, ignore.strand = TRUE)

# If BED windows are non-overlapping, each mutation should hit at most one window.
# If not, you should clean the BED first.
mut_dt[, mut_id := .I]

hit_dt <- data.table(
  mut_id = queryHits(hits),
  win_id = mcols(win_gr)$win_id[subjectHits(hits)]
)

mut_dt <- merge(mut_dt, hit_dt, by = "mut_id")

###############################################################################
# Count mutations per window × SBS192 type
###############################################################################

count_dt <- mut_dt[, .N, by = .(win_id, mut_type)]

count_wide <- dcast(
  count_dt,
  win_id ~ mut_type,
  value.var = "N",
  fill = 0
)

###############################################################################
# Count trinucleotide opportunities per replication window
###############################################################################

# Extend each window by 1 bp on each side.
# This gives trinucleotide contexts centered on bases inside the original window.
win_ext_gr <- GRanges(
  seqnames = seqnames(win_gr),
  ranges = IRanges(start = start(win_gr) - 1L, end = end(win_gr) + 1L)
)

# You may need to trim if windows touch chromosome boundaries
# Requires seqinfo from FASTA if available
seqinfo(win_ext_gr) <- seqinfo(fa)[seqlevels(win_ext_gr)]
win_ext_gr <- trim(win_ext_gr)

seqs <- getSeq(fa, win_ext_gr)

# Counts all overlapping 3-mers in each extended sequence
opp_mat <- oligonucleotideFrequency(
  seqs,
  width = 3,
  step = 1,
  as.array = FALSE
)

opp_dt <- as.data.table(opp_mat)
opp_dt[, win_id := bed_dt$win_id]

###############################################################################
# Construct all possible SBS192 mutation IDs
###############################################################################

contexts <- colnames(opp_mat)

mut_map <- CJ(
  context = contexts,
  alt = c("A", "C", "G", "T")
)

mut_map <- mut_map[substr(context, 2, 2) != alt]
mut_map[, mut_type := paste0(context, "_", alt)]

###############################################################################
# Merge counts + opportunities and compute frequencies
###############################################################################

# Start with bed table
ann <- copy(bed_dt)

# Merge mutation counts
ann <- merge(ann, count_wide, by = "win_id", all.x = TRUE)

# Fill missing mutation counts with 0
count_cols <- setdiff(names(count_wide), "win_id")
for (cc in count_cols) {
  set(ann, which(is.na(ann[[cc]])), cc, 0)
}

# Add opportunity columns with suffix _opp
opp_long <- melt(
  opp_dt,
  id.vars = "win_id",
  variable.name = "context",
  value.name = "opportunity"
)

# Convert opportunities to wide with context_opp names
opp_wide <- dcast(
  opp_long,
  win_id ~ context,
  value.var = "opportunity",
  fill = 0
)

setnames(
  opp_wide,
  old = setdiff(names(opp_wide), "win_id"),
  new = paste0(setdiff(names(opp_wide), "win_id"), "_opp")
)

ann <- merge(ann, opp_wide, by = "win_id", all.x = TRUE)

# Make frequency columns for every SBS192 type
for (i in seq_len(nrow(mut_map))) {
  
  mt  <- mut_map$mut_type[i]
  ctx <- mut_map$context[i]
  opp_col <- paste0(ctx, "_opp")
  
  if (!mt %in% names(ann)) {
    ann[, (mt) := 0]
  }
  
  freq_col <- paste0(mt, "_freq")
  
  ann[, (freq_col) := fifelse(
    get(opp_col) > 0,
    get(mt) / get(opp_col),
    NA_real_
  )]
}