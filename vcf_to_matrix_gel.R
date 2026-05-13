##################################################
# VCF to matrix form
##################################################
## plan
# POLE/POLD1 VCFs
# SigProfilerMatrixGenerator
# Windowed SBS192 counts matrix (10–20kb bins, same resolution as BED)
# normalise to frequencies per window
# Overlap with replication direction BED (existing GRanges code)
#
#Classify windows: leftward / rightward / ambiguous
#
#For each SBS192 channel pair (pyrimidine / RC):
#  log2FC(leftward counts / rightward counts)  ← asymmetry index
#
# Compare to TOPMebaseline:
#  Is the asymmetry stronger in POLE? Which channels?
#  Are the POLE-specific channels (SBS10a/10b contexts) the most asymmetric?
#
# vcf <- readVcf("patient1.vcf.gz", genome = "hg38")
# Keep SNVs only
# snvs <- vcf[isSNV(vcf)]
# Keep PASS only (if FILTER field is populated)
# snvs <- snvs[fixed(snvs)$FILTER == "PASS"]
#
# library(MutationalPatterns)
# library(BSgenome.Hsapiens.UCSC.hg38)  # or hg19
# ref_genome <- BSgenome.Hsapiens.UCSC.hg38
# Read all VCFs at once — one column per patient
# vcf_files <- list.files("input_vcfs/", pattern = "\\.vcf(\\.gz)?$", full.names = TRUE)
# sample_names <- sub("\\.vcf.*", "", basename(vcf_files))
# grl <- read_vcfs_as_granges(vcf_files, sample_names, ref_genome)
# SBS96 matrix
# mat96 <- mut_matrix(grl, ref_genome)
# SBS192 — strand-aware, this is the key one for you
# "transcript" uses gene strand; "replication" uses rep strand if you supply it
# mat192 <- mut_matrix_stranded(grl, ref_genome, regions = NULL, mode = "replication")  # or "transcription"
#
# Get per-mutation GRanges with channel annotation from MutationalPatterns
# grl is already a GRangesList — one GRanges per patient
# Each range has metadata columns including the SBS channel
# Combine all patients, tag with sample ID
#all_muts <- unlist(grl)
#all_muts$sample <- rep(sample_names, lengths(grl))
# Your windows (same as mat_dt in your current code)
#windows_gr <- GRanges(
#  seqnames = mat_dt$chr,
#  ranges   = IRanges(mat_dt$start + 1L, mat_dt$end),
#  mat_id   = mat_dt$mat_id
#)
# Intersect mutations with windows
#hits <- findOverlaps(all_muts, windows_gr, ignore.strand = TRUE)
# Build count table: window × SBS192 channel
#mut_dt <- data.table(
#  mat_id  = mcols(windows_gr)$mat_id[subjectHits(hits)],
#  channel = all_muts$channel[queryHits(hits)]   # SBS192 channel label
#)
# Pivot to wide: rows = windows, columns = channels (like your mat_dt)
#windowed_counts <- dcast(mut_dt, mat_id ~ channel, 
#                         fun.aggregate = length, value.var = "channel", fill = 0L)


library(GenomicRanges)
library(data.table)
library(VariantAnnotation)
library(MutationalPatterns)
library(BSgenome.Hsapiens.UCSC.hg38)

ref <- BSgenome.Hsapiens.UCSC.hg38

# --- Load and combine all patient VCFs ---
vcf_files   <- list.files("vcfs/", pattern = "\\.vcf(\\.gz)?$", full.names = TRUE)
sample_names <- sub("\\.vcf.*", "", basename(vcf_files))

grl <- read_vcfs_as_granges(vcf_files, sample_names, ref)

# Add sample ID to each mutation before unlisting
for (i in seq_along(grl)) mcols(grl[[i]])$sample <- sample_names[i]
all_muts_gr <- unlist(grl, use.names = FALSE)

# --- Load replication direction BED (your existing bed_clean after GRanges cleaning) ---
# bed_clean already has isLeft, isRight columns from your code
rep_gr <- GRanges(
  seqnames = bed_clean$chr,
  ranges   = IRanges(bed_clean$start + 1L, bed_clean$end),
  isLeft   = bed_clean$isLeft,
  isRight  = bed_clean$isRight,
  timing   = bed_clean$timing
)

# --- Intersect each mutation with replication direction ---
hits <- findOverlaps(all_muts_gr, rep_gr, ignore.strand = TRUE)

# Keep only mutations that fall in unambiguous windows
mut_rep <- data.table(
  mut_idx = queryHits(hits),
  isLeft  = mcols(rep_gr)$isLeft[subjectHits(hits)],
  isRight = mcols(rep_gr)$isRight[subjectHits(hits)],
  timing  = mcols(rep_gr)$timing[subjectHits(hits)]
)

# If a mutation overlaps multiple windows, take majority direction
# (rare at mutation level but possible at window boundaries)
mut_rep <- mut_rep[, .(
  isLeft  = as.integer(mean(isLeft)  > 0.5),
  isRight = as.integer(mean(isRight) > 0.5),
  timing  = mean(timing)
), by = mut_idx]

# Assign direction label
mut_rep[, rep_class := fcase(
  isLeft == 1L & isRight == 0L, "leftward",
  isRight == 1L & isLeft == 0L, "rightward",
  default = "ambiguous"
)]

# Attach back to mutation GRanges metadata
mcols(all_muts_gr)$rep_class <- NA_character_
mcols(all_muts_gr)$rep_class[mut_rep$mut_idx] <- mut_rep$rep_class
mcols(all_muts_gr)$timing[mut_rep$mut_idx]    <- mut_rep$timing


# --- Get trinucleotide context per mutation ---
# MutationalPatterns attaches this when building the mut_matrix
# but we can also get it per-mutation:
mut_context <- get_mut_context(all_muts_gr, ref)  
# returns e.g. "T[C>A]G" per mutation

mcols(all_muts_gr)$channel <- mut_context

# --- Build callable bases per patient per rep-class region ---
# You need mosdepth or GATK callable BEDs for this
# If you have them:
callable_files <- list.files("callable/", pattern = "\\.bed$", full.names = TRUE)

callable_list <- lapply(callable_files, function(f) {
  dt <- fread(f, col.names = c("chr", "start", "end"))
  GRanges(dt$chr, IRanges(dt$start + 1L, dt$end))
})
names(callable_list) <- sample_names

# Total callable bases per patient intersected with each rep class
get_callable_bp <- function(callable_gr, rep_gr, rep_class_label) {
  class_gr <- rep_gr[mcols(rep_gr)$isLeft == (rep_class_label == "leftward") &
                       mcols(rep_gr)$isRight == (rep_class_label == "rightward")]
  sum(width(intersect(callable_gr, class_gr)))
}

callable_bp <- rbindlist(lapply(sample_names, function(s) {
  data.table(
    sample    = s,
    leftward  = get_callable_bp(callable_list[[s]], rep_gr, "leftward"),
    rightward = get_callable_bp(callable_list[[s]], rep_gr, "rightward")
  )
}))

total_callable <- data.table(
  rep_class = c("leftward", "rightward"),
  callable_bp = c(sum(callable_bp$leftward), sum(callable_bp$rightward))
)

# --- Count reference trinucleotide opportunities in left vs right regions ---
left_gr  <- rep_gr[mcols(rep_gr)$isLeft == 1L & mcols(rep_gr)$isRight == 0L]
right_gr <- rep_gr[mcols(rep_gr)$isRight == 1L & mcols(rep_gr)$isLeft == 0L]

trinuc_left  <- count_trinucleotides(left_gr,  ref)  # 64-element named vector
trinuc_right <- count_trinucleotides(right_gr, ref)


# --- Build mutation count table ---
mut_dt <- as.data.table(mcols(all_muts_gr))[
  rep_class %in% c("leftward", "rightward")
]

# Raw counts: SBS192 channel x rep_class
count_mat <- dcast(mut_dt, channel ~ rep_class, fun.aggregate = length, fill = 0L)

# --- Normalise by trinucleotide opportunity ---
# Each SBS192 channel has a reference trinucleotide (the central 3 bases before mutation)
# Map channel -> trinucleotide
count_mat[, trinuc := substr(channel, 1, 3)]  # adjust parsing to your channel format

# Merge with opportunity counts
opp_dt <- data.table(
  trinuc    = names(trinuc_left),
  opp_left  = as.numeric(trinuc_left),
  opp_right = as.numeric(trinuc_right)
)
count_mat <- merge(count_mat, opp_dt, by = "trinuc")

# Also merge with callable bp for depth correction
count_mat[, opp_left_adj  := opp_left  * (total_callable[rep_class=="leftward",  callable_bp] / 1e9)]
count_mat[, opp_right_adj := opp_right * (total_callable[rep_class=="rightward", callable_bp] / 1e9)]

# Final normalised rate: mutations per trinucleotide opportunity per callable base
count_mat[, rate_left  := leftward  / opp_left_adj]
count_mat[, rate_right := rightward / opp_right_adj]

# Asymmetry index (same as your existing log2fc metric)
eps <- 1e-6
count_mat[, log2fc_asym := log2((rate_left + eps) / (rate_right + eps))]

