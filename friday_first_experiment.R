###############################################################################
# EXPERIMENT 1 | WEEK 1
# Leading vs Lagging strand mutational asymmetry
# Germline low-frequency VAFs 
# TOPMed SBS192 matrix, 10kb, aggregate across patients
# Gnomad SBS192 matrix, 100kb, aggregate across patients
# 
# all is expressed in SBS192 relative to the plus strand on hg38
#
# rightward replicating + strand is lagging template (isRightReplicating == 1)
# TCG_T on rightward is CGA_A reverse complemented
# TCG_T rightward = lagging template
# TCG_T leftward = leading template
# CGA_A rightward = leading template
# CGA_A leftward = lagging template
# leftward replicating + strand is the leading template (isLeftReplicating == 1)
# ambiguous = in the bed file that I got given = both are 0 or both are 
# 
# to compare leading vs lagging, split windows by replication class and compare
# channels across the two groups
#
# the - strand rates are encoded as rc partner channels, we are comparing across
# window classes (not within), composition differences average out and within 
# window rc is not needed
#
# for any genuine strand-dependent process, if channel X is enriched on
# leading, its rc partner must be enriched on lagging (same physical mutation,
# seen from the other strand)
###############################################################################

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

bioc_pkgs <- c("GenomicRanges", "IRanges")
for (p in bioc_pkgs) {
  if (!require(p, character.only = TRUE, quietly = TRUE))
    BiocManager::install(p)
  library(p, character.only = TRUE)
}

###############################################################################
# SETTINGS
###############################################################################
matrix_file      <- "TOPMed_10kb.txt"
bed_file         <- "final.hg38.bed"
output_dir       <- "strand_specific_output_topmed_hg38"
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# pseudocount added before log2 to avoid log(0)
PSEUDOCOUNT <- 1e-6

###############################################################################
# HELPER FUNCTIONS
###############################################################################

# all chromosome names "chrN" format
normalize_chr <- function(x) {
  x <- as.character(x)
  paste0("chr", sub("^chr", "", x, ignore.case = TRUE))
}

# complement a single DNA character
rev_comp_base <- function(x) chartr("ACGT", "TGCA", x)

# reverse complement a DNA string
rev_comp_seq <- function(x) {
  vapply(x, function(s) {
    chars <- strsplit(s, "", fixed = TRUE)[[1]]
    paste(rev(rev_comp_base(chars)), collapse = "")
  }, character(1))
}

# parse SBS192 column names in "ACT_G" (trinucleotide_alt) format.
# returns a data.table with three columns per label:
# orig_label = the label as it appears in the matrix (e.g. "ACT_G")
# sbs96 = the COSMIC SBS96 canonical form (pyrimidine ref: C or T)
# rc_label = the reverse-complement partner label (e.g. "AGT_C")
# used for rc QC, not matrix manipulation
parse_sbs_label <- function(labels) {
  labels <- as.character(labels)
  rbindlist(lapply(labels, function(lbl) {
    lbl <- trimws(lbl)
    # expect "XYZ_W" where Y is the mutated base as that is their sbs frq data
    if (!grepl("^[ACGT]{3}_[ACGT]$", lbl))
      stop(sprintf("Unrecognised SBS label: '%s'", lbl))
    tri     <- sub("_.*$", "", lbl)
    alt     <- sub("^.*_",  "", lbl)
    ref     <- substr(tri, 2, 2)
    tri_rc  <- rev_comp_seq(tri)
    alt_rc  <- rev_comp_base(alt)
    rc_label <- paste0(tri_rc, "_", alt_rc)
    # canonical SBS96: ref must be C or T (pyrimidine convention)
    if (ref %in% c("A", "G")) {
      sbs96 <- paste0(substr(tri_rc,1,1), "[", substr(tri_rc,2,2), ">",
                      alt_rc, "]", substr(tri_rc,3,3))
    } else {
      sbs96 <- paste0(substr(tri,1,1), "[", ref, ">", alt, "]", substr(tri,3,3))
    }
    data.table(orig_label = lbl, sbs96 = sbs96, rc_label = rc_label)
  }))
}

# build an index map: for each of the 192 columns, store the position of
# its rc partner in the column list, used for anticorrelation QC (about 0.992 for topmed file)
make_rc_map <- function(sbs_cols) {
  meta   <- parse_sbs_label(sbs_cols)
  rc_idx <- match(meta$rc_label, meta$orig_label)
  if (any(is.na(rc_idx)))
    stop("RC partner missing for some SBS columns. Check column names.")
  list(meta = meta, rc_idx = rc_idx)
}

# write a named numeric vector as a two-column tsv (channel, value)
write_vec <- function(x, file)
  fwrite(data.table(channel = names(x), value = as.numeric(x)), file, sep = "\t")

###############################################################################
# READ & STANDARDISE INPUT FILES
###############################################################################

mat_dt <- fread(matrix_file)
bed_dt <- fread(bed_file)

mat_dt[, chr        := normalize_chr(chr)]
bed_dt[, chromosome := normalize_chr(chromosome)]

mat_lengths <- mat_dt[, .(mat_max_end = max(end)), by = chr]
bed_lengths <- bed_dt[, .(chr = chromosome, bed_max_end = max(end)),
                      by = chromosome][, chromosome := NULL]

length_check <- merge(mat_lengths, bed_lengths, by = "chr", all = TRUE)
length_check[, diff := mat_max_end - bed_max_end]

# Order by absolute difference -- a build mismatch shows up as large diffs
print(length_check[order(-abs(diff))])

mat_dt[, `:=`(
  chr    = normalize_chr(chr),
  start  = start,
  end    = end,
  mat_id = .I                   
)]

# bed file sanity checks
stopifnot(all(c("chromosome", "start", "end",
                "isLeftReplicating", "isRightReplicating") %in% names(bed_dt)))

bed_dt <- bed_dt[, .(
  chr    = normalize_chr(chromosome),
  start  = start,
  end    = end,
  isLeft  = isLeftReplicating,
  isRight = isRightReplicating,
  timing = replicationTiming
)]

# report bed composition as a sanity check before any overlap
message(sprintf(
  "BED intervals -- leftward only: %d | rightward only: %d | both: %d | neither: %d",
  sum(bed_dt$isLeft == 1 & bed_dt$isRight == 0),
  sum(bed_dt$isRight == 1 & bed_dt$isLeft == 0),
  sum(bed_dt$isLeft == 1 & bed_dt$isRight == 1),
  sum(bed_dt$isLeft == 0 & bed_dt$isRight == 0)
))

###############################################################################
# if the bed file for replication direction is lifted over, there are overlapping 
# regions which have both left and right replication at the same time, so these
# must be excluded

bed_dt <- copy(bed_dt)
bed_dt[, row_id := .I]

bed_gr <- GRanges(
  bed_dt$chr,
  IRanges(bed_dt$start + 1L, bed_dt$end)
)

left_idx  <- which(bed_dt$isLeft  == 1L & bed_dt$isRight == 0L)
right_idx <- which(bed_dt$isRight == 1L & bed_dt$isLeft  == 0L)

# find any overlap between left and right intervals
hits_lr <- findOverlaps(
  bed_gr[left_idx],
  bed_gr[right_idx],
  ignore.strand = TRUE,
  minoverlap = 1L
)

# rows to remove: both sides of every conflicting pair
drop_idx <- unique(c(
  left_idx[queryHits(hits_lr)],
  right_idx[subjectHits(hits_lr)]
))

bed_clean <- bed_dt[-drop_idx]

message(sprintf(
  "BED intervals -- leftward only: %d | rightward only: %d | both: %d | neither: %d",
  sum(bed_clean$isLeft == 1 & bed_clean$isRight == 0),
  sum(bed_clean$isRight == 1 & bed_clean$isLeft == 0),
  sum(bed_clean$isLeft == 1 & bed_clean$isRight == 1),
  sum(bed_clean$isLeft == 0 & bed_clean$isRight == 0)
))

###############################################################################
# overlap of windows: assign replication direction to each topmed window
#
# bed format is 0-based half-open [start, end).
# granges is 1-based closed [start, end].
# conversion syntax: granges_start = bed_start + 1, granges_end = bed_end.
###############################################################################

mat_gr <- GRanges(mat_dt$chr, IRanges(mat_dt$start + 1L, mat_dt$end)) 
bed_gr <- GRanges(bed_clean$chr, IRanges(bed_clean$start + 1L, bed_clean$end)) 

mcols(bed_gr)$isLeft <- bed_clean$isLeft 
mcols(bed_gr)$isRight <- bed_clean$isRight 
mcols(bed_gr)$timing <- bed_clean$timing 

hits <- findOverlaps(mat_gr, bed_gr, ignore.strand = TRUE) 
if (length(hits) == 0) stop("No overlaps found. Check chromosome names and coordinate systems.")

q <- queryHits(hits)
s <- subjectHits(hits)

ann <- data.table(
  mat_id = q,
  ov_bp  = width(pintersect(ranges(mat_gr)[q], ranges(bed_gr)[s])),
  isLeft  = mcols(bed_gr)$isLeft[s],
  isRight = mcols(bed_gr)$isRight[s],
  timing  = mcols(bed_gr)$timing[s]
)[
  , .(
    bp_left  = sum(ov_bp[isLeft == 1 & isRight == 0], na.rm = TRUE),
    bp_right = sum(ov_bp[isRight == 1 & isLeft == 0], na.rm = TRUE),
    timing_w = weighted.mean(timing, w = ov_bp, na.rm = TRUE)
  ),
  by = mat_id
][
  mat_dt, on = "mat_id"
]

ann[is.na(bp_left),  bp_left  := 0]
ann[is.na(bp_right), bp_right := 0]

ann[, win_bp := end - start]
ann[, classified_frac := (bp_left + bp_right) / win_bp]
ann[, left_cov  := bp_left  / win_bp]
ann[, right_cov := bp_right / win_bp]

# we are asking only full windows to be overlapped
ann[, rep_class := fcase(
  left_cov  >= 1 & right_cov < 1, "leftward",
  right_cov >= 1 & left_cov  < 1, "rightward",
  default = "ambiguous"
)]

# report and save
message("\nWindow classification summary:")
print(ann[, .N, by = rep_class])
message(sprintf("Median classified_frac: %.3f", median(ann$classified_frac)))

fwrite(ann[, .(chr, start, end, timing_w, bp_left, bp_right,
               classified_frac, left_cov, right_cov, rep_class)],
       file.path(output_dir, "window_annotations.tsv"), sep = "\t")

a <- ann[ann$bp_right != 0 & ann$bp_left != 0]
a # should be empty
###############################################################################
# EXTRACT SBS192 MUTATION RATE MATRIX
###############################################################################

# identify which columns are SBS mutation channels (format: "ACG_T")
meta_cols <- c("mat_id", "chr", "start", "end", "bp_left", "bp_right",
               "win_bp", "classified_frac", "left_cov", "right_cov", "rep_class")

sbs_cols <- setdiff(names(ann), meta_cols)
sbs_cols <- sbs_cols[grepl("^[ACGT]{3}_[ACGT]$", sbs_cols)]

if (length(sbs_cols) == 0)
  stop("No SBS192 columns detected. Check column names in the matrix file.")
message(sprintf("SBS192 columns detected: %d (expected: 192)", length(sbs_cols)))

# Extract as a plain numeric matrix; rows = genomic windows, cols = mutation channels
mat192 <- as.matrix(ann[, ..sbs_cols])
rownames(mat192) <- paste0(ann$chr, ":", ann$start, "-", ann$end)

###############################################################################
# ROW-SUM DIAGNOSTIC
# SBS192 channels are mutation rates in this case (not counts), so they do not sum
# to 1. channels should be on a comparable scale (order 1e-3 to 1e-2).
# bimodal distribution in row sums indicates low-quality windows that
# passed the coverage filter but have unusual mutation densities
###############################################################################

row_sums <- rowSums(mat192, na.rm = TRUE)
message(sprintf("Row sum stats: min=%.4f  median=%.4f  max=%.4f",
                min(row_sums), median(row_sums), max(row_sums)))

pdf(file.path(output_dir, "diagnostic_rowsums.pdf"), width = 7, height = 5)
print(
  ggplot(data.frame(total = row_sums), aes(x = total)) +
    geom_histogram(bins = 60, fill = "steelblue", colour = "white") +
    labs(title = "Row sums of SBS192 channels per 10 kb window",
         subtitle = "These are rates, not probabilities -- no requirement to sum to 1",
         x = "Sum of all 192 mutation rates", y = "Number of windows")
)
dev.off()

###############################################################################
# leading vs lagging window indices based on the origins
#
# leftward windows  (+ strand = leading template) --> idx_leading
# rightward windows (+ strand = lagging template) --> idx_lagging
#
# windows below assign_threshold are excluded from both groups
###############################################################################
# of the ann matrix, we assign everything that is going to the right as lagging
# and everything going to the left as leading (which makes sense as we only have the + strand)
idx_leading <- which(ann$rep_class == "leftward"  & ann$classified_frac >= 1)
idx_lagging <- which(ann$rep_class == "rightward" & ann$classified_frac >= 1)
message(sprintf("Leading windows (leftward forks):  %d", length(idx_leading)))
message(sprintf("Lagging windows (rightward forks): %d", length(idx_lagging)))

# if we lost a lot of windows by filtering, we will know that if they are fewer than 100
if (length(idx_leading) < 100 || length(idx_lagging) < 100)
  warning("Very few windows in one class -- check BED file coverage and assign_threshold.")

###############################################################################
# compute profiles and log2fc on the sbs192 matrix level
#
# leading_192[j] = mean rate of mutation channel j across all leading windows
# lagging_192[j] = mean rate of mutation channel j across all lagging windows
# log2fc_192[j]  = log2(leading / lagging) for channel j
#
# positive log2fc = channel j is more frequent when the + strand is the
# leading template (leftward forks)
# negative log2fc = more frequent when the + strand is the
# lagging template (rightward forks)
###############################################################################

leading_192 <- colMeans(mat192[idx_leading, ])
lagging_192 <- colMeans(mat192[idx_lagging, ])

log2fc_192 <- log2((leading_192 + PSEUDOCOUNT) / (lagging_192 + PSEUDOCOUNT))

###############################################################################
# rc anticorrelation sanity check
#
# for any real strand-dependent signal, channel X and its RC partner
# should show OPPOSITE enrichment patterns (one up on leading, other up
# on lagging). This means log2fc[X] and log2fc[RC(X)] should be
# strongly anti-correlated (Spearman rho close to -1).
#
# rho is near 0: no strand signal, or window annotation is wrong
# rho is near +1: something is inverted (strand assignment or RC map)
###############################################################################

rc_map <- make_rc_map(sbs_cols)

rc_anticor <- cor(
  log2fc_192[rc_map$meta$orig_label],
  log2fc_192[rc_map$meta$rc_label],
  method = "spearman"
)
message(sprintf(
  "\nRC anticorrelation QC: Spearman rho = %.3f  (expect strongly negative for real signal)",
  rc_anticor))

###############################################################################
# parse mutation labels for downstream use
#
# tri = the 3-mer (e.g. "ACG")
# ref = middle base of tri = the mutated base (e.g. "C")
# alt = the new base after mutation (e.g. "T")
###############################################################################

tri <- sub("_.*$", "", names(log2fc_192))   # e.g. "ACG"
ref <- substr(tri, 2, 2)                    # e.g. "C"  (the reference base)
alt <- sub("^.*_",  "", names(log2fc_192))  # e.g. "T"  (the alternate base)

# indices for specific mutation types (on the + strand as recorded in mat192)
# C>T on + strand: ref == "C", alt == "T"
ct_idx <- which(ref == "C" & alt == "T")
# G>A on + strand: ref == "G", alt == "A"
# G>A on the + strand is C>T on the - strand (same physical event,
# opposite strand). comparing ct_idx vs ga_idx enrichment between leading
# and lagging windows is therefore a direct measure of C>T strand asymmetry.
ga_idx <- which(ref == "G" & alt == "A")

# C>G on + strand
cg_idx <- which(ref == "C" & alt == "G")
# G>C on + strand (= C>G on - strand)
gc_idx <- which(ref == "G" & alt == "C")

###############################################################################
# plot a = full SBS192 log2FC
# shows the strand asymmetry for every one of the 192 channels.
###############################################################################

df_fc <- data.frame(
  mutation = names(log2fc_192),
  log2fc   = as.numeric(log2fc_192),
  ref      = ref,
  alt      = alt,
  class    = paste0(ref, ">", alt)
)

pdf(file.path(output_dir, "log2FC_SBS192_full.pdf"), width = 18, height = 6)
print(
  ggplot(df_fc, aes(x = seq_along(log2fc), y = log2fc,
                    colour = class)) +
    geom_hline(yintercept = 0, colour = "black", linewidth = 0.5) +
    geom_point(size = 0.8, alpha = 0.7) +
    labs(title = sprintf("SBS192 strand asymmetry (RC anticorrelation = %.3f)", rc_anticor),
         subtitle = "Positive = enriched on leading (leftward) windows; negative = enriched on lagging",
         x = "SBS192 channel index", y = "log2(leading / lagging)",
         colour = "Mutation class") +
    theme_bw(base_size = 11)
)
dev.off()

###############################################################################
# plot b:  mean log2FC per mutation class (12 strand-specific classes)
# this summarises the direction and magnitude of asymmetry per mutation type.
###############################################################################

group_means <- tapply(df_fc$log2fc, df_fc$class, mean, na.rm = TRUE)

class_colours <- c(
  "A>C"="#4E79A7","A>G"="#F28E2B","A>T"="#E15759",
  "C>A"="#76B7B2","C>G"="#59A14F","C>T"="#EDC948",
  "G>A"="#B07AA1","G>C"="#FF9DA7","G>T"="#9C755F",
  "T>A"="#BAB0AC","T>C"="#D37295","T>G"="#A0CBE8"
)

pdf(file.path(output_dir, "mean_log2FC_by_class_SBS192.pdf"), width = 9, height = 5)
par(mar = c(6, 5, 4, 2))
bp <- barplot(group_means,
              col    = class_colours[names(group_means)],
              las    = 2,
              ylab   = "Mean log2(leading / lagging)",
              main   = "Mean strand asymmetry by mutation class (SBS192)")
abline(h = 0, lwd = 2)
dev.off()

###############################################################################
# plot c: C>T vs G>A summary
#
# C>T (ct_idx) = C>T on the + strand = C>T on the leading template
#                when in a leftward window
# G>A (ga_idx) = G>A on the + strand = C>T on the - strand = C>T on the
#                lagging template when in a leftward window
#
# if transcription-coupled repair or replication asymmetry is active,
# these two groups will have opposite log2FC values.
###############################################################################

ct_mean <- mean(log2fc_192[ct_idx], na.rm = TRUE)
ga_mean <- mean(log2fc_192[ga_idx], na.rm = TRUE)

pdf(file.path(output_dir, "CtoT_vs_GtoA_summary.pdf"), width = 5, height = 5)
barplot(c("C>T\n(+ strand)" = ct_mean, "G>A\n(+ strand)\n= C>T on - strand" = ga_mean),
        col  = c("firebrick", "steelblue"),
        ylab = "Mean log2(leading / lagging)",
        main = "C>T strand asymmetry\n(should be opposite if signal is real)")
abline(h = 0, lwd = 2)
dev.off()

###############################################################################
# what Marketa asked: CpG paired violin plots
#
# for each CpG trinucleotide context, plot the per-window mutation rate of
# the forward mutation (e.g. ACG_T = C>T in ACG context) vs its rc partner
# (e.g. CGT_A = T>A in CGT context, which is C>T on the - strand reading
# the same CpG) split by replication direction class.
#
# what to look for:
#   - for a spontaneous CpG deamination signal (process 8 in the paper),
#     ACG_T should be similarly elevated in both leftward and rightward windows
#     (strand-independent process -- same on both templates).
#   - for a replication-asymmetric signal, the two distributions will differ
#     between leftward and rightward windows.
###############################################################################

# CpG C>T contexts and their rc partners (C>T on the - strand)
cpg_pairs <- list(
  c("ACG_T", "CGT_A"),   # A[C>T]G  =  C[G>A]T on - strand
  c("CCG_T", "CGG_A"),   # C[C>T]G  =  C[G>A]G on - strand
  c("GCG_T", "CGC_A"),   # G[C>T]G  =  C[G>A]C on - strand
  c("TCG_T", "CGA_A")    # T[C>T]G  =  C[G>A]A on - strand
)

for (pair in cpg_pairs) {
  ctx1 <- pair[1]   # forward: C>T on + strand in CpG context
  ctx2 <- pair[2]   # reverse: G>A on + strand = C>T on - strand in CpG context
  
  # columns exist before plotting
  if (!all(c(ctx1, ctx2) %in% names(ann))) {
    warning(sprintf("Columns %s or %s not found -- skipping", ctx1, ctx2))
    next
  }
  
  df <- ann[rep_class %in% c("leftward", "rightward"),
            .(rep_class, v1 = get(ctx1), v2 = get(ctx2))]
  
  df_long <- melt(df,
                  id.vars     = "rep_class",
                  measure.vars = c("v1", "v2"),
                  variable.name = "strand_view",
                  value.name    = "rate")
  
  # label  v1 = + strand (forward), v2 = - strand (reverse complement)
  df_long[, strand_view := fifelse(strand_view == "v1",
                                   paste0(ctx1, "\n(+ strand)"),
                                   paste0(ctx2, "\n(- strand RC)"))]
  
  p <- ggplot(df_long, aes(x = interaction(strand_view, rep_class),
                           y = rate, fill = rep_class)) +
    geom_violin(trim = FALSE, alpha = 0.4) +
    geom_boxplot(width = 0.15, outlier.size = 0.3, alpha = 0.8) +
    scale_fill_manual(values = c("leftward" = "firebrick", "rightward" = "steelblue"),
                      labels = c("leftward" = "Leading (leftward fork)",
                                 "rightward" = "Lagging (rightward fork)")) +
    labs(title   = sprintf("CpG C>T asymmetry: %s vs %s", ctx1, ctx2),
         subtitle = "Leftward = + strand is leading template; Rightward = + strand is lagging template",
         x = "Mutation context × replication direction",
         y = "Per-window mutation rate",
         fill = "Replication class") +
    theme_bw(base_size = 11) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  
  ggsave(file.path(output_dir, paste0(ctx1, "_vs_", ctx2, "_violin.pdf")),
         p, width = 8, height = 5)
}

###############################################################################
# Paper-like attempt per chromosome CpG asymmetry vs replication direction track
#
# For each window, compute the CpG C>T asymmetry:
#   log2(sum of CpG C>T rates on + strand / sum of CpG G>A rates on + strand)
#
# CpG C>T on + strand (idx_cpg_t): mutation where a C in a CpG dinucleotide
#   changes to T, recorded relative to the + strand.
# CpG G>A on + strand (idx_cpg_a): mutation where a G changes to A in a GpC
#   context on the + strand -- this IS C>T on the OPPOSITE strand at a CpG site.
#
# If spontaneous deamination (process 8) dominates and is strand-symmetric,
# this ratio should be near 0 everywhere.
# If there is a replication strand bias, the ratio should correlate with
# rep_signal (right_frac - left_frac), which captures the dominant fork direction.
#
# The rep_signal variable:
#   right_frac - left_frac > 0 => mostly rightward (lagging on + strand)
#   right_frac - left_frac < 0 => mostly leftward  (leading on + strand)
###############################################################################

# CpG trinucleotides (C is the middle base, next base is G)
cpg_tris <- c("ACG", "CCG", "GCG", "TCG")
rc_tris  <- vapply(cpg_tris, rev_comp_seq, character(1))  # CGT, CGG, CGC, CGA

# Forward: C>T at CpG on + strand (C in middle of ACG/CCG/GCG/TCG, alt = T)
idx_cpg_t <- which(tri %in% cpg_tris & alt == "T")

# Reverse: G>A on + strand in GpC context = C>T on - strand at CpG
# (G is the middle base of CGT/CGG/CGC/CGA, alt = A)
idx_cpg_a <- which(tri %in% rc_tris & alt == "A")

if (length(idx_cpg_t) != 4 || length(idx_cpg_a) != 4)
  stop(sprintf("CpG context matching failed: found %d forward and %d reverse contexts",
               length(idx_cpg_t), length(idx_cpg_a)))

message("CpG C>T forward contexts: ", paste(colnames(mat192)[idx_cpg_t], collapse=", "))
message("CpG G>A reverse contexts: ", paste(colnames(mat192)[idx_cpg_a], collapse=", "))

# Per-window CpG asymmetry: log2(C>T on + / G>A on +)
# Positive => C>T more frequent on + strand => if in leading window, more
#             deamination on the leading template
ann$cpg_asym_log <- log2(
  (rowSums(mat192[, idx_cpg_t, drop = FALSE], na.rm = TRUE) + PSEUDOCOUNT) /
    (rowSums(mat192[, idx_cpg_a, drop = FALSE], na.rm = TRUE) + PSEUDOCOUNT)
)
ann$cpg_asym_log[!is.finite(ann$cpg_asym_log)] <- NA

# Replication direction signal: positive = rightward dominant, negative = leftward
# This is a continuous version of rep_class for correlation analysis
ann$rep_signal <- ann$right_cov - ann$left_cov

ann$rep_timing <- ann$timing_w
###############################################################################
# CpG asymmetry vs replication direction per chromosome
###############################################################################

all_asym <- numeric(0)
all_rep  <- numeric(0)
all_timing <- numeric(0)

cor_results <- data.table(chr = character(), spearman = numeric(), n = integer())

pdf(file.path(output_dir, "CpG_asymmetry_vs_replication.pdf"), width = 14, height = 4)

for (chr_i in paste0("chr", c(1:22, "X"))) {
  
  sub <- ann[chr == chr_i &
               classified_frac >= 0.6 &
               !is.na(cpg_asym_log)]
  
  if (nrow(sub) < 50) next
  
  sub <- sub[order(start)]
  x_vals <- sub$start
  
  # Clip extreme values before smoothing to avoid a few outlier windows
  # dominating the smoothed line
  y_asym <- pmax(pmin(sub$cpg_asym_log, 3), -3)
  
  # k=200 gives ~2 Mb smoothing scale at 10kb windows; adjust if needed
  k <- min(200, floor(nrow(sub) / 2))
  smooth_asym <- zoo::rollmean(y_asym,          k = k, fill = NA)
  rep_smooth  <- zoo::rollmean(sub$rep_signal,  k = k, fill = NA)
  time_smooth  <- zoo::rollmean(sub$rep_timing,  k = k, fill = NA)
  
  # Scale rep_signal to the same visual range as the asymmetry for overlay
  rep_range  <- max(abs(rep_smooth), na.rm = TRUE)
  asym_range <- max(abs(smooth_asym), na.rm = TRUE)
  rep_scaled <- if (rep_range > 0)
    rep_smooth * (asym_range / rep_range) else rep_smooth
  
  time_range <- max(abs(time_smooth), na.rm = TRUE)
  time_scaled <- if (time_range > 0) 
    time_smooth * (asym_range / time_range) else time_smooth
  
  plot(x_vals, smooth_asym,
       type = "l", col = "darkgreen", lwd = 2,
       ylim = range(c(smooth_asym, rep_scaled, time_scaled), na.rm = TRUE) * 1.2,
       xlab = paste0(chr_i, " coordinate (bp)"),
       ylab = "log2(CpG C>T / CpG G>A)",
       main = paste0(chr_i, ": CpG mutation asymmetry vs replication direction"))
  
  lines(x_vals, rep_scaled, col = "red", lwd = 2, lty = 2)
  lines(x_vals, time_scaled, col = "grey", lwd = 2, lty = 3)
  
  abline(h = 0, lty = 3)
  legend("topright",
         legend = c("CpG C>T asymmetry (log2)", "Replication direction (scaled)", "Replication timing (scaled)"),
         col = c("darkgreen", "red", "grey"), lwd = 2, lty = c(1, 2, 3), bty = "n")
  
  # Accumulate for global correlation
  valid <- is.finite(smooth_asym) & is.finite(rep_smooth)
  if (sum(valid) >= 10) {
    corr <- cor(smooth_asym[valid], rep_smooth[valid], method = "spearman")
    cor_results <- rbind(cor_results,
                         data.table(chr = chr_i, spearman = corr, n = sum(valid)))
    all_asym <- c(all_asym, smooth_asym[valid])
    all_rep  <- c(all_rep,  rep_smooth[valid])
    all_timing <- c(all_timing, time_smooth[valid])
  }
}

dev.off()

# FIX: global correlation now correctly uses values from ALL chromosomes
global_corr <- cor(all_asym, all_rep, method = "spearman")
message(sprintf("\nGlobal Spearman correlation (CpG asymmetry vs replication): %.4f",
                global_corr))
print(cor_results)

fwrite(cor_results, file.path(output_dir, "CpG_asymmetry_correlation_by_chr.tsv"), sep = "\t")

###############################################################################
# SAVE CORE OUTPUT VECTORS
###############################################################################

write_vec(leading_192, file.path(output_dir, "leading_192.tsv"))
write_vec(lagging_192, file.path(output_dir, "lagging_192.tsv"))
write_vec(log2fc_192,  file.path(output_dir, "log2FC_192.tsv"))

message("\nAll done. Results saved to: ", output_dir)

