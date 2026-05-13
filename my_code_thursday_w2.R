###############################################################################
# EXPERIMENT 1 | WEEK 2
# Leading vs Lagging strand mutational asymmetry
# Germline low-frequency VAFs 
# TOPMed SBS192 matrix, 10kb, aggregate across patients, hg38
# Gnomad SBS192 matrix, 100kb, aggregate across patients, hg19
# bed file of directionalities is 20kb bins, hg19
# no need to liftover for gnomad, but there is issue with overlapping
# liftover for topmed, I only take windows that match perfectly
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
          "ggplot2", "data.table", "zoo", "ggpubr")
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

# write a named numeric vector as a two-column tsv (channel, value)
write_vec <- function(x, file)
  fwrite(data.table(channel = names(x), value = as.numeric(x)), file, sep = "\t")

###############################################################################
# READ & STANDARDISE INPUT FILES
###############################################################################
# read files
mat_dt <- fread(matrix_file)
bed_dt <- fread(bed_file)
# normalise chrom label
mat_dt[, chr        := normalize_chr(chr)]
bed_dt[, chromosome := normalize_chr(chromosome)]
# check chromosome lengths and the differences of the lengths by files
mat_lengths <- mat_dt[, .(mat_max_end = max(end)), by = chr]
bed_lengths <- bed_dt[, .(chr = chromosome, bed_max_end = max(end)),
                      by = chromosome][, chromosome := NULL]

length_check <- merge(mat_lengths, bed_lengths, by = "chr", all = TRUE)
length_check[, diff := mat_max_end - bed_max_end]
# order by absolute difference, if very large then bad
print(length_check[order(-abs(diff))])

# build the mat_dt fully with start, end, id
mat_dt[, `:=`(
  chr    = normalize_chr(chr),
  start  = start,
  end    = end,
  mat_id = .I                   
)]

# bed file sanity checks
stopifnot(all(c("chromosome", "start", "end",
                "isLeftReplicating", "isRightReplicating") %in% names(bed_dt)))

# build bed_dt fully with start, end, id
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
  sum(bed_dt$isLeft == 1 & bed_dt$isRight == 1), # nothing is both but we are not checking intervals yet
  sum(bed_dt$isLeft == 0 & bed_dt$isRight == 0)
))

###############################################################################
# if the bed file for replication direction is lifted over, there are overlapping 
# regions which have both left and right replication at the same time, so these
# must be excluded
# granges uses its own s4 syntax, not base r nor data table not tidyverse
# copy my lifted over hg38 bed and the row indices
bed_dt <- copy(bed_dt)
bed_dt[, row_id := .I]
# makes ranges out of my start and end, different indexing than bed file
bed_gr <- GRanges(
  bed_dt$chr,
  IRanges(bed_dt$start + 1L, bed_dt$end)
)
# left indices are all rows where we only have left and vice versa, indices are unique per all rows
left_idx  <- which(bed_dt$isLeft  == 1L & bed_dt$isRight == 0L)
right_idx <- which(bed_dt$isRight == 1L & bed_dt$isLeft  == 0L)

# find any overlap between left and right intervals
# overlap has to be at least 1 nt so quite a lot of elimination here for spurious 
# liftover blocks, only matching by indices and I eliminate both of them since unsure
# whether they are right or left, so both rows eliminated
hits_lr <- findOverlaps(
  bed_gr[left_idx],
  bed_gr[right_idx],
  ignore.strand = TRUE,
  minoverlap = 1L
)
# this just tells the indices of the rows to remove and which match
data.table(
  left_row  = left_idx[queryHits(hits_lr)],
  right_row = right_idx[subjectHits(hits_lr)]
)

# rows to remove: both sides of every conflicting pair
drop_idx <- unique(c(
  left_idx[queryHits(hits_lr)],
  right_idx[subjectHits(hits_lr)]
))

# new clean bed to compare with the previous bed intervals
bed_clean <- bed_dt[-drop_idx]

message(sprintf(
  "BED intervals -- leftward only: %d | rightward only: %d | both: %d | neither: %d",
  sum(bed_clean$isLeft == 1 & bed_clean$isRight == 0),
  sum(bed_clean$isRight == 1 & bed_clean$isLeft == 0),
  sum(bed_clean$isLeft == 1 & bed_clean$isRight == 1),
  sum(bed_clean$isLeft == 0 & bed_clean$isRight == 0)
))

# now the lifted-over bed file is cleaned of all windows

###############################################################################
# overlap of windows: assign replication direction to each topmed window
#
# bed format is 0-based half-open [start, end).
# granges is 1-based closed [start, end].
# conversion syntax: granges_start = bed_start + 1, granges_end = bed_end.
###############################################################################

# here I am constructing absolute overlaps of my matrix and bed windows
# again making ranges and putting left, right and timing also there
mat_gr <- GRanges(mat_dt$chr, IRanges(mat_dt$start + 1L, mat_dt$end)) 
bed_gr <- GRanges(bed_clean$chr, IRanges(bed_clean$start + 1L, bed_clean$end)) 

mcols(bed_gr)$isLeft <- bed_clean$isLeft 
mcols(bed_gr)$isRight <- bed_clean$isRight 
mcols(bed_gr)$timing <- bed_clean$timing 

# finding overlaps
hits <- findOverlaps(mat_gr, bed_gr, ignore.strand = TRUE) 
if (length(hits) == 0) stop("No overlaps found. Check chromosome names and coordinate systems.")

# query and subject definition
q <- queryHits(hits)
s <- subjectHits(hits)

# for ann which will be the overlapped thing
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
  c("TCG_T", "CGA_A"), 
  c("ACA_T", "TGT_A"),   
  c("CCA_T", "TGG_A"),   
  c("GCA_T", "TGC_A"),   
  c("TCA_T", "TGA_A"),
  c("ACT_T", "AGT_A"),   
  c("CCT_T", "AGG_A"),   
  c("GCT_T", "AGC_A"),  
  c("TCT_T", "AGA_A"),
  c("ACC_T", "GGT_A"),  
  c("CCC_T", "GGG_A"),   
  c("GCC_T", "GGC_A"),   
  c("TCC_T", "GGA_A"),
  c("ATG_C", "CAT_G"),   
  c("CTG_C", "CAG_G"),   
  c("GTG_C", "CAC_G"),   
  c("TTG_C", "CAA_G"),
  c("ATA_C", "TAT_G"),
  c("CTA_C", "TAG_G"),   
  c("GTA_C", "TAC_G"),   
  c("TTA_C", "TAA_G"),
  c("ATT_C", "AAT_G"),   
  c("CTT_C", "AAG_G"),  
  c("GTT_C", "AAC_G"),  
  c("TTT_C", "AAA_G"),
  c("ATC_C", "CAT_G"),   
  c("CTC_C", "CAG_G"),   
  c("GTC_C", "CAC_G"),   
  c("TTC_C", "CAA_G")
)


### ad p values on the plots among all combos - Lead vs lead should be ns and vice versa
for (pair in cpg_pairs) {
  ctx1 <- pair[1]
  ctx2 <- pair[2]
  
  if (!all(c(ctx1, ctx2) %in% names(ann))) {
    warning(sprintf("Columns %s or %s not found -- skipping", ctx1, ctx2))
    next
  }
  
  plus_label  <- paste0(ctx1, "\n(+ strand)")
  minus_label <- paste0(ctx2, "\n(- strand RC)")
  
  df <- ann[rep_class %in% c("leftward", "rightward"),
            .(rep_class, v1 = get(ctx1), v2 = get(ctx2))]
  
  df_long <- melt(
    df,
    id.vars = "rep_class",
    measure.vars = c("v1", "v2"),
    variable.name = "strand_view",
    value.name = "rate"
  )
  
  df_long[, strand_view := fifelse(
    strand_view == "v1",
    plus_label,
    minus_label
  )]
  
  # Biological role:
  # leftward:  + = leading, - = lagging
  # rightward: + = lagging, - = leading
  df_long[, template_role := fcase(
    rep_class == "leftward"  & strand_view == plus_label,  "leading",
    rep_class == "leftward"  & strand_view == minus_label, "lagging",
    rep_class == "rightward" & strand_view == plus_label,  "lagging",
    rep_class == "rightward" & strand_view == minus_label, "leading"
  )]
  
  # Keep x-axis order fixed and biologically readable
  df_long[, group := factor(
    paste(strand_view, rep_class, sep = "_"),
    levels = c(
      paste0(minus_label, "_leftward"),
      paste0(plus_label,  "_leftward"),
      paste0(minus_label, "_rightward"),
      paste0(plus_label,  "_rightward")
    )
  )]
  
  g1 <- paste0(plus_label,  "_leftward")
  g2 <- paste0(minus_label, "_leftward")
  g3 <- paste0(plus_label,  "_rightward")
  g4 <- paste0(minus_label, "_rightward")
  
  # one-vs-all comparisons using the first group as reference
  ref_group <- g1
  comparisons <- lapply(setdiff(c(g1, g2, g3, g4), ref_group), function(g) c(ref_group, g))
  
  p <- ggplot(df_long, aes(x = group, y = rate, fill = template_role)) +
    geom_violin(trim = FALSE, alpha = 0.4) +
    geom_boxplot(width = 0.15, outlier.size = 0.3, alpha = 0.8) +
    scale_fill_manual(values = c("leading" = "firebrick",
                                 "lagging" = "steelblue")) +
    stat_compare_means(
      comparisons = comparisons,
      method = "wilcox.test",
      label = "p.signif",
      hide.ns = FALSE
    ) +
    labs(
      title = sprintf("CpG C>T asymmetry: %s vs %s", ctx1, ctx2),
      subtitle = "Colors: leading (red) vs lagging (blue)",
      x = "Mutation context × replication direction",
      y = "Per-window mutation rate",
      fill = "Template role"
    ) +
    theme_bw(base_size = 11) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  
  ggsave(file.path(output_dir, paste0(ctx1, "_vs_", ctx2, "_violin.pdf")),
         p, width = 8, height = 6)
}


for (pair in cpg_pairs) {
  
  ctx1 <- pair[1]
  ctx2 <- pair[2]
  
  if (!all(c(ctx1, ctx2) %in% names(ann))) {
    warning(sprintf("Skipping %s vs %s (missing columns)", ctx1, ctx2))
    next
  }
  
  df <- ann[rep_class %in% c("leftward", "rightward"),
            .(rep_class,
              plus  = get(ctx1),
              minus = get(ctx2))]
  
  # correct mapping
  df[, leading := fifelse(rep_class == "leftward", plus, minus)]
  df[, lagging := fifelse(rep_class == "leftward", minus, plus)]
  
  # remove empty windows
  df <- df[(leading + lagging) > 0]
  
  # reshape for plotting (preserves pairing)
  df_long <- melt(df,
                  id.vars = c("rep_class"),
                  measure.vars = c("leading", "lagging"),
                  variable.name = "template",
                  value.name = "rate")
  
  # paired test (correct)
  # stat compare means on ggplot
  test <- wilcox.test(df$leading, df$lagging, paired = TRUE)
  pval <- test$p.value
  p_label <- ifelse(pval < 2.2e-16, "p < 2.2e-16", sprintf("p = %.2e", pval))
  
  # corrected plot
  p <- ggplot(df_long, aes(x = template, y = rate, fill = template)) +
    geom_violin(trim = FALSE, alpha = 0.4) +
    geom_boxplot(width = 0.15, outlier.size = 0.3) +
    scale_fill_manual(values = c("leading" = "firebrick",
                                 "lagging" = "steelblue")) +
    labs(
      title = sprintf("Leading vs lagging mutation rates: %s vs %s", ctx1, ctx2),
      subtitle = sprintf("Paired Wilcoxon %s", p_label),
      x = "Template strand",
      y = "Per-window mutation rate"
    ) +
    theme_bw(base_size = 11)
  
  ggsave(
    file.path(output_dir, paste0(ctx1, "_vs_", ctx2, "_leading_vs_lagging_violin.pdf")),
    p,
    width = 6,
    height = 5
  )
  
  message(sprintf("%s vs %s: p = %.3e", ctx1, ctx2, pval))
}

for (pair in cpg_pairs) {
  ctx1 <- pair[1]
  ctx2 <- pair[2]
  
  df <- ann[rep_class %in% c("leftward", "rightward"),
            .(rep_class,
              plus  = get(ctx1),
              minus = get(ctx2))]
  
  df <- df[!is.na(plus) & !is.na(minus) & (plus + minus) > 0]
  
  # Raw strand asymmetry: + vs - encoded channel
  df[, plus_share := plus / (plus + minus)]
  
  # Replication-aware asymmetry
  df[, leading := fifelse(rep_class == "leftward", plus, minus)]
  df[, lagging := fifelse(rep_class == "leftward", minus, plus)]
  df[, leading_share := leading / (leading + lagging)]
  
  print(ctx1)
  print(df[, .(
    n = .N,
    median_plus_share = median(plus_share),
    median_leading_share = median(leading_share),
    pct_plus_gt_minus = mean(plus > minus),
    pct_leading_gt_lagging = mean(leading > lagging)
  ), by = rep_class])
}

df[, log_plus_minus := log2((plus + 1e-12) / (minus + 1e-12))]

df[, .(
  median_log_plus_minus = median(log_plus_minus),
  pct_plus_gt_minus = mean(plus > minus)
), by = rep_class]

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
###############################################################################
###############################################################################
# CpG trinucleotides (C is the middle base, next base is G)
colnames_mat <- colnames(mat192)

# Split into context and alt
tri <- sub("_.*", "", colnames_mat) 
alt <- sub(".*_", "", colnames_mat)  

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
# Per-window CpG asymmetry: log2(C>T on + / G>A on +)
# Positive => C>T more frequent on + strand
# data-based psudocount rather
ann$cpg_asym_log <- log2(
  (rowSums(mat192[, idx_cpg_t, drop = FALSE]) + PSEUDOCOUNT) /
    (rowSums(mat192[, idx_cpg_a, drop = FALSE]) + PSEUDOCOUNT)
)

# this fixes the mistake in scaling
#ann$cpg_asym_log[!is.finite(ann$cpg_asym_log)] <- NA

# Replication direction signal:
# positive = rightward dominant
# negative = leftward dominant
ann$rep_signal <- ann$right_cov - ann$left_cov

# Raw replication timing
ann$rep_timing <- ann$timing_w

###############################################################################
# CpG asymmetry vs replication direction and replication timing derivative
###############################################################################

all_asym <- numeric(0)
all_rep  <- numeric(0)
all_timing_diff <- numeric(0)

cor_results <- data.table(chr = character(), spearman = numeric(), n = integer())

pdf(file.path(output_dir, "CpG_asymmetry_vs_replication.pdf"), width = 14, height = 4)

for (chr_i in paste0("chr", c(1:22, "X"))) {
  
  sub <- ann[
    chr == chr_i &
      #classified_frac >= 1 &
      !is.na(cpg_asym_log) &
      !is.na(rep_timing)
  ]
  
  if (nrow(sub) < 50) next
  
  # Important: order windows before taking diff()
  sub <- sub[order(start)]
  x_vals <- sub$start
  
  # Timing derivative / discrete timing change between adjacent windows
  # First window has no previous value; set to 0 for plotting
  sub[, rep_timing_diff := c(0, diff(rep_timing))]
  
  # If you want a derivative per bp instead of simple difference, use:
  # sub[, rep_timing_diff := c(0, diff(rep_timing) / diff(start))]
  
  # Clip extreme CpG asymmetry values before smoothing
  y_asym <- pmax(pmin(sub$cpg_asym_log, 3), -3)
  
  # k=200 gives ~2 Mb smoothing scale at 10kb windows
  k <- min(200, floor(nrow(sub) / 2))
  
  smooth_asym <- zoo::rollmean(y_asym,                 k = k, fill = NA)
  rep_smooth  <- zoo::rollmean(sub$rep_signal,         k = k, fill = NA)
  time_smooth <- zoo::rollmean(sub$rep_timing_diff,    k = k, fill = NA)
  
  # Scale rep_signal to same visual range as asymmetry
  rep_range  <- max(abs(rep_smooth), na.rm = TRUE)
  asym_range <- max(abs(smooth_asym), na.rm = TRUE)
  
  rep_scaled <- if (is.finite(rep_range) && rep_range > 0) {
    rep_smooth * (asym_range / rep_range)
  } else {
    rep_smooth
  }
  
  # Scale timing derivative to same visual range as asymmetry
  time_range <- max(abs(time_smooth), na.rm = TRUE)
  
  time_scaled <- if (is.finite(time_range) && time_range > 0) {
    time_smooth * (asym_range / time_range)
  } else {
    time_smooth
  }
  
  plot(
    x_vals, smooth_asym,
    type = "l",
    col = "darkgreen",
    lwd = 2,
    ylim = range(c(smooth_asym, rep_scaled, time_scaled), na.rm = TRUE) * 1.2,
    xlab = paste0(chr_i, " coordinate (bp)"),
    ylab = "scaled tracks",
    main = paste0(chr_i, ": CpG mutation asymmetry vs replication direction")
  )
  
  lines(x_vals, rep_scaled, col = "red", lwd = 2, lty = 2)
  lines(x_vals, time_scaled, col = "grey40", lwd = 2, lty = 3)
  
  abline(h = 0, lty = 3)
  
  legend(
    "topright",
    legend = c(
      "CpG C>T asymmetry (log2)",
      "Replication direction (scaled)",
      "Replication timing diff (scaled)"
    ),
    col = c("darkgreen", "red", "grey40"),
    lwd = 2,
    lty = c(1, 2, 3),
    bty = "n"
  )
  
  # Accumulate for global correlation - check that (without that)
  valid <- is.finite(smooth_asym) & is.finite(rep_smooth)
  
  if (sum(valid) >= 10) {
    corr <- cor(smooth_asym[valid], rep_smooth[valid], method = "spearman")
    
    cor_results <- rbind(
      cor_results,
      data.table(chr = chr_i, spearman = corr, n = sum(valid))
    )
    
    all_asym <- c(all_asym, smooth_asym[valid])
    all_rep  <- c(all_rep,  rep_smooth[valid])
    all_timing_diff <- c(all_timing_diff, time_smooth[valid])
  }
}

dev.off()

# Global correlation: CpG asymmetry vs replication direction
global_corr <- cor(all_asym, all_rep, method = "spearman")

message(sprintf(
  "\nGlobal Spearman correlation (CpG asymmetry vs replication direction): %.4f",
  global_corr
))

print(cor_results)

fwrite(
  cor_results,
  file.path(output_dir, "CpG_asymmetry_correlation_by_chr.tsv"),
  sep = "\t"
)

###############################################################################
# NOTES
###############################################################################
#
# take mat_192 and order it depending on RC
# then do a A-B columns of which some will be R and some L replicating
# or take log2fc(A/B) whether CpG to T are more asymmetrical than other types


# ----------------------------
# 1) Reverse-complement helper
# ----------------------------
rev_comp_sbs <- function(x) {
  comp <- c(A = "T", T = "A", C = "G", G = "C")
  tri <- sub("_.*$", "", x)
  alt <- sub("^.*_", "", x)
  
  tri_rc <- paste(comp[rev(strsplit(tri, "", fixed = TRUE)[[1]])], collapse = "")
  alt_rc <- comp[[alt]]
  
  paste0(tri_rc, "_", alt_rc)
}

# ----------------------------
# 2) Build the 96 pyrimidine/RC pairs
# ----------------------------
make_pair_table <- function(mat192) {
  cols <- colnames(mat192)
  tri  <- sub("_.*$", "", cols)
  
  # canonical pyrimidine-side channels: middle base C or T
  pyr_cols <- cols[substr(tri, 2, 2) %in% c("C", "T")]
  rc_cols  <- vapply(pyr_cols, rev_comp_sbs, character(1))
  
  pair_tbl <- data.table(
    pyr_col = pyr_cols,
    rc_col  = rc_cols
  )
  
  pair_tbl <- pair_tbl[rc_col %chin% cols]
  pair_tbl[, subclass := sub("^([ACGT]{3}).*$", "\\1", pyr_col)]
  pair_tbl[, pair_label := paste0(pyr_col, " / ", rc_col)]
  
  pair_tbl
}

pair_tbl <- make_pair_table(mat192)
if (nrow(pair_tbl) != 96L) {
  warning(sprintf("Expected 96 pairs, found %d. Check matrix column names.", nrow(pair_tbl)))
}

eps <- 1e-6

# ----------------------------
# 3) Build a long table with both metrics
# ----------------------------
pair_long <- rbindlist(lapply(seq_len(nrow(pair_tbl)), function(i) {
  pcol <- pair_tbl$pyr_col[i]
  rcol <- pair_tbl$rc_col[i]
  subclass <- pair_tbl$subclass[i]
  pair_label <- pair_tbl$pair_label[i]
  
  dt <- ann[
    rep_class %in% c("leftward", "rightward"),
    .(
      rep_class,
      pair = pair_label,
      subclass = subclass,
      pyr = get(pcol),
      rc  = get(rcol)
    )
  ]
  
  dt <- dt[!is.na(pyr) & !is.na(rc)]
  
  # Metric 1: subtraction
  dt[, diff := pyr - rc]
  
  # Metric 2: log2 fold-change
  dt[, log2fc := log2((pyr + eps) / (rc + eps))]
  
  dt
}), use.names = TRUE, fill = TRUE)

# ----------------------------
# 4) Order pairs within subclass by absolute effect size
# ----------------------------
pair_order <- pair_long[, .(
  med_abs = median(abs(log2fc), na.rm = TRUE)
), by = .(subclass, pair)][order(subclass, -med_abs)]$pair

pair_long[, pair := factor(pair, levels = pair_order)]

# ----------------------------
# 5) Generic plotting function
# ----------------------------
plot_signed_metric <- function(dt, metric, file, title, ylab) {
  dt <- copy(dt)
  
  # per-pair median sign + p-value vs 0
  sum_dt <- dt[, .(
    med = median(get(metric), na.rm = TRUE),
    p   = wilcox.test(get(metric), mu = 0, exact = FALSE)$p.value
  ), by = .(subclass, rep_class, pair)]
  
  # Compute a small positive offset above 0 based on actual data range
  y_offset <- quantile(abs(dt[[metric]]), 0.8, na.rm = TRUE) * 0.3
  sum_dt[, y := y_offset]
  
  sum_dt[, p_adj := p.adjust(p, method = "BH")]
  
  sum_dt[, sig_color := fifelse(
    p_adj < 0.05 & med > 0, "above 0",
    fifelse(p_adj < 0.05 & med < 0, "below 0", "ns")
  )]
  
  sum_dt[, p_lab := fifelse(
    p_adj < 2.2e-16,
    sprintf("med=%.2f\npval<2e-16", med),
    sprintf("med=%.2f\npval=%.1e", med, p_adj)
  )]
  
  # join sign/color back for fill color
  dt <- merge(
    dt,
    sum_dt[, .(subclass, rep_class, pair, sig_color)],
    by = c("subclass", "rep_class", "pair"),
    all.x = TRUE
  )
  
  # symmetric display range for readability
  lim <- quantile(abs(dt[[metric]]), 0.85, na.rm = TRUE) #use 0.85 for log2fc, 0.9 for subtraction
  if (!is.finite(lim) || lim == 0) {
    lim <- max(abs(dt[[metric]]), na.rm = TRUE)
  }
  
  p <- ggplot(dt, aes(x = pair, y = get(metric), fill = sig_color)) +
    geom_violin(trim = TRUE, color = NA, alpha = 0.45) +
    geom_boxplot(width = 0.12, outlier.shape = NA, alpha = 0.85) +
    geom_hline(yintercept = 0, linetype = "dashed") +
    geom_text(
      data = sum_dt,
      aes(x = pair, y = y, label = p_lab),
      inherit.aes = FALSE,
      angle = 90,
      hjust = 0,        # left-align text so it reads upward from anchor
      vjust = -0.2,     # small gap from the anchor point
      nudge_x = 0.55,   # shift slightly right of violin center
      size = 2.1
    ) +
    facet_grid(rep_class ~ subclass, scales = "free_x", space = "free_x") +
    coord_cartesian(ylim = c(-lim, lim * 1.15)) +
    scale_fill_manual(
      values = c(
        "above 0" = "firebrick",
        "below 0" = "steelblue",
        "ns"      = "grey70"
      )
    ) +
    labs(
      title = title,
      subtitle = "Rows = leftward/rightward; columns = CpG subclass; dashed line = 0",
      x = "",
      y = ylab,
      fill = "Median sign"
    ) +
    theme_bw(base_size = 10) +
    theme(
      axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
      panel.spacing = unit(0.7, "lines"),
      strip.text.x = element_text(size = 9),
      strip.text.y = element_text(size = 9)
    )
  
  ggsave(
    file.path(output_dir, file),
    p,
    width = 32,
    height = 6,
    device = cairo_pdf
  )
  
  p
}

# ----------------------------
# 6) Subtraction plot
# ----------------------------
plot_signed_metric(
  pair_long,
  metric = "diff",
  file = "CpG_96pairs_subtraction_by_direction_signed.pdf",
  title = "CpG pairwise subtraction across 96 SBS channels",
  ylab = "pyrimidine-side frequency - reverse-complement frequency"
)

# ----------------------------
# 7) Log2 fold-change plot
# ----------------------------
plot_signed_metric(
  pair_long,
  metric = "log2fc",
  file = "CpG_96pairs_log2fc_by_direction_signed.pdf",
  title = "CpG pairwise log2 fold-change across 96 SBS channels",
  ylab = "log2(pyrimidine-side / reverse-complement)"
)
