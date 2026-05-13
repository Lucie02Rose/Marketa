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
matrix_file      <- "gnomAD_100kb.txt"
bed_file         <- "per_base_territories_20kb_line_numbers.bed"
output_dir       <- "strand_specific_output_gnomad_hg19"
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
bed_gr <- GRanges(bed_dt$chr, IRanges(bed_dt$start + 1L, bed_dt$end)) 

mcols(bed_gr)$isLeft <- bed_dt$isLeft 
mcols(bed_gr)$isRight <- bed_dt$isRight 
mcols(bed_gr)$timing <- bed_dt$timing 

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
# leading vs lagging window indices based on the origins
#
# leftward windows  (+ strand = leading template) --> idx_leading
# rightward windows (+ strand = lagging template) --> idx_lagging
#
# windows below assign_threshold are excluded from both groups
###############################################################################
# of the ann matrix, we assign everything that is going to the right as lagging
# and everything going to the left as leading (which makes sense as we only have the + strand)
idx_leading <- which(ann$rep_class == "leftward")
idx_lagging <- which(ann$rep_class == "rightward")
message(sprintf("Leading windows (leftward forks):  %d", length(idx_leading)))
message(sprintf("Lagging windows (rightward forks): %d", length(idx_lagging)))

# if we lost a lot of windows by filtering, we will know that if they are fewer than 100
if (length(idx_leading) < 100 || length(idx_lagging) < 100)
  warning("Very few windows in one class -- check BED file coverage and assign_threshold.")


##############################################################################

leading_192 <- colMeans(mat192[idx_leading, ])
lagging_192 <- colMeans(mat192[idx_lagging, ])

log2fc_192 <- log2((leading_192 + PSEUDOCOUNT) / (lagging_192 + PSEUDOCOUNT))
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
  c("TCG_T", "CGA_A"), 
  c("ACA_T", "TGT_A"),   # A[C>T]G  =  C[G>A]T on - strand
  c("CCA_T", "TGG_A"),   # C[C>T]G  =  C[G>A]G on - strand
  c("GCA_T", "TGC_A"),   # G[C>T]G  =  C[G>A]C on - strand
  c("TCA_T", "TGA_A"),
  c("ACT_T", "AGT_A"),   # A[C>T]G  =  C[G>A]T on - strand
  c("CCT_T", "AGG_A"),   # C[C>T]G  =  C[G>A]G on - strand
  c("GCT_T", "AGC_A"),   # G[C>T]G  =  C[G>A]C on - strand
  c("TCT_T", "AGA_A"),
  c("ACC_T", "GGT_A"),   # A[C>T]G  =  C[G>A]T on - strand
  c("CCC_T", "GGG_A"),   # C[C>T]G  =  C[G>A]G on - strand
  c("GCC_T", "GGC_A"),   # G[C>T]G  =  C[G>A]C on - strand
  c("TCC_T", "GGA_A"),
  c("ATG_C", "CAT_G"),   # A[C>T]G  =  C[G>A]T on - strand
  c("CTG_C", "CAG_G"),   # C[C>T]G  =  C[G>A]G on - strand
  c("GTG_C", "CAC_G"),   # G[C>T]G  =  C[G>A]C on - strand
  c("TTG_C", "CAA_G"),
  c("ATA_C", "TAT_G"),   # A[C>T]G  =  C[G>A]T on - strand
  c("CTA_C", "TAG_G"),   # C[C>T]G  =  C[G>A]G on - strand
  c("GTA_C", "TAC_G"),   # G[C>T]G  =  C[G>A]C on - strand
  c("TTA_C", "TAA_G"),
  c("ATT_C", "AAT_G"),   # A[C>T]G  =  C[G>A]T on - strand
  c("CTT_C", "AAG_G"),   # C[C>T]G  =  C[G>A]G on - strand
  c("GTT_C", "AAC_G"),   # G[C>T]G  =  C[G>A]C on - strand
  c("TTT_C", "AAA_G"),
  c("ATC_C", "CAT_G"),   # A[C>T]G  =  C[G>A]T on - strand
  c("CTC_C", "CAG_G"),   # C[C>T]G  =  C[G>A]G on - strand
  c("GTC_C", "CAC_G"),   # G[C>T]G  =  C[G>A]C on - strand
  c("TTC_C", "CAA_G")# T[C>T]G  =  C[G>A]A on - strand
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
      subtitle = sprintf("Paired Wilcoxon p = %.2e", pval),
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

eps <- 1e-12

for (pair in cpg_pairs) {
  
  ctx1 <- pair[1]
  ctx2 <- pair[2]
  
  if (!all(c(ctx1, ctx2, "rep_class") %in% names(ann))) {
    warning(sprintf("Skipping %s vs %s: missing columns", ctx1, ctx2))
    next
  }
  
  df <- ann[
    rep_class %in% c("leftward", "rightward"),
    .(
      rep_class,
      plus  = get(ctx1),
      minus = get(ctx2)
    )
  ]
  
  df <- df[!is.na(plus) & !is.na(minus) & (plus + minus) > 0]
  
  # leftward:  + = leading, - = lagging
  # rightward: + = lagging, - = leading
  df[, leading := fifelse(rep_class == "leftward", plus,  minus)]
  df[, lagging := fifelse(rep_class == "leftward", minus, plus)]
  
  # within-window normalization
  df[, pair_mean := (leading + lagging) / 2]
  df <- df[pair_mean > 0]
  
  df[, leading_rel := leading / pair_mean]
  df[, lagging_rel := lagging / pair_mean]
  
  # useful effect summaries
  df[, log2_lead_lag := log2((leading + eps) / (lagging + eps))]
  med_log2 <- median(df$log2_lead_lag, na.rm = TRUE)
  med_leading_rel <- median(df$leading_rel, na.rm = TRUE)
  med_lagging_rel <- median(df$lagging_rel, na.rm = TRUE)
  
  wt <- wilcox.test(df$log2_lead_lag, mu = 0, exact = FALSE)
  pval <- wt$p.value
  p_label <- ifelse(pval < 2.2e-16, "p < 2.2e-16", sprintf("p = %.2e", pval))
  
  df_long <- melt(
    df,
    measure.vars = c("leading_rel", "lagging_rel"),
    variable.name = "template",
    value.name = "relative_rate"
  )
  
  df_long[, template := fifelse(
    template == "leading_rel",
    "leading",
    "lagging"
  )]
  
  p <- ggplot(df_long, aes(x = template, y = relative_rate, fill = template)) +
    geom_violin(trim = FALSE, alpha = 0.45) +
    geom_boxplot(width = 0.15, outlier.size = 0.25, alpha = 0.85) +
    geom_hline(yintercept = 1, linetype = "dashed") +
    scale_fill_manual(
      values = c("leading" = "firebrick", "lagging" = "steelblue")
    ) +
    labs(
      title = sprintf("Composition-aware leading vs lagging: %s vs %s", ctx1, ctx2),
      subtitle = sprintf(
        "%s; median log2(Ld/Lg) = %.3f; median relative rates: leading %.3f, lagging %.3f",
        p_label, med_log2, med_leading_rel, med_lagging_rel
      ),
      x = "Template strand",
      y = "Rate relative to within-window pair mean",
      fill = "Template"
    ) +
    theme_bw(base_size = 11)
  
  ggsave(
    file.path(output_dir, paste0(ctx1, "_vs_", ctx2, "_leading_lagging_relative_rate.pdf")),
    p,
    width = 7,
    height = 5
  )
  
  message(sprintf(
    "%s vs %s: median log2(Ld/Lg) = %.3f; median leading_rel = %.3f; median lagging_rel = %.3f; %s; n = %d",
    ctx1, ctx2, med_log2, med_leading_rel, med_lagging_rel, p_label, nrow(df)
  ))
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
ann$cpg_asym_log <- log2(
  (rowSums(mat192[, idx_cpg_t, drop = FALSE], na.rm = TRUE) + PSEUDOCOUNT) /
    (rowSums(mat192[, idx_cpg_a, drop = FALSE], na.rm = TRUE) + PSEUDOCOUNT)
)

ann$cpg_asym_log[!is.finite(ann$cpg_asym_log)] <- NA

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
      classified_frac >= 0.6 &
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
    ylab = "CpG asymmetry / scaled tracks",
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
  
  # Accumulate for global correlation
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
