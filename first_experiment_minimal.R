############################################################
# CLEAN STRAND-ASYMMETRY PIPELINE
############################################################

library(data.table)
library(GenomicRanges)
library(IRanges)
library(zoo)

############################################################
# SETTINGS
############################################################

matrix_file <- "TOPMed_10kb.txt"
bed_file    <- "final.hg38.bed"
output_dir  <- "strand_specific_output"
assign_threshold <- 0.5
PSEUDOCOUNT <- 1e-6

dir.create(output_dir, showWarnings = FALSE)

############################################################
# LOAD DATA
############################################################

mat_dt <- fread(matrix_file)
bed_dt <- fread(bed_file)

# Normalize chromosome names
normalize_chr <- function(x) {
  paste0("chr", sub("^chr", "", x))
}

mat_dt[, chr := normalize_chr(chr)]
bed_dt[, chromosome := normalize_chr(chromosome)]
setnames(bed_dt, "chromosome", "chr")

############################################################
# BUILD GRanges
############################################################
mat_dt[, mat_id := .I]
mat_gr <- GRanges(
  seqnames = mat_dt$chr,
  ranges   = IRanges(mat_dt$start + 1, mat_dt$end)
)

bed_gr <- GRanges(
  seqnames = bed_dt$chr,
  ranges   = IRanges(bed_dt$start + 1, bed_dt$end)
)

mcols(bed_gr)$left  <- bed_dt$isLeftReplicating
mcols(bed_gr)$right <- bed_dt$isRightReplicating

############################################################
# OVERLAP + CLASSIFY WINDOWS
############################################################

hits <- findOverlaps(mat_gr, bed_gr)
qh <- queryHits(hits)
sh <- subjectHits(hits)

ov_bp <- width(pintersect(ranges(mat_gr)[qh], ranges(bed_gr)[sh]))

hit_dt <- data.table(
  mat_id = qh,
  bp = ov_bp,
  left  = mcols(bed_gr)$left[sh],
  right = mcols(bed_gr)$right[sh]
)

ov_sum <- hit_dt[, .(
  bp_left  = sum(bp[left == 1 & right == 0]),
  bp_right = sum(bp[right == 1 & left == 0])
), by = mat_id]

ann <- merge(mat_dt, ov_sum,
             by = "mat_id", all.x = TRUE)

ann[is.na(bp_left),  bp_left := 0]
ann[is.na(bp_right), bp_right := 0]

ann[, win_bp := end - start]
ann[, classified_frac := (bp_left + bp_right) / win_bp]

ann[, rep_class := fifelse(
  classified_frac < assign_threshold, "ambiguous",
  fifelse(bp_right > bp_left, "rightward",
          fifelse(bp_left > bp_right, "leftward", "ambiguous"))
)]

############################################################
# EXTRACT SBS192 MATRIX
############################################################

meta_cols <- c("chr","start","end","bp_left","bp_right","win_bp","classified_frac","rep_class")
sbs_cols <- setdiff(names(mat_dt), meta_cols)

mat192 <- as.matrix(mat_dt[, ..sbs_cols])

############################################################
# BUILD RC MATRIX
############################################################

rev_comp_base <- function(x) chartr("ACGT","TGCA",x)

rev_comp_seq <- function(s) {
  paste(rev(rev_comp_base(strsplit(s,"")[[1]])), collapse = "")
}

make_rc_map <- function(cols) {
  tri <- sub("_.*","",cols)
  alt <- sub(".*_","",cols)
  rc  <- paste0(rev_comp_seq(tri), "_", rev_comp_base(alt))
  match(rc, cols)
}

rc_idx <- make_rc_map(colnames(mat192))
mat192_rc <- mat192[, rc_idx]

############################################################
# === A. REPLICATION ASYMMETRY ===
############################################################

# Align to template strand
mat192_template <- mat192
mat192_template[ann$rep_class == "rightward", ] <-
  mat192_rc[ann$rep_class == "rightward", ]

# Collapse AFTER alignment
collapse_192_to_96 <- function(mat) {
  tri <- substr(colnames(mat),1,3)
  ref <- substr(tri,2,2)
  alt <- sub(".*_","",colnames(mat))
  
  keep <- ref %in% c("C","T")
  
  mat2 <- mat[, keep]
  colnames(mat2) <- paste0(substr(tri[keep],1,1),
                           "[",ref[keep],">",alt[keep],"]",
                           substr(tri[keep],3,3))
  mat2
}

mat96_template <- collapse_192_to_96(mat192_template)

# define windows
idx_leading <- which(ann$rep_class == "leftward")
idx_lagging <- which(ann$rep_class == "rightward")

############################################################
# BOXPLOTS (FULL CONTEXT, NO AGGREGATION)
############################################################

pdf(file.path(output_dir,"boxplots_SBS96_template.pdf"),
    width=18,height=12)

par(mfrow=c(8,12), mar=c(3,2,2,1))

for(i in seq_along(colnames(mat96_template))){
  
  lead_vals <- mat96_template[idx_leading,i]
  lag_vals  <- mat96_template[idx_lagging,i]
  
  lead_vals <- lead_vals[is.finite(lead_vals)]
  lag_vals  <- lag_vals[is.finite(lag_vals)]
  
  if(length(lead_vals)<10 || length(lag_vals)<10){
    plot.new(); next
  }
  
  boxplot(list(L=lead_vals,G=lag_vals),
          main=colnames(mat96_template)[i],
          col=c("red","blue"),
          outline=FALSE, xaxt="n", yaxt="n")
}

dev.off()

############################################################
# STATISTICS
############################################################

results <- lapply(seq_along(colnames(mat96_template)), function(i){
  
  lead_vals <- mat96_template[idx_leading,i]
  lag_vals  <- mat96_template[idx_lagging,i]
  
  lead_vals <- lead_vals[is.finite(lead_vals)]
  lag_vals  <- lag_vals[is.finite(lag_vals)]
  
  if(length(lead_vals)<10 || length(lag_vals)<10) return(NULL)
  
  test <- wilcox.test(lead_vals, lag_vals)
  
  data.frame(
    context = colnames(mat96_template)[i],
    p = test$p.value,
    effect = median(lead_vals) - median(lag_vals)
  )
})

results <- rbindlist(results)
results[, padj := p.adjust(p,"BH")]

fwrite(results, file.path(output_dir,"strand_stats.csv"))

############################################################
# === B. PAPER-STYLE CpG ASYMMETRY ===
############################################################

idx_tpg <- grep("CG_T$", colnames(mat192))
idx_cpa <- grep("CG_A$", colnames(mat192))

tpg_vals <- rowMeans(mat192[, idx_tpg])
cpa_vals <- rowMeans(mat192[, idx_cpa])

asym <- (tpg_vals - cpa_vals) / (tpg_vals + cpa_vals + PSEUDOCOUNT)
asym[!is.finite(asym)] <- NA

############################################################
# CHROMOSOME PLOT (paper-style)
############################################################

pdf(file.path(output_dir,"CpG_asymmetry_by_chr.pdf"),
    width=12,height=4)

for(chr_i in unique(ann$chr)){
  
  idx <- which(ann$chr == chr_i & !is.na(asym))
  if(length(idx)<50) next
  
  ord <- idx[order(ann$start[idx])]
  x <- ann$start[ord]
  y <- asym[ord]
  
  smooth <- rollmean(y,200,fill=NA)
  
  plot(x,y,pch=16,cex=0.2,col="grey70",
       main=chr_i, ylab="CpG asymmetry", xlab="position")
  lines(x,smooth,col="darkgreen",lwd=2)
  abline(h=0,lty=2)
}

dev.off()

############################################################
message("DONE")
