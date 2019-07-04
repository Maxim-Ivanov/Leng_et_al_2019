# This script calculates read-through (RT) distance using plaNET-Seq, TSS-Seq and DR-Seq data;
# Technically, it uses the sliding window approach and a statistical model which takes into account the variability of plaNET-Seq tag density in both the gene of interest and the intergenic intervals;
# The RT tail was considered to extend from PAS to the position of the sliding window where the likelihood of observing at most this tag count under the "genic" model 
# for the first time drops below the likelihood of observing at least this tag count under the "intergenic" model;
# This file also contains the code to reproduce Fig. 5D (boxplot of RT distances in WT and Y732F samples);

library(rtracklayer)
library(SummarizedExperiment)
library(CAGEfightR)
library(edgeR)
library(ggplot2)

r_dir <- "." # change to the folder where you saved the custom functions from https://github.com/Maxim-Ivanov/Leng_et_al_2019
scripts <- c("batchReadTrackData.R", "getOverlappingScores.R", "parallelOverlapType.R", "computeRTtails.R")
for (script in scripts) { source(file.path(r_dir, script)) }

genes_araport_adj <- readRDS("genes_araport_adj.RDS") # this file was returned by 03-Adjustment_of_Araport11_gene_boundaries.R
genes_npcd <- genes_araport_adj[mcols(genes_araport_adj)$tx_type == "mRNA" & seqnames(genes_araport_adj) %in% 1:5]

## Load TSS clusters called from TSS-Seq data on hen2-2 samples (see https://github.com/Maxim-Ivanov/Kindgren_et_al_2019/10-Call_TSS_and_PAS_clusters.R):
tss <- rowRanges(readRDS("TSS_clusters_hen2-2.RDS"))

# Load PAS clusters called from DR-Seq data on WT samples (see https://github.com/Maxim-Ivanov/Kindgren_et_al_2019/10-Call_TSS_and_PAS_clusters.R):
pas <- rowRanges(readRDS("PAS_clusters_WT.RDS"))

# Import plaNET data (raw counts):
planet_dir <- "." # change to the directory containing merged plaNET-Seq Bedgraph files obtained from https://github.com/Maxim-Ivanov/Leng_et_al_2019/02-Postprocessing_of_plaNET-Seq_data.R
planet_files <- list.files(planet_dir, pattern = "merged_fw_rev.bedgraph.gz$")
planet_data <- batchReadTrackData(planet_files, dir = planet_dir, format = "bedGraph", seqinfo = seqinfo(genes_araport_adj))
names(planet_data) <- paste0("plaNET_", sub("_merged_fw_rev.bedgraph.gz", "", names(planet_data)))

# Quantify plaNET-Seq signal on nuclear protein-coding genes:
cov <- getOverlappingScores(genes_npcd, planet_data, value = "count_matrix")
fpkm <- round(cpm(cov) / width(genes_npcd) * 1000, 3)
colnames(fpkm) <- paste0(colnames(fpkm), "_FPKM")
mcols(genes_npcd) <- cbind(mcols(genes_npcd), as.data.frame(cov), as.data.frame(fpkm))

# Subset genes to the expressed ones (WT FPKM >= 5):
genes_expr <- genes_npcd[mcols(genes_npcd)$WT_new_FPKM >= 5] # n = 14382

# Try to find the relevant PAS clusters (if multiple PAS per gene, choose the one with the highest score):
gene_end_ext <- resize(resize(genes_expr, 200, "end"), 400, "start")
genes_wo_pas <- genes_expr[gene_end_ext %outside% pas]
hits <- findOverlaps(gene_end_ext, pas)
genes_with_pas <- genes_expr[unique(queryHits(hits))]
pas_par <- pas[subjectHits(hits)]
highest_score <- as.logical(unlist(tapply(score(pas_par), queryHits(hits), function(x) { x == max(x) }, simplify = FALSE)))
best_pas <- pas_par[highest_score]

# Find the gap between PAS of the gene of interest and TSS of the nearest downstream gene;
mcols(genes_wo_pas)$gap_start <- unname(granges(resize(genes_wo_pas, 1, "end")))
mcols(genes_with_pas)$gap_start <- unname(GRanges(seqnames = seqnames(best_pas), ranges = mcols(best_pas)$thick, strand = strand(best_pas)))
genes_expr <- sort(c(genes_with_pas, genes_wo_pas))

# Exclude genes with gap_start overlapping a TSS cluster:
genes_expr <- genes_expr[mcols(genes_expr)$gap_start %outside% tss] # n = 14363

# Find the nearest downstream TSS-Seq cluster: 
idx <- precede(mcols(genes_expr)$gap_start, tss)
not_na <- !is.na(idx)
genes_expr <- genes_expr[not_na]
tss_par <- tss[idx[not_na]]
gap <- pgap(mcols(genes_expr)$gap_start, tss_par)
mcols(genes_expr)$gap <- unname(gap)

# Skip genes with too short gaps ( < 1 Kb):
genes_chosen <- genes_expr[width(mcols(genes_expr)$gap) >= 1000] # n = 10806

# Make 100 bp windows with 10 bp offset along the whole nuclear genome:
tiles <- tileGenome(seqinfo(genes_araport_adj), tilewidth = 10, cut.last.tile.in.chrom = TRUE)
win_fw <- suppressWarnings(trim(resize(tiles, 100)))
win_fw <- win_fw[width(win_fw) == 100 & seqnames(win_fw) %in% 1:5]
strand(win_fw) <- "+"
win_rev <- win_fw
strand(win_rev) <- "-"

# Quantify plaNET-Seq tags on all windows:
chunk_size <- 1000000
chunk_starts <- seq(1, length(win_fw), by = chunk_size)
out_fw <- vector("list", length(chunk_starts))
out_rev <- vector("list", length(chunk_starts))
for (i in seq_along(chunk_starts)) {
  chunk_start <- chunk_starts[[i]]
  chunk_end <- chunk_start + chunk_size - 1
  if (chunk_end > length(win_fw)) { chunk_end <- length(win_fw) }
  win_fw_chunk <- win_fw[chunk_start:chunk_end]
  out_fw[[i]] <- getOverlappingScores(win_fw_chunk, planet_data)
  gc()
  win_rev_chunk <- win_rev[chunk_start:chunk_end]
  out_rev[[i]] <- getOverlappingScores(win_rev_chunk, planet_data)
  gc()
}
win_fw_cov <- do.call(c, out_fw)
win_rev_cov <- do.call(c, out_rev)
win_cov <- sort(c(win_fw_cov, win_rev_cov))

# Compute intergenic (control) lambda:
grohmm_all <- readRDS("GroHMM_transcripts_merged_across_samples.RDS") # all transcripts called by groHMM (see https://github.com/Maxim-Ivanov/Kindgren_et_al_2019/06-groHMM_pipeline.R)
all_tu <- c(granges(genes_araport_adj), granges(grohmm_all))
all_tu_ext <- sort(reduce(suppressWarnings(trim(resize(all_tu, width(all_tu) + 2000, "center"))))) # all TU extended by 1Kb each side
whole_genome <- GRanges(seqnames = rep(1:5, 2), ranges = IRanges(start = 1, end = seqlengths(seqinfo(genes_araport_adj))[rep(1:5, 2)]), strand = rep(c("+", "-"), each = 5))
non_genic <- setdiff(whole_genome, all_tu_ext) # 89 Mb of untranscribed sequences
win_nongenic <- subsetByOverlaps(win_cov, non_genic, type = "within")

# Compute the RT tails:
for (i in seq_along(planet_data)) {
  name <- names(planet_data)[[i]]
  message(name); flush.console()
  mean_nongenic_cov <- mean(mcols(win_nongenic)[, name])
  genes_chosen <- computeRTtails(genes = genes_chosen, cov_mcol = name, win_cov = win_cov, lambda_ctrl = mean_nongenic_cov * 5)
}

# Save the results:
saveRDS(genes_chosen, "Genes_with_RT_tails.RDS")

# Export RT tails as BED files:
for (i in seq_along(planet_data)) {
  name <- names(planet_data)[[i]]
  export(mcols(genes_chosen)[, paste0("RT_", name)], paste0("RT_", name, ".bed"), format = "BED")
}

# Fig. 5D (boxplot of RT distances in WT vs Y732F samples):
df1 <- data.frame("Readthrough" = width(mcols(genes_chosen)$RT_WT), "Sample" = "WT")
df2 <- data.frame("Readthrough" = width(mcols(genes_chosen)$RT_nrpb2-Y732F), "Sample" = "Y732F")
df <- rbind(df1, df2)
filename <- "RT distances (WT vs Y732F)"
p <- ggplot(df, aes(x = Sample, y = Readthrough)) + geom_boxplot(outlier.shape = NA) + ggtitle(filename) + ylim(0, 1500)
for (ext in c("png", "pdf")) {
  suppressMessages(ggsave(paste(filename, ext, sep = "."), plot = p, width = 6, height = 8, units = "in"))
}

# Median RT distances:
tapply(df$Readthrough, df$Sample, median) # 534 bp for WT, 649 bp for Y732F

# Statistical significance of the differences between samples:
wilcox.test(width(mcols(genes_chosen)$RT_WT), width(mcols(genes_chosen)$RT_nrpb2-Y732F))$p.value # 9.918692e-62

