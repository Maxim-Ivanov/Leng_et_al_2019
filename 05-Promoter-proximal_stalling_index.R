# The script below was used to compute the promoter-proximal stalling index;
# This file also contains the code to reproduce Fig. 3D (boxplot of stalling indexes in WT vs Y732F);

# The methodology for calculating promoter-proximal RNAPII stalling index was inspired by Zhu 2018 (PMID 30374093):
# We expect to observe the TSS stalling event somewhere in the interval [TSS-100, TSS+300];
# A sliding window of fixed width (100 bp) moves in 10 bp steps along this interval;
# PlaNET-Seq coverage in WT sample is calculated for each position of the window;
# Position with the strongest signal is considered as the most representative position for TSS stalling in given gene;
# (if multiple positions had the same highest signal, then the position closest to [TSS+100] was chosen);
# Finally, stalling index is calculated for each sample as the plaNET-Seq FPKM at the TSS stalling position divided by FPKM in the coding region [TSS+300, PAS-100];

library(SummarizedExperiment)
library(GenomicFeatures)
library(edgeR)
library(rtracklayer)
library(ggplot2)
library(reshape2)

r_dir <- "." # change to the folder where you saved the custom functions from https://github.com/Maxim-Ivanov/Leng_et_al_2019
scripts <- c("getOverlappingScores.R", "windowsFromGRanges.R", "findBestWindow.R", "normalizeGR.R")
for (script in scripts) { source(file.path(r_dir, script)) }

genes_araport_adj <- readRDS("genes_araport_adj.RDS") # see https://github.com/Maxim-Ivanov/Leng_et_al_2019/03-Adjustment_of_Araport11_gene_boundaries.R
genes_npcd <- genes_araport_adj[mcols(genes_araport_adj)$tx_type == "mRNA" & seqnames(genes_araport_adj) %in% 1:5]

txdb_araport <- makeTxDbFromGFF("Araport11_GFF3_genes_transposons.201606.gff.gz") # https://www.arabidopsis.org/download_files/Genes/Araport11_genome_release/Araport11_GFF3_genes_transposons.201606.gff.gz
ebg_araport <- exonsBy(txdb_araport, by = "gene") # exons grouped by gene
seqinfo(ebg_araport, new2old = as.integer(c(1:5, 7, 6))) <- seqinfo(genes_araport_adj)

library(BSgenome.Athaliana.TAIR.TAIR9)
bsgen <- BSgenome.Athaliana.TAIR.TAIR9
seqlevels(bsgen) <- seqlevels(genes_araport_adj)

# Load plaNET-Seq data (merged biological replicates):
planet_dir <- "." # change to the directory containing merged plaNET-Seq Bedgraph files obtained from https://github.com/Maxim-Ivanov/Leng_et_al_2019/02-Postprocessing_of_plaNET-Seq_data.R
planet_files <- list.files(planet_dir, pattern = "merged_fw_rev.bedgraph.gz$")
planet_data <- batchReadTrackData(planet_files, dir = planet_dir, format = "bedGraph", seqinfo = seqinfo(genes_araport_adj))
names(planet_data) <- paste0("plaNET_", sub("_merged_fw_rev.bedgraph.gz", "", names(planet_data)))
planet_data <- endoapply(planet_data, normalizeGR, by = genes_npcd) # normalize track to 1M tags in nuclear protein-coding genes

# Calculate FPKM transcription of whole genes in the wild type sample:
cov <- as.numeric(getOverlappingScores(genes_araport_adj, planet_data[1], value = "count_matrix"))
mcols(genes_araport_adj)$fpkm_wt <- round(cov / width(genes_araport_adj) * 1000, 5)

# Subset Araport11 genes to width >= 1 Kb, WT FPKM >= 10 and no overlap with other genes:
genes_ext <- suppressWarnings(trim(resize(granges(genes_araport_adj), width(genes_araport_adj) + 100, "end")))
good <- width(genes_araport_adj) >= 1000 & mcols(genes_araport_adj)$fpkm_wt >= 10 & countOverlaps(genes_ext, genes_araport_adj) == 1
genes_good <- genes_araport_adj[good] # n = 6863

# Find the position of TSS stalling peak:
prom <- trim(resize(resize(granges(genes_good), 300, "start"), 400, "end"))
win <- windowsFromGRanges(prom, window_width = 100, window_offset = 10)
wcov <- getOverlappingScores(win, planet_data, value = "count_matrix")
mcols(win)$cov <- round(wcov[, 1], 5)
win_c <- resize(win, 1, "center")
prior <- resize(resize(granges(genes_good), 100, "start"), 1, "end") # [TSS+100] was taken as the best prior estimate for the TSS stalling position
prior_par <- prior[mcols(win)$ID]
dist <- start(win_c) - start(prior_par)
mcols(win)$dist <- ifelse(strand(win_c) == "+", dist, -dist) # negative distance: window is upstream from the gene TSS
mcols(win)$idx <- 1:length(win)
best_win_idx <- by(mcols(win)[, c("cov", "dist", "idx")], INDICES = mcols(win)$ID, FUN = findBestWindow, simplify = FALSE) %>% unlist()
w1 <- win[best_win_idx]

# PlaNET-Seq coverage of the TSS stalling peak:
cov1 <- getOverlappingScores(w1, planet_data, value = "count_matrix")

# PlaNET-Seq coverage in coding regions of respective genes:
w2 <- resize(granges(genes_good), width(genes_good) - 300, "end")
w2 <- resize(w2, width(w2) - 100, "start")
cov2 <- getOverlappingScores(w2, planet_data, value = "count_matrix")
cov2_norm <- cov2 / width(w2) * 100 # normalize coverage to 100 bp of gene width

# Compute TSS stalling index:
tss_si <- cov1 / cov2_norm
tss_si[is.nan(tss_si)] <- 0
tss_si[is.infinite(tss_si)] <- 0
colnames(tss_si) <- sub("plaNET", "TSS_SI", colnames(tss_si))
tss_si <- round(tss_si, 3)

# Fig. 3D (boxplot of stalling indexes in WT vs Y732F):
df <- melt(as.data.frame(tss_si))
ttl <- "Promoter-proximal stalling index in WT and Y732F"
p <- ggplot(df, aes(x = variable, y = value)) + geom_boxplot(outlier.colour = NA) + coord_cartesian(ylim = c(0, 6)) + ylab("RNAPII stalling index") + ggtitle(ttl)
for (ext in c("png", "pdf")) {
  suppressMessages(ggsave(paste(ttl, ext, sep = "."), plot = p, width = 6, height = 8, units = "in"))
}

# Statistical significance of the difference:
wilcox.test(tss_si[, 1], tss_si[, 2])$p.value # 1.127924e-268
