# This code finds tandem and convergent gene pairs in Araport11;
# It also allows to reproduce Fig. 5E-F (plaNET-Seq metagenes of tandem and convergent pairs) and Fig. S5A-B (boxplot of plaNET-Seq signal in the second half of the gap between genes in a pair);

library(GenomicRanges)
library(ggplot2)
library (RColorBrewer)
library(reshape2)

r_dir <- "." # change to the folder where you saved the custom functions from https://github.com/Maxim-Ivanov/Leng_et_al_2019
scripts <- c("getOverlappingScores.R", "metageneMatrix.R", "drawMetagenePlot.R")
for (script in scripts) { source(file.path(r_dir, script)) }

# Load adjusted Araport11 genes:
genes_araport_adj <- readRDS("genes_araport_adj.RDS") # see https://github.com/Maxim-Ivanov/Leng_et_al_2019/03-Adjustment_of_Araport11_gene_boundaries
# Flip them to the opposite strand:
genes_flipped <- genes_araport_adj
strand(genes_flipped) <- ifelse(strand(genes_flipped) == "+", "-", "+")
# Extract nuclear protein-coding genes (for normalization of plaNET-Seq tracks):
genes_npcd <- genes_araport_adj[mcols(genes_araport_adj)$tx_type == "mRNA" & seqnames(genes_araport_adj) %in% 1:5] # nuclear protein-coding genes
genes_npcd_ext <- suppressWarnings(trim(resize(genes_npcd, width(genes_npcd) + 100, "end"))) # extend by 100 bp upstream
genes_npcd_ext <- suppressWarnings(trim(resize(genes_npcd_ext, width(genes_npcd_ext) + 500, "start")))  # extend by 500 bp downstream to capture pA peaks in plaNET-Seq data
genes_npcd_m <- reduce(genes_npcd_ext)

# Load and normalize plaNET-Seq data (merged biological replicates):
planet_dir <- "." # change to the directory containing merged plaNET-Seq Bedgraph files obtained from https://github.com/Maxim-Ivanov/Leng_et_al_2019/02-Postprocessing_of_plaNET-Seq_data.R
planet_files <- list.files(planet_dir, pattern = "merged_fw_rev.bedgraph.gz$")
planet_data <- batchReadTrackData(planet_files, dir = planet_dir, format = "bedGraph", seqinfo = seqinfo(genes_araport_adj))
names(planet_data) <- paste0("plaNET_", sub("_merged_fw_rev.bedgraph.gz", "", names(planet_data)))
planet_data <- endoapply(planet_data, normalizeGR, by = genes_npcd_m) # normalize track to 1M tags in nuclear protein-coding genes


##### PART 1: TANDEM GENE PAIRS #####

# Find gaps between tTI candidate gene pairs:
df <- data.frame("idx1" = 1:length(genes_araport_adj), "idx2" = precede(genes_araport_adj, genes_araport_adj))
df <- df[!is.na(df$idx2), ]
par1 <- genes_araport_adj[df$idx1]
par2 <- genes_araport_adj[df$idx2]
gap <- pgap(par1, par2)
mcols(gap)$g1 <- granges(par1)
mcols(gap)$g1_name <- mcols(par1)$gene_id
mcols(gap)$g2 <- granges(par2)
mcols(gap)$g2_name <- mcols(par2)$gene_id
# Choose good tTI pairs:
good <- width(gap) >= 500 & width(gap) <= 1500 & countOverlaps(gap, genes_araport_adj) == 0 & 
  mcols(par1)$tx_type == "mRNA" & mcols(par2)$tx_type == "mRNA" & seqnames(par1) %in% 1:5
gap <- gap[good]

### Draw metagene plot anchored at both upstream gene PAS and downstream gene TTS (Fig. 5E):
# Make 200 bp windows upstream and downstream from the gap:
win_up <- resize(mcols(gap)$g1, 200, "end")
win_down <- resize(mcols(gap)$g2, 200, "start")
# Compute metagene matrices:
ml1 <- lapply(planet_data, metageneMatrix, intervals = win_up, skip.zeros = FALSE, skip.top.obs = TRUE, scaling = FALSE)
ml2 <- lapply(planet_data, metageneMatrix, intervals = gap, skip.zeros = FALSE, skip.top.obs = TRUE, scaling = TRUE, matrix.length = 300)
ml3 <- lapply(planet_data, metageneMatrix, intervals = win_down, skip.zeros = FALSE, skip.top.obs = TRUE, scaling = FALSE)
# Combine matrices for the same sample:
ml <- mapply(function(x, y, z) { return(cbind(x, y, z)) }, ml1, ml2, ml3, SIMPLIFY = FALSE)
ttl <- "Metagene plot (plaNET-Seq WT vs Y732F) of tTI pairs"
drawMetagenePlot(ml, x.axis = seq(-199, 500), vline = c(0, 300), hline = 0, title = paste0(ttl, " (n=", nrow(ml[[1]]), ") "), 
                 xlabel = "Gap intervals 0.5-1.5Kb scaled to 300 bins (with 200 bp flanks)", plotPDF = TRUE, 
                 custom.colors = brewer.pal(n = 3, "Dark2"), width = 8, height = 8, units = "in")

# Assess statistical significance of the difference between A2 and A6 in the last 50% of the gap (Fig. S5A):
w05 <- resize(gap, width(gap) * 0.5, "end")
cm <- getOverlappingScores(intervals = w05, signal_grl = planet_data, value = "count_matrix") * 5
df <- melt(as.data.frame(cm[rowSums(cm) > 0, ]))
p <- ggplot(df, aes(x = variable, y = value)) + geom_boxplot(outlier.colour = NA) + coord_cartesian(ylim = c(0, 12))  + ggtitle("tTI gene pairs")
for (ext in c("png", "pdf")) { ggsave(paste("Boxplot of A2 vs A6 coverage in the gap between tTI genes (the last 50 percent)", ext, sep = "."), plot = p, 
                                      width = 7, height = 7, units = "in") }
wilcox.test(cm[, 1], cm[, 2])$p.value # 1.69587e-43

##### PART 2: CONVERGENT GENE PAIRS #####

### Find convergent gene pairs:
# Make a matrix of indexes of convergent gene pairs:
idx1 <- 1:length(genes_araport_adj)
idx2 <- precede(genes_araport_adj, genes_flipped)
mat <- cbind(idx1, idx2)
mat <- mat[!is.na(mat[, 2]), ]
# Remove duplicated pairs of indexes (x-y and y-x):
test <- mat[, 1] < mat[, 2]
c1 <- ifelse(test, mat[, 1], mat[, 2])
c2 <- ifelse(test, mat[, 2], mat[, 1])
pasted <- paste(c1, c2, sep = ":")
dupl <- duplicated(pasted)
mat2 <- mat[!dupl, ]
# Calculate distances between gene pairs:
par1 <- genes_araport_adj[mat2[, 1]]
par2 <- genes_araport_adj[mat2[, 2]]
gap <- pgap(par1, par2, ignore.strand = TRUE)
# Choose good convergent pairs:
good <- width(gap) >= 500 & width(gap) <= 1500 & countOverlaps(gap, genes_araport_adj) == 0 & 
  mcols(par1)$tx_type == "mRNA" & mcols(par2)$tx_type == "mRNA" & seqnames(gap) %in% 1:5
gap <- gap[good]

### Draw sense/antisense metagene plot anchored at both PAS (Fig. 5F):
# Make 200 bp windows upstream and downstream from the gap between PAS:
win_up <- flank(gap, 200)
win_down <- flank(gap, 200, start = FALSE)
# Compute metagene matrices (sense strand):
mls1 <- lapply(planet_data, metageneMatrix, intervals = win_up, skip.zeros = FALSE, skip.top.obs = TRUE, scaling = FALSE)
mls2 <- lapply(planet_data, metageneMatrix, intervals = gap, skip.zeros = FALSE, skip.top.obs = TRUE, scaling = TRUE, matrix.length = 300)
mls3 <- lapply(planet_data, metageneMatrix, intervals = win_down, skip.zeros = FALSE, skip.top.obs = TRUE, scaling = FALSE)
mls <- mapply(function(x, y, z) { return(cbind(x, y, z)) }, mls1, mls2, mls3, SIMPLIFY = FALSE)
# Compute metagene matrices for the antisense strand:
mla1 <- lapply(planet_data, metageneMatrix, intervals = win_up, skip.zeros = FALSE, skip.top.obs = TRUE, scaling = FALSE, antisenseMode = TRUE)
mla2 <- lapply(planet_data, metageneMatrix, intervals = gap, skip.zeros = FALSE, skip.top.obs = TRUE, scaling = TRUE, matrix.length = 300, antisenseMode = TRUE)
mla3 <- lapply(planet_data, metageneMatrix, intervals = win_down, skip.zeros = FALSE, skip.top.obs = TRUE, scaling = FALSE, antisenseMode = TRUE)
mla <- mapply(function(x, y, z) { return(cbind(x, y, z)) }, mla1, mla2, mla3, SIMPLIFY = FALSE)
mla <- lapply(mla, function(mat) { return(-mat) }) # convert antisense values to negative
# Combine sense and antisense data:
names(mla) <- paste0(names(mla), "_as")
ml <- c(mls, mla)
len <- length(planet_data)
ml <- ml[as.integer(mapply(c, seq(1, len), seq(len + 1, len * 2)))] # interleave the lists of sense and antisense matrices
my_colors <- suppressWarnings(rep(brewer.pal(n = len, "Dark2"), each = 2))
ttl <- "Metagene plot (plaNET-Seq A2 vs A6) of convergent pairs"
drawMetagenePlot(ml, x.axis = seq(-199, 500), vline = c(0, 300), hline = 0, title = paste0(ttl, " (n=", nrow(ml[[1]]), ") "),
                 xlabel = "Gap intervals 0.5-1.5Kb scaled to 300 bins (with 200 bp flanks)", plotPDF = TRUE, custom.colors = my_colors, width = 8, height = 8, units = "in")

# Boxplot of A2 vs A6 coverage in the last 50% of gap (on each side) (Fig. S5B):
w1_05 <- resize(gap, width(gap) * 0.5, "end")
w2_05 <- resize(gap, width(gap) * 0.5, "start")
strand(w2_05) <- ifelse(strand(w2_05) == "+", "-", "+")
w05 <- c(w1_05, w2_05)
cm <- getOverlappingScores(intervals = w05, signal_grl = planet_data, value = "count_matrix") * 5
df <- melt(as.data.frame(cm[rowSums(cm) > 0, ]))
p <- ggplot(df, aes(x = variable, y = value)) + geom_boxplot(outlier.colour = NA) + coord_cartesian(ylim = c(0, 15))  + ggtitle("tTI gene pairs")
for (ext in c("png", "pdf")) { ggsave(paste("Boxplot of A2 vs A6 coverage in the gap between Convergent genes (the last 50 percent)", ext, sep = "."), 
                                      plot = p, width = 7, height = 7, units = "in") }
wilcox.test(cm[, 1], cm[, 2])$p.value # 7.09949e-14
