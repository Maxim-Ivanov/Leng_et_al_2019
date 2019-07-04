# This script allows to reproduce the metagene plots of plaNET-Seq signal in WT vs Y732F (Fig. 3C, 4B-D, 5C, S3C, S4A-B);

library(rtracklayer)
library(ggplot2)

# Load adjusted Araport11 genes:
genes_araport_adj <- readRDS("genes_araport_adj.RDS") # see https://github.com/Maxim-Ivanov/Leng_et_al_2019/03-Adjustment_of_Araport11_gene_boundaries

# Extract nuclear protein-coding genes (for normalization of plaNET-Seq tracks):
genes_npcd <- genes_araport_adj[mcols(genes_araport_adj)$tx_type == "mRNA" & seqnames(genes_araport_adj) %in% 1:5] # nuclear protein-coding genes
genes_npcd_ext <- suppressWarnings(trim(resize(genes_npcd, width(genes_npcd) + 100, "end"))) # extend by 100 bp upstream
genes_npcd_ext <- suppressWarnings(trim(resize(genes_npcd_ext, width(genes_npcd_ext) + 500, "start")))  # extend by 500 bp downstream to capture pA peaks in plaNET-Seq data
genes_npcd_m <- reduce(genes_npcd_ext) # merge overlapping intervals

# Load the original Araport11 annotation:
txdb_araport <- makeTxDbFromGFF("Araport11_GFF3_genes_transposons.201606.gff.gz") # https://www.arabidopsis.org/download_files/Genes/Araport11_genome_release/Araport11_GFF3_genes_transposons.201606.gff.gz
ebg_araport <- exonsBy(txdb_araport, by = "gene") # exons grouped by gene
seqinfo(ebg_araport, new2old = as.integer(c(1:5, 7, 6))) <- seqinfo(genes_araport_adj)
exons <- unlist(ebg_araport)
ebg_npcd <- ebg_araport[names(ebg_araport) %in% names(genes_npcd)] # exons in nuclear protein-coding genes

# Load custom functions:
r_dir <- "." # change to the folder where you saved the custom functions from https://github.com/Maxim-Ivanov/Leng_et_al_2019
scripts <- c("batchReadTrackData.R", "metageneMatrix.R", "drawMetagenePlot.R", "normalizeGR.R", "removeFirstAndLastExons.R", "findFirstNucleosome.R")
for (script in scripts) { source(file.path(r_dir, script)) }

# Load and normalize plaNET-Seq data (merged biological replicates):
planet_dir <- "." # change to the directory containing merged plaNET-Seq Bedgraph files obtained from https://github.com/Maxim-Ivanov/Leng_et_al_2019/02-Postprocessing_of_plaNET-Seq_data.R
planet_files <- list.files(planet_dir, pattern = "merged_fw_rev.bedgraph.gz$")
planet_data <- batchReadTrackData(planet_files, dir = planet_dir, format = "bedGraph", seqinfo = seqinfo(genes_araport_adj))
names(planet_data) <- paste0("plaNET_", sub("_merged_fw_rev.bedgraph.gz", "", names(planet_data)))
planet_data <- endoapply(planet_data, normalizeGR, by = genes_npcd_m) # normalize track to 1M tags in nuclear protein-coding genes

# Load nucleosome positions (PlantDHS):
nps <- read.table("Ath_leaf_NPS.gff", sep = "\t", header = FALSE) # download from http://plantdhs.org/static/download/Ath_leaf_NPS.gff.gz
nps <- GRanges(seqnames=sub("Chr", "", nps$V1), ranges=IRanges(nps$V4, end=nps$V5), seqinfo=seqinfo(genes_araport_adj))


### Make genomic intervals ###

# [TSS-500bp, TSS+500bp]:
tss_500 <- suppressWarnings(trim(resize(resize(genes_npcd, 0, "start"), 1000, "center")))
tss_500_filt <- tss_500[width(tss_500) == 1000 & countOverlaps(tss_500, genes_araport_adj) == 1] # skip windows which overlap with other annotated genes
# [PAS-500bp, PAS+500bp]:
pas_500 <- suppressWarnings(trim(resize(resize(genes_npcd, 0, "end"), 1000, "center")))
pas_500_filt <- pas_500[width(pas_500) == 1000 & countOverlaps(pas_500, genes_araport_adj) == 1]
# Internal exons (without the first and the last continuous exonic intervals):
int_exons <- removeFirstAndLastExons(ebg_npcd)[[3]]
int_exons_filt <- int_exons[countOverlaps(int_exons, genes_araport_adj) == 1 & countOverlaps(int_exons, exons) == 1]
int_exons_50_300 <- int_exons_filt[width(int_exons_filt) >= 50 & width(int_exons_filt) <= 300] # subset exons to width from 50 bp to 300 bp
int_exons_50_300 <- resize(int_exons_50_300, width(int_exons_50_300) - 10, "center") # trim exons by 5 bp each side to avoid possible edge effects
# Introns:
ibg_araport <- psetdiff(unlist(range(ebg_npcd)), ebg_npcd) # intronic intervals grouped by gene
introns <- unlist(ibg_araport)
introns_filt <- introns[countOverlaps(introns, genes_araport_adj) == 1]
introns_50_300 <- introns_filt[width(introns_filt) >= 50 & width(introns_filt) <= 300] # subset introns to width from 50 bp to 300 bp
introns_50_300 <- resize(introns_50_300, width(introns_50_300) - 10, "center") # trim introns by 5 bp each side
# Centers of first nucleosomes:
nps_first <- findFirstNucleosome(nps, genes_npcd) # for each gene, find the first nucleosome (within 500 bp downstream from the annotated TSS)
nps_first_500 <- suppressWarnings(trim(resize(nps_first, 1000, "center"))) # make 1 Kb windows around centers of the first nucleosomes
nps_first_500_filt <- nps_first_500[width(nps_first_500) == 1000 & countOverlaps(nps_first_500, genes_araport_adjust) == 1]
# Whole genes (500 bp to 5 Kb):
whole_genes_5kb <- genes_npcd[countOverlaps(genes_npcd, genes_araport_adj) == 1 & width(genes_npcd) <= 5000 & width(genes_npcd) >= 500]

### Draw metagene plots of plaNET-Seq signal ###

all_intervals <- list(list(tss_500_filt, "TSS 500 bp", TRUE, 200, "Window (1 Kb) centered at TSS (5 bp bins)", TRUE), # Fig. S3C
                      list(pas_500_filt, "PAS 500 bp", TRUE, 200, "Window (1 Kb) centered at PAS (5 bp bins)", TRUE), # Fig. 5C
                      list(nps_first_500_filt, "First nucl 500 bp", TRUE, 200, "Window (1 Kb) centered at first nucleosome (5 bp bins)", TRUE), # Fig. 3C
                      list(int_exons_50_300, "Exons 50bp to 300bp", TRUE, 100,  "Exons 50-300 bp scaled to 100 bins", FALSE), # Fig. 4C
                      list(int_exons_50_300, "Exons 50bp to 300bp unscaled", FALSE, 150, "Exons 50-300 bp unscaled (trimmed to 150 bp) anchored at start", FALSE), # Fig. S4A
                      list(introns_50_300, "Introns 50bp to 300bp", TRUE, 100, "Introns 50-300 bp scaled to 100 bins", FALSE), # Fig. 4D
                      list(introns_50_300, "Introns 50bp to 300bp unscaled", FALSE, 100, "Introns 50-300 bp unscaled (trimmed to 100 bp) anchored at start", FALSE), # Fig. S4B
                      list(whole_genes_5kb, "Whole genes 5Kb", TRUE, 500, "Whole genes 0.5-5 Kb scaled (500 bins)", FALSE)) # Fig. 4B

for (i in seq_along(all_data)) {
  curr_data <- all_data[[i]]
  data <- eval(parse(text = curr_data[[1]])) # convert string to variable name
  ttl2 <- curr_data[[2]] # part 2 of the title
  message(ttl2); flush.console()
}

for (i in seq_along(all_intervals)) {
  curr_int <- all_intervals[[i]]
  win <- curr_int[[1]] # current genomic windows
  ttl <- curr_int[[2]] # plot title
  sc <- curr_int[[3]] # scaling of genomic windows allowed (TRUE/FALSE)
  mlen <- curr_int[[4]] # ncol of the metagene matrix
  xlab <- curr_int[[5]] # label for the X axis
  vl <- curr_int[[6]] # draw vertical line at 0 (TRUE/FALSE)
  matlist <- lapply(planet_data, metageneMatrix, intervals = win, skip.zeros = FALSE, scaling = sc, matrix.length = mlen)
  if (isTRUE(vl)) {
    ml <- ncol(matlist[[1]])
    x.axis <- seq(-(ml/2-1), ml/2)
    vline <- 0
  } else {
    x.axis <- FALSE
    vline <- FALSE
  }
  drawMetagenePlot(matlist, x.axis = x.axis, vline = vline, title = paste0(ttl, " (n=", nrow(matlist[[1]]), ")"), xlabel = xlab,
                   ylim = c(0, NA), width = 8, height = 8, units = "in")
}
