library(Rsamtools)
library(GenomicAlignments)
library(GenomicFeatures)
library(DESeq2)
library(DEXSeq)

genes_araport_adj <- readRDS("genes_araport_adj.RDS") # see https://github.com/Maxim-Ivanov/Leng_et_al_2019/03-Adjustment_of_Araport11_gene_boundaries.R
genes_npcd <- genes_araport_adj[mcols(genes_araport_adj)$tx_type == "mRNA" & seqnames(genes_araport_adj) %in% 1:5]

txdb_araport <- makeTxDbFromGFF("Araport11_GFF3_genes_transposons.201606.gff.gz") # https://www.arabidopsis.org/download_files/Genes/Araport11_genome_release/Araport11_GFF3_genes_transposons.201606.gff.gz
ebg <- exonsBy(txdb_araport, by = "gene")
seqinfo(ebg, new2old = as.integer(c(1:5, 7, 6))) <- seqinfo(genes_araport_adj)
ebg_npcd <- ebg[names(ebg) %in% names(genes_npcd)]
ebg_npcd_m <- reduce(ebg_npcd) # non-overlapping exonic intervals within each gene

##### PART 1: call DE genes on RNA-Seq data using DESeq2 #####

# Load BAM files from unstranded RNA-Seq:
bamdir <- "." # change to the directory with BAM files produced by https://github.com/Maxim-Ivanov/Leng_et_al_2019/06-Alignment_of_RNA-Seq_data.sh
bamfiles <- BamFileList(list.files(path = bamdir, pattern="bam$", full.names = TRUE), yieldSize = 2000000)
# Count RNA-Seq reads on exons grouped by genes:
se_genes <- summarizeOverlaps(features = ebg_npcd_m, reads = bamfiles, inter.feature = FALSE, singleEnd = FALSE, ignore.strand = TRUE)
colData(se_genes)$genotype <- factor(rep(c("WT", "Y732F"), each = 2), levels = c("WT", "Y732F"))
# Skip features with low coverage:
se_genes <- se_genes[rowSums(assay(se_genes)) >= 10, ] # n = 22928
# Modify SE containers to DESeq input format:
dds_genes <- DESeqDataSet(se_genes, design = ~ genotype)
# Run DESeq2:
dds_genes <- DESeq(dds_genes)
# Process the results:
df <- as.data.frame(results(dds_genes, alpha = 0.05))
input_genes <- unlist(range(rowRanges(se_genes)))
mcols(input_genes)$gene_id <- names(input_genes)
names(input_genes) <- NULL
mcols(input_genes)$padj <- df$padj
mcols(input_genes)$lfc <- df$log2FoldChange
de_idx <- which(mcols(input_genes)$padj <= 0.05 & abs(mcols(input_genes)$lfc) >= 1)
de_lfc <- mcols(input_genes)$lfc[de_idx]
mcols(input_genes)$de <- "No"
mcols(input_genes)$de[de_idx] <- ifelse(de_lfc > 0, "Up", "Down")
mcols(input_genes) <- cbind(mcols(input_genes), as.data.frame(assay(dds_genes)))
names(mcols(input_genes)) <- c("padj", "lfc", "de", "WT_rep1", "WT_rep2", "Y732F_rep1", "Y732F_rep2")
# Numbers of DE genes:
table(mcols(input_genes)$de) # 1005 down-, 1646 upregulated
# Export the results (Table S2):
write.table(as.data.frame(input_genes[mcols(input_genes)$de %in% c("Up", "Down")]), "DE_genes.txt", sep = "\t", quote = FALSE, row.names = FALSE)

##### PART 2: call DE introns on RNA-Seq data using DEXSeq #####

# Prepare exons for DEXSeq:
exonic <- exonicParts(txdb3, linked.to.single.gene.only = TRUE) # extract disjoint exonic intervals
seqinfo(exonic) <- seqinfo(genes3)
exonic <- exonic[mcols(exonic)$gene_id %in% names(genes3_npcd)] # subset to nuclear protein-coding
mcols(exonic) <- mcols(exonic)[, c("gene_id", "tx_name")] # skip unnecessary columns
exonic_grl <- split(exonic, mcols(exonic)$gene_id)
exon_counts <- lengths(exonic_grl) # enumerate exons
exon_enum_fw <- as.vector(unlist(lapply(exon_counts, function(x) { seq(1, x) })))
exon_enum_rev <- as.vector(unlist(lapply(exon_counts, function(x) { seq(x, 1) })))
exons <- unlist(exonic_grl)
mcols(exons)$type <- "exon"
mcols(exons)$id <- paste(mcols(exons)$gene_id, "exon", ifelse(strand(exons) == "+", exon_enum_fw, exon_enum_rev), sep = "_")
# Do the same for introns:
intronic <- intronicParts(txdb3, linked.to.single.gene.only = TRUE)
seqinfo(intronic) <- seqinfo(genes3)
intronic <- intronic[mcols(intronic)$gene_id %in% names(genes3_npcd)]
intronic <- intronic[width(intronic) >= 25]
intronic <- resize(intronic, width(intronic) - 10, "center") # shrink introns by 5 bp each side to avoid edge effects
mcols(intronic) <- mcols(intronic)[, c("gene_id", "tx_name")]
intronic_grl <- split(intronic, mcols(intronic)$gene_id)
intron_counts <- lengths(intronic_grl) # enumerate exons
intron_enum_fw <- as.vector(unlist(lapply(intron_counts, function(x) { seq(1, x) })))
intron_enum_rev <- as.vector(unlist(lapply(intron_counts, function(x) { seq(x, 1) })))
introns <- unlist(intronic_grl)
mcols(introns)$type <- "intron"
mcols(introns)$id <- paste(mcols(introns)$gene_id, "intron", ifelse(strand(introns) == "+", intron_enum_fw, intron_enum_rev), sep = "_")
# Interleave exons and introns:
ex_intr <- c(exons, introns)
ex_intr <- ex_intr[order(mcols(ex_intr)$gene_id, start(ex_intr))] # sort by gene name, then by coordinates
ex_intr_grl <- split(ex_intr, mcols(ex_intr)$gene_id)
# Add exonic_part column (required by DEXSeqDataSetFromSE()):
mcols(ex_intr)$exonic_part <- unlist(lapply(lengths(ex_intr_grl), function(x) { seq(1, x) }))
names(ex_intr) <- NULL
# Count reads on exons and introns:
se_ex_intr <- summarizeOverlaps(features = ex_intr, reads = bamfiles, inter.feature = FALSE, singleEnd = FALSE, ignore.strand = TRUE)
colData(se_ex_intr)$genotype <- factor(rep(c("WT", "Y732F"), each = 2), levels = c("WT", "Y732F"))
# Skip features with low coverage:
table(mcols(se_ex_intr)$type) # 196885 exons, 130188 introns
se_ex_intr <- se_ex_intr[rowSums(assay(se_ex_intr)) >= 10, ]
table(mcols(se_ex_intr)$type) # 149365 exons, 39008 introns
# Import to DEXSeq:
dxd <- DEXSeqDataSetFromSE(se_ex_intr, design = ~ sample + exon + genotype:exon)
# Run DEXSeq:
dxr <- DEXSeq(dxd, fitExpToVar = "genotype")
# Process the results:
input_ex_intr <- rowRanges(se_ex_intr)
mcols(input_ex_intr) <- mcols(input_ex_intr)[, c("gene_id", "type", "id")]
de_idx2 <- which(dxr$padj <= 0.05 & abs(dxr$log2fold_A6_A2) >= 1)
de_lfc2 <- dxr$log2fold_A6_A2[de_idx2]
input_ex_intr$padj <- dxr$padj
input_ex_intr$lfc <- dxr$log2fold_A6_A2
input_ex_intr$de <- "No"
input_ex_intr$de[de_idx2] <- ifelse(de_lfc2 > 0, "Up", "Down")
mcols(input_ex_intr) <- cbind(mcols(input_ex_intr), as.data.frame(assay(se_ex_intr)))
names(mcols(input_ex_intr)) <- sub("_sorted.bam", "", names(mcols(input_ex_intr)))
# Numbers of DE exons and introns:
table(mcols(input_ex_intr)$type, mcols(input_ex_intr)$de) # exons: 481 down, 405 up; introns: 1334 down, 183 up
