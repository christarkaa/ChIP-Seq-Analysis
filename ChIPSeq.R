# Load necessary libraries
library(GenomicRanges)
library(rtracklayer)
library(Gviz)
library(biomaRt)
library(BSgenome.Hsapiens.UCSC.hg38)

# Load data
setwd("/Users/christophertarkaa/ChIPSeq/Peak_Calling")
peaks <- import.bed("H128-0-H3K27ac_peaks.bed")

# Create annotation track for the peaks
peaks.track <- AnnotationTrack(peaks, genome="hg38", name='H3K27ac Peaks', shape='box', fill='blue3', size=2)

# Define the chromosome range for visualization
chromosomes <- unique(as.character(seqnames(peaks)))
from <- min(start(peaks))
to <- max(end(peaks))

# Use biomaRt to get gene annotations
mart <- useMart("ensembl", dataset="hsapiens_gene_ensembl")
genes <- getBM(attributes=c("chromosome_name", "start_position", "end_position", "strand", "ensembl_gene_id", "external_gene_name"),
               filters="chromosome_name",
               values=gsub("chr", "", chromosomes),
               mart=mart)

# Ensure chromosome names match the format in the BSgenome object
genes$chromosome_name <- paste0("chr", genes$chromosome_name)

# Create a GRanges object for the genes
gene_ranges <- GRanges(seqnames = genes$chromosome_name,
                       ranges = IRanges(start = genes$start_position, end = genes$end_position),
                       strand = ifelse(genes$strand == 1, "+", "-"),
                       gene = genes$ensembl_gene_id,
                       symbol = genes$external_gene_name)

# Create a gene annotation track
gene.track <- GeneRegionTrack(gene_ranges, genome="hg38", name="Genes", transcriptAnnotation="symbol")

# Plot the tracks in a smaller region to ensure fit
small_region_from <- 629000
small_region_to <- 635000
plotTracks(list(peaks.track, gene.track), chromosome=chromosomes[1], from=small_region_from, to=small_region_to, transcriptAnnotation="symbol")

# Save the plot
png("annotation_tracks.png")
plotTracks(list(peaks.track, gene.track), chromosome=chromosomes[1], from=small_region_from, to=small_region_to, transcriptAnnotation="symbol")
dev.off()

# Extend regions around TSS to ±1000 bp and create tiles
genes$TSS <- ifelse(genes$strand == 1, genes$start_position, genes$end_position)
tiles <- sapply(1:nrow(genes), function(i) {
  if (genes$strand[i] == 1) {
    genes$TSS[i] + seq(-1000, 900, length.out=20)
  } else {
    genes$TSS[i] + seq(900, -1000, length.out=20)
  }
})
tiles <- GRanges(tilename=paste(rep(genes$ensembl_gene_id, each=20), 1:20, sep="_"),
                 seqnames=Rle(rep(genes$chromosome_name, each=20)),
                 ranges=IRanges(start=as.vector(tiles), width=100),
                 strand=Rle(rep("*", length(as.vector(tiles)))), seqinfo=seqinfo(BSgenome.Hsapiens.UCSC.hg38))

# Count reads mapping to each tile
H3K27ac_p <- countOverlaps(tiles, peaks)
H3K27ac_p_matrix <- matrix(H3K27ac_p, nrow=nrow(genes), ncol=20, byrow=TRUE)

# Open a new graphical device with larger dimensions
png("heatmap.png", width = 800, height = 600)

# Define color palette
colors <- colorRampPalette(c('white','red','gray','black'))(100)

# Set up a simple layout and margins for the plots
layout(matrix(c(1, 2), nrow = 1), widths = c(1, 2))
par(mar = c(5, 5, 2, 2))  # Conservative margins

# Plot the color scale
image(seq(0, max(H3K27ac_p_matrix), length.out = 100), 1, matrix(seq(0, max(H3K27ac_p_matrix), length.out = 100), 100, 1),
      col = colors, xlab = 'Number of reads', ylab = '', main = 'Color Scale', axes = FALSE)
axis(1)

# Plot the heatmap
par(mar = c(5, 5, 2, 2))  # Conservative margins
image(x = seq(-1000, 1000, length.out = 20), y = 1:nrow(H3K27ac_p_matrix),
      z = t(H3K27ac_p_matrix[order(rowSums(H3K27ac_p_matrix)), ]), col = colors,
      xlab = 'Distance from TSS (bp)', ylab = 'Promoters', lwd = 2)
abline(v = 0, lwd = 1, col = 'gray')

# Close the graphical device
dev.off()

# Open a new graphical device for the average profile
png("average_profile.png", width = 800, height = 600)

# Plot average profile
plot(x = seq(-1000, 1000, length.out = 20), y = colMeans(H3K27ac_p_matrix), type = 'b', pch = 19, col = 'red4', lwd = 2,
     ylab = 'Mean tag count', xlab = 'Distance from TSS (bp)')
abline(h = seq(1, 100, by = 5), v = seq(-1000, 1000, length.out = 20), lwd = 0.25, col = 'gray')
box(col = 'black', lwd = 2)

# Close the graphical device
dev.off()


