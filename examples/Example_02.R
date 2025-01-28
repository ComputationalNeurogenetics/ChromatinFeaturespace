library(ChromatinFeaturespace)
library(Signac)
library(parallel)
library(circlize)
library(qs)
library(fpc)
library(ComplexHeatmap)
library(valr)
library(common)
library(readr)
library(magrittr)
library(tidyverse)
library(stringr)
library(Seurat)
library(Repitools)
library(BSgenome.Hsapiens.UCSC.hg19)

source("examples/local_settings.R") # just local paths and core count

# Data used in this example are from 10x webpage
# wget https://cf.10xgenomics.com/samples/cell-atac/1.0.1/atac_v1_pbmc_10k/atac_v1_pbmc_10k_singlecell.csv
# wget https://cf.10xgenomics.com/samples/cell-atac/1.0.1/atac_v1_pbmc_10k/atac_v1_pbmc_10k_fragments.tsv.gz
# wget https://cf.10xgenomics.com/samples/cell-atac/1.0.1/atac_v1_pbmc_10k/atac_v1_pbmc_10k_fragments.tsv.gz.tbi

# ProcessReplicates uses data from a single sample origin

# Sample origin #1
p<-paste(data_path,"PBMC_example_data",sep="")
pbmc <- ProcessReplicates(p,sample.names=c("_0"), is.cell.name="is__cell_barcode", genome="hg19", genome.obj=BSgenome.Hsapiens.UCSC.hg19, buffer_length=512L, used.chromosomes = c("chr1", "chr2","chr3","chr4","chr5","chr6","chr7","chr8","chr9","chr10","chr11","chr12","chr13","chr14","chr15","chr16","chr17","chr18","chr19","chr20","chr21","chr22"))

# Use object with merged replicates (nor merging in this example)
pbmc.merged <- pbmc$merged

# Gene annotations from EnsDb
annotations <- readRDS("metadata/Ensembel_Hsapiens75_annotation.Rds")

# Change to UCSC style
seqlevelsStyle(annotations) <- 'UCSC'
genome(annotations) <- "hg19"

# Add the gene information to the object
Signac::Annotation(pbmc.merged) <- annotations

# Calculate and plot nucleosome signal
pbmc.merged <- NucleosomeSignal(object = pbmc.merged)
pbmc.merged$nucleosome_group <- ifelse(pbmc.merged$nucleosome_signal > 4, 'NS > 4', 'NS < 4')
FragmentHistogram(object = pbmc.merged, group.by = 'nucleosome_group', region = 'chr1-1-10000000')

# TSS Enrichment
pbmc.merged <- TSSEnrichment(pbmc.merged, fast = FALSE)

# Plot QC data
pbmc.merged$pct_reads_in_peaks <- pbmc.merged$peak_region_fragments / pbmc.merged$passed_filters * 100
pbmc.merged$blacklist_ratio <- pbmc.merged$blacklist_region_fragments / pbmc.merged$peak_region_fragments

VlnPlot(
  object = pbmc.merged,
  features = c('pct_reads_in_peaks', 'peak_region_fragments',
               'TSS.enrichment', 'blacklist_ratio', 'nucleosome_signal'),
  pt.size = 0.1,
  ncol = 5
)

# Filter based on QC
pbmc.merged.filt <- subset(
  x = pbmc.merged,
  subset = peak_region_fragments < 25000 &
    pct_reads_in_peaks > 50 &
    blacklist_ratio < 0.02 &
    nucleosome_signal < 3 &
    TSS.enrichment > 1
)

pbmc.merged.filt <- FrequencyFilter(pbmc.merged.filt, min.features = 0.01)

temp <- pbmc.merged.filt[["cladepeaks_"]]
Key(temp) <- "cladepeakscount_"
pbmc.merged.filt[["cladepeakscount"]] <- temp

DefaultAssay(pbmc.merged.filt) <- 'cladepeakscount'
pbmc.merged.filt <- BinarizeCounts(pbmc.merged.filt, assay = "cladepeakscount")
pbmc.merged.filt <- RunTFIDF(pbmc.merged.filt)
pbmc.merged.filt <- FindTopFeatures(pbmc.merged.filt)
pbmc.merged.filt <- RunSVD(object = pbmc.merged.filt)

#qsave(pbmc.merged.filt,"./examples/pbmc.merged.filt.qs",nthreads = cores)
#pbmc.merged.filt <- qread("./examples/pbmc.merged.filt.qs",nthreads = cores)

# Measure depth correlation with reduced dimension components
DepthCor(pbmc.merged.filt)

# By default SVD component 1 is dropped due to the correlation with read depth
hcd_02 <- RightSVD_dist(pbmc.merged.filt)

heights.to.test <- seq(4,10,by=1)
clade.stats_02 <- CladeStats(pbmc.merged.filt,hcd_02,heights=heights.to.test, cores=cores)

PlotCladeStats(clade.stats = clade.stats_02)

# h is decided based on the graph
h=7
PlotCladeTree(hcd=hcd_02, cut.height = h)

cells_tree_cut_02 = cutree(hcd_02$hclust.cells, h=h)
lsi_cells_02 = dplyr::tibble(barcode = Cells(pbmc.merged.filt), cells_tree_cut = cells_tree_cut_02)

pbmc.merged.filt <- AddMetaData(pbmc.merged.filt, lsi_cells_02$cells_tree_cut, "clade")
clade.peaks_02_bedpe <- CallPeaks(pbmc.merged.filt, macs2.path = macs2_path, group.by = "clade", combine.peaks = TRUE, name = paste("Example_02_per_clade_peaks_merged", sep=""), format="BEDPE")
hist(width(clade.peaks_02_bedpe), xlab = "Peak length (bp)", ylab = "Count", main="Histogram of peak lengths")
#qsave(clade.peaks_02_bedpe,file="examples/clade.peaks_02_bedpe.qs", nthreads = cores)
#qsave(pbmc.merged.filt, file="examples/pbmc.merged.filt.qs", nthreads = cores)
#clade.peaks_02_bedpe <- qread("examples/clade.peaks_02_bedpe.qs")

#### Code for example plots

five.k.features <- StringToGRanges(rownames(pbmc.merged))
five.k.features.filt <- StringToGRanges(rownames(pbmc.merged.filt))

# Create alternating colors
num_elements <- length(five.k.features.filt)
colors <- rep(c("steelblue", "orange"), length.out = num_elements)

# Add the "color" metadata column
mcols(five.k.features.filt)$color <- colors

circos.initializeWithIdeogram(plotType = c("labels", "axis"),species="hg19",ideogram.height = convert_height(5, "mm"),axis.labels.cex = 0.8*par("cex"), labels.cex = 1.2*par("cex"),chromosome.index="chr6")
circos.genomicIdeogram(species="hg19")
circos.genomicDensity(Repitools::annoGR2DF(five.k.features), col="yellow",track.height = 0.075)
circos.genomicDensity(Repitools::annoGR2DF(five.k.features.filt), col="orange",track.height = 0.15)
circos.genomicDensity(Repitools::annoGR2DF(clade.peaks_02_bedpe), col="steelblue",track.height = 0.15)
circos.clear()

region.of.interest <- "chr12-68010000-68012897"
DefaultAssay(pbmc.merged.filt) <- "cladepeakscount"
coverage_track <- CoveragePlot(pbmc.merged.filt,
                               region = region.of.interest,
                               annotation = F,
                               peaks = F,
                               group.by = "clade",
                               downsample.rate = 1,
                               label.size=8)

fivek.filt.track <- PeakPlot(object = pbmc.merged.filt, region=region.of.interest,peaks=five.k.features.filt,color=mcols(five.k.features.filt)$color)
fivek.filt.track$labels$y <- "5k filt"

clade.track <- PeakPlot(object = pbmc.merged.filt, region=region.of.interest,peaks=clade.peaks_02_bedpe,color = "forestgreen")
clade.track$labels$y <- "clade opt."

theme_set(theme_minimal(base_size = 10))

combined.p <- CombineTracks(
  plotlist = list(coverage_track,fivek.filt.track,clade.track),
  heights = c(.7,0.1,0.1),
  widths = c(10)
)
ggsave(plot=combined.p, file="PBMC_chr12_LYZ.pdf")
