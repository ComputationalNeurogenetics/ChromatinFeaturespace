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
library(BSgenome.Mmusculus.UCSC.mm10)
library(Seurat)
library(Repitools)
library(dendextend)

source("R/local_settings.R") # just local paths and core count

# ProcessReplicates uses data from a single sample origin

# Example 01 files can be downloaded as tar archive
# wget http://tegex.helsinki.fi/chromatin_featurespace/example_01.tar

# Sample origin #1
p<-paste(data_path,"e14vmb1_data/scATAC/E14VMB_0",sep="")
E14VMB <- ProcessReplicates(p,sample.names=c("_0"),is.cell.name="is__cell_barcode",buffer_length=512L)

# Sample origin #2
p<-paste(data_path,"e14di_data/scATAC/E14DI1_0",sep="")
E14DI <- ProcessReplicates(p,sample.names=c("_0"),is.cell.name="is__cell_barcode",buffer_length=512L)

# In the case of two or more sample origins
neuronal.merged <- merge(x = E14VMB$merged, y = c(E14DI$merged), add.cell.ids = c("_e14vmb","_e14di"))

# Gene annotations from EnsDb
annotations <- readRDS("./metadata/Ensembel_Mmusculus79_annotation.Rds")

# Change to UCSC style
seqlevelsStyle(annotations) <- 'UCSC'
genome(annotations) <- "mm10"

# Add the gene information to the object
Signac::Annotation(neuronal.merged) <- annotations

# Calculate and plot nucleosome signal
neuronal.merged <- NucleosomeSignal(object = neuronal.merged)
neuronal.merged$nucleosome_group <- ifelse(neuronal.merged$nucleosome_signal > 4, 'NS > 4', 'NS < 4')
FragmentHistogram(object = neuronal.merged, group.by = 'nucleosome_group', region = 'chr1-1-10000000')

# TSS Enrichment
neuronal.merged <- TSSEnrichment(neuronal.merged, fast = FALSE)

# Plot QC data
neuronal.merged$pct_reads_in_peaks <- neuronal.merged$peak_region_fragments / neuronal.merged$passed_filters * 100
neuronal.merged$blacklist_ratio <- neuronal.merged$blacklist_region_fragments / neuronal.merged$peak_region_fragments

VlnPlot(
  object = neuronal.merged,
  features = c('pct_reads_in_peaks', 'peak_region_fragments',
               'TSS.enrichment', 'blacklist_ratio', 'nucleosome_signal'),
  pt.size = 0.1,
  ncol = 5
)

# Filter based on QC
neuronal.merged.filt <- subset(
  x = neuronal.merged,
  subset = peak_region_fragments > 2500 &
    peak_region_fragments < 45000 &
    pct_reads_in_peaks > 55 &
    blacklist_ratio < 0.015 &
    nucleosome_signal < 3 &
    TSS.enrichment > 4.5
)

neuronal.merged.filt <- FrequencyFilter(neuronal.merged.filt)

temp <- neuronal.merged.filt[["cladepeaks_"]]
Key(temp) <- "cladepeakscount_"
neuronal.merged.filt[["cladepeakscount"]] <- temp

DefaultAssay(neuronal.merged.filt) <- 'clade_peaks'
neuronal.merged.filt <- BinarizeCounts(neuronal.merged.filt, assay = "clade_peaks")
neuronal.merged.filt <- RunTFIDF(neuronal.merged.filt)
neuronal.merged.filt <- FindTopFeatures(neuronal.merged.filt)
neuronal.merged.filt <- RunSVD(object = neuronal.merged.filt)

# Measure depth correlation with reduced dimension components
DepthCor(neuronal.merged.filt)

# By default SVD component 1 is dropped due to the correlation with read depth
hcd_01 <- RightSVD_dist(neuronal.merged.filt)

heights.to.test <- seq(4,10,by=1)
clade.stats_01 <- CladeStats(neuronal.merged.filt,hcd_01,heights=heights.to.test, cores=cores)

PlotCladeStats(clade.stats = clade.stats_01)
# h is decided manually based on the graph
h=8
PlotCladeTree(hcd=hcd_01, cut.height = h)

#qsave(neuronal.merged.filt,"examples/Example_01_filt.qs", nthreads = 6)
#neuronal.merged.filt <- qread("examples/Example_01_filt.qs", nthreads = 6)

cells_tree_cut_01 = cutree(hcd_01$hclust.cells, h=h)
lsi_cells_01 = dplyr::tibble(barcode = Cells(neuronal.merged.filt), cells_tree_cut = cells_tree_cut_01)

neuronal.merged.filt <- AddMetaData(neuronal.merged.filt, lsi_cells_01$cells_tree_cut, "clade")
clade.peaks_01_bedpe <- Signac::CallPeaks(neuronal.merged.filt, macs2.path = macs2_path, effective.genome.size = 1.87e9, group.by = "clade", combine.peaks = TRUE, name = paste("Example_01_per_clade_peaks_merged", sep=""), format="BEDPE")
hist(width(clade.peaks_01_bedpe), xlab = "Peak length (bp)", ylab = "Count", main="Histogram of peak lengths")
#qsave(clade.peaks_01_bedpe,file="examples/clade.peaks_01_bedpe.qs")
#qsave(neuronal.merged.filt, file="examples/neuronal.merged.filt.qs")
#clade.peaks_01_bedpe <- qread("examples/clade.peaks_01_bedpe.qs")

#### Code for example plots

sample.clade <- table(str_extract(string = lsi_cells_01$barcode, pattern = "_.*_") %>% str_replace_all(pattern = "_*", replacement = ""), lsi_cells_01$cells_tree_cut)
sample.clade.prop <- apply(sample.clade,2,function(x){x/sum(x)})

sample.prop.long <- rownames_to_column(as.data.frame(sample.clade.prop),var = "origin") %>% as_tibble() %>% pivot_longer(names_to = "branch", values_to = "prop", cols=2:ncol(sample.clade.prop)+1)
ggpubr::ggbarplot(sample.prop.long, x="branch", y="prop", fill="origin", color="origin", label=FALSE,position = position_dodge(0.9))


five.k.features <- StringToGRanges(rownames(neuronal.merged))
five.k.features.filt <- StringToGRanges(rownames(neuronal.merged.filt))

# Create alternating colors
num_elements <- length(five.k.features.filt)
colors <- rep(c("steelblue", "orange"), length.out = num_elements)

# Add the "color" metadata column
mcols(five.k.features.filt)$color <- colors

circos.initializeWithIdeogram(plotType = c("labels", "axis"),species="mm10",ideogram.height = convert_height(5, "mm"),axis.labels.cex = 0.8*par("cex"), labels.cex = 1.2*par("cex"),chromosome.index="chr6")
circos.genomicIdeogram(species="mm10")
circos.genomicDensity(Repitools::annoGR2DF(five.k.features), col="yellow",track.height = 0.075)
circos.genomicDensity(Repitools::annoGR2DF(five.k.features.filt), col="orange",track.height = 0.15)
circos.genomicDensity(Repitools::annoGR2DF(clade.peaks_01_bedpe), col="steelblue",track.height = 0.15)
circos.clear()


region.of.interest <- "chr6-88188000-88200000"
DefaultAssay(neuronal.merged.filt) <- "clade_peaks"
coverage_track <- CoveragePlot(neuronal.merged.filt,
                               region = region.of.interest,
                               annotation = F,
                               peaks = F,
                               group.by = "clade",
                               downsample.rate = 1,
                               label.size=8)

fivek.filt.track <- PeakPlot(object = neuronal.merged.filt, region=region.of.interest,peaks=five.k.features.filt, color=mcols(five.k.features.filt)$color)
fivek.filt.track$labels$y <- "5k filt"

clade.track <- PeakPlot(object = neuronal.merged.filt,region=region.of.interest,peaks=clade.peaks_01_bedpe, color = "forestgreen")
clade.track$labels$y <- "clade opt."

theme_set(theme_minimal(base_size = 10))

combined.p <- CombineTracks(
  plotlist = list(coverage_track,fivek.filt.track,clade.track),
  heights = c(.7,0.1,0.1),
  widths = c(10)
)

ggsave(plot=combined.p, file="Example_01_chr6.pdf")
