#' @import ggplot2 dplyr Seurat Signac dendextend fpc proxy rlang magrittr GenomeInfoDb
#' @importFrom stats hclust as.hclust as.dendrogram
#' @importFrom graphics abline par
#' @importFrom utils read.csv
#' @inheritParams documentation_helpers
NULL

#' CladeStats
#'
#'
#' @param data.obj Seurat object containing scATACseq data.
#' @param rightSVD.cells List containing right singular values and clustering information.
#' @param heights Numeric vector of heights for dendrogram cuts.
#' @param cores Integer, number of CPU cores for parallel processing.
#' @return A list containing cluster statistics calculated for each height.
#' @examples
#' \dontrun{
#' heights <- seq(0.01, 1, by=0.01)
#' clade.stats <- CladeStats(data.obj, RightSVD_dist(data.obj), heights, cores)
#' print(clade.stats)
#' }
#' @export
CladeStats <- function(data.obj, rightSVD.cells, heights, cores){
  hclust.obj <- rightSVD.cells$hclust.cells
  cell.dist <- rightSVD.cells$cell.dist
  if (cores>1){
    cluster.stats.list <- parallel::mclapply(heights,function(h){
      cells_tree_cut = dendextend::cutree(hclust.obj, h=h)
      lsi_cells = dplyr::tibble(barcode = Cells(data.obj), cells_tree_cut = cells_tree_cut)
      tmp.stats <- fpc::cluster.stats(d=cell.dist, clustering=lsi_cells$cells_tree_cut)
      return(tmp.stats)
    },mc.cores=cores)
  } else if (cores==1){
    cluster.stats.list <- lapply(heights,function(h){
      cells_tree_cut = dendextend::cutree(hclust.obj, h=h)
      lsi_cells = dplyr::tibble(barcode = Cells(data.obj), cells_tree_cut = cells_tree_cut)
      tmp.stats <- fpc::cluster.stats(d=cell.dist, clustering=lsi_cells$cells_tree_cut)
      return(tmp.stats)
    })
  }
  return(cluster.stats.list)
}


#' FrequencyFilter
#'
#' Removes cells with too few or too many features and features found in too few or too many cells.
#' @param data.object A \code{Seurat} object including a "cladepeaks_" \code{\link[Signac]{ChromatinAssay}}.
#' @param min.cells A threshold for the smallest amount of cells per feature. Given as a percentage of all cells in "data.object".
#' @param max.cells A threshold for the largest amount of cells per feature. Given as a percentage of all cells in "data.object".
#' @param min.features A threshold for the smallest amount of features per cell. Given as a percentage of all features in "data.object".
#' @param max.features A threshold for the largest amount of features per cell. Given as a percentage of all features in "data.object".
#' @return Returns a \code{Seurat} object.
#' @examples
#' \dontrun{
#' data.obj <- CreateSeuratObject(matrix(rnorm(100), nrow=10))
#' filtered_obj <- FrequencyFilter(data.obj, min.cells=0.025, max.cells=0.975,\
#'  min.features=0.025, max.features=0.975)
#' print(filtered_obj)
#' }
#' @export
FrequencyFilter <- function(data.object, min.cells=0.025, max.cells=0.975, min.features=0.025, max.features=0.975){
  min.cells.thr <- round(length(unique(Signac::Cells(data.object)))*min.cells, digits = 0)
  max.cells.thr <- round(length(unique(Signac::Cells(data.object)))*max.cells, digits = 0)

  min.features.thr <- round(nrow(data.object)*min.features, digits = 0)
  max.features.thr <- round(nrow(data.object)*max.features, digits = 0)

  cell.per.feature <- SeuratObject::rowSums(Seurat::GetAssayData(data.object[["cladepeaks_"]])>0)
  feature.per.cell <- SeuratObject::colSums(Seurat::GetAssayData(data.object[["cladepeaks_"]])>0)

  feature.filter.index <- which(cell.per.feature >= min.cells.thr & cell.per.feature <= max.cells.thr)
  cell.filter.index <- which(feature.per.cell >= min.features.thr & feature.per.cell <= max.features.thr)

  features.to.keep <- rownames(data.object)[feature.filter.index]
  cells.to.keep <- Signac::Cells(data.object)[cell.filter.index]

  data.object <- subset(data.object, cells = cells.to.keep, features = features.to.keep)
  return(data.object)
}

#' CreateChromAssay
#'
#' Creates a \code{\link[Signac]{ChromatinAssay}} object for every replicate on a list.
#' @param data.list A list of lists for each replicate. The sublists need to include list components "feature.matrix" and "fragmentObject".
#' @param genome UCSC name of the used genome.
#' @return Returns a list
#' @details Part of the wrapper function \code{\link{ProcessReplicates}}. Based on the Signac function \code{\link[Signac]{CreateChromatinAssay}}.
#' @examples
#' \dontrun{
#' data.list <- list(list(feature.matrix=matrix(rpois(100, lambda=5), nrow=10),\
#'  fragmentObject="fragment_path"))
#' chrom_assay <- CreateChromAssay(data.list, genome="mm10")
#' print(chrom_assay)
#' }
#' @export
CreateChromAssay <- function(data.list, genome="mm10"){
  data.list <- lapply(data.list, function(p){
    p$chr.assay <- Signac::CreateChromatinAssay(counts = p$feature.matrix, fragments = list(p$fragmentObject), genome = GenomeInfoDb::Seqinfo(genome=genome), ranges = Signac::StringToGRanges(rownames(p$feature.matrix)))
    return(p)
  })
  return(data.list)
}

#' CreateFeatureMatrix
#'
#' Constructs a bin x cell matrix for each replicate on a list.
#' @param data.list A list of lists for each replicate. The sublists need to have list components "fragmentObject" and "genome.lengths".
#' @param binsize Size of the genome bins to use
#' @param process_n Number of regions to load into memory at a time, per thread. Processing more regions at once can be faster but uses more memory.
#' @return Returns a list
#' @details Part of the wrapper function \code{\link{ProcessReplicates}}. Based on the Signac function \code{\link[Signac]{GenomeBinMatrix}}.
#' @examples
#' \dontrun{
#' data.list <- list(list(fragmentObject="fragment_path", genome.lengths=c(chr1=100000)))
#' feature_matrix <- CreateFeatureMatrix(data.list, binsize=5000, process_n=5000)
#' print(feature_matrix)
#' }
#' @export
CreateFeatureMatrix <- function(data.list, binsize=5000, process_n=5000){
  data.list <- lapply(data.list, function(p){
    p$feature.matrix <- Signac::GenomeBinMatrix(fragments = list(p$fragmentObject), binsize = binsize, genome = p$genome.lengths, process_n = process_n)
    return(p)
  })
  return(data.list)
}

#' CreateFragmentObjects
#'
#' Creates a \code{Fragment} object for every replicate on a list.
#' @param data.list A list with list components "filtered.fragment.name" (a path to a fragment file) and "true.cell.barcodes" (the barcodes of the fragments).
#' @param validate.fragments Check that expected cells are present in the fragment file.
#' @return Returns a list
#' @details Part of the wrapper function \code{\link{ProcessReplicates}}. Based on the Signac function \code{\link[Signac]{CreateFragmentObject}}.
#' @examples
#' \dontrun{
#' data.list <- list(
#'   filtered.fragment.name = "path/to/fragment.tsv.gz",
#'   true.cell.barcodes = c("cell1", "cell2")
#' )
#' fragment_objects <- CreateFragmentObjects(data.list, validate.fragments = TRUE)
#' print(fragment_objects)
#' }
#' @export
CreateFragmentObjects <- function(data.list, validate.fragments = FALSE){
  data.list <- lapply(data.list, function(p){
    p$fragmentObject <- Signac::CreateFragmentObject(path = p$filtered.fragment.name, cells = p$true.cell.barcodes, validate.fragments = validate.fragments)
    return(p)
  })
  return(data.list)
}

#' CreateSeuratObj
#'
#' Creates a \code{Seurat} object for each replicate on a list from the data on that list.
#' @param data.list A list with list components "barcode.metadata", "chr.assay" and "true.cell.barcodes".
#' @return Returns a list
#' @details Part of the wrapper function \code{\link{ProcessReplicates}}. Based on the SeuratObject function \code{\link[SeuratObject]{CreateSeuratObject}}.
#' @examples
#' \dontrun{
#' data.list <- list(
#'   barcode.metadata = data.frame(barcode = c("cell1", "cell2"), count = c(100, 200)),
#'   chr.assay = list(),
#'   true.cell.barcodes = c("cell1", "cell2")
#' )
#' seurat_obj <- CreateSeuratObj(data.list)
#' print(seurat_obj)
#' }
#' @export
CreateSeuratObj <- function(data.list){
  data.files <- lapply(data.list, function(p){
    tmp.metadata <- as.data.frame(dplyr::filter(p$barcode.metadata, barcode %in% p$true.cell.barcodes))
    rownames(tmp.metadata) <- tmp.metadata$barcode

    p$seurat.obj <- Seurat::CreateSeuratObject(
      counts =  p$chr.assay,
      assay = 'cladepeaks_',
      project = 'ATAC',
      meta.data = tmp.metadata
    )

    return(p)
  })
  return(data.files)
}


#' FilterFragments
#'
#' Takes a list of replicates each with a fragments file path and a list component "true.cell.barcodes". Removes fragments from the fragment file if their barcodes are not included in "true.cell.barcodes". The output is a gzip-compressed and indexed file.
#' @param data.list A list of replicates each with a fragments file path and a list component "true.cell.barcodes".
#' @param force If TRUE recreates an output file even if it already exists.
#' @param buffer_length Size of buffer to be read from the fragment file. This must be longer than the longest line in the file.
#' @return Returns a list with a path to a fragment file.
#' @details Part of the wrapper function \code{\link{ProcessReplicates}}. Based on the Signac function \code{\link[Signac]{FilterCells}}.
#' @examples
#' \dontrun{
#' filtered_fragments <- FilterFragments(fragments)
#' print(filtered_fragments)
#' }
#' @export
FilterFragments <- function(data.list, force=FALSE,buffer_length=256L){
  data.files <- lapply(data.list, function(p){
    p$filtered.fragment.name <- paste(dirname(p$fragment.file),"/",stringr::str_replace(basename(p$fragment.file),pattern = ".tsv.gz",replacement = ".filtered.tsv.gz"),sep="")
    if (!file.exists(p$filtered.fragment.name) | force){
      Signac::FilterCells(
        fragments = p$fragment.file,
        cells = p$true.cell.barcodes,
        outfile = p$filtered.fragment.name,
        buffer_length = buffer_length,
        verbose = TRUE
      )
    } else {
      print("Filtered fragments file exists, recreate with force=TRUE")
    }
    return(p)
  })
  return(data.files)
}


#' FindDataFiles
#'
#' Finds the relevant scATAC data and barcode files and writes their paths in a list. This is done for each replicate, so the result is a list of lists.
#' @param data.paths Path to the folder with scATAC data including fragments and barcodes.
#' @param sample.names Names for each replicate.
#' @param barcode.file The file name or pattern of the file with cell barcodes.
#' @return Returns a list
#' @details Part of the wrapper function \code{\link{ProcessReplicates}}.
#' @examples
#' \dontrun{
#' data.paths <- c("path/to/sample1", "path/to/sample2")
#' sample.names <- c("sample1", "sample2")
#' barcode.file <- "singlecell.csv"
#' data_files <- FindDataFiles(data.paths, sample.names, barcode.file)
#' print(data_files)
#' }
#' @export
FindDataFiles <- function(data.paths, sample.names, barcode.file = "*singlecell.csv"){

  data.files <- lapply(data.paths, function(p){
    barcode.file <- common::file.find(path=p, pattern = barcode.file)
    fragment.file <- common::file.find(path=p, pattern = "*fragments.tsv.gz")
    return(list(barcode.file=barcode.file, fragment.file=fragment.file))
  })

  names(data.files) <- sample.names

  return(data.files)
}

#' MergeReplicates
#'
#' Merges the \code{Seurat} objects (one per replicate) on a list and adds the result as a new list component.
#' @param data.list A list of lists for each replicate. The sublists need to have a list component with the name "seurat.obj".
#' @param sample.names Names for each replicate.
#' @return Returns a list
#' @details Part of the wrapper function \code{\link{ProcessReplicates}}.
#' @examples
#' \dontrun{
#' data.list <- list(
#'   list(seurat.obj = CreateSeuratObject(matrix(rnorm(100), nrow=10))),
#'   list(seurat.obj = CreateSeuratObject(matrix(rnorm(100), nrow=10)))
#' )
#' sample.names <- c("replicate1", "replicate2")
#' merged_data <- MergeReplicates(data.list, sample.names)
#' print(merged_data)
#' }
#' @export
MergeReplicates <- function(data.list, sample.names){
  if (length(data.list)>1){
    data.list$merged <- merge(x = data.list[[1]]$seurat.obj, y = lapply(2:length(data.list),function(d){data.list[[d]]$seurat.obj}), add.cell.ids = sample.names)
  } else {
    data.list$merged <- data.list[[1]]$seurat.obj
  }
  return(data.list)
}

#' ProcessReplicates
#'
#' A wrapper function that creates a list of processed genomic and scATAC information for each replicate, thus the output is a list of lists.
#' @param data.paths Path to the folder with scATAC data including fragments and barcodes.
#' @param sample.names Names for each replicate.
#' @param barcode.file The file name or pattern of the file with cell barcodes.
#' @param is.cell.name A column in the barcode file that tells if the barcode in question is a cell barcode (1) or not (0).
#' @param force If TRUE recreates an output file for filtered fragments even if it already exists.
#' @param validate.fragments Check that expected cells are present in the fragment file.
#' @param genome UCSC name of the used genome.
#' @param genome.obj An object containing sequence information.
#' @param used.chromosomes The chromosomes used in the analysis.
#' @param binsize Size of the genome bins to use
#' @param process_n Number of regions to load into memory at a time, per thread. Processing more regions at once can be faster but uses more memory.
#' @param buffer_length Size of buffer to be read from the fragment file. This must be longer than the longest line in the file.
#' @return Returns a list
#' @details Uses FeatureSpace functions \code{\link{FindDataFiles}}, \code{\link{ReadBarcodeMetadata}}, \code{\link{FilterFragments}}, \code{\link{CreateFragmentObjects}}, \code{\link{SetUsedChromosomes}}, \code{\link{CreateFeatureMatrix}}, \code{\link{CreateChromAssay}}, \code{\link{CreateSeuratObj}} and \code{\link{MergeReplicates}}.
#' @examples
#' \dontrun{
#' p <- paste(data_path, "e14vmb1_data/scATAC/E14VMB_0", sep="")
#' E14VMB <- ProcessReplicates(p, sample.names=c("_0"),\
#'  is.cell.name="is__cell_barcode", buffer_length=512L)
#' print(E14VMB)
#' }
#' @export
ProcessReplicates <- function(data.paths,sample.names, barcode.file="*singlecell.csv", is.cell.name="is__cell", force=FALSE, validate.fragments=FALSE, genome="mm10", genome.obj=BSgenome.Mmusculus.UCSC.mm10, used.chromosomes = c("chr1", "chr2","chr3","chr4","chr5","chr6","chr7","chr8","chr9","chr10","chr11","chr12","chr13","chr14","chr15","chr16","chr17","chr18","chr19"),binsize=5000, process_n=5000,buffer_length=256L){
  data.files <- FindDataFiles(p,sample.names =sample.names, barcode.file=barcode.file)
  data.files <- ReadBarcodeMetadata(data.files, is.cell.name = is.cell.name)
  data.files <- FilterFragments(data.files, force=force, buffer_length=buffer_length)
  data.files <- CreateFragmentObjects(data.files, validate.fragments = validate.fragments)
  data.files <- SetUsedChromosomes(data.files, genome.obj=genome.obj, used.chromosomes = used.chromosomes)
  data.files <- CreateFeatureMatrix(data.files, binsize=binsize, process_n=process_n)
  data.files <- CreateChromAssay(data.files, genome=genome)
  data.files <- CreateSeuratObj(data.files)
  data.files <- MergeReplicates(data.files, sample.names = sample.names)
  return(data.files)
}

#' ReadBarcodeMetadata
#'
#' Adds list components to a list including a path to a barcode file. The components are the barcode file as a tibble ("barcode.metadata") and a character vector including only cell barcodes ("true.cell.barcodes").
#' @param data.list A list including a path to a file with barcode information for each replicate.
#' @param is.cell.name A column in the barcode file that tells if the barcode in question is a cell barcode (1) or not (0).
#' @return Returns a list
#' @details Part of the wrapper function \code{\link{ProcessReplicates}}.
#' @examples
#' \dontrun{
#' barcode_metadata <- ReadBarcodeMetadata("path/to/barcode_metadata.csv")
#' }
#' @export
ReadBarcodeMetadata <- function(data.list, is.cell.name="is__cell_barcode"){
  data.list <- lapply(data.list, function(p){
    p$barcode.metadata <- readr::read_csv(file = p$barcode.file, show_col_types = FALSE)
    p$true.cell.barcodes <- dplyr::filter(p$barcode.metadata, !!rlang::sym(is.cell.name) == 1) %>% pull(barcode)
    return(p)
  })
  return(data.list)
}

#' RightSVD_dist
#'
#' @param data.object A Seurat object containing chromatin accessibility data.
#' @param drop.components Numeric or integer vector specifying the SVD components to exclude from the analysis. Default is 1.
#' @return A list containing:
#' \item{hclust.cells}{Hierarchical clustering object of cells.}
#' \item{cell.dist}{A distance matrix of cells computed using cosine distance.}
#' @examples
#' \dontrun{
#' rightSVD <- RightSVD_dist(E14VMB)
#' print(rightSVD$hclust.cells)
#' print(rightSVD$cell.dist)
#' }
#' @export
RightSVD_dist <- function(data.object, drop.components=1){
  SVD.d <- data.object@reductions$lsi@misc$d
  SVD.u <- data.object@reductions$lsi@misc$u
  num.dimensions <- ncol(data.object@reductions$lsi)
  d_diagtsne = matrix(0, nrow=num.dimensions, ncol=num.dimensions)
  diag(d_diagtsne) = SVD.d
  if (!any(is.na(drop.components))){
    d_diagtsne[drop.components,drop.components] = 0
  }
  SVDtsne_vd = t(d_diagtsne %*% t(SVD.u))
  cell.dist <- proxy::dist(SVDtsne_vd, method="cosine")
  hclust_cells <- stats::hclust(cell.dist, method="ward.D2")
  return(list(hclust.cells=hclust_cells, cell.dist=cell.dist))
}

#' SetUsedChromosomes
#'
#' Creates a list component (for every replicate) with chromosome sequence lengths for every selected chromosome.
#' @param data.list A list of 1 or more list components.
#' @param genome.obj An object containing sequence information.
#' @param used.chromosomes The chromosomes used in the analysis.
#' @return Returns a list
#' @details Part of the wrapper function \code{\link{ProcessReplicates}}.
#' @examples
#' \dontrun{
#' SetUsedChromosomes(c("chr1", "chr2", "chr3"))
#' }
#' @export
SetUsedChromosomes <- function(data.list, genome.obj=BSgenome.Mmusculus.UCSC.mm10, used.chromosomes = c("chr1", "chr2","chr3","chr4","chr5","chr6","chr7","chr8","chr9","chr10","chr11","chr12","chr13","chr14","chr15","chr16","chr17","chr18","chr19")){
  data.files <- lapply(data.list, function(p){
    p$genome.lengths  <- GenomeInfoDb::seqlengths(genome.obj)[names(GenomeInfoDb::seqlengths(genome.obj)) %in% used.chromosomes]
    return(p)
  })
  return(data.files)
}

#' PlotCladeTree
#'
#' Draws a tree of cell clades based on cut height.
#' @param hcd Out put from RightSVD_dist()
#' @param cut.height The height under which all branches are combined into a cell clade.
#' @return Draws a plot.
#' @examples
#' \dontrun{
#' hclust.obj <- hclust(dist(matrix(rnorm(100), nrow=10)))
#' PlotCladeTree(hclust.obj)
#' }
#' @export
PlotCladeTree <- function(hcd, cut.height){
  par(mai=c(2,.5,.25,.25),cex=.75)
  hcd.cut <- cut(stats::as.dendrogram(hcd$hclust.cells), h = cut.height)
  hcd.upper<-hcd.cut$upper
  hcd.upper.clusters<-dendextend::cutree(hcd.upper,h=cut.height)
  hcd.upper<-dendextend::reindex_dend(hcd.upper)
  # Find out why this needs to be done, what is not reset in the object after cut but is reset after conversion back and worth
  hcd.upper<-stats::as.dendrogram(as.hclust(hcd.upper))
  n.term.nodes <- sapply(hcd.cut$lower,dendextend::count_terminal_nodes)
  leaf.labels <- paste("Clade #",hcd.upper.clusters,": ",n.term.nodes," cells",sep="")
  hcd.upper<-dendextend::set_labels(hcd.upper,leaf.labels)
  #color_branches(hcd.upper,clusters=hcd.upper.clusters,col=color.vector(hcd.upper.clusters))
  hcd.upper %>% dendextend::set("branches_lwd",5) %>% dendextend::set("labels_cex",1.2) %>% plot(center=T)
  abline(h = h, lty = 2, col="red", lwd=2)
}

#' PlotCladeStats
#'
#' @param clade.stats A list containing clustering statistics including Dunn index and silhouette width.
#' @return A ggplot2 object visualizing clustering statistics across different cut heights.
#' @examples
#' \dontrun{
#' heights <- seq(0.01, 1, by=0.01)
#' clade.stats <- CladeStats(E14VMB, RightSVD_dist(E14VMB), heights, cores = 2)
#' PlotCladeStats(clade.stats)
#' }
#' @export
PlotCladeStats <- function(clade.stats){
  avg.silhouette <- sapply(clade.stats, function(x){x$avg.silwidth})
  dunn2 <- sapply(clade.stats, function(x){x$dunn2})
  clusterint.stats <- tibble(clus.param=heights.to.test, dunn2=dunn2, avg.silhouette=avg.silhouette)
  ggplot(clusterint.stats, aes(x=clus.param, group=1)) + geom_line(aes(y=dunn2, color="Dunn2")) + theme_minimal() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + scale_y_continuous("Dunn2", sec.axis = sec_axis(~ (. / 2.5), name = "Avg. Silhouette")) + geom_line(aes(y=avg.silhouette*2.5, color="Avg. Silhouette")) + xlab("Clade cutting height") + scale_x_continuous(breaks = seq(4, 10, by = 0.5))  + theme(legend.title=element_blank(), text = element_text(size=16))
}
