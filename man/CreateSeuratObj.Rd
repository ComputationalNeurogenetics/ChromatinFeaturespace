% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/FeatureSpace.R
\name{CreateSeuratObj}
\alias{CreateSeuratObj}
\title{CreateSeuratObj}
\usage{
CreateSeuratObj(data.list)
}
\arguments{
\item{data.list}{A list with list components "barcode.metadata", "chr.assay" and "true.cell.barcodes".}
}
\value{
Returns a list
}
\description{
Creates a \code{Seurat} object for each replicate on a list from the data on that list.
}
\details{
Part of the wrapper function \code{\link{ProcessReplicates}}. Based on the SeuratObject function \code{\link[SeuratObject]{CreateSeuratObject}}.
}
\examples{
\dontrun{
data.list <- list(
  barcode.metadata = data.frame(barcode = c("cell1", "cell2"), count = c(100, 200)),
  chr.assay = list(),
  true.cell.barcodes = c("cell1", "cell2")
)
seurat_obj <- CreateSeuratObj(data.list)
print(seurat_obj)
}
}
