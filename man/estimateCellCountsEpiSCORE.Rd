\name{estimateCellCountsEpiSCORE}
\alias{\S4method{estimateCellCountsEpiSCORE}{NetworkDataCompanion}}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Run the EpiSCORE algorithm to estimate cell type proportions.}
\description{
This function applies the `constAvBetaTSS` and `wRPC` functions from the EpiSCORE R package within the TCGA data structure. The `wRPC` parameters used are the defaults: `useW=TRUE`, `wth=0.4`, and `maxit=100`.
}
\usage{
\S4method{estimateCellCountsEpiSCORE}{NetworkDataCompanion}(methylation_betas, tissue, array = "450k")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{methylation_betas}{A data frame of methylation beta values, with CGs in rows and samples in columns. The first column must be "probeID" and contain the Illumina probeIDs matching the specified array or a subset thereof.}
  \item{tissue}{Tissue type. Must be one of the tissues with a reference available in EpiSCORE. Acceptable values are "Bladder","Brain","Breast","Colon","Heart","Kidney","Liver","Lung","OE", "Pancreas_6ct","Pancreas_9ct","Prostate","Skin".}
  \item{array}{Methylation array identifier. Acceptable values are "450k" or "850k" (EPIC).}
}

\value{A data frame containing samples in rows and estimated cell type proportions in columns. The first two columns are the TCGA barcode and the TCGA UUID.
}

\references{
Teschendorff, A.E., Zhu, T., Breeze, C.E. et al. EPISCORE: cell type deconvolution of bulk tissue DNA methylomes from single-cell RNA-Seq data. Genome Biol 21, 221 (2020). \url{https://doi.org/10.1186/s13059-020-02126-9}
}

\author{
Kate Hoff Shutta (kshutta@hsph.harvard.edu)
}
