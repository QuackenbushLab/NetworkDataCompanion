\name{filterPurity}
\alias{\S4method{filterPurity}{NetworkDataCompanion}}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Filter samples based on tumor purity.
}
\description{
This function filters a character vector of TCGA barcodes for tumor purity based on a particular method and threshold.
}
\usage{
\S4method{filterPurity}{NetworkDataCompanion}(TCGA_barcodes, method="ESTIMATE", threshold=.6)
}
\arguments{
  \item{TCGA_barcodes}{Character vector of TCGA barcodes that the user wishes to filter based on tumor purity.}
  \item{method}{One of "ESTIMATE","ABSOLUTE", "LUMP", "IHC", or "CPE". Default is "ESTIMATE".}
  \item{threshold}{Threshold for purity-based filtering. Samples with a purity below \code{threshold} will be filtered out.}
}

\details{
Describe the method options.
}
\value{
Integer vector of indices indicating which samples in \code{TCGA_barcodes} should be kept.
}

\references{
This code is based on the \code{TCGAPurityFiltering} package found at \url{https://github.com/pmandros/TCGAPurityFiltering}.}

\author{
Panagiotis Mandros 
}

