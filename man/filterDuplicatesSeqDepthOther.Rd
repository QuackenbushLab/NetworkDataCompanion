\name{filterDuplicatesSeqDepthOther}
\alias{\S4method{filterDuplicatesSeqDepthOther}{NetworkDataCompanion}}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{A version of \code{filterDuplicatesSeqDepth} to handle the case when sequencing depth is not available.
}
\description{
This function takes a random duplicate if no info is available on sequencing depth for all vials.
}

\usage{
\S4method{filterDuplicatesSeqDepthOther}{NetworkDataCompanion}(expression_count_matrix, tcga_barcodes)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{expression_count_matrix}{Matrix of count data from RNA-seq experiment, with genes in rows and samples in columns.}
  \item{tcga_barcodes}{List of TCGA barcodes for filtering.}
%%     ~~Describe \code{x} here~~
}

\value{
Integer vector of indices to KEEP in \code{tcga_barcodes}.
}

\author{
Jonas Fischer (jfischer@hsph.harvard.edu)
}
