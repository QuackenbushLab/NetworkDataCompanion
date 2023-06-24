\name{filterSampleType}
\alias{filterSampleType}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Filter samples based on sample type.
}
\description{This function filters samples based on sample types. Some examples are: "Primary Tumor", "Solid Tissue Normal", "Primary Blood Derived Cancer - Peripheral Blood".
}
\usage{
  filterSampleType(TCGA_barcodes, types_of_samples, rds_info)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{TCGA_barcodes}{A character vector of TCGA barcodes.}
  \item{types_of_samples}{A character vector representing the types of samples to select.}
  \item{rds_info}{A data frame extracted from a \code{RangedSummarizedExperiment} object containing TCGA metadata on samples.}
}

\details{
  Candidate values for \code{types_of_samples}: "Primary Tumor", "Solid Tissue Normal", "Primary Blood Derived Cancer - Peripheral Blood". If unsure, call this function with an invalid argument. The error message will list available sample types in rds_info.
}
\value{Named list of containing "index", an integer vector of indices in \code{TCGA_barcodes} to keep, and "type", a character vector of sample type corresponding to each index.
}

\author{
  Jonas Fischer (jfischer@hsph.harvard.edu)
}
