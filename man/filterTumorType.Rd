\name{filterTumorType}
\alias{\S4method{filterTumorType}{NetworkDataCompanion}}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Filter samples based on tumor type.
}
\description{This function filters samples based on tumor type. Some examples are: "Primary Tumor", "Solid Tissue Normal", "Primary Blood Derived Cancer - Peripheral Blood". This function is particularly useful for excluding normal samples from analyses.
}
\usage{
\S4method{filterTumorType}{NetworkDataCompanion}(TCGA_barcodes, type_of_tumor, rds_info)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{TCGA_barcodes}{A character vector of TCGA barcodes.}
  \item{type_of_tumor}{A string representing the type of tumor to select. Currently, only a single tumor type is supported.}
  \item{rds_info}{A data frame extracted from a \code{RangedSummarizedExperiment} object containing expression data using the \code{extractSampleAndGeneInfo} function.}
}

\details{
Candidate values for \code{type_of_tumor}: "Primary Tumor", "Solid Tissue Normal", "Primary Blood Derived Cancer - Peripheral Blood", etc. There are other options that show up in TCGA that are not listed here. Make sure it is an exact match - check spaces, case, etc.
}
\value{Integer vector of indices in \code{TCGA_barcodes} to keep.
}

\author{
Jonas Fischer
}
