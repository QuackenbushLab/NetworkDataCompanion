\name{filterSampleType}
\alias{filterSampleType}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Filter samples based on sample type.
}
\description{This function filters samples based on TCGA sample types.
}
\usage{
  filterSampleType(TCGA_barcodes, types_of_samples)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{TCGA_barcodes}{A character vector of TCGA barcodes.}
  \item{types_of_samples}{A character vector representing the types of samples to select.}
}

\details{
  Candidate values for \code{types_of_samples} are characters of the form "01","02", etc. that correspond to TCGA sample type codes. Valid arguments for \code{types_of_samples} can be found on the GDC website: \url{https://gdc.cancer.gov/resources-tcga-users/tcga-code-tables/sample-type-codes} and a table mapping sample types to values can be found by using \code{NetSciDataCompanion::getSampleTypeMap()}.
}

\value{Named list of containing "index", an integer vector of indices in \code{TCGA_barcodes} to keep, and "type", a character vector of sample type corresponding to each index.
}

\author{
  Jonas Fischer (jfischer@hsph.harvard.edu)
}
