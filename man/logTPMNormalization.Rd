\name{logTPMNormalization}
\alias{\S4method{logTPMNormalization}{NetworkDataCompanion}}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Function for within-array normalization of a RangedSummarizedExperiment object with log transcripts per million (TPM) normalization.
}
\description{
Returns a named list with raw counts (useful for duplicate filtering based on sequencing depth, see \code{?filterDuplicatesSeqDepth}), TPM, (useful for TPM-based filtering, see \code{?filterGenesByTPM}),
and the actual log TPM. A pseudocount of 1 is added to each TPM value for this function, so returned "log TPM" values actually correspond to log(TPM + 1).
}
\usage{
\S4method{logTPMNormalization}{NetworkDataCompanion}(expression_rds_obj)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{expression_rds_obj}{A \code{RangedSummarizedExperiment} object.}
}

\value{
\item{counts}{A data frame of RNA sequencing counts matching the row and column ordering of \code{expression_rds_obj}.}
\item{TPM}{A data frame of TPM matching the row and column ordering of \code{expression_rds_obj}.}
\item{logTPM}{A data frame of log-transformed TPM with pseudocounts (i.e., log(TPM + 1)) matching the row and column ordering of \code{expression_rds_obj}.}
}
\references{
RangedSummarizedExperiment documentation
}
\author{
Jonas Fischer (jfischer@hsph.harvard.edu)
}

