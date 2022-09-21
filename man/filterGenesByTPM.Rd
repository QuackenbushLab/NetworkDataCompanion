\name{filterGenesByTPM}
\alias{filterGenesByTPM}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Filter genes based on minimum expression level (TPM) across samples.
}
\description{
Filter all genes which have at least \code{tpm_threshold} TPM scores in at least \code{sample_fraction} of samples.
}
\usage{
filterGenesByTPM(expression_tpm_matrix, tpm_threshold, sample_fraction)
}
\arguments{
  \item{expression_tpm_matrix}{A data frame extracted from a \code{RangedSummarizedExperiment} object containing expression data using the \code{extractSampleAndGeneInfo} function.}
  \item{tpm_threshold}{Numeric > 0. Genes with TPM below this values in more than \code{sample_fraction} of the data will be excluded from the analysis.}
  \item{sample_fraction}{Numeric in [0,1]. Genes with TPM below \code{tpm_threshold} in more than this fraction of the data will be excluded from the analysis.}
}

\value{Integer vector indexing the rows of \code{expression_tpm_matrix} that correspond to genes that should be kept.}

\author{
Jonas Fischer (jfischer@hsph.harvard.edu)
}
