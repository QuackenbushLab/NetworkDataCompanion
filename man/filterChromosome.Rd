\name{filterChromosome}
\alias{filterChromosome}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Filter for genes in a particular chromosome or chromosomes.
}
\description{
This function filters for genes in a particular chromosome or chromosomes.
}
\usage{
filterTumorType(rds_gene_info, chroms)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{rds_gene_info}{A data frame extracted from a \code{RangedSummarizedExperiment} object containing expression data using the \code{extractSampleAndGeneInfo} function.}
  \item{chroms}{A character vector of chromosomes. Must exactly match chromosomes in the \code{seqname} attribute of \code{rds_gene_info}.}
}

\value{Integer vector of the row indices (genes) in \code{rds_gene_info} to keep.
}

\author{
Jonas Fischer (jfischer@hsph.harvard.edu)
}
