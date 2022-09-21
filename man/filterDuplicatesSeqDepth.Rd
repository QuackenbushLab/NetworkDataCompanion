\name{filterDuplicatesSeqDepth}
\alias{filterDuplicatesSeqDepth}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
A function to filter duplicates based on RNA sequencing depth.
}
\description{
This function filters out duplicates based on RNA-seq, keeping the samples with maximum read depth. Returns indices of samples to KEEP.
}
\usage{
filterDuplicatesSeqDepth(expression_count_matrix)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{expression_count_matrix}{Matrix of count data from RNA-seq experiment, with genes in rows and samples in columns.}
}

\value{Integer vector of indices to keep, corresponding to columns of \code{expression_count_matrix}.
}

\author{
Jonas Fischer(jfischer@hsph.harvard.edu)
}
