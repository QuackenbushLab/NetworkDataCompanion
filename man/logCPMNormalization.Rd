\name{logCPMNormalization}
\alias{\S4method{logCPMNormalization}{NetworkDataCompanion}}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Function to compute CPM values from raw counts.
}
\description{
This function computes CPM values from raw expression counts using the edgeR package as a backend.
}
\usage{
\S4method{logCPMNormalization}{NetworkDataCompanion}(exp_count_mat)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{exp_count_mat}{Matrix or data.frame of raw expression counts.}
}
\value{
\item{counts}{The original count matrix passed as argument to this function.}
\item{CPM}{CPM transformed values of the same shape as the count matrix.}
\item{logCPM}{log(CPM + 1) transformed values of the same shape as the count matrix.}
}
\author{
Jonas Fischer
}
\references{
\url{https://bioconductor.org/packages/release/bioc/html/edgeR.html}
}
\seealso{
\code{\link{edgeR}}.
}
