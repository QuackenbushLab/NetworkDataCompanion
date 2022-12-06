\name{mapProbesToGenes}
\alias{mapProbesToGenes}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Maps input probe IDs to gene TSS within a certain range upstream and downstream.
}
\description{
}
\usage{
mapProbesToGenes(probelist, rangeUp, rangeDown, localManifestPath=NA)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{probelist}{A character vector of Illumina array probes, e.g., c("cg03636183","cg19859270").}
\item{rangeUp}{The number of base pairs upstream to search for a TSS. Must be a non-negative number.}
\item{rangeDown}{The number of base pairs downstream to search for a TSS. Must be a non-negative number.}
\item{localManifestPath}{If you wish to use a manifest file other than the Illumina manifest found at \url{https://zhouserver.research.chop.edu/InfiniumAnnotation/20210615/HM450/HM450.hg38.manifest.gencode.v36.tsv.gz}, you can pass a path to that file here. It should be formatted in the same way as the Illumina manifest.}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
A matrix with four columns: probeID, geneName, ensemblID, distToTSS. When a probe maps to more than one TSS within the upstream and downstream parameters provided, the geneName, ensemblID, and distToTSS columns wil contain lists of genes separated by a semicolon (";"). Ordering of the lists matches between the three columns.
}
\references{
%% ~put references to the literature/web site here ~
https://zwdzwd.github.io/InfiniumAnnotation
}
\author{
%%  ~~who you are~~
Kate Hoff Shutta (kshutta@hsph.harvard.edu)
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{}
