\name{convertBetaToM}
\alias{\S4method{convertBetaToM}{NetworkDataCompanion}}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Convert methylation beta values to M-values.}
\description{
This function uses the typical logit base 2 transformation to convert from methylation beta values (in the [0,1] range) to m-values (on the real line). The formula is m = log2(beta/(1-beta)).
}
\usage{
\S4method{convertBetaToM}{NetworkDataCompanion}(methylation_betas)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{methylation_betas}{A numeric vector of values in the range [0,1].
  }
}

\value{
A numeric vector of m-values corresponding to the converted values of \code{methylation_betas}.
}

\references{
\url{https://bmcbioinformatics.biomedcentral.com/articles/10.1186/1471-2105-11-587}
}
\author{
Kate Hoff Shutta (kshutta@hsph.harvard.edu)
}
