\name{filterBarcodesIntersection}
\alias{filterBarcodesIntersection}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Convenience wrapper function for \code{mapBarcodeToBarcode} that applies the function directly to two data frames.
}
\description{
This function returns a list of the two argument data frames, intersected, and the second frame ordered to match the first. NOTE: Ordering is done based on columns, which are expected to be named by TCGA barcodes.
}
\usage{
filterBarcodesIntersection(exp1, exp2)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{exp1}{A matrix or dataframe with TCGA barcodes as the \code{colnames} attribute.}
  \item{exp2}{A matrix or dataframe with TCGA barcodes as the \code{colnames} attribute.}
}

\details{
No additional details at this time.
}
\value{
\item{mappedExp1}{A data frame filtered to include only columns with TCGA barcodes that are in both \code{colnames(exp1)} and \code{colnames(exp2)}.}
\item{mappedExp2}{A data frame filtered to include only columns with TCGA barcodes that are in both \code{colnames(exp1)} and \code{colnames(exp2)}. IMPORTANT: The columns of \code{mappedExp2} are ordered to match the column ordering of \code{mappedExp1}.}
}

\author{
Jonas Fischer (jfischer@hsph.harvard.edu)
}

