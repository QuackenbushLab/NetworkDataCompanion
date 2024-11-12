\name{mapBarcodeToBarcode}
\alias{\S4method{mapBarcodeToBarcode}{NetworkDataCompanion}}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Helper function for mapping two sets of TCGA barcodes to each other.
}
\description{
There are 4 different pieces of information returned in a named list that are all useful depending on the context in which they are used.

\code{is_inter1} is an indicator (boolean) vector of the same length as \code{bc1} that indicates which elements of \code{bc1} are present in \code{bc2}.

\code{idcs1} indicates where to find each barcode of \code{bc1} in \code{bc2}, returning NA if there is no match. That is, \code{idcs1[i] != NA}, then \code{bc1[i] == bc2[idcs1[i]]}.

The same information is provided for \code{bc2}.
}
\usage{
\S4method{mapBarcodeToBarcode}{NetworkDataCompanion}(bc1,bc2)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{bc1}{Character vector of barcodes in the first set.}
  \item{bc2}{Character vector of barcodes in the second set.}
}
\value{
\item{is_inter1}{Boolean vector of the same length as \code{bc1} that indicates which elements of \code{bc1} are present in \code{bc2}.}
\item{idcs1}{Integer vector of the same length as \code{bc1} that indicates where to find each barcode of \code{bc1} in \code{bc2}, returning \code{NA} if there is no match. That is, \code{idcs1[i] != NA}, then \code{bc1[i] == bc2[idcs1[i]]}}
\item{is_inter2}{Boolean vector of the same length as \code{bc2} that indicates which elements of \code{bc2} are present in \code{bc1}.}
\item{idcs2}{Integer vector of the same length as \code{bc2} that indicates where to find each barcode of \code{bc2} in \code{bc1}, returning \code{NA} if there is no match. That is, \code{idcs2[i] != NA}, then \code{bc2[i] == bc1[idcs2[i]]}.}
}
\author{
Jonas Fischer (jfischer@hsph.harvard.edu)
}
\note{For example, if you want to map experiment 1 onto experiment 2, keeping only the information for samples that are present in both, and reordering the first experiment to match the samples of the second, you can do:

\code{exp1[,is_inter1] # this will remove samples that are not in experiment 2)}

\code{exp2[,idcs1[is_inter1]] # this will remove samples that are not in exp1 and reorder to match exp1}
}

