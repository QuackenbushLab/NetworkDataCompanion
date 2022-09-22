\name{geneNameToENSG}
\alias{geneNameToENSG}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Convert from gene name to Ensembl stable id.
}
\description{Given an input character vector of gene names, this function converts them to Ensembl stable IDs. Note from \url{https://useast.ensembl.org/Help/Faq?id=488}: "An Ensembl stable ID consists of five parts: ENS(species)(object type)(identifier).(version)."}

\usage{
geneNameToENSG(gene_names, version = FALSE)
}
\arguments{
  \item{gene_names}{Character vector of gene names.}
  \item{version}{Boolean; retrieve Ensembl version along with Ensembl identifier.}
}


\value{Character vector of Ensembl stable IDs.
}
\references{
\url{https://useast.ensembl.org/index.html}
}
\author{
Panagiotis Mandros (pmandros@hsph.harvard.edu)}

