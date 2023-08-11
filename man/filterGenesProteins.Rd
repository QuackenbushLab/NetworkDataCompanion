\name{filterGenesProteins}
\alias{filterGenesProteins}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Filtering protein coding genes through an rds object.
}
\description{
Filtering protein coding genes through an rds object.
}
\usage{
filterGenesProteins(rds_gene_info)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{rds_gene_info}{rds info object of genes, usually extracted from row information from recount retrieved rds expression objects.}
}
\details{
%%  ~~ If necessary, more details than the description above ~~

}
\value{
 Array of indices that correspond to the protein coding genes in the rds gene info table.
}

\author{
Jonas Fischer (jfischer@hsph.harvard.edu)
}
\references{
}
