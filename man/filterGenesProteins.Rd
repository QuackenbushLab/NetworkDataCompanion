\name{filterGenesProteins}
\alias{\S4method{filterGenesProteins}{NetworkDataCompanion}}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Filtering protein coding genes through an rds object.
}
\description{
Filtering protein coding genes through an rds object.
}
\usage{
\S4method{filterGenesProteins}{NetworkDataCompanion}(rds_gene_info)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{rds_gene_info}{rds info object of genes, usually extracted from row information from recount retrieved rds expression objects.}
}
\value{
 Array of indices that correspond to the protein coding genes in the rds gene info table.
}

\author{
Jonas Fischer (jfischer@hsph.harvard.edu)
}
