\name{getGeneInfo}
\alias{\S4method{getGeneInfo}{NetworkDataCompanion}}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Retrieve a variety of gene information based on gene name or Ensemble stable ID.
}
\description{ This function uses the \code{gene_mapping} attribute of the NetworkDataCompanion object to provide information on \code{seqid}, \code{source},\code{start},\code{end},\code{strand}, \code{gene_id}, \code{gene_name}, \code{gene_type}, and \code{gene_id_no_ver}.
}
\usage{
\S4method{getGeneInfo}{NetworkDataCompanion}(gene_names_or_ids)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{gene_names_or_ids}{A character vector of gene names or Ensembl stable IDs.
}
}
\details{
This function will determine the input type based on the presence of the string "ENSG".
}
\value{
A data frame with rows representing genes and columns representing gene attributes (e.g., source, start, end.)
}

\author{
Panagiotis Mandros (pmandros@hsph.harvard.edu)
}

