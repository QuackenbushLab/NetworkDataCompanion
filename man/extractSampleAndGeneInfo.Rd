\name{extractSampleAndGeneInfo}
\alias{\S4method{extractSampleAndGeneInfo}{NetworkDataCompanion}}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Extracts experiment-specific information and metadata from ranged summarized experiment object.
}
\description{
Extracts experiment-specific information and metadata from ranged summarized experiment object.
}
\usage{
\S4method{extractSampleAndGeneInfo}{NetworkDataCompanion}(expression_rds_obj)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{expression_rds_obj}{A ranged summarized experiment object}
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
\item{rds_sample_info}{metadata about the samples (columns)}
\item{rds_gene_info}{metadata about the genes (rows)}
}
\author{
Jonas Fischer 
}

%% ~Make other sections like Warning with \section{Warning }{....} ~
