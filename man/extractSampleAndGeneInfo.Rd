\name{extractSampleAndGeneInfo}
\alias{extractSampleAndGeneInfo}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Extracts experiment-specific information and metadata from ranged summarized experiment object.
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
extractSampleAndGeneInfo(expression_rds_obj)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{expression_rds_obj}{A ranged summarized experiment object}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
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
\references{
%% ~put references to the literature/web site here ~
}
\author{
Jonas Fischer (jfischer@hsph.harvard.edu)
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
}
