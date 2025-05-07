\name{probeToMeanPromoterMethylation}
\alias{\S4method{probeToMeanPromoterMethylation}{NetworkDataCompanion}}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{probeToMeanPromoterMethylation}
\description{Calculates the average promoter methylation within a certain window around the transcription start site (TSS), as defined by the input probe_gene_map.}
\usage{\S4method{probeToMeanPromoterMethylation}{NetworkDataCompanion}(methylation_betas, probe_gene_map, genesOfInterest)}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{methylation_betas}{A data frame of methylation beta values, with CGs in rows and samples in columns. The first column must be "probeID" and contain the Illumina probeIDs matching the probe_gene_map argument or a subset thereof.}
  \item{probe_gene_map}{Output from mapProbesToGenes, or otherwise a bespoke matrix with four columns: probeID, geneName, ensemblID, distToTSS.}
  \item{genesOfInterest}{Character vector of gene names for which mean promoter methylation should be calculated. 
  Note that each entry in the vector should contain only one gene and that gene should not be repeated in the entry; e.g., \code{c("ASMTL","CRLF2")} is a correct input while \code{c("ASMTL;ASMTL", "CRLF2")} will fail to find probes mapping to ASMTL. 
  This could happen, for example, if you pull a column of genes from a resource that combines splice variants with a semicolon.}
}

\value{Matrix of samples in rows and genes in columns. row.names stores sample names and colnames stores gene names.
}

\author{
Kate Hoff Shutta
}
