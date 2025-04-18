% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/xegaGaDecode.R
\name{xegaGaGeneMapPerm}
\alias{xegaGaGeneMapPerm}
\title{Map the bit strings of a binary gene to a permutation.}
\usage{
xegaGaGeneMapPerm(gene, penv)
}
\arguments{
\item{gene}{A binary gene (the genotype).}

\item{penv}{A problem environment.}
}
\value{
A permutation (the decoded gene (the phenotype))
}
\description{
\code{xegaGaGeneMapPerm()} maps 
             the bit strings of a binary string 
             to a permutation of integers.
             Example: Traveling Salesman Problem (TSP).
}
\examples{
gene<-xegaGaInitGene(lFxegaGaGene)
xegaGaGeneMapPerm(gene$gene1, lFxegaGaGene$penv)

}
\seealso{
Other Decoder: 
\code{\link{xegaGaDecodeGene}()},
\code{\link{xegaGaGeneMap}()},
\code{\link{xegaGaGeneMapGray}()},
\code{\link{xegaGaGeneMapIdentity}()}
}
\concept{Decoder}
