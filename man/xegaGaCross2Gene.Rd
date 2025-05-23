% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/xegaGaCrossover.R
\name{xegaGaCross2Gene}
\alias{xegaGaCross2Gene}
\title{One point crossover of 2 genes.}
\usage{
xegaGaCross2Gene(gg1, gg2, lF)
}
\arguments{
\item{gg1}{A binary gene.}

\item{gg2}{A binary gene.}

\item{lF}{The local configuration of the genetic algorithm.}
}
\value{
A list of 2 binary genes.
}
\description{
\code{xegaGaCross2Gene()} randomly determines a cut point.
             It combines the bits before the cut point of the first gene
             with the bits after the cut point from the second gene (kid 1).
             It combines the bits before the cut point of the second gene
             with the bits after the cut point from the first gene (kid 2).
             It returns 2  genes.
}
\examples{
gene1<-xegaGaInitGene(lFxegaGaGene)
gene2<-xegaGaInitGene(lFxegaGaGene)
xegaGaDecodeGene(gene1, lFxegaGaGene)
xegaGaDecodeGene(gene2, lFxegaGaGene)
newgenes<-xegaGaCross2Gene(gene1, gene2, lFxegaGaGene)
xegaGaDecodeGene(newgenes[[1]], lFxegaGaGene)
xegaGaDecodeGene(newgenes[[2]], lFxegaGaGene)
}
\seealso{
Other Crossover (Returns 2 Kids): 
\code{\link{xegaGaUCross2Gene}()},
\code{\link{xegaGaUPCross2Gene}()}
}
\concept{Crossover (Returns 2 Kids)}
