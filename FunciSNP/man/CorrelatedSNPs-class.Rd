\name{CorrelatedSNPs-class}
\Rdversion{1.1}
\docType{class}
\alias{CorrelatedSNPs-class}
\alias{class:CorrelatedSNPs}
\alias{alt.allele<-,CorrelatedSNPs-method}
\alias{alt.allele,CorrelatedSNPs-method}
\alias{chr<-,CorrelatedSNPs-method}
\alias{chr,CorrelatedSNPs-method}
\alias{overlapping.features<-}
\alias{overlapping.features<-,CorrelatedSNPs-method}
\alias{overlapping.features}
\alias{overlapping.features,CorrelatedSNPs-method}
\alias{pop.genotype<-}
\alias{pop.genotype<-,CorrelatedSNPs-method}
\alias{pop.genotype}
\alias{pop.genotype,CorrelatedSNPs-method}
\alias{position<-,CorrelatedSNPs-method}
\alias{position,CorrelatedSNPs-method}
\alias{ref.allele<-,CorrelatedSNPs-method}
\alias{ref.allele,CorrelatedSNPs-method}
\alias{snpid<-,CorrelatedSNPs-method}
\alias{snpid,CorrelatedSNPs-method}

\title{Class \code{"CorrelatedSNPs"}}
\description{
Class for CorrelatedSNPs
}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("CorrelatedSNPs", ...)}.
}
\usage{
\S4method{alt.allele}{CorrelatedSNPs}(x)
\S4method{alt.allele}{CorrelatedSNPs}(x) <- value
\S4method{chr}{CorrelatedSNPs}(x)
\S4method{chr}{CorrelatedSNPs}(x) <- value
\S4method{overlapping.features}{CorrelatedSNPs}(x)
\S4method{overlapping.features}{CorrelatedSNPs}(x) <- value
\S4method{pop.genotype}{CorrelatedSNPs}(x, population)
\S4method{pop.genotype}{CorrelatedSNPs}(x) <- value
\S4method{position}{CorrelatedSNPs}(x)
\S4method{position}{CorrelatedSNPs}(x) <- value
\S4method{ref.allele}{CorrelatedSNPs}(x)
\S4method{ref.allele}{CorrelatedSNPs}(x) <- value
\S4method{snpid}{CorrelatedSNPs}(x)
\S4method{snpid}{CorrelatedSNPs}(x) <- value
}
\arguments{
  \item{x}{
    The CorrelatedSNPs object from/on which to get/set snp information.
  }
  \item{value}{
    Examine the section "Slots" for the values taken by each method.
  }
  \item{population}{
    Population for which to set/get genotype information, either "AFR", "AMR", "ASN", "EUR" or "ALL".
  }
}
\section{Slots}{
  \describe{
    \item{\code{chromosome}:}{Object of class \code{"integer"} ~~ }
    \item{\code{position}:}{Object of class \code{"integer"} ~~ }
    \item{\code{snpid}:}{Object of class \code{"character"} ~~ }
    \item{\code{ref.allele}:}{Object of class \code{"character"} ~~ }
    \item{\code{alt.allele}:}{Object of class \code{"character"} ~~ }
    \item{\code{overlapping.features}:}{Object of class \code{"GRanges"} ~~ }
    \item{\code{genotype}:}{Object of class \code{"CorrGeno"} ~~ }
  }
}
\section{Methods}{
  \describe{
    \item{alt.allele<-}{\code{signature(x = "CorrelatedSNPs")}: ... }
    \item{alt.allele}{\code{signature(x = "CorrelatedSNPs")}: ... }
    \item{chr<-}{\code{signature(x = "CorrelatedSNPs")}: ... }
    \item{chr}{\code{signature(x = "CorrelatedSNPs")}: ... }
    \item{overlapping.features<-}{\code{signature(x = "CorrelatedSNPs")}: ... }
    \item{overlapping.features}{\code{signature(x = "CorrelatedSNPs")}: ... }
    \item{pop.genotype<-}{\code{signature(x = "CorrelatedSNPs")}: ... }
    \item{pop.genotype}{\code{signature(x = "CorrelatedSNPs")}: ... }
    \item{position<-}{\code{signature(x = "CorrelatedSNPs")}: ... }
    \item{position}{\code{signature(x = "CorrelatedSNPs")}: ... }
    \item{ref.allele<-}{\code{signature(x = "CorrelatedSNPs")}: ... }
    \item{ref.allele}{\code{signature(x = "CorrelatedSNPs")}: ... }
    \item{snpid<-}{\code{signature(x = "CorrelatedSNPs")}: ... }
    \item{snpid}{\code{signature(x = "CorrelatedSNPs")}: ... }
	 }
}
\references{
SG. Coetzee, SK. Rhie, BP. Berman, GA. Coetzee and H. Noushmehr, FunciSNP: An R/Bioconductor Tool Integrating Functional Non-coding Datasets with Genetic Association Studies to Identify Candidate Regulatory SNPs., Nucleic Acids Research, In press, 2012 (doi:10.1093/nar/gks542).
}
\author{
Simon G. Coetzee (maintainer: scoetzee@gmail.com); Houtan Noushmehr, PhD (houtan@usp.br)
}
\note{
NA
}


\seealso{
\code{\link{getFSNPs}}, \code{\link{FunciSNPplot}}, \code{\link{FunciSNPAnnotateSummary}}, \code{\link{FunciSNPtable}}, \code{\link{FunciSNPbed}}
}
\examples{
showClass("CorrelatedSNPs")
}
\keyword{classes}
