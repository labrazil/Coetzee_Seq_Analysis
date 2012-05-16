\name{CorrelatedSNPs-class}
\Rdversion{1.1}
\docType{class}
\alias{CorrelatedSNPs-class}
\alias{class:CorrelatedSNPs}
\alias{alt.allele<-,CorrelatedSNPs-method}
\alias{alt.allele,CorrelatedSNPs-method}
\alias{chr<-,CorrelatedSNPs-method}
\alias{chr,CorrelatedSNPs-method}
\alias{overlapping.features<-,CorrelatedSNPs-method}
\alias{overlapping.features,CorrelatedSNPs-method}
\alias{pop.genotype<-,CorrelatedSNPs-method}
\alias{pop.genotype,CorrelatedSNPs-method}
\alias{position<-,CorrelatedSNPs-method}
\alias{position,CorrelatedSNPs-method}
\alias{ref.allele<-,CorrelatedSNPs-method}
\alias{ref.allele,CorrelatedSNPs-method}
\alias{snpid<-,CorrelatedSNPs-method}
\alias{snpid,CorrelatedSNPs-method}
\alias{overlapping.features.snps-methods}
\alias{pop.genotype-methods}

\title{Class \code{"CorrelatedSNPs"}}
\description{
Class for CorrelatedSNPs
}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("CorrelatedSNPs", ...)}.
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
Coetzee SG et al. submitted for review. 2012
}
\author{
Simon G. Coetzee; Houtan Noushmehr, PhD
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
