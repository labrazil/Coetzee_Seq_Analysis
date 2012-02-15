\name{TagSNP-class}
\Rdversion{1.1}
\docType{class}
\alias{TagSNP-class}
\alias{AFR.overlapping.snps.geno,TagSNP-method}
\alias{ALL.overlapping.snps.geno,TagSNP-method}
\alias{alt.allele<-,TagSNP-method}
\alias{alt.allele,TagSNP-method}
\alias{AMR.overlapping.snps.geno,TagSNP-method}
\alias{ASN.overlapping.snps.geno,TagSNP-method}
\alias{chr<-,TagSNP-method}
\alias{chr,TagSNP-method}
\alias{correlated.snps<-,TagSNP-method}
\alias{correlated.snps,TagSNP-method}
\alias{D.prime.corrsnps<-,TagSNP-method}
\alias{D.prime.corrsnps,TagSNP-method}
\alias{EUR.overlapping.snps.geno,TagSNP-method}
\alias{genotype<-,TagSNP-method}
\alias{genotype,TagSNP-method}
\alias{overlapping.features<-,TagSNP-method}
\alias{overlapping.features,TagSNP-method}
\alias{population<-,TagSNP-method}
\alias{population,TagSNP-method}
\alias{position<-,TagSNP-method}
\alias{position,TagSNP-method}
\alias{ref.allele<-,TagSNP-method}
\alias{ref.allele,TagSNP-method}
\alias{R.squared.corrsnps<-,TagSNP-method}
\alias{R.squared.corrsnps,TagSNP-method}
\alias{show,TagSNP-method}
\alias{snpid<-,TagSNP-method}
\alias{snpid,TagSNP-method}

\title{Class \code{"TagSNP"}}
\description{
hello%%  ~~ A concise (1-5 lines) description of what the class is. ~~
}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("TagSNP", ...)}.
%%  ~~ describe objects here ~~ 
}
\section{Slots}{
  \describe{
    \item{\code{chromosome}:}{Object of class \code{"integer"} ~~ }
    \item{\code{position}:}{Object of class \code{"integer"} ~~ }
    \item{\code{snpid}:}{Object of class \code{"character"} ~~ }
    \item{\code{population}:}{Object of class \code{"character"} ~~ }
    \item{\code{ref.allele}:}{Object of class \code{"character"} ~~ }
    \item{\code{alt.allele}:}{Object of class \code{"character"} ~~ }
    \item{\code{overlapping.features}:}{Object of class \code{"GRanges"} ~~ }
    \item{\code{genotype}:}{Object of class \code{"SnpMatrix"} ~~ }
    \item{\code{R.squared.corrsnps}:}{Object of class \code{"dgCMatrix"} ~~ }
    \item{\code{D.prime.corrsnps}:}{Object of class \code{"dgCMatrix"} ~~ }
    \item{\code{correlated.snps}:}{Object of class \code{"CorrelatedSNPs"} ~~ }
  }
}
\section{Methods}{
  \describe{
    \item{AFR.overlapping.snps.geno}{\code{signature(object = "TagSNP")}: ... }
    \item{ALL.overlapping.snps.geno}{\code{signature(object = "TagSNP")}: ... }
    \item{alt.allele<-}{\code{signature(x = "TagSNP")}: ... }
    \item{alt.allele}{\code{signature(x = "TagSNP")}: ... }
    \item{AMR.overlapping.snps.geno}{\code{signature(object = "TagSNP")}: ... }
    \item{ASN.overlapping.snps.geno}{\code{signature(object = "TagSNP")}: ... }
    \item{chr<-}{\code{signature(x = "TagSNP")}: ... }
    \item{chr}{\code{signature(x = "TagSNP")}: ... }
    \item{correlated.snps<-}{\code{signature(x = "TagSNP")}: ... }
    \item{correlated.snps}{\code{signature(x = "TagSNP")}: ... }
    \item{D.prime.corrsnps<-}{\code{signature(x = "TagSNP")}: ... }
    \item{D.prime.corrsnps}{\code{signature(x = "TagSNP")}: ... }
    \item{EUR.overlapping.snps.geno}{\code{signature(object = "TagSNP")}: ... }
    \item{genotype<-}{\code{signature(x = "TagSNP")}: ... }
    \item{genotype}{\code{signature(x = "TagSNP")}: ... }
    \item{overlapping.features<-}{\code{signature(x = "TagSNP")}: ... }
    \item{overlapping.features}{\code{signature(x = "TagSNP")}: ... }
    \item{population<-}{\code{signature(x = "TagSNP")}: ... }
    \item{population}{\code{signature(x = "TagSNP")}: ... }
    \item{position<-}{\code{signature(x = "TagSNP")}: ... }
    \item{position}{\code{signature(x = "TagSNP")}: ... }
    \item{ref.allele<-}{\code{signature(x = "TagSNP")}: ... }
    \item{ref.allele}{\code{signature(x = "TagSNP")}: ... }
    \item{R.squared.corrsnps<-}{\code{signature(x = "TagSNP")}: ... }
    \item{R.squared.corrsnps}{\code{signature(x = "TagSNP")}: ... }
    \item{show}{\code{signature(object = "TagSNP")}: ... }
    \item{snpid<-}{\code{signature(x = "TagSNP")}: ... }
    \item{snpid}{\code{signature(x = "TagSNP")}: ... }
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
showClass("TagSNP")
}
\keyword{classes}
