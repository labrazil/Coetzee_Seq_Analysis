\name{tagsnp-class}
\Rdversion{1.1}
\docType{class}
\alias{tagsnp-class}
\alias{alt.allele<-,tagsnp-method}
\alias{alt.allele,tagsnp-method}
\alias{chromosome<-,tagsnp-method}
\alias{chromosome,tagsnp-method}
\alias{correlated.snps<-,tagsnp-method}
\alias{correlated.snps,tagsnp-method}
\alias{D.prime.corrsnps<-,tagsnp-method}
\alias{D.prime.corrsnps,tagsnp-method}
\alias{genotype<-,tagsnp-method}
\alias{genotype,tagsnp-method}
\alias{overlapping.features<-,tagsnp-method}
\alias{overlapping.features,tagsnp-method}
\alias{population<-,tagsnp-method}
\alias{population,tagsnp-method}
\alias{position<-,tagsnp-method}
\alias{position,tagsnp-method}
\alias{ref.allele<-,tagsnp-method}
\alias{ref.allele,tagsnp-method}
\alias{R.squared.corrsnps<-,tagsnp-method}
\alias{R.squared.corrsnps,tagsnp-method}
\alias{snpid<-,tagsnp-method}
\alias{snpid,tagsnp-method}
\alias{summary,tagsnp-method}

\title{titlehereClass \code{"tagsnp"}}
\description{
%%  ~~ A concise (1-5 lines) description of what the class is. ~~
}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("tagsnp", ...)}.
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
    \item{\code{correlated.snps}:}{Object of class \code{"corrsnps"} ~~ }
  }
}
\section{Methods}{
  \describe{
    \item{alt.allele<-}{\code{signature(x = "tagsnp")}: ... }
    \item{alt.allele}{\code{signature(x = "tagsnp")}: ... }
    \item{chromosome<-}{\code{signature(x = "tagsnp")}: ... }
    \item{chromosome}{\code{signature(x = "tagsnp")}: ... }
    \item{correlated.snps<-}{\code{signature(x = "tagsnp")}: ... }
    \item{correlated.snps}{\code{signature(x = "tagsnp")}: ... }
    \item{D.prime.corrsnps<-}{\code{signature(x = "tagsnp")}: ... }
    \item{D.prime.corrsnps}{\code{signature(x = "tagsnp")}: ... }
    \item{genotype<-}{\code{signature(x = "tagsnp")}: ... }
    \item{genotype}{\code{signature(x = "tagsnp")}: ... }
    \item{overlapping.features<-}{\code{signature(x = "tagsnp")}: ... }
    \item{overlapping.features}{\code{signature(x = "tagsnp")}: ... }
    \item{population<-}{\code{signature(x = "tagsnp")}: ... }
    \item{population}{\code{signature(x = "tagsnp")}: ... }
    \item{position<-}{\code{signature(x = "tagsnp")}: ... }
    \item{position}{\code{signature(x = "tagsnp")}: ... }
    \item{ref.allele<-}{\code{signature(x = "tagsnp")}: ... }
    \item{ref.allele}{\code{signature(x = "tagsnp")}: ... }
    \item{R.squared.corrsnps<-}{\code{signature(x = "tagsnp")}: ... }
    \item{R.squared.corrsnps}{\code{signature(x = "tagsnp")}: ... }
    \item{snpid<-}{\code{signature(x = "tagsnp")}: ... }
    \item{snpid}{\code{signature(x = "tagsnp")}: ... }
    \item{summary}{\code{signature(object = "tagsnp")}: ... }
	 }
}
\references{
%%  ~~put references to the literature/web site here~~
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%%  ~~objects to See Also as \code{\link{~~fun~~}}, ~~~
%%  ~~or \code{\linkS4class{CLASSNAME}} for links to other classes ~~~
}
\examples{
showClass("tagsnp")
}
\keyword{classes}
