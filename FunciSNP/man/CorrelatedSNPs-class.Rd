\name{CorrelatedSNPs-class}
\Rdversion{1.1}
\docType{class}
\alias{CorrelatedSNPs-class}
\alias{AFR.D.prime<-,CorrelatedSNPs-method}
\alias{AFR.D.prime,CorrelatedSNPs-method}
\alias{AFR.p.value<-,CorrelatedSNPs-method}
\alias{AFR.p.value,CorrelatedSNPs-method}
\alias{AFR.R.squared<-,CorrelatedSNPs-method}
\alias{AFR.R.squared,CorrelatedSNPs-method}
\alias{ALL.D.prime<-,CorrelatedSNPs-method}
\alias{ALL.D.prime,CorrelatedSNPs-method}
\alias{ALL.p.value<-,CorrelatedSNPs-method}
\alias{ALL.p.value,CorrelatedSNPs-method}
\alias{ALL.R.squared<-,CorrelatedSNPs-method}
\alias{ALL.R.squared,CorrelatedSNPs-method}
\alias{alt.allele<-,CorrelatedSNPs-method}
\alias{alt.allele,CorrelatedSNPs-method}
\alias{AMR.D.prime<-,CorrelatedSNPs-method}
\alias{AMR.D.prime,CorrelatedSNPs-method}
\alias{AMR.p.value<-,CorrelatedSNPs-method}
\alias{AMR.p.value,CorrelatedSNPs-method}
\alias{AMR.R.squared<-,CorrelatedSNPs-method}
\alias{AMR.R.squared,CorrelatedSNPs-method}
\alias{ASN.D.prime<-,CorrelatedSNPs-method}
\alias{ASN.D.prime,CorrelatedSNPs-method}
\alias{ASN.p.value<-,CorrelatedSNPs-method}
\alias{ASN.p.value,CorrelatedSNPs-method}
\alias{ASN.R.squared<-,CorrelatedSNPs-method}
\alias{ASN.R.squared,CorrelatedSNPs-method}
\alias{chr<-,CorrelatedSNPs-method}
\alias{chr,CorrelatedSNPs-method}
\alias{EUR.D.prime<-,CorrelatedSNPs-method}
\alias{EUR.D.prime,CorrelatedSNPs-method}
\alias{EUR.p.value<-,CorrelatedSNPs-method}
\alias{EUR.p.value,CorrelatedSNPs-method}
\alias{EUR.R.squared<-,CorrelatedSNPs-method}
\alias{EUR.R.squared,CorrelatedSNPs-method}
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

\title{Class \code{"CorrelatedSNPs"}}
\description{
ccc%%  ~~ A concise (1-5 lines) description of what the class is. ~~
}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("CorrelatedSNPs", ...)}.
%%  ~~ describe objects here ~~ 
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
    \item{\code{ALL.R.squared}:}{Object of class \code{"matrix"} ~~ }
    \item{\code{AFR.R.squared}:}{Object of class \code{"matrix"} ~~ }
    \item{\code{AMR.R.squared}:}{Object of class \code{"matrix"} ~~ }
    \item{\code{ASN.R.squared}:}{Object of class \code{"matrix"} ~~ }
    \item{\code{EUR.R.squared}:}{Object of class \code{"matrix"} ~~ }
    \item{\code{ALL.D.prime}:}{Object of class \code{"matrix"} ~~ }
    \item{\code{AFR.D.prime}:}{Object of class \code{"matrix"} ~~ }
    \item{\code{AMR.D.prime}:}{Object of class \code{"matrix"} ~~ }
    \item{\code{ASN.D.prime}:}{Object of class \code{"matrix"} ~~ }
    \item{\code{EUR.D.prime}:}{Object of class \code{"matrix"} ~~ }
    \item{\code{ALL.p.value}:}{Object of class \code{"list"} ~~ }
    \item{\code{AFR.p.value}:}{Object of class \code{"list"} ~~ }
    \item{\code{AMR.p.value}:}{Object of class \code{"list"} ~~ }
    \item{\code{ASN.p.value}:}{Object of class \code{"list"} ~~ }
    \item{\code{EUR.p.value}:}{Object of class \code{"list"} ~~ }
  }
}
\section{Methods}{
  \describe{
    \item{AFR.D.prime<-}{\code{signature(x = "CorrelatedSNPs")}: ... }
    \item{AFR.D.prime}{\code{signature(x = "CorrelatedSNPs")}: ... }
    \item{AFR.p.value<-}{\code{signature(x = "CorrelatedSNPs")}: ... }
    \item{AFR.p.value}{\code{signature(x = "CorrelatedSNPs")}: ... }
    \item{AFR.R.squared<-}{\code{signature(x = "CorrelatedSNPs")}: ... }
    \item{AFR.R.squared}{\code{signature(x = "CorrelatedSNPs")}: ... }
    \item{ALL.D.prime<-}{\code{signature(x = "CorrelatedSNPs")}: ... }
    \item{ALL.D.prime}{\code{signature(x = "CorrelatedSNPs")}: ... }
    \item{ALL.p.value<-}{\code{signature(x = "CorrelatedSNPs")}: ... }
    \item{ALL.p.value}{\code{signature(x = "CorrelatedSNPs")}: ... }
    \item{ALL.R.squared<-}{\code{signature(x = "CorrelatedSNPs")}: ... }
    \item{ALL.R.squared}{\code{signature(x = "CorrelatedSNPs")}: ... }
    \item{alt.allele<-}{\code{signature(x = "CorrelatedSNPs")}: ... }
    \item{alt.allele}{\code{signature(x = "CorrelatedSNPs")}: ... }
    \item{AMR.D.prime<-}{\code{signature(x = "CorrelatedSNPs")}: ... }
    \item{AMR.D.prime}{\code{signature(x = "CorrelatedSNPs")}: ... }
    \item{AMR.p.value<-}{\code{signature(x = "CorrelatedSNPs")}: ... }
    \item{AMR.p.value}{\code{signature(x = "CorrelatedSNPs")}: ... }
    \item{AMR.R.squared<-}{\code{signature(x = "CorrelatedSNPs")}: ... }
    \item{AMR.R.squared}{\code{signature(x = "CorrelatedSNPs")}: ... }
    \item{ASN.D.prime<-}{\code{signature(x = "CorrelatedSNPs")}: ... }
    \item{ASN.D.prime}{\code{signature(x = "CorrelatedSNPs")}: ... }
    \item{ASN.p.value<-}{\code{signature(x = "CorrelatedSNPs")}: ... }
    \item{ASN.p.value}{\code{signature(x = "CorrelatedSNPs")}: ... }
    \item{ASN.R.squared<-}{\code{signature(x = "CorrelatedSNPs")}: ... }
    \item{ASN.R.squared}{\code{signature(x = "CorrelatedSNPs")}: ... }
    \item{chr<-}{\code{signature(x = "CorrelatedSNPs")}: ... }
    \item{chr}{\code{signature(x = "CorrelatedSNPs")}: ... }
    \item{EUR.D.prime<-}{\code{signature(x = "CorrelatedSNPs")}: ... }
    \item{EUR.D.prime}{\code{signature(x = "CorrelatedSNPs")}: ... }
    \item{EUR.p.value<-}{\code{signature(x = "CorrelatedSNPs")}: ... }
    \item{EUR.p.value}{\code{signature(x = "CorrelatedSNPs")}: ... }
    \item{EUR.R.squared<-}{\code{signature(x = "CorrelatedSNPs")}: ... }
    \item{EUR.R.squared}{\code{signature(x = "CorrelatedSNPs")}: ... }
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
showClass("CorrelatedSNPs")
}
\keyword{classes}
