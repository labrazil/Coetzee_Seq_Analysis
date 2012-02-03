\name{corrsnps-class}
\Rdversion{1.1}
\docType{class}
\alias{corrsnps-class}
\alias{AFR.D.prime<-,corrsnps-method}
\alias{AFR.D.prime,corrsnps-method}
\alias{AFR.geno<-,corrsnps-method}
\alias{AFR.geno,corrsnps-method}
\alias{AFR.overlapping.snps.geno<-,corrsnps-method}
\alias{AFR.overlapping.snps.geno,corrsnps-method}
\alias{AFR.p.value<-,corrsnps-method}
\alias{AFR.p.value,corrsnps-method}
\alias{AFR.R.squared<-,corrsnps-method}
\alias{AFR.R.squared,corrsnps-method}
\alias{ALL.D.prime<-,corrsnps-method}
\alias{ALL.D.prime,corrsnps-method}
\alias{ALL.geno<-,corrsnps-method}
\alias{ALL.geno,corrsnps-method}
\alias{ALL.overlapping.snps.geno<-,corrsnps-method}
\alias{ALL.overlapping.snps.geno,corrsnps-method}
\alias{ALL.p.value<-,corrsnps-method}
\alias{ALL.p.value,corrsnps-method}
\alias{ALL.R.squared<-,corrsnps-method}
\alias{ALL.R.squared,corrsnps-method}
\alias{alt.allele<-,corrsnps-method}
\alias{alt.allele,corrsnps-method}
\alias{AMR.D.prime<-,corrsnps-method}
\alias{AMR.D.prime,corrsnps-method}
\alias{AMR.geno<-,corrsnps-method}
\alias{AMR.geno,corrsnps-method}
\alias{AMR.overlapping.snps.geno<-,corrsnps-method}
\alias{AMR.overlapping.snps.geno,corrsnps-method}
\alias{AMR.p.value<-,corrsnps-method}
\alias{AMR.p.value,corrsnps-method}
\alias{AMR.R.squared<-,corrsnps-method}
\alias{AMR.R.squared,corrsnps-method}
\alias{ASN.D.prime<-,corrsnps-method}
\alias{ASN.D.prime,corrsnps-method}
\alias{ASN.geno<-,corrsnps-method}
\alias{ASN.geno,corrsnps-method}
\alias{ASN.overlapping.snps.geno<-,corrsnps-method}
\alias{ASN.overlapping.snps.geno,corrsnps-method}
\alias{ASN.p.value<-,corrsnps-method}
\alias{ASN.p.value,corrsnps-method}
\alias{ASN.R.squared<-,corrsnps-method}
\alias{ASN.R.squared,corrsnps-method}
\alias{chromosome<-,corrsnps-method}
\alias{chromosome,corrsnps-method}
\alias{EUR.D.prime<-,corrsnps-method}
\alias{EUR.D.prime,corrsnps-method}
\alias{EUR.geno<-,corrsnps-method}
\alias{EUR.geno,corrsnps-method}
\alias{EUR.overlapping.snps.geno<-,corrsnps-method}
\alias{EUR.overlapping.snps.geno,corrsnps-method}
\alias{EUR.p.value<-,corrsnps-method}
\alias{EUR.p.value,corrsnps-method}
\alias{EUR.R.squared<-,corrsnps-method}
\alias{EUR.R.squared,corrsnps-method}
\alias{overlapping.features<-,corrsnps-method}
\alias{overlapping.features,corrsnps-method}
\alias{position<-,corrsnps-method}
\alias{position,corrsnps-method}
\alias{ref.allele<-,corrsnps-method}
\alias{ref.allele,corrsnps-method}
\alias{snpid<-,corrsnps-method}
\alias{snpid,corrsnps-method}

\title{titlehereClass \code{"corrsnps"}}
\description{
%%  ~~ A concise (1-5 lines) description of what the class is. ~~
}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("corrsnps", ...)}.
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
    \item{\code{ALL.geno}:}{Object of class \code{"SnpMatrix"} ~~ }
    \item{\code{AFR.geno}:}{Object of class \code{"SnpMatrix"} ~~ }
    \item{\code{AMR.geno}:}{Object of class \code{"SnpMatrix"} ~~ }
    \item{\code{ASN.geno}:}{Object of class \code{"SnpMatrix"} ~~ }
    \item{\code{EUR.geno}:}{Object of class \code{"SnpMatrix"} ~~ }
    \item{\code{ALL.overlapping.snps.geno}:}{Object of class \code{"SnpMatrix"} ~~ }
    \item{\code{AFR.overlapping.snps.geno}:}{Object of class \code{"SnpMatrix"} ~~ }
    \item{\code{AMR.overlapping.snps.geno}:}{Object of class \code{"SnpMatrix"} ~~ }
    \item{\code{ASN.overlapping.snps.geno}:}{Object of class \code{"SnpMatrix"} ~~ }
    \item{\code{EUR.overlapping.snps.geno}:}{Object of class \code{"SnpMatrix"} ~~ }
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
    \item{AFR.D.prime<-}{\code{signature(x = "corrsnps")}: ... }
    \item{AFR.D.prime}{\code{signature(x = "corrsnps")}: ... }
    \item{AFR.geno<-}{\code{signature(x = "corrsnps")}: ... }
    \item{AFR.geno}{\code{signature(x = "corrsnps")}: ... }
    \item{AFR.overlapping.snps.geno<-}{\code{signature(x = "corrsnps")}: ... }
    \item{AFR.overlapping.snps.geno}{\code{signature(x = "corrsnps")}: ... }
    \item{AFR.p.value<-}{\code{signature(x = "corrsnps")}: ... }
    \item{AFR.p.value}{\code{signature(x = "corrsnps")}: ... }
    \item{AFR.R.squared<-}{\code{signature(x = "corrsnps")}: ... }
    \item{AFR.R.squared}{\code{signature(x = "corrsnps")}: ... }
    \item{ALL.D.prime<-}{\code{signature(x = "corrsnps")}: ... }
    \item{ALL.D.prime}{\code{signature(x = "corrsnps")}: ... }
    \item{ALL.geno<-}{\code{signature(x = "corrsnps")}: ... }
    \item{ALL.geno}{\code{signature(x = "corrsnps")}: ... }
    \item{ALL.overlapping.snps.geno<-}{\code{signature(x = "corrsnps")}: ... }
    \item{ALL.overlapping.snps.geno}{\code{signature(x = "corrsnps")}: ... }
    \item{ALL.p.value<-}{\code{signature(x = "corrsnps")}: ... }
    \item{ALL.p.value}{\code{signature(x = "corrsnps")}: ... }
    \item{ALL.R.squared<-}{\code{signature(x = "corrsnps")}: ... }
    \item{ALL.R.squared}{\code{signature(x = "corrsnps")}: ... }
    \item{alt.allele<-}{\code{signature(x = "corrsnps")}: ... }
    \item{alt.allele}{\code{signature(x = "corrsnps")}: ... }
    \item{AMR.D.prime<-}{\code{signature(x = "corrsnps")}: ... }
    \item{AMR.D.prime}{\code{signature(x = "corrsnps")}: ... }
    \item{AMR.geno<-}{\code{signature(x = "corrsnps")}: ... }
    \item{AMR.geno}{\code{signature(x = "corrsnps")}: ... }
    \item{AMR.overlapping.snps.geno<-}{\code{signature(x = "corrsnps")}: ... }
    \item{AMR.overlapping.snps.geno}{\code{signature(x = "corrsnps")}: ... }
    \item{AMR.p.value<-}{\code{signature(x = "corrsnps")}: ... }
    \item{AMR.p.value}{\code{signature(x = "corrsnps")}: ... }
    \item{AMR.R.squared<-}{\code{signature(x = "corrsnps")}: ... }
    \item{AMR.R.squared}{\code{signature(x = "corrsnps")}: ... }
    \item{ASN.D.prime<-}{\code{signature(x = "corrsnps")}: ... }
    \item{ASN.D.prime}{\code{signature(x = "corrsnps")}: ... }
    \item{ASN.geno<-}{\code{signature(x = "corrsnps")}: ... }
    \item{ASN.geno}{\code{signature(x = "corrsnps")}: ... }
    \item{ASN.overlapping.snps.geno<-}{\code{signature(x = "corrsnps")}: ... }
    \item{ASN.overlapping.snps.geno}{\code{signature(x = "corrsnps")}: ... }
    \item{ASN.p.value<-}{\code{signature(x = "corrsnps")}: ... }
    \item{ASN.p.value}{\code{signature(x = "corrsnps")}: ... }
    \item{ASN.R.squared<-}{\code{signature(x = "corrsnps")}: ... }
    \item{ASN.R.squared}{\code{signature(x = "corrsnps")}: ... }
    \item{chromosome<-}{\code{signature(x = "corrsnps")}: ... }
    \item{chromosome}{\code{signature(x = "corrsnps")}: ... }
    \item{EUR.D.prime<-}{\code{signature(x = "corrsnps")}: ... }
    \item{EUR.D.prime}{\code{signature(x = "corrsnps")}: ... }
    \item{EUR.geno<-}{\code{signature(x = "corrsnps")}: ... }
    \item{EUR.geno}{\code{signature(x = "corrsnps")}: ... }
    \item{EUR.overlapping.snps.geno<-}{\code{signature(x = "corrsnps")}: ... }
    \item{EUR.overlapping.snps.geno}{\code{signature(x = "corrsnps")}: ... }
    \item{EUR.p.value<-}{\code{signature(x = "corrsnps")}: ... }
    \item{EUR.p.value}{\code{signature(x = "corrsnps")}: ... }
    \item{EUR.R.squared<-}{\code{signature(x = "corrsnps")}: ... }
    \item{EUR.R.squared}{\code{signature(x = "corrsnps")}: ... }
    \item{overlapping.features<-}{\code{signature(x = "corrsnps")}: ... }
    \item{overlapping.features}{\code{signature(x = "corrsnps")}: ... }
    \item{position<-}{\code{signature(x = "corrsnps")}: ... }
    \item{position}{\code{signature(x = "corrsnps")}: ... }
    \item{ref.allele<-}{\code{signature(x = "corrsnps")}: ... }
    \item{ref.allele}{\code{signature(x = "corrsnps")}: ... }
    \item{snpid<-}{\code{signature(x = "corrsnps")}: ... }
    \item{snpid}{\code{signature(x = "corrsnps")}: ... }
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
showClass("corrsnps")
}
\keyword{classes}
