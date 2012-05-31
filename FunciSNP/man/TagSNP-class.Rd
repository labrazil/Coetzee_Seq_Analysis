\name{TagSNP-class}
\Rdversion{1.1}
\docType{class}

\alias{TagSNP}
\alias{TagSNP-class}
\alias{class:TagSNP}
\alias{alt.allele<-}
\alias{alt.allele<-,TagSNP-method}
\alias{alt.allele}
\alias{alt.allele,TagSNP-method}
\alias{chr<-,TagSNP-method}
\alias{chr<-}
\alias{chr,TagSNP-method}
\alias{chr}
\alias{correlated.snps<-,TagSNP-method}
\alias{correlated.snps<-}
\alias{correlated.snps,TagSNP-method}
\alias{correlated.snps}
\alias{genotype<-,TagSNP-method}
\alias{genotype<-}
\alias{genotype,TagSNP-method}
\alias{genotype}
\alias{population<-,TagSNP-method}
\alias{population<-}
\alias{population,TagSNP-method}
\alias{population}
\alias{position<-,TagSNP-method}
\alias{position<-}
\alias{position,TagSNP-method}
\alias{position}
\alias{ref.allele<-,TagSNP-method}
\alias{ref.allele<-}
\alias{ref.allele,TagSNP-method}
\alias{ref.allele}
\alias{show,TagSNP-method}
\alias{snpid<-,TagSNP-method}
\alias{snpid<-}
\alias{snpid,TagSNP-method}
\alias{snpid}
\alias{yafsnp.rsq<-,TagSNP-method}
\alias{yafsnp.rsq<-}
\alias{yafsnp.rsq,TagSNP-method}
\alias{yafsnp.rsq}
\alias{yafsnp.dprime<-,TagSNP-method}
\alias{yafsnp.dprime<-}
\alias{yafsnp.dprime,TagSNP-method}
\alias{yafsnp.dprime}
\alias{yafsnp.pvalue<-,TagSNP-method}
\alias{yafsnp.pvalue<-}
\alias{yafsnp.pvalue,TagSNP-method}
\alias{yafsnp.pvalue}
\alias{AFR.overlapping.snps.geno,TagSNP-method}
\alias{AFR.overlapping.snps.geno}
\alias{AMR.overlapping.snps.geno,TagSNP-method}
\alias{AMR.overlapping.snps.geno}
\alias{ASN.overlapping.snps.geno,TagSNP-method}
\alias{ASN.overlapping.snps.geno}
\alias{EUR.overlapping.snps.geno,TagSNP-method}
\alias{EUR.overlapping.snps.geno}


\title{Class \code{"TagSNP"}}
\description{
  In the code snippets below, \code{x} is a \code{TagSNP} object.
  for the usage of \code{alt.allele}, \code{chr}, \code{position}, \code{ref.allele}, and \code{snpid} for the object \code{CorrelatedSNPs}, see \code{?CorrelatedSNPs-class}
}
\usage{
\S4method{alt.allele}{TagSNP}(x)
\S4method{alt.allele}{TagSNP}(x) <- value
\S4method{chr}{TagSNP}(x)
\S4method{chr}{TagSNP}(x) <- value
\S4method{correlated.snps}{TagSNP}(x)
\S4method{correlated.snps}{TagSNP}(x) <- value
\S4method{genotype}{TagSNP}(x)
\S4method{genotype}{TagSNP}(x) <- value
\S4method{population}{TagSNP}(x)
\S4method{population}{TagSNP}(x) <- value
\S4method{position}{TagSNP}(x)
\S4method{position}{TagSNP}(x) <- value
\S4method{ref.allele}{TagSNP}(x)
\S4method{ref.allele}{TagSNP}(x) <- value
\S4method{snpid}{TagSNP}(x)
\S4method{snpid}{TagSNP}(x) <- value
\S4method{yafsnp.rsq}{TagSNP}(x)
\S4method{yafsnp.rsq}{TagSNP}(x) <- value
\S4method{yafsnp.dprime}{TagSNP}(x)
\S4method{yafsnp.dprime}{TagSNP}(x) <- value
\S4method{yafsnp.pvalue}{TagSNP}(x)
\S4method{yafsnp.pvalue}{TagSNP}(x) <- value
\S4method{AFR.overlapping.snps.geno}{TagSNP}(x)
\S4method{AMR.overlapping.snps.geno}{TagSNP}(x)
\S4method{ASN.overlapping.snps.geno}{TagSNP}(x)
\S4method{EUR.overlapping.snps.geno}{TagSNP}(x)
}
\arguments{
  \item{x}{
    The TagSNP object from/on which to get/set snp information.
  }
  \item{value}{
    Examine the section "Slots" for the values taken by each method.
  }
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
    \item{\code{genotype}:}{Object of class \code{"SnpMatrix"} ~~ }
    \item{\code{yafsnp.rsq}:}{Object of class \code{"dgCMatrix"} ~~ }
    \item{\code{yafsnp.dprime}:}{Object of class \code{"dgCMatrix"} ~~ }
    \item{\code{yafsnp.pvalue}:}{Object of class \code{"list"} ~~ }
    \item{\code{correlated.snps}:}{Object of class \code{"CorrelatedSNPs"} ~~ }
  }
}
\section{Methods}{
  \describe{
    \item{alt.allele<-}{\code{signature(x = "TagSNP")}: ... }
    \item{alt.allele}{\code{signature(x = "TagSNP")}: ... }
    \item{chr<-}{\code{signature(x = "TagSNP")}: ... }
    \item{chr}{\code{signature(x = "TagSNP")}: ... }
    \item{correlated.snps<-}{\code{signature(x = "TagSNP")}: ... }
    \item{correlated.snps}{\code{signature(x = "TagSNP")}: ... }
    \item{genotype<-}{\code{signature(x = "TagSNP")}: ... }
    \item{genotype}{\code{signature(x = "TagSNP")}: ... }
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
    \item{yafsnp.rsq<-}{\code{signature(x = "TagSNP")}: ... }
    \item{yafsnp.rsq}{\code{signature(x = "TagSNP")}: ... }
    \item{yafsnp.dprime<-}{\code{signature(x = "TagSNP")}: ... }
    \item{yafsnp.dprime}{\code{signature(x = "TagSNP")}: ... }
    \item{yafsnp.pvalue<-}{\code{signature(x = "TagSNP")}: ... }
    \item{yafsnp.pvalue}{\code{signature(x = "TagSNP")}: ... }
    \item{AFR.overlapping.snps.geno}{\code{signature(x = "TagSNP")}: ... }
    \item{AMR.overlapping.snps.geno}{\code{signature(x = "TagSNP")}: ... }
    \item{ASN.overlapping.snps.geno}{\code{signature(x = "TagSNP")}: ... }
    \item{EUR.overlapping.snps.geno}{\code{signature(x = "TagSNP")}: ... }
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
showClass("TagSNP")
}
\keyword{classes}
