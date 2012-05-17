\name{TagSNP-class}
\Rdversion{1.1}
\docType{class}

\alias{TagSNP-class}
\alias{class:TagSNP}
\alias{alt.allele-methods}
\alias{alt.allele<-,TagSNP-method}
\alias{alt.allele,TagSNP-method}
\alias{chr-methods}
\alias{chr<-,TagSNP-method}
\alias{chr,TagSNP-method}
\alias{correlated.snps-methods}
\alias{correlated.snps<-,TagSNP-method}
\alias{correlated.snps,TagSNP-method}
\alias{genotype-methods}
\alias{genotype<-,TagSNP-method}
\alias{genotype,TagSNP-method}
\alias{population-methods}
\alias{population<-,TagSNP-method}
\alias{population,TagSNP-method}
\alias{position-methods}
\alias{position<-,TagSNP-method}
\alias{position,TagSNP-method}
\alias{ref.allele-methods}
\alias{ref.allele<-,TagSNP-method}
\alias{ref.allele,TagSNP-method}
\alias{show,TagSNP-method}
\alias{snpid-methods}
\alias{snpid<-,TagSNP-method}
\alias{snpid,TagSNP-method}
\alias{yafsnp.rsq-methods}
\alias{yafsnp.rsq<-,TagSNP-method}
\alias{yafsnp.rsq,TagSNP-method}
\alias{yafsnp.dprime-methods}
\alias{yafsnp.dprime<-,TagSNP-method}
\alias{yafsnp.dprime,TagSNP-method}
\alias{yafsnp.pvalue-methods}
\alias{yafsnp.pvalue<-,TagSNP-method}
\alias{yafsnp.pvalue,TagSNP-method}
\alias{AFR.overlapping.snps.geno-methods}
\alias{AFR.overlapping.snps.geno,TagSNP-method}
\alias{AMR.overlapping.snps.geno-methods}
\alias{AMR.overlapping.snps.geno,TagSNP-method}
\alias{ASN.overlapping.snps.geno-methods}
\alias{ASN.overlapping.snps.geno,TagSNP-method}
\alias{EUR.overlapping.snps.geno-methods}
\alias{EUR.overlapping.snps.geno,TagSNP-method}


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
SG. Coetzee, SK. Rhie, BP. Berman, GA. Coetzee and H. Noushmehr, \emph{FunciSNP: An R/Bioconductor Tool Integrating Functional Non-coding Datasets with Genetic Association Studies to Identify Candidate Regulatory SNPs}, Nucleic Acids Research, In press, 2012 (doi:10.1093/nar/GKS542).
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
