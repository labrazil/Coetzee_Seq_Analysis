\name{lincRNA.hg19.rda}
\alias{lincRNA.hg19.rda}
\alias{lincRNA}
\docType{data}
\title{Known long intergentic non coding RNA}
\description{
    This data is taken from the UCSC table browser. Data retrieved on March of 2012. This data was created to be used in the vignette for the \pkg{FunciSNP} pacakge.

}
\usage{lincRNA}

\seealso{
	\code{\link[FunciSNP.data]{FunciSNP.builtin.features}}
	\code{\link[FunciSNP.data]{refseqgenes.rda}}
}

\references{
SG. Coetzee, SK. Rhie, BP. Berman, GA. Coetzee and H. Noushmehr, FunciSNP: An
R/Bioconductor Tool Integrating Functional Non-coding Datasets with Genetic
Association Studies to Identify Candidate Regulatory SNPs., Nucleic Acids
Research, In press, 2012 (doi:10.1093/nar/gks542).
}

\examples{
lincRNA
dim(lincRNA)
head(lincRNA)
}
\keyword{datasets}
