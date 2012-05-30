\name{refseqgenes.rda}
\alias{refseqgenes.rda}
\alias{refseqgenes}
\docType{data}
\title{Known transcripts location}
\description{
    This data is taken from the UCSC table browser. Data retrieved on March of 2012. This data was created to be used in the vignette for the \pkg{FunciSNP} pacakge.

}
\usage{refseqgenes.rda}

\seealso{
	\code{\link[FunciSNP.data]{lincRNA.hg19.rda}}
	\code{\link[FunciSNP.data]{FunciSNP.builtin.features}}
}

\references{
SG. Coetzee, SK. Rhie, BP. Berman, GA. Coetzee and H. Noushmehr, FunciSNP: An
R/Bioconductor Tool Integrating Functional Non-coding Datasets with Genetic
Association Studies to Identify Candidate Regulatory SNPs., Nucleic Acids
Research, In press, 2012 (doi:10.1093/nar/gks542).
}

\examples{
refseqgenes
dim(refseqgenes)
head(refseqgenes)
}
\keyword{datasets}
