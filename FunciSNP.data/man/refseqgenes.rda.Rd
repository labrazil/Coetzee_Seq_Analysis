\name{refseqgenes.rda}
\alias{refseqgenes.rda}
\docType{data}
\title{Known transcripts location}
\description{
    This data is taken from the UCSC table browser. Data retrieved on March of 2012. This data was created to be used in the vignette for the \pkg{FunciSNP} pacakge.

}
\usage{data(refseqgenes.rda)}

\seealso{
	\code{\link[FunciSNP.data]{EncodeDnaseI_only.known.bed}}
	\code{\link[FunciSNP.data]{EncodeDnaseI_withCTCF.known.bed}}
	\code{\link[FunciSNP.data]{CTCF_only.known.bed}}
	\code{\link[FunciSNP.data]{EncodeFaire.known.bed}}
	\code{\link[FunciSNP.data]{knownGene.Promoters.known.bed}}
	\code{\link[FunciSNP.data]{lincRNA.hg19.rda}}
	\code{\link[FunciSNP.data]{refseqgenes.rda}}
}
\examples{
data(refseqgenes)
dim(refseqgenes)
head(refseqgenes)
}
\keyword{datasets}
