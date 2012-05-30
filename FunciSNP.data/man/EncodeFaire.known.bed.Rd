\name{EncodeFaire.known.bed}
\alias{EncodeFaire.known.bed}
\docType{data}
\title{Known Encode FAIRE cluster sites extracted from ENCODE db}
\description{
    This data is taken from the UCSC Encode website. Data retrieved on March of 2012. This data was created to be used in the vignette for the \pkg{FunciSNP} pacakge.

}
\usage{data(EncodeFaire.known.bed)}

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
data(EncodeFaire.known.bed)
dim(EncodeFaire.known.bed)
head(EncodeFaire.known.bed)
}
\keyword{datasets}
