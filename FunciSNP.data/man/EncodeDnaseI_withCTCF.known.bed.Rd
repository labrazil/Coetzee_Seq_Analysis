\name{EncodeDnaseI_withCTCF.known.bed}
\alias{EncodeDnaseI_withCTCF.known.bed}
\docType{data}
\title{Known Encode DNAseI cluster sites with overlapping CTCF sites extracted from ENCODE db}
\description{
    This data is taken from the UCSC Encode website. Data retrieved on March of 2012. These CTCF sites overlap a known Encode DNAseI site. This data was created to be used in the vignette for the \pkg{FunciSNP} pacakge.

}
\usage{data(EncodeDnaseI_withCTCF.known.bed)}

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
data(EncodeDnaseI_withCTCF.known.bed)
dim(EncodeDnaseI_withCTCF.known.bed)
head(EncodeDnaseI_withCTCF.known.bed)
}
\keyword{datasets}
