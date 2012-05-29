\name{EncodeDnaseI_only.known.bed}
\alias{EncodeDnaseI_only.known.bed}
\docType{data}
\title{Known Encode DNAseI cluster sites extracted from ENCODE db}
\description{
    This data is taken from the UCSC Encode website. Data retrieved on March of 2012. These sites do not overlap a known CTCF site. To find Encode DNAseI sites overlapping known CTCF, see EncodeDnaseI_withCTCF.known.bed. This data was created to be used in the vignette for the \pkg{FunciSNP} pacakge.

}
\usage{data(EncodeDnaseI_only.known.bed)}

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
data(EncodeDnaseI_only.known.bed)
dim(EncodeDnaseI_only.known.bed)
head(EncodeDnaseI_only.known.bed)
}
\keyword{datasets}
