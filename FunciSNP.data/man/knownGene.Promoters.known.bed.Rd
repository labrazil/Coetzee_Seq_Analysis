\name{knownGene.Promoters.known.bed}
\alias{knownGene.Promoters.known.bed}
\docType{data}
\title{Known Gene Promoters}
\description{
    This data is taken from the UCSC table browser. We defined promoters to be -1000/+100 bp from the known Transcription start sites. Data retrieved on March of 2012. This data was created to be used in the vignette for the \pkg{FunciSNP} pacakge.

}
\usage{data(knownGene.Promoters.known.bed)}

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
data(knownGene.Promoters.known.bed)
dim(knownGene.Promoters.known.bed)
head(knownGene.Promoters.known.bed)
}
\keyword{datasets}
