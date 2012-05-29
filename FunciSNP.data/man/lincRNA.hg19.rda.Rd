\name{lincRNA.hg19.rda}
\alias{lincRNA.hg19.rda}
\docType{data}
\title{Known long intergentic non coding RNA}
\description{
    This data is taken from the UCSC table browser. Data retrieved on March of 2012. This data was created to be used in the vignette for the \pkg{FunciSNP} pacakge.

}
\usage{data(lincRNA.hg19.rda)}

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
data(lincRNA.hg19)
dim(lincRNA.hg19)
head(lincRNA.hg19)
}
\keyword{datasets}
