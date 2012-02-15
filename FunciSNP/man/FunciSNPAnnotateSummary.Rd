\name{FunciSNPAnnotateSummary}
\alias{FunciSNPAnnotateSummary}

\title{
Genomic Annotation of Func-y-SNPs.
}
\description{
This will annotate all identified Func-y-SNP for it's distance to the nearest known TSS, whether it overlapps a known exon, intron, 5'UTR, 3'UTR, promoter, lincRNA or in gene desert (intergentic) regions.
}
\usage{
FunciSNPAnnotateSummary(snp.list)
}

\arguments{
  \item{snp.list}{
          a FunciSNP object: \code{snp.list} represents the FunciSNP object output from FunciSNP. See \code{\link{getFSNPs}}.
}
}
\details{
All known genomic features (exon, intron, 5'UTR, 3'UTR, promoter, lincRNA or in gene desert (intergentic)) are used to annotate the newly identified Func-y-SNP. Information described in this data.frame is used for all summary plots, table, and bed file generations.
}
\value{
data.frame with rows for each correlated SNP.
}
\references{
	Coetzee SG et al. submitted for review. 2012
}
\author{
	Simon Coetzee, Houtan Noushmehr
}
\note{
NA
}



\seealso{
	\code{\link{getFSNPs}}, \code{\link{FunciSNPplot}}, \code{\link{FunciSNPAnnotateSummary}}, \code{\link{FunciSNPtable}}, \code{\link{FunciSNPbed}}
}
\examples{
data(glioma);
gl <- FunciSNPAnnotateSummary(glioma);
dim(gl)
head(gl)
names(gl)
}
\keyword{ annotation }
