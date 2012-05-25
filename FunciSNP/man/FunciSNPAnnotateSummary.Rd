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
	SG. Coetzee, SK. Rhie, BP. Berman, GA. Coetzee and H. Noushmehr, FunciSNP: An R/Bioconductor Tool Integrating Functional Non-coding Datasets with Genetic Association Studies to Identify Candidate Regulatory SNPs., Nucleic Acids Research, In press, 2012 (doi:10.1093/nar/gks542).
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
data(glioma);
gl <- FunciSNPAnnotateSummary(glioma);
dim(gl)
head(gl)
names(gl)
}
\keyword{ annotation }
