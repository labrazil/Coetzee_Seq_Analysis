\name{FunciSNPtable}
\alias{FunciSNPtable}

\title{
Will output a summary report from FunciSNP at specified Rsquare cut-offs.
}
\description{
Using a specified Rsquare value (0-1) to subset the data, a table is generated which summarizes the total number of YAFSNPs, associated tagSNPs, and number of overlapping biofeatures.
}
\usage{
FunciSNPtable(dat, rsq, geneSum = FALSE)
}

\arguments{
  \item{dat}{
            FunciSNP data.frame: \code{dat} is a data.frame object from FunciSNPAnnotateSummary. Need to run \code{\link{FunciSNPAnnotateSummary}} first.
}
  \item{rsq}{
            an interger (0-1): \code{rsq} is the Rsquared cutoff used to subset.
}
  \item{geneSum}{
            logical: \code{geneSum} is set to FALSE. Setting to TRUE will output a list of Gene names which are nearest to the YAFSNP.
}
}
\details{
Using a specified Rsquare value (0-1) to subset the data, a table is generated which summarizes the total number of YAFSNPs, associated tagSNPs, and number of overlapping biofeatures. This will provide user a first look at the total number of available YAFSNP at a particular Rsquare cutoff.  If geneSum is set to TRUE, a list of gene names is reported instead.
}
\value{
Standard output which summarizes FunciSNP results.
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
FunciSNPtable(gl, rsq=0.5);
FunciSNPtable(gl, rsq=0.5, geneSum=TRUE);
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
