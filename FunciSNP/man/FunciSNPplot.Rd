\name{FunciSNPplot}
\alias{FunciSNPplot}

\title{
FunciSNPplot to visualize Func-y-SNP summary.
}
\description{
FunciSNPplot is a function developed to plot various types of plots to summarize and assist end-user in making informed discoveries of FunciSNP results. Plots can be stored in a folder for future reference.
}
\usage{
FunciSNPplot(dat, rsq = 0, split = FALSE, splitbysnp = FALSE, tagSummary = FALSE, heatmap = FALSE, genomicSum = FALSE, save = FALSE)
}

\arguments{
  \item{dat}{
             FunciSNP data.frame: \code{dat} is a data.frame object from FunciSNPAnnotateSummary. Need to run \code{\link{FunciSNPAnnotateSummary}} first.
}
  \item{rsq}{
            an interger (0-1): \code{rsq} is the Rsquared cutoff used to subset.
}
  \item{split}{
		logical: \code{split} will generate distribution plot of all Correlated SNPs by Rsquare values.  
}
  \item{splitbysnp}{
		logical: \code{splitbysnp} is similar to \code{split} but instead split the distribution by tagSNP.
}
  \item{tagSummary}{
		logical: \code{tagSummary} here~~
}
  \item{heatmap}{
		logical: \code{heatmap} correlation heatmap to visualize the number of correlated SNPs at each tagSNP overlapping each biological feature.
}
  \item{genomicSum}{
		logical: \code{genomicSum} Stacked bar chart summarizing all correlated SNPs for each of the identified genomie features (exon, intron, 5'UTR, 3'UTR, promoter, lincRNA or in gene desert (intergentic))
}
  \item{save}{
		logical: \code{save} to save outputs to folder. Set at getwd(), in folder 'FunciSNP.VERSION/plots
}
}
\details{

}
\value{
All outputs are plots in standard x11().
}
\references{
	Coetzee SG et al. submitted for review. 2012
}
\author{
	Simon Coetzee, Houtan Noushmehr
}
\note{

}



\seealso{
	\code{\link{FunciSNP}}, \code{\link{FunciSNPplot}}, \code{\link{FunciSNPAnnotateSummary}}, \code{\link{FunciSNPtable}}, \code{\link{FunciSNPbed}}
}
\examples{

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
