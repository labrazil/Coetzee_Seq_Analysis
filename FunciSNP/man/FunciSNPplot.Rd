\name{FunciSNPplot}
\alias{FunciSNPplot}

\title{
FunciSNPplot to visualize Func-y-SNP summary.
}
\description{
FunciSNPplot is a function developed to plot various types of plots to summarize and assist end-user in making informed discoveries of FunciSNP results. Plots can be stored in a folder for future reference.
}
\usage{
FunciSNPplot(dat, rsq = 0, split = FALSE, splitbysnp = FALSE, tagSummary = FALSE, heatmap = FALSE, heatmap.key = FALSE,genomicSum = FALSE, save = FALSE, pathplot=getwd(), text.size=10, save.width=7, save.height=7)
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
		logical: \code{tagSummary} Will output two plots per biofeature. The first one is a scatter plot showing the relationship between Rsquare and Distance to tagSNP for each Func-y-SNP. The second plot is a histogram distribution of number of correlated SNPs at each Rsquare value. Each set of plot is further divided by tagSNP. Best if used with \code{rsq} value.
}
  \item{heatmap}{
		logical: \code{heatmap} correlation heatmap to visualize the number of correlated SNPs at each tagSNP overlapping each biological feature. Most informative if used with a \code{rsq} value.
}
  \item{heatmap.key}{
		logical: \code{heatmap.key} Places the count of each cell on the heatmap.
}
  \item{genomicSum}{
		logical: \code{genomicSum} Stacked bar chart summarizing all correlated SNPs for each of the identified genomie features (exon, intron, 5'UTR, 3'UTR, promoter, lincRNA or in gene desert (intergentic)). Most informative if used with a \code{rsq} value.
}
  \item{save}{
		logical: \code{save} to save outputs to folder. Set at getwd(), in folder 'FunciSNP.VERSION/plots
}
  \item{pathplot}{
           a character: \code{pathplot} is the path to the folder where to save the plots. Default to getwd() or current working directory.
}
  \item{text.size}{
    Text size passed to graphing functions
  }
  \item{save.width}{
    Width of saved images in inches
  }
  \item{save.height}{
    Height of saved images in inches
  }
}
\details{
NA
}
\value{
Plots are generated either in X11 or in specified folder.
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
data(glioma)
gl <- FunciSNPAnnotateSummary(glioma)
FunciSNPplot(gl)
FunciSNPplot(gl, rsq=0, genomicSum=TRUE, save=FALSE)
FunciSNPplot(gl, rsq=0.5, genomicSum=TRUE, save=FALSE)
# DO NOT RUN
#FunciSNPplot(gl, tagSummary=TRUE, rsq=0.5)
#
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
