\name{FunciSNPsummaryOverlaps}
\alias{FunciSNPsummaryOverlaps}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Summarizes FunciSNP results by overlaps}
\description{
This will summarize the total number of YAFSNPs identified by the number of biofeatures per tagSNP}
\usage{
FunciSNPsummaryOverlaps(dat, rsq = 0)
}

\arguments{
  \item{dat}{FunciSNP data object}
  \item{rsq}{Rsquare value used to filter the summary report}

}
\details{
}
\value{
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
FunciSNPsummaryOveralps(glioma, rsq = 0.2);
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1}
\keyword{ ~kwd2}% __ONLY ONE__ keyword per line
