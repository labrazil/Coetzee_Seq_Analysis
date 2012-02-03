\name{FunciSNPtable}
\alias{FunciSNPtable}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{titlehere
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
FunciSNPtable(R.2, summary.snp.list)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{R.2}{
%%     ~~Describe \code{R.2} here~~
}
  \item{summary.snp.list}{
%%     ~~Describe \code{summary.snp.list} here~~
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (R.2, summary.snp.list) 
{
    total.tagSNPs <- length(unique(summary.snp.list[, "tag.snp.id"]))
    total.1kSNPs <- length(unique(summary.snp.list[, "corr.snp.id"]))
    total.feature <- length(unique(subset(summary.snp.list, R.squared > 
        R.2)[, "bio.feature"]))
    total.tagSNPs.cutoff <- length(unique(subset(summary.snp.list, 
        R.squared > R.2)[, "tag.snp.id"]))
    total.1kSNPs.cutoff <- length(unique(subset(summary.snp.list, 
        R.squared > R.2)[, "corr.snp.id"]))
    total.feature.cutoff <- length(unique(subset(summary.snp.list, 
        R.squared > R.2)[, "bio.feature"]))
    total.summary.snp.list <- matrix(c(total.tagSNPs, total.1kSNPs, 
        total.feature, total.tagSNPs.cutoff, total.1kSNPs.cutoff, 
        total.feature.cutoff), nrow = 3, ncol = 2, byrow = FALSE, 
        dimnames = list(c("tagSNPs", "1kSNPs", "bio.features"), 
            c("Total", paste("R.squared.cuff.", R.2, sep = ""))))
    total.summary.snp.list <- as.data.frame(total.summary.snp.list)
    total.summary.snp.list$Percent <- round((total.summary.snp.list[, 
        2]/total.summary.snp.list[, 1]) * 100, 2)
    return(total.summary.snp.list)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
