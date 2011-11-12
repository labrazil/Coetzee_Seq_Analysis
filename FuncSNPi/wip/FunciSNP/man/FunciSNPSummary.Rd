\name{FunciSNPSummary}
\alias{FunciSNPSummary}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
FunciSNPSummary(R.2, dat)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{R.2}{
%%     ~~Describe \code{R.2} here~~
}
  \item{dat}{
%%     ~~Describe \code{dat} here~~
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
function (R.2, dat) 
{
    total.tagSNPs <- length(unique(dat[, "snp.names.chosen"]))
    total.1kSNPs <- length(unique(dat[, "snp.names.cor"]))
    total.feature <- length(unique(subset(dat, R.squared > R.2)[, 
        "feature"]))
    total.tagSNPs.cutoff <- length(unique(subset(dat, R.squared > 
        R.2)[, "snp.names.chosen"]))
    total.1kSNPs.cutoff <- length(unique(subset(dat, R.squared > 
        R.2)[, "snp.names.cor"]))
    total.feature.cutoff <- length(unique(subset(dat, R.squared > 
        R.2)[, "feature"]))
    total.dat <- matrix(c(total.tagSNPs, total.1kSNPs, total.feature, 
        total.tagSNPs.cutoff, total.1kSNPs.cutoff, total.feature.cutoff), 
        nrow = 3, ncol = 2, byrow = FALSE, dimnames = list(c("tagSNPs", 
            "1kSNPs", "bio.features"), c("Total", paste("R.squared.cuff.", 
            R.2, sep = ""))))
    total.dat <- as.data.frame(total.dat)
    total.dat$Percent <- round((total.dat[, 2]/total.dat[, 1]) * 
        100, 2)
    return(total.dat)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
