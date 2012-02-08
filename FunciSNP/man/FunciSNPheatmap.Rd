\name{FunciSNPheatmap}
\alias{FunciSNPheatmap}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{titlehere
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
FunciSNPheatmap(R.2, summary.snp.list)
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
    require("gplots")
    require("matlab")
    try(dir.create(path = "funcisnp_results/plots", showWarnings = FALSE, 
        recursive = TRUE), silent = TRUE)
    cat("Heatmap generation will begin")
    cat(" using R square cut off of ", R.2, "\n")
    all.s <- (table(subset(summary.snp.list, R.squared >= R.2)[, 
        "bio.feature"], subset(summary.snp.list, R.squared >= 
        R.2)[, "tag.snp.id"]))
    all.s <- as.matrix(all.s)
    png(filename = "funcisnp_results/plots/heatmap.png", bg = "white")
    heatmap.2(all.s, na.rm = TRUE, scale = "none", col = jet.colors(max(all.s, 
        na.rm = T)/2), key = T, symkey = FALSE, density.info = "none", 
        trace = "none", Rowv = F, Colv = T, cexRow = 1, cexCol = 1, 
        keysize = 1, dendrogram = c("col"), main = paste("Rsquare values >= ", 
            R.2))
    dev.off()
    cat("\n\nHeatmap plot is complete, see 'funcisnp_results/plots/' folder in ", 
        getwd(), "\n\n")
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
