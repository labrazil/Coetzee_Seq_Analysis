\name{FunciSNPplot2}
\alias{FunciSNPplot2}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{titlehere
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
FunciSNPplot2(R.2, summary.snp.list)
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
    try(dir.create(path = "funcisnp_results/plots", showWarnings = FALSE, 
        recursive = TRUE), silent = TRUE)
    require("ggplot2")
    cat("\n\nPlotting will begin ")
    cat("using R square cut off of", R.2, "\n")
    theme_white()
    all.s <- try(subset(summary.snp.list, R.squared >= R.2), 
        silent = TRUE)
    all.ss <- try(subset(summary.snp.list, R.squared < R.2), 
        silent = TRUE)
    try(all.s$r2 <- c("Yes"), silent = TRUE)
    try(all.ss$r2 <- c("No"), silent = TRUE)
    if (nrow(all.s) > 0 && nrow(all.ss) > 0) {
        all <- try(rbind(all.s, all.ss), silent = TRUE)
    }
    else {
        if (nrow(all.s) > 0 && nrow(all.ss) <= 0) {
            all <- all.s
        }
        if (nrow(all.s) <= 0 && nrow(all.ss) > 0) {
            all <- all.ss
        }
        if (nrow(all.s) <= 1 && nrow(all.ss) <= 0) {
            return()
        }
    }
    for (i in 1:length(summary(as.factor(all[, "bio.feature"])))) {
        bio <- names(summary(as.factor(all[, "bio.feature"])))
        tmp <- subset(all, bio.feature == bio[i])
        ggplot(tmp, aes(x = R.squared, fill = factor(r2))) + 
            geom_histogram() + geom_vline(xintercept = R.2, linetype = 2) + 
            scale_x_continuous("Rsquare Values (0-1)", limits = c(0, 
                1)) + scale_y_continuous("Total # of Surrogate SNPs associated with riskSNP") + 
            scale_fill_manual(values = c(Yes = "Red", No = "Black")) + 
            opts(legend.position = "none", axis.text.y = theme_text(), 
                axis.text.x = theme_text(angle = 90), title = paste("riskSNP\nOverlapping: ", 
                  bio[i], sep = "")) + facet_wrap(~tag.snp.id)
        ggplot(tmp, aes(x = R.squared, y = dist, colour = r2, 
            size = factor(r2))) + geom_point() + geom_vline(xintercept = R.2, 
            linetype = 2) + scale_x_continuous("Rsquare Values (0-1)", 
            limits = c(0, 1)) + scale_y_continuous("Distance to Surrogate SNPs associated with riskSNP (bp)", 
            formatter = "comma") + scale_colour_manual(values = c(Yes = "Red", 
            No = "Black")) + scale_size_manual(values = c(Yes = 2, 
            No = 1)) + opts(legend.position = "none", axis.text.y = theme_text(), 
            axis.text.x = theme_text(angle = 90), title = paste("Distance between riskSNP\nand Surrogate SNP\nOverlapping: ", 
                bio[i], sep = "")) + facet_wrap(~tag.snp.id)
        cat("Finished plotting ", i, "/", length(bio), "\n")
    }
    cat("\n\nPlots are done, see 'funcisnp_results/plots/' folder in ", 
        getwd(), "\n\n")
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
