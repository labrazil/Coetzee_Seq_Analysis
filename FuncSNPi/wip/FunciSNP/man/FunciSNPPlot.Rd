\name{FunciSNPPlot}
\alias{FunciSNPPlot}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
FunciSNPPlot(R.2, dat)
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
    try(dir.create(path = "funcisnp_results/plots", showWarnings = FALSE, 
        recursive = TRUE), silent = TRUE)
    require("ggplot2")
    cat("\n\nPlotting will begin ")
    cat("using R square cut off of", R.2, "\n")
    theme_white <- function() {
        theme_update(plot.background = theme_blank(), panel.background = theme_rect(colour = "black", 
            size = 1), axis.text.x = theme_text(colour = "black", 
            vjust = 1, size = 12), axis.text.y = theme_text(colour = "black", 
            hjust = 1, size = 12), axis.title.x = theme_text(colour = "black", 
            face = "bold", size = 12), axis.title.y = theme_text(colour = "black", 
            face = "bold", angle = 90, size = 12))
    }
    theme_white()
    all.s <- try(subset(dat, R.squared >= R.2), silent = TRUE)
    all.ss <- try(subset(dat, R.squared < R.2), silent = TRUE)
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
    for (i in 1:length(summary(as.factor(all[, "feature"])))) {
        bio <- names(summary(as.factor(all[, "feature"])))
        tmp <- subset(all, feature == bio[i])
        p.all <- ggplot(tmp, aes(x = R.squared, fill = factor(r2)))
        p.all + geom_histogram() + geom_vline(xintercept = R.2, 
            linetype = 2) + scale_x_continuous("Rsquare Values (0-1)", 
            limits = c(0, 1)) + scale_y_continuous("Total # of Surrogate SNPs associated with riskSNP") + 
            scale_fill_manual(values = c(Yes = "Red", No = "Black")) + 
            opts(legend.position = "none", axis.text.y = theme_text(), 
                axis.text.x = theme_text(angle = 90), title = paste("riskSNP\nOverlapping: ", 
                  bio[i], sep = "")) + facet_wrap(~snp.names.chosen)
        ggsave(file = paste("funcisnp_results/plots/", bio[i], 
            "_R2summary_riskSNP.pdf", sep = ""))
        p.all.d <- ggplot(tmp, aes(x = R.squared, y = dist, colour = r2, 
            size = factor(r2)))
        p.all.d + geom_point() + geom_vline(xintercept = R.2, 
            linetype = 2) + scale_x_continuous("Rsquare Values (0-1)", 
            limits = c(0, 1)) + scale_y_continuous("Distance to Surrogate SNPs associated with riskSNP (bp)", 
            formatter = "comma") + scale_colour_manual(values = c(Yes = "Red", 
            No = "Black")) + scale_size_manual(values = c(Yes = 2, 
            No = 1)) + opts(legend.position = "none", axis.text.y = theme_text(), 
            axis.text.x = theme_text(angle = 90), title = paste("Distance between riskSNP\nand Surrogate SNP\nOverlapping: ", 
                bio[i], sep = "")) + facet_wrap(~snp.names.chosen)
        ggsave(file = paste("funcisnp_results/plots/", bio[i], 
            "_R2vsDist_riskSNP.pdf", sep = ""))
        cat("Finished plotting ", i, "/", length(bio), "\n")
    }
    cat("\n\nMaking plots is finished, see 'funcisnp_results/plots/' folder in ", 
        getwd(), "\n\n")
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
