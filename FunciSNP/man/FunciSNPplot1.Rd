\name{FunciSNPplot1}
\alias{FunciSNPplot1}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{titlehere
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
FunciSNPplot1(summary.snp.list, splitbysnp = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{summary.snp.list}{
%%     ~~Describe \code{summary.snp.list} here~~
}
  \item{splitbysnp}{
%%     ~~Describe \code{splitbysnp} here~~
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
function (summary.snp.list, splitbysnp = TRUE) 
{
    if (splitbysnp == TRUE) {
        theme_white()
        ggplot(summary.snp.list, aes(x = R.squared)) + geom_histogram(binwidth = 0.05) + 
            geom_vline(xintercept = 0.5, linetype = 2) + scale_x_continuous("Correlated SNPs R squared to Tag SNP (0-1)") + 
            scale_y_continuous("Total # of Correlated SNPs associated with tagSNP") + 
            opts(legend.position = "none", axis.text.y = theme_text(), 
                axis.text.x = theme_text(angle = 90), title = "Distribution of correlated SNPs for each tagSNP\nat Rsquared values") + 
            facet_wrap(chromosome ~ tag.snp.id)
    }
    else {
        tt <- count(df = summary.snp.list, vars = "R.squared")
        tt <- na.omit(tt)
        ht <- range(tt[, "freq"])[2] + 10
        plot(tt, xlim = c(0, 1), ylim = c(0, ht), pch = "*", 
            main = paste("Distribution of Correlated SNPs at Rsquared values\n(n=", 
                length(na.omit(summary.snp.list[, "R.squared"])), 
                "; ", round(c(length(na.omit(summary.snp.list[, 
                  "R.squared"]))/dim(summary.snp.list)[1]) * 
                  100, digits = 1), "\%)", sep = ""), xlab = "Correlated SNPs R squared to Tag SNP (0-1)", 
            ylab = "Count of Correlated SNP")
        abline(v = 0.1, lty = 2, col = "red")
        abline(v = 0.2, lty = 2, col = "red")
        abline(v = 0.3, lty = 2, col = "red")
        abline(v = 0.4, lty = 2, col = "red")
        abline(v = 0.5, lty = 2, col = "black")
        abline(v = 0.6, lty = 2, col = "green")
        abline(v = 0.7, lty = 2, col = "green")
        abline(v = 0.8, lty = 2, col = "green")
        abline(v = 0.9, lty = 2, col = "green")
        abline(h = ht - 7, col = "black", lty = 2)
        text(0.05, ht - 5, as.character(dim(subset(summary.snp.list, 
            R.squared < 0.1))[1]))
        text(0.15, ht - 5, as.character(dim(subset(summary.snp.list, 
            R.squared >= 0.1 & R.squared < 0.2))[1]))
        text(0.25, ht - 5, as.character(dim(subset(summary.snp.list, 
            R.squared >= 0.2 & R.squared < 0.3))[1]))
        text(0.35, ht - 5, as.character(dim(subset(summary.snp.list, 
            R.squared >= 0.3 & R.squared < 0.4))[1]))
        text(0.45, ht - 5, as.character(dim(subset(summary.snp.list, 
            R.squared >= 0.4 & R.squared < 0.5))[1]))
        text(0.55, ht - 5, as.character(dim(subset(summary.snp.list, 
            R.squared >= 0.5 & R.squared < 0.6))[1]))
        text(0.65, ht - 5, as.character(dim(subset(summary.snp.list, 
            R.squared >= 0.6 & R.squared < 0.7))[1]))
        text(0.75, ht - 5, as.character(dim(subset(summary.snp.list, 
            R.squared >= 0.7 & R.squared < 0.8))[1]))
        text(0.85, ht - 5, as.character(dim(subset(summary.snp.list, 
            R.squared >= 0.8 & R.squared < 0.9))[1]))
        text(0.95, ht - 5, as.character(dim(subset(summary.snp.list, 
            R.squared >= 0.9))[1]))
    }
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
