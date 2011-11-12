\name{FilterByFeatures}
\alias{FilterByFeatures}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
FilterByFeatures(features.file, snp.names.chosen, ethno.chosen, variants.data)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{features.file}{
%%     ~~Describe \code{features.file} here~~
}
  \item{snp.names.chosen}{
%%     ~~Describe \code{snp.names.chosen} here~~
}
  \item{ethno.chosen}{
%%     ~~Describe \code{ethno.chosen} here~~
}
  \item{variants.data}{
%%     ~~Describe \code{variants.data} here~~
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
function (features.file, snp.names.chosen, ethno.chosen, variants.data) 
{
    close.snp.ranges <- GRanges(seqnames = paste("chr", as.character(variants.data$CHROM), 
        sep = ""), ranges = (IRanges(start = as.integer(variants.data$POS), 
        width = 1)))
    names(close.snp.ranges) <- variants.data$ID
    features.file.interval <<- import(features.file, asRangedData = FALSE)
    snps.included <- !is.na(match(close.snp.ranges, features.file.interval) > 
        0)
    if (sum(snps.included) >= 1) {
        snps.included <- subset(variants.data, snps.included)
        snp.geno <- data.frame(t(snps.included[, c(6:length(snps.included))]))
        snp.geno <- as.data.frame(snp.geno[, colSums(is.na(snp.geno)) < 
            nrow(snp.geno)])
        colnames(snp.geno) <- snps.included$ID
        if (snp.names.chosen \%in\% colnames(snp.geno)) {
            cat(snp.names.chosen, ", the risk snp, already overlaps with the feature", 
                features.file, "\n")
        }
        else {
            temp <- data.frame(t(subset(variants.data, ID == 
                snp.names.chosen))[6:length(variants.data), ])
            colnames(temp) <- snp.names.chosen
            snp.geno <- cbind(snp.geno, temp)
            rm(temp)
        }
        for (i in 1:ncol(snp.geno)) {
            snp.geno[, i] <- genotype(snp.geno[, i])
        }
        close.snp.ranges <<- close.snp.ranges
        snp.geno
    }
    else {
        NULL
        close.snp.ranges <<- close.snp.ranges
    }
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
