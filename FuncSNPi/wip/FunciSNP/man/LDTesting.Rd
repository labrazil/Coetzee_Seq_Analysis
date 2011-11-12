\name{LDTesting}
\alias{LDTesting}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
LDTesting(snps, snp.names.chosen, ethno.chosen, snp.ld.frame, bf.count, snp.regions, var.data)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{snps}{
%%     ~~Describe \code{snps} here~~
}
  \item{snp.names.chosen}{
%%     ~~Describe \code{snp.names.chosen} here~~
}
  \item{ethno.chosen}{
%%     ~~Describe \code{ethno.chosen} here~~
}
  \item{snp.ld.frame}{
%%     ~~Describe \code{snp.ld.frame} here~~
}
  \item{bf.count}{
%%     ~~Describe \code{bf.count} here~~
}
  \item{snp.regions}{
%%     ~~Describe \code{snp.regions} here~~
}
  \item{var.data}{
%%     ~~Describe \code{var.data} here~~
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
function (snps, snp.names.chosen, ethno.chosen, snp.ld.frame, 
    bf.count, snp.regions, var.data) 
{
    try(snp.ld <- LD(snps), silent = TRUE)
    if (exists("snp.ld")) {
        if (snp.names.chosen \%in\% colnames(snp.ld$D)) {
            d <- na.omit(as.data.frame(snp.ld$D[, c(snp.names.chosen)]))
            colnames(d) <- "d"
            D.prime <- na.omit(as.data.frame(snp.ld$"D'"[, c(snp.names.chosen)]))
            colnames(D.prime) <- "D.prime"
            r <- na.omit(as.data.frame(snp.ld$r[, c(snp.names.chosen)]))
            colnames(r) <- "r"
            R.squared <- na.omit(as.data.frame(snp.ld$"R^2"[, 
                c(snp.names.chosen)]))
            colnames(R.squared) <- "R.squared"
            n <- na.omit(as.data.frame(snp.ld$n[, c(snp.names.chosen)]))
            colnames(n) <- "n"
            chi.squared <- na.omit(as.data.frame(snp.ld$"X^2"[, 
                c(snp.names.chosen)]))
            colnames(chi.squared) <- "chi.squared"
            p.val <- na.omit(as.data.frame(snp.ld$"P-value"[, 
                c(snp.names.chosen)]))
            colnames(p.val) <- "p.val"
            snp.names.cor <- dimnames(p.val)[[1]]
            if (identical(snp.ld.frame, NULL)) {
                snp.ld.frame <- data.frame(snp.names.chosen, 
                  snp.names.cor, d, D.prime, r, R.squared, n, 
                  chi.squared, p.val)
                snp.ld.frame$feature <- unlist(strsplit(bio.features[bf.count], 
                  "\\.bed"))
                snp.ld.frame$loc.feature <- paste(seqnames(features.file.interval[c(as.matrix(findOverlaps(features.file.interval, 
                  close.snp.ranges[snp.names.cor]))[, 1])]), 
                  ":", start(features.file.interval[c(as.matrix(findOverlaps(features.file.interval, 
                    close.snp.ranges[snp.names.cor]))[, 1])]), 
                  "-", end(features.file.interval[c(as.matrix(findOverlaps(features.file.interval, 
                    close.snp.ranges[snp.names.cor]))[, 1])]), 
                  sep = "")
                snp.ld.frame$eth.risk <- snp.regions$snp.ethno[snp.regions$snp.name == 
                  snp.names.chosen][1]
                snp.ld.frame$eth.chosen <- ethno.chosen
                snp.ld.frame$chr <- snp.regions$snp.chromosome[snp.regions$snp.name == 
                  snp.names.chosen][1]
                snp.ld.frame$loc.risk <- subset(var.data, ID == 
                  snp.names.chosen)$POS
                snp.ld.frame$loc.cor <- subset(var.data, var.data$ID \%in\% 
                  snp.names.cor)$POS
                snp.ld.frame$dist <- abs(snp.ld.frame$loc.risk - 
                  snp.ld.frame$loc.cor)
            }
            else {
                tmp.snp.ld.frame <- data.frame(snp.names.chosen, 
                  snp.names.cor, d, D.prime, r, R.squared, n, 
                  chi.squared, p.val)
                tmp.snp.ld.frame$feature <- unlist(strsplit(bio.features[bf.count], 
                  "\\.bed"))
                tmp.snp.ld.frame$loc.feature <- paste(seqnames(features.file.interval[c(as.matrix(findOverlaps(features.file.interval, 
                  close.snp.ranges[snp.names.cor]))[, 1])]), 
                  ":", start(features.file.interval[c(as.matrix(findOverlaps(features.file.interval, 
                    close.snp.ranges[snp.names.cor]))[, 1])]), 
                  "-", end(features.file.interval[c(as.matrix(findOverlaps(features.file.interval, 
                    close.snp.ranges[snp.names.cor]))[, 1])]), 
                  sep = "")
                tmp.snp.ld.frame$eth.risk <- snp.regions$snp.ethno[snp.regions$snp.name == 
                  snp.names.chosen][1]
                tmp.snp.ld.frame$eth.chosen <- ethno.chosen
                tmp.snp.ld.frame$chr <- snp.regions$snp.chromosome[snp.regions$snp.name == 
                  snp.names.chosen][1]
                tmp.snp.ld.frame$loc.risk <- subset(var.data, 
                  ID == snp.names.chosen)$POS
                tmp.snp.ld.frame$loc.cor <- subset(var.data, 
                  var.data$ID \%in\% snp.names.cor)$POS
                tmp.snp.ld.frame$dist <- abs(tmp.snp.ld.frame$loc.risk - 
                  tmp.snp.ld.frame$loc.cor)
                snp.ld.frame <- rbind(snp.ld.frame, tmp.snp.ld.frame)
            }
            rownames(snp.ld.frame) <- NULL
            snp.ld.frame
        }
        else {
            snp.ld.frame
        }
    }
    else {
        snp.ld.frame
    }
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
