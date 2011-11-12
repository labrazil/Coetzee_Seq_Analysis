\name{ReadRegionsFile}
\alias{ReadRegionsFile}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
ReadRegionsFile(regions.file)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{regions.file}{
%%     ~~Describe \code{regions.file} here~~
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
function (regions.file) 
{
    snp.regions <- read.table(regions.file)
    snp.region.split <- unlist(strsplit(as.vector(snp.regions[, 
        1]), ":"))
    snp.chromosome <- grep("^\\d{0,1}[0-9X-Y]$", snp.region.split, 
        value = TRUE, perl = TRUE)
    m <- regexpr("^\\d{1,}-", snp.region.split)
    snp.region.start <- substr(snp.region.split, m, m + attr(m, 
        "match.length") - 1)
    snp.region.start <- as.vector(na.omit(as.numeric(gsub("-$", 
        "", snp.region.start))))
    m <- regexpr("-\\d{1,}$", snp.region.split)
    snp.region.end <- substr(snp.region.split, m, m + attr(m, 
        "match.length") - 1)
    snp.region.end <- as.vector(na.omit(as.numeric(gsub("^-", 
        "", snp.region.end))))
    snp.name <- snp.regions[, 2]
    snp.ethno <- snp.regions[, 3]
    snp.regions <- data.frame(snp.chromosome, snp.region.start, 
        snp.region.end, snp.name, snp.ethno, stringsAsFactors = FALSE)
    snp.regions
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
