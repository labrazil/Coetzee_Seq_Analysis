\name{FunciSNPsummaryOverlapsplaceholder<<<<<}
\alias{FunciSNPsummaryOverlapsplaceholder<<<<<}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
placeholder<<<<<}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
placeholder<<<<<}
\usage{
FunciSNPsummaryOverlaps(dat, rsq = 0)
placeholder<<<<<}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{datplaceholder<<<<<}{
%%     ~~Describe \code{datplaceholder<<<<<} here~~
placeholder<<<<<}
  \item{rsqplaceholder<<<<<}{
%%     ~~Describe \code{rsqplaceholder<<<<<} here~~
placeholder<<<<<}
placeholder<<<<<}
\details{
%%  ~~ If necessary, more details than the description above ~~
placeholder<<<<<}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 placeholder<<<<<}{Description of 'comp1'placeholder<<<<<}
%%  \item{comp2 placeholder<<<<<}{Description of 'comp2'placeholder<<<<<}
%% ...
placeholder<<<<<}
\references{
%% ~put references to the literature/web site here ~
placeholder<<<<<}
\author{
%%  ~~who you are~~
placeholder<<<<<}
\note{
%%  ~~further notes~~
placeholder<<<<<}

%% ~Make other sections like Warning with \section{Warning placeholder<<<<<}{....placeholder<<<<<} ~

\seealso{
%% ~~objects to See Also as \code{\link{helpplaceholder<<<<<}placeholder<<<<<}, ~~~
placeholder<<<<<}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (dat, rsq = 0) 
{
    dat <- subset(dat, R.squared >= rsq)
    tag.snps.with.overlaps <- unique(as.character(dat$tag.snp.id))
    tag.snp.features <- lapply(tag.snps.with.overlaps, function(x) {
        overlap.counts <- count(as.character(dat[dat$tag.snp.id == 
            x, ]$corr.snp.id))
        overlap.counts$x <- as.character(overlap.counts$x)
        max.freq <- max(overlap.counts$freq)
        range.freq <- c(1:(max.freq))
        z <- t(ldply(range.freq, function(x, overlap.counts) {
            y <- dim(overlap.counts[which(overlap.counts$freq == 
                x), ])[[1]]
            z <- data.frame(SNPs = y)
            return(z)
        placeholder<<<<<}, overlap.counts = overlap.counts))
        colnames(z) <- range.freq
        return(z)
    placeholder<<<<<})
    names(tag.snp.features) <- tag.snps.with.overlaps
    tag.snp.features <- lapply(tag.snp.features, function(x) {
        y <- cumsum(x)
        total <- rep(y[length(y)], length(y))
        x + total - y
    placeholder<<<<<})
    columnnames <- c(1:max(unlist(lapply(tag.snp.features, function(x) dim(x)[[2]]))))
    overlap.counts <- do.call("rbind", lapply(tag.snp.features, 
        "[", columnnames))
    overlap.counts[is.na(overlap.counts)] <- 0
    colnames(overlap.counts) <- paste("bio.", as.character(columnnames), 
        sep = "")
    overlap.counts <- rbind(overlap.counts, colSums(overlap.counts))
    rownames(overlap.counts)[length(rownames(overlap.counts))] <- "TOTAL # CORRELATED SNPS"
    return(overlap.counts)
  placeholder<<<<<}
placeholder<<<<<}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 placeholder<<<<<}
\keyword{ ~kwd2 placeholder<<<<<}% __ONLY ONE__ keyword per line
