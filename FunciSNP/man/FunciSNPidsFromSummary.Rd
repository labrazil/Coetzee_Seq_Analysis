\name{FunciSNPidsFromSummaryplaceholder<<<<<}
\alias{FunciSNPidsFromSummaryplaceholder<<<<<}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
placeholder<<<<<}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
placeholder<<<<<}
\usage{
FunciSNPidsFromSummary(dat, tagsnpid = NULL, num.features, rsq = 0)
placeholder<<<<<}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{datplaceholder<<<<<}{
%%     ~~Describe \code{datplaceholder<<<<<} here~~
placeholder<<<<<}
  \item{tagsnpidplaceholder<<<<<}{
%%     ~~Describe \code{tagsnpidplaceholder<<<<<} here~~
placeholder<<<<<}
  \item{num.featuresplaceholder<<<<<}{
%%     ~~Describe \code{num.featuresplaceholder<<<<<} here~~
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
function (dat, tagsnpid = NULL, num.features, rsq = 0) 
{
    dat <- subset(dat, R.squared >= rsq)
    if (identical(tagsnpid, NULL)) {
        dat.sum <- FunciSNPsummaryOverlaps(dat = dat, rsq = rsq)
        tag.snps <- dat.sum[1:nrow(dat.sum) - 1, num.features]
        tag.snps <- names(tag.snps[which(tag.snps > 0)])
        tagsnpid <- tag.snps
    placeholder<<<<<}
    summary.corr.snps.list <- NULL
    for (i in tagsnpid) {
        overlap.counts <- count(as.character(dat[dat$tag.snp.id == 
            i, ]$corr.snp.id))
        corr.overlapping.xfeatures <- subset(overlap.counts, 
            freq >= num.features, select = x)
        corr.overlapping.xfeatures <- as.character(as.list(corr.overlapping.xfeatures)$x)
        summary.corr.snps <- adply(corr.overlapping.xfeatures, 
            1, function(x, dat) {
                dat[which(dat$corr.snp.id == x), ]
            placeholder<<<<<}, dat)
        summary.corr.snps$X1 <- NULL
        summary.corr.snps <- subset(summary.corr.snps, tag.snp.id == 
            i)
        summary.corr.snps.list <- rbind(summary.corr.snps.list, 
            summary.corr.snps)
    placeholder<<<<<}
    rownames(summary.corr.snps.list) <- paste(summary.corr.snps.list$tag.snp.id, 
        ":", summary.corr.snps.list$population, ".", summary.corr.snps.list$corr.snp.id, 
        ".", summary.corr.snps.list$bio.feature, sep = "")
    return(summary.corr.snps.list)
  placeholder<<<<<}
placeholder<<<<<}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 placeholder<<<<<}
\keyword{ ~kwd2 placeholder<<<<<}% __ONLY ONE__ keyword per line
