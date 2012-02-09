\name{FunciSNPtableplaceholder<<<<<}
\alias{FunciSNPtableplaceholder<<<<<}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
placeholder<<<<<}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
placeholder<<<<<}
\usage{
FunciSNPtable(dat, rsq, geneSum = FALSE)
placeholder<<<<<}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{datplaceholder<<<<<}{
%%     ~~Describe \code{datplaceholder<<<<<} here~~
placeholder<<<<<}
  \item{rsqplaceholder<<<<<}{
%%     ~~Describe \code{rsqplaceholder<<<<<} here~~
placeholder<<<<<}
  \item{geneSumplaceholder<<<<<}{
%%     ~~Describe \code{geneSumplaceholder<<<<<} here~~
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
function (dat, rsq, geneSum = FALSE) 
{
    if (!geneSum) {
        total.tagSNPs <- length(unique(dat[, "tag.snp.id"]))
        total.1kSNPs <- length(unique(dat[, "corr.snp.id"]))
        total.feature <- length(unique(subset(dat, R.squared > 
            rsq)[, "bio.feature"]))
        total.tagSNPs.cutoff <- length(unique(subset(dat, R.squared > 
            rsq)[, "tag.snp.id"]))
        total.1kSNPs.cutoff <- length(unique(subset(dat, R.squared > 
            rsq)[, "corr.snp.id"]))
        total.feature.cutoff <- length(unique(subset(dat, R.squared > 
            rsq)[, "bio.feature"]))
        total.summary.snp.list <- matrix(c(total.tagSNPs, total.1kSNPs, 
            total.feature, total.tagSNPs.cutoff, total.1kSNPs.cutoff, 
            total.feature.cutoff), nrow = 3, ncol = 2, byrow = FALSE, 
            dimnames = list(c("tagSNPs", "1kSNPs", "bio.features"), 
                c("Total", paste("R.squared.cuff.", rsq, sep = ""))))
        total.summary.snp.list <- as.data.frame(total.summary.snp.list)
        total.summary.snp.list$Percent <- round((total.summary.snp.list[, 
            2]/total.summary.snp.list[, 1]) * 100, 2)
        return(total.summary.snp.list)
    placeholder<<<<<}
    else {
        dat.s <- subset(dat, R.squared > rsq)
        y <- colSums(table(dat.s$corr.snp.id, dat.s$nearest.TSS.GeneSymbol))
        y[y == 0] <- NA
        y <- na.omit(y)
        y <- as.matrix(names(y))
        y <- as.data.frame(y)
        dimnames(y)[[2]] <- "Gene_Names"
        return(y)
    placeholder<<<<<}
  placeholder<<<<<}
placeholder<<<<<}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 placeholder<<<<<}
\keyword{ ~kwd2 placeholder<<<<<}% __ONLY ONE__ keyword per line
