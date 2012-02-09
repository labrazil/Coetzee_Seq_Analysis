\name{FunciSNPbedplaceholder<<<<<}
\alias{FunciSNPbedplaceholder<<<<<}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
placeholder<<<<<}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
placeholder<<<<<}
\usage{
FunciSNPbed(dat, rsq, path = getwd(), filename = NULL)
placeholder<<<<<}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{datplaceholder<<<<<}{
%%     ~~Describe \code{datplaceholder<<<<<} here~~
placeholder<<<<<}
  \item{rsqplaceholder<<<<<}{
%%     ~~Describe \code{rsqplaceholder<<<<<} here~~
placeholder<<<<<}
  \item{pathplaceholder<<<<<}{
%%     ~~Describe \code{pathplaceholder<<<<<} here~~
placeholder<<<<<}
  \item{filenameplaceholder<<<<<}{
%%     ~~Describe \code{filenameplaceholder<<<<<} here~~
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
function (dat, rsq, path = getwd(), filename = NULL) 
{
    if (identical(filename, NULL)) {
        filename <- paste("FunciSNP_results_rsq.", rsq, ".bed", 
            sep = "")
    placeholder<<<<<}
    else {
        filename <- filename
    placeholder<<<<<}
    d.s <- subset(dat, R.squared > rsq)
    d.cor <- d.s[which(!(duplicated(d.s[, "corr.snp.id"]))), 
        ]
    d.tag <- d.s[which(!(duplicated(d.s[, "tag.snp.id"]))), ]
    d.tag <- (d.tag[, c(1, 8, 8, 7, 10, 8, 8)])
    d.cor <- (d.cor[, c(1, 6, 6, 5, 10, 6, 6)])
    d.cor$chromosome <- paste("chr", d.cor$chromosome, sep = "")
    d.tag$chromosome <- paste("chr", d.tag$chromosome, sep = "")
    d.cor[, 3] <- d.cor[, 3] + 1
    d.tag[, 3] <- d.tag[, 3] + 1
    d.cor[, 7] <- d.cor[, 7] + 1
    d.tag[, 7] <- d.tag[, 7] + 1
    d.cor$strand <- "+"
    d.tag$strand <- "+"
    d.cor$color <- "255,0,0"
    d.tag$color <- "0,0,0"
    d.cor <- d.cor[, c(1:5, 8, 6:7, 9)]
    d.tag <- d.tag[, c(1:5, 8, 6:7, 9)]
    dimnames(d.cor)[[2]] <- c("chr", "snp.pos.s", "snp.pos.e", 
        "snp.id", "rsquare", "strand", "snp.pos.s", "snp.pos.e", 
        "color")
    dimnames(d.tag)[[2]] <- dimnames(d.cor)[[2]]
    y <- rbind(d.tag, d.cor)
    y$rsquare <- round(y$rsquare, digits = 4)
    con <- file(paste(path, filename, sep = "/"), open = "wt")
    writeLines(paste("browser position chr", d.s[1, 1], ":", 
        if (d.s[1, 6] < d.s[1, 8]) {
            paste(d.s[1, 6] - 500, "-", d.s[1, 8] + 500, sep = "")
        placeholder<<<<<}
        else {
            paste(d.s[1, 8] - 500, "-", d.s[1, 6] + 500, sep = "")
        placeholder<<<<<}, "\ntrack name=\"FunciSNP_results\" description=\"FunciSNP ver. ", 
        package.version("FunciSNP"), " Results\" visibility=3 itemRgb=\"On\"", 
        sep = ""), con)
    write.table(y, row.names = F, col.names = F, sep = "\t", 
        file = con, quote = F)
    close(con)
    message("####\nBed file \"", filename, "\" created successfully.\n(See folder: \"", 
        path, "\")")
    cat("Total corSNP (RED): ", dim(d.cor)[1], "\nTotal tagSNP (BLK): ", 
        dim(d.tag)[1], "\n")
    message("\nTo view results, submit bed file as a\n  custom track in UCSC Genome Browser (genome.ucsc.edu), \n\nNow have fun with your new Func-y SNP", 
        if (dim(d.cor)[1] > 1) {
            "s"
        placeholder<<<<<}
        else {
            ""
        placeholder<<<<<}, "!!\n####")
  placeholder<<<<<}
placeholder<<<<<}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 placeholder<<<<<}
\keyword{ ~kwd2 placeholder<<<<<}% __ONLY ONE__ keyword per line
