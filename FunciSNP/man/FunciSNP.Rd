\name{FunciSNPplaceholder<<<<<}
\alias{FunciSNPplaceholder<<<<<}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
placeholder<<<<<}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
placeholder<<<<<}
\usage{
FunciSNP(snp.regions.file, bio.features.loc = NULL, bio.features.TSS = TRUE, par.threads = detectCores()/2, verbose = par.threads < 2, method.p = "BH", reduce.by.features = TRUE, search.window = 2e+05)
placeholder<<<<<}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{snp.regions.fileplaceholder<<<<<}{
%%     ~~Describe \code{snp.regions.fileplaceholder<<<<<} here~~
placeholder<<<<<}
  \item{bio.features.locplaceholder<<<<<}{
%%     ~~Describe \code{bio.features.locplaceholder<<<<<} here~~
placeholder<<<<<}
  \item{bio.features.TSSplaceholder<<<<<}{
%%     ~~Describe \code{bio.features.TSSplaceholder<<<<<} here~~
placeholder<<<<<}
  \item{par.threadsplaceholder<<<<<}{
%%     ~~Describe \code{par.threadsplaceholder<<<<<} here~~
placeholder<<<<<}
  \item{verboseplaceholder<<<<<}{
%%     ~~Describe \code{verboseplaceholder<<<<<} here~~
placeholder<<<<<}
  \item{method.pplaceholder<<<<<}{
%%     ~~Describe \code{method.pplaceholder<<<<<} here~~
placeholder<<<<<}
  \item{reduce.by.featuresplaceholder<<<<<}{
%%     ~~Describe \code{reduce.by.featuresplaceholder<<<<<} here~~
placeholder<<<<<}
  \item{search.windowplaceholder<<<<<}{
%%     ~~Describe \code{search.windowplaceholder<<<<<} here~~
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
function (snp.regions.file, bio.features.loc = NULL, bio.features.TSS = TRUE, 
    par.threads = detectCores()/2, verbose = par.threads < 2, 
    method.p = "BH", reduce.by.features = TRUE, search.window = 2e+05) 
{
    message("\n", "| | _  |  _  _ __  _    _|_ _ \n", "|^|(/_ | (_ (_)|||(/_    |_(_)\n", 
        "\n", " __             __    _ \n", "|_    __  _  o (_ |\\||_)\n", 
        "|  |_|| |(_  | __)| ||  \n", "\nVersion: ", package.version("FunciSNP"), 
        "\n", "System: ", if (Sys.info()[["sysname"]] != "Windows") {
            Sys.info()[["sysname"]]
        placeholder<<<<<}
        else {
            paste(Sys.info()[["sysname"]], " :-( ; parallel code not in effect, reverting to 1 core", 
                sep = "")
        placeholder<<<<<}, "\n", "::args used::", "\n", "          verbose:                          ", 
        verbose, "\n", "          cores in use:                     ", 
        par.threads, "\n", "          snp.regions.file:                 ", 
        as.character(snp.regions.file), "\n", "          p-value adjustment by:            ", 
        method.p)
    if (identical(bio.features.loc, NULL)) {
        message("          Bio Features:                     no biofeatures selected")
    placeholder<<<<<}
    else {
        message("          Bio Features:                     ", 
            if (bio.features.TSS) {
                (length(list.files(bio.features.loc, pattern = "*.bed$", 
                  full.names = FALSE))) + 1
            placeholder<<<<<}
            else {
                length(list.files(bio.features.loc, pattern = "*.bed$", 
                  full.names = FALSE))
            placeholder<<<<<}, ": ", gsub(".bed$", ", ", list.files(bio.features.loc, 
                pattern = "*.bed$", full.names = FALSE)), if (bio.features.TSS) 
                "knownGene.TSS.hg19")
    placeholder<<<<<}
    primary.server <- sample(c("ncbi", "ebi"), size = 1)
    snp.region <- ReadRegionsFile(snp.regions.file, search.window)
    message("          Number of TagSNPs Interrogated:   ", nrow(snp.region), 
        " representing ", length(unique(snp.region$snp.name)), 
        " unique tag SNPs")
    options(mc.cores = par.threads)
    populations <- CreatePopulations(primary.server)
    if (identical(bio.features.loc, NULL)) {
        bio.features.file <- NULL
        if (bio.features.TSS) 
            bio.features.file <- list.files(system.file("extdata", 
                package = "FunciSNP"), pattern = "knownGene.TSS.hg19.bed", 
                full.names = TRUE)
    placeholder<<<<<}
    else {
        bio.features.file <- list.files(bio.features.loc, pattern = "*.bed$", 
            full.names = TRUE)
        if (bio.features.TSS) {
            bio.features.file.TSS <- list.files(system.file("extdata", 
                package = "FunciSNP"), pattern = "knownGene.TSS.hg19.bed", 
                full.names = TRUE)
            bio.features.file <- c(bio.features.file, bio.features.file.TSS)
        placeholder<<<<<}
    placeholder<<<<<}
    tag.snp.names <- paste(snp.region$snp.name, ":", snp.region$snp.ethno, 
        sep = "")
    snp.list <- lapply(tag.snp.names, CreateTagSNP)
    names(snp.list) <- tag.snp.names
    snp.list <- mclapply(tag.snp.names, CreateCorrelatedSNPs, 
        snp.list = snp.list, primary.server = primary.server, 
        snp.region = snp.region, bio.features.file = bio.features.file, 
        populations = populations, method.p = method.p, verbose = verbose, 
        reduce.by.features = reduce.by.features, window.size = search.window, 
        par.threads = par.threads)
    names(snp.list) <- tag.snp.names
    snp.list <- TSList(snp.list)
    return(snp.list)
  placeholder<<<<<}
placeholder<<<<<}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 placeholder<<<<<}
\keyword{ ~kwd2 placeholder<<<<<}% __ONLY ONE__ keyword per line
