\name{FunciSNP}
\alias{FunciSNP}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
FunciSNP(ethno = c("AFR", "ASN", "EUR", "AMR", "ALL"), bio.features.loc, snp.regions.file, R.squared.cutoff = 0, primary.server, make.plots = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{ethno}{
%%     ~~Describe \code{ethno} here~~
}
  \item{bio.features.loc}{
%%     ~~Describe \code{bio.features.loc} here~~
}
  \item{snp.regions.file}{
%%     ~~Describe \code{snp.regions.file} here~~
}
  \item{R.squared.cutoff}{
%%     ~~Describe \code{R.squared.cutoff} here~~
}
  \item{primary.server}{
%%     ~~Describe \code{primary.server} here~~
}
  \item{make.plots}{
%%     ~~Describe \code{make.plots} here~~
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
function (ethno = c("AFR", "ASN", "EUR", "AMR", "ALL"), bio.features.loc, 
    snp.regions.file, R.squared.cutoff = 0, primary.server, make.plots = TRUE) 
{
    cat("\n        ####################################\n        ##                                ##\n        ##      Welcome to FunciSNP       ##\n        ##                                ##\n        ####################################\n        args used:\n\n        snp.regions.file:   ", 
        as.character(snp.regions.file), "\n", "       ethno:            ", 
        ethno, "\n", "       bio.features.loc: ", list.files(bio.features.loc, 
            pattern = "*.bed$", full.names = FALSE), "\n", "       Make Plots?:      ", 
        make.plots, "\n")
    if (primary.server == "ebi" || primary.server == "ncbi") {
        cat("you have selected", primary.server, "as your primary server \n")
    }
    else {
        stop("please select ebi or ncbi as your primary server")
    }
    if (R.squared.cutoff < 1 || R.squared.cutoff > 0) {
        cat("you have selected", R.squared.cutoff, "as your prefered R squared cut off used for summary reports and plots \n")
    }
    else {
        stop("please select an R squared cutoff between 0-1")
    }
    cat("\n ### WARNING:: existing vcf.gz files will be removed from", 
        getwd(), "/funcisnp/", "\n ### and new ones will be placed there, press\n ### <Esc> or <Ctrl-C> to cancel if you do not want this to happen after the\n ### 10 second countdown:", 
        "\n")
    timer <- 10
    pb <- txtProgressBar(min = 0, max = timer, style = 2)
    for (i in 1:timer) {
        Sys.sleep(1)
        setTxtProgressBar(pb, i)
    }
    close(pb)
    try(dir.create(path = "funcisnp_results", showWarnings = FALSE), 
        silent = TRUE)
    try(dir.create(path = "funcisnp_results/tmp", showWarnings = FALSE), 
        silent = TRUE)
    try(dir.create(path = "funcisnp_results/log", showWarnings = FALSE), 
        silent = TRUE)
    try(dir.create(path = "funcisnp_results/tables", showWarnings = FALSE), 
        silent = TRUE)
    snp.region <- ReadRegionsFile(snp.regions.file)
    snp.names <- snp.region$snp.name
    system("find funcisnp_results/ -iname \\*vcf.gz\\* | xargs rm", 
        ignore.stdout = TRUE, ignore.stderr = TRUE)
    onek.genome.server <- ServerCheck(primary.server)
    manifest.read <- read.delim(paste(onek.genome.server, "ftp/release/20110521/phase1_integrated_calls.20101123.ALL.panel", 
        sep = ""), sep = "\t", header = FALSE)
    bio.features.file <<- list.files(bio.features.loc, pattern = "*.bed$", 
        full.names = TRUE)
    bio.features <<- list.files(bio.features.loc, pattern = "*.bed$", 
        full.names = FALSE)
    snp.ld.frame <- NULL
    snp.names.count <- 0
    for (j in as.character(snp.names)) {
        snp.names.count <- snp.names.count + 1
        cat(paste("\n#####\n", "Pulling In Variants File for SNP: ", 
            j, "\n", sep = ""))
        for (i in as.character(ethno)) {
            variants.data <- PullInVariants(i, j, snp.region, 
                primary.server, manifest.read)
            cat(paste("subsetting Variants File for racial/ethnic group: ", 
                i, "\n", sep = ""))
            bio.features.count <- 0
            for (h in bio.features.file) {
                bio.features.count <- bio.features.count + 1
                cat(paste("time: ", format(Sys.time(), "\%H:\%M:\%S"), 
                  "\tRisk SNP: ", j, ", ", i, " (", snp.names.count, 
                  "/", length(snp.names), ")\n>>>>> (", bio.features.count, 
                  "/", length(bio.features), ")", "\t", bio.features[bio.features.count], 
                  "\n", sep = ""))
                snp.geno <- FilterByFeatures(h, j, i, variants.data)
                snp.ld.frame <- LDTesting(snp.geno, j, i, snp.ld.frame, 
                  bio.features.count, snp.region, variants.data)
                try(save(snp.ld.frame, file = "funcisnp_results/snp_table.Rda"), 
                  silent = TRUE)
                try(write.table(snp.ld.frame, file = "funcisnp_results/snp_table.txt", 
                  sep = "\t", quote = FALSE, row.names = FALSE), 
                  silent = TRUE)
            }
        }
        file.remove(paste("funcisnp_results/", j, ".vcf.gz", 
            sep = ""))
        file.remove(paste("funcisnp_results/", j, ".vcf.gz.tbi", 
            sep = ""))
        summary(snp.ld.frame)
    }
    write.table(snp.ld.frame, file = "funcisnp_results/snp_table.txt", 
        sep = "\t", quote = FALSE, row.names = FALSE)
    if (make.plots == TRUE) {
        snp.ld.frame <<- snp.ld.frame
        tmp <- FunciSNPSummary(R.squared.cutoff, snp.ld.frame)
        write.table(tmp, file = "funcisnp_results/tables/summary.txt", 
            sep = "\t", row.names = T, col.names = T, quote = F)
        rm(tmp)
        FunciSNPPlot(R.squared.cutoff, snp.ld.frame)
        FunciSNPHeatmap(R.squared.cutoff, snp.ld.frame)
        cat("##########################################\nFunciSNP is complete!\nSee 'funcisnp_results/tables/' and 'funcisnp_results/plots/' folder for overall results\nSee 'funcisnp_results/snp_table' for complete results.\n##########################################\n")
    }
    else {
        tmp <- FunciSNPSummary(R.squared.cutoff, snp.ld.frame)
        write.table(tmp, file = "funcisnp_results/tables/summary.txt", 
            sep = "\t", row.names = T, col.names = T, quote = F)
        rm(tmp)
        cat("##########################################\nFunciSNP is complete!\nSee 'funcisnp_results/tables/' folder for overall results\nSee 'funcisnp_results/snp_table' for complete results.\n##########################################")
    }
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
