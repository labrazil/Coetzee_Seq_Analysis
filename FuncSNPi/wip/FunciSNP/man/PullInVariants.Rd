\name{PullInVariants}
\alias{PullInVariants}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
PullInVariants(ethno.chosen, snp.names.chosen, snp.regions, primary.server, manifest)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{ethno.chosen}{
%%     ~~Describe \code{ethno.chosen} here~~
}
  \item{snp.names.chosen}{
%%     ~~Describe \code{snp.names.chosen} here~~
}
  \item{snp.regions}{
%%     ~~Describe \code{snp.regions} here~~
}
  \item{primary.server}{
%%     ~~Describe \code{primary.server} here~~
}
  \item{manifest}{
%%     ~~Describe \code{manifest} here~~
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
function (ethno.chosen, snp.names.chosen, snp.regions, primary.server, 
    manifest) 
{
    ifelse(ethno.chosen != "ALL", ethno.sample.set <- as.character(subset(manifest, 
        V3 == ethno.chosen)[, 1]), ethno.sample.set <- as.character(manifest[, 
        1]))
    intermediate.vcf <- paste("funcisnp_results/", snp.names.chosen, 
        ".vcf", sep = "")
    variants.file <- paste(intermediate.vcf, ".gz", sep = "")
    if (file.exists(variants.file) == FALSE) {
        onek.genome.server <- ServerCheck(primary.server)
        system(paste("tabix -hf ", onek.genome.server, "/ftp/release/20101123/interim_phase1_release/ALL.chr", 
            snp.regions$snp.chromosome[snp.regions$snp.name == 
                snp.names.chosen][1], ".phase1.projectConsensus.genotypes.vcf.gz", 
            " ", snp.regions$snp.chromosome[snp.regions$snp.name == 
                snp.names.chosen][1], ":", snp.regions$snp.region.start[snp.regions$snp.name == 
                snp.names.chosen][1], "-", snp.regions$snp.region.end[snp.regions$snp.name == 
                snp.names.chosen][1], " > ", intermediate.vcf, 
            " && ", "bgzip -f ", intermediate.vcf, " && ", "tabix -hf ", 
            intermediate.vcf, ".gz", sep = ""))
        subset.variants.file <<- unpackVcf(scanVcf(variants.file))[[1]]
    }
    ifelse(ethno.chosen != "ALL", genotype.data <- subset.variants.file$GENO$GT[, 
        c(ethno.sample.set)], genotype.data <- subset.variants.file$GENO$GT)
    variants.data <- data.frame(subset.variants.file$CHROM, subset.variants.file$POS, 
        subset.variants.file$ID, subset.variants.file$REF, subset.variants.file$ALT, 
        genotype.data, stringsAsFactors = FALSE)
    core.names <- c("CHROM", "POS", "ID", "REF", "ALT")
    colnames(variants.data)[1:5] <- core.names
    variants.data.length <- length(variants.data)
    variants.data[, 6:variants.data.length] <- ifelse(variants.data[, 
        6:variants.data.length] == "0|0", paste(variants.data$REF, 
        "/", variants.data$REF, sep = ""), ifelse(variants.data[, 
        6:variants.data.length] == "0|1", paste(variants.data$REF, 
        "/", variants.data$ALT, sep = ""), ifelse(variants.data[, 
        6:variants.data.length] == "1|0", paste(variants.data$ALT, 
        "/", variants.data$REF, sep = ""), ifelse(variants.data[, 
        6:variants.data.length] == "1|1", paste(variants.data$ALT, 
        "/", variants.data$ALT, sep = ""), NA))))
    variants.data$ID <- ifelse(variants.data$ID == ".", variants.data$POS, 
        variants.data$ID)
    variants.data
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
