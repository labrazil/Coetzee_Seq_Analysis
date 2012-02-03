\name{AnnotateSummary}
\alias{AnnotateSummary}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{titlehere
  AnnotateSummary
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
AnnotateSummary(snp.list, verbose = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{snp.list}{
%%     ~~Describe \code{snp.list} here~~
}
  \item{verbose}{
%%     ~~Describe \code{verbose} here~~
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
function (snp.list, verbose = TRUE) 
{
    if (identical(snp.list, NULL)) {
        return(NULL)
    }
    else {
        summary.snp.list <- lapply(snp.list, SNPSummary)
        summary.snp.list <- summary.snp.list[!sapply(summary.snp.list, 
            is.null)]
        summary.snp.list <- IRanges::unlist(GRangesList(summary.snp.list))
        names(summary.snp.list) <- paste(names(summary.snp.list), 
            elementMetadata(summary.snp.list)[, "feature"], sep = ".")
        summary.snp.list <- as.data.frame(summary.snp.list)
        summary.snp.list$width <- NULL
        summary.snp.list$strand <- NULL
        colnames(summary.snp.list) <- c("chromosome", "bio.feature.start", 
            "bio.feature.end", "bio.feature", "corr.snp.id", 
            "corr.snp.position", "tag.snp.id", "tag.snp.position", 
            "D.prime", "R.squared", "p.value", "distance.from.tag", 
            "population.count", "population")
        gr.corr.snp.loc <- GRanges(seqnames = summary.snp.list$chromosome, 
            ranges = IRanges(start = summary.snp.list$corr.snp.position, 
                width = 1), snpid = rownames(summary.snp.list))
        x <- strsplit(as.character(summary.snp.list$chromosome), 
            split = "chr")
        summary.snp.list$chromosome <- sapply(x, "[", 2)
        rd.corr.snp.loc <- RangedData(space = summary.snp.list$chromosome, 
            ranges = IRanges(start = summary.snp.list$corr.snp.position, 
                width = 1), snpid = rownames(summary.snp.list))
        data(TSS.human.GRCh37)
        data(lincRNA.hg19)
        nearest.RNA <- annotatePeakInBatch(myPeakList = rd.corr.snp.loc, 
            AnnotationData = lincRNA, output = "nearestStart")
        summary.snp.list$nearest.lincRNA.ID <- nearest.RNA$feature
        summary.snp.list$nearest.lincRNA.ID <- as.factor(summary.snp.list$nearest.lincRNA.ID)
        summary.snp.list$nearest.lincRNA.distancetoFeature <- nearest.RNA$distancetoFeature
        summary.snp.list$nearest.lincRNA.coverage <- nearest.RNA$insideFeature
        summary.snp.list$nearest.lincRNA.coverage <- as.factor(summary.snp.list$nearest.lincRNA.coverage)
        nearest.TSS <- annotatePeakInBatch(myPeakList = rd.corr.snp.loc, 
            AnnotationData = TSS.human.GRCh37, output = "nearestStart")
        nearest.TSS <- addGeneIDs(nearest.TSS, "org.Hs.eg.db", 
            IDs2Add = c("symbol", "refseq"), silence = TRUE)
        summary.snp.list$nearest.TSS.GeneSymbol <- nearest.TSS$symbol
        summary.snp.list$nearest.TSS.refseq <- nearest.TSS$refseq
        summary.snp.list$nearest.TSS.ensembl <- nearest.TSS$feature
        summary.snp.list$nearest.TSS.ensembl <- as.factor(summary.snp.list$nearest.TSS.ensembl)
        summary.snp.list$nearest.TSS.coverage <- nearest.TSS$insideFeature
        summary.snp.list$nearest.TSS.coverage <- as.factor(summary.snp.list$nearest.TSS.coverage)
        summary.snp.list$nearest.TSS.distancetoFeature <- nearest.TSS$distancetoFeature
        gf.overlaps <- locateVariants(gr.corr.snp.loc, TxDb.Hsapiens.UCSC.hg19.knownGene)
        genomic.feature <- as.character(gf.overlaps$Location)
        queryRow <- (gf.overlaps$queryHits)
        ddd <- (cbind(queryRow, genomic.feature))
        summary.snp.list$Promoter <- "NO"
        summary.snp.list$utr5 <- "NO"
        summary.snp.list$Exon <- "NO"
        summary.snp.list$Intron <- "NO"
        summary.snp.list$utr3 <- "NO"
        summary.snp.list$Intergenic <- "NO"
        promoter.state <- subset(summary.snp.list, (nearest.TSS.distancetoFeature < 
            100) & (nearest.TSS.distancetoFeature > -1000))
        if (dim(promoter.state)[1] > 0) 
            summary.snp.list[rownames(promoter.state), ]$Promoter <- "YES"
        summary.snp.list$Promoter <- as.factor(summary.snp.list$Promoter)
        utr5.rows <- as.numeric(subset(ddd, genomic.feature == 
            "5'UTR")[, 1])
        if (isTRUE(length(unique(utr5.rows)) > 0)) 
            summary.snp.list[utr5.rows, "utr5"] <- "YES"
        summary.snp.list$utr5 <- as.factor(summary.snp.list$utr5)
        exon.rows <- as.numeric(subset(ddd, genomic.feature == 
            "coding")[, 1])
        if (isTRUE(length(unique(exon.rows)) > 0)) 
            summary.snp.list[exon.rows, "Exon"] <- "YES"
        summary.snp.list$Exon <- as.factor(summary.snp.list$Exon)
        intron.rows <- as.numeric(subset(ddd, genomic.feature == 
            "intron")[, 1])
        if (isTRUE(length(unique(intron.rows)) > 0)) 
            summary.snp.list[unique(intron.rows), "Intron"] <- "YES"
        summary.snp.list$Intron <- as.factor(summary.snp.list$Intron)
        utr3.rows <- as.numeric(subset(ddd, genomic.feature == 
            "3'UTR")[, 1])
        if (isTRUE(length(unique(utr3.rows)) > 0)) 
            summary.snp.list[utr3.rows, "utr3"] <- "YES"
        summary.snp.list$utr3 <- as.factor(summary.snp.list$utr3)
        intergenic.rows <- as.numeric(subset(ddd, genomic.feature == 
            "intergenic")[, 1])
        if (isTRUE(length(unique(intergenic.rows)) > 0)) 
            summary.snp.list[intergenic.rows, "Intergenic"] <- "YES"
        summary.snp.list$Intergenic <- as.factor(summary.snp.list$Intergenic)
        promoter.intergenic.rows <- dimnames(subset(summary.snp.list, 
            Intergenic == "YES" & Promoter == "YES"))[[1]]
        if (isTRUE(length(promoter.intergenic.rows) > 0)) 
            summary.snp.list[promoter.intergenic.rows, "Intergenic"] <- "NO"
        return(summary.snp.list)
    }
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
