\name{FunciSNPplotplaceholder<<<<<}
\alias{FunciSNPplotplaceholder<<<<<}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
placeholder<<<<<}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
placeholder<<<<<}
\usage{
FunciSNPplot(dat, rsq = 0, split = FALSE, splitbysnp = FALSE, tagSummary = FALSE, heatmap = FALSE, genomicSum = FALSE, save = FALSE)
placeholder<<<<<}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{datplaceholder<<<<<}{
%%     ~~Describe \code{datplaceholder<<<<<} here~~
placeholder<<<<<}
  \item{rsqplaceholder<<<<<}{
%%     ~~Describe \code{rsqplaceholder<<<<<} here~~
placeholder<<<<<}
  \item{splitplaceholder<<<<<}{
%%     ~~Describe \code{splitplaceholder<<<<<} here~~
placeholder<<<<<}
  \item{splitbysnpplaceholder<<<<<}{
%%     ~~Describe \code{splitbysnpplaceholder<<<<<} here~~
placeholder<<<<<}
  \item{tagSummaryplaceholder<<<<<}{
%%     ~~Describe \code{tagSummaryplaceholder<<<<<} here~~
placeholder<<<<<}
  \item{heatmapplaceholder<<<<<}{
%%     ~~Describe \code{heatmapplaceholder<<<<<} here~~
placeholder<<<<<}
  \item{genomicSumplaceholder<<<<<}{
%%     ~~Describe \code{genomicSumplaceholder<<<<<} here~~
placeholder<<<<<}
  \item{saveplaceholder<<<<<}{
%%     ~~Describe \code{saveplaceholder<<<<<} here~~
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
function (dat, rsq = 0, split = FALSE, splitbysnp = FALSE, tagSummary = FALSE, 
    heatmap = FALSE, genomicSum = FALSE, save = FALSE) 
{
    if (sum(c(split, tagSummary, heatmap, genomicSum)) == 0) {
        split = TRUE
    placeholder<<<<<}
    require(ggplot2)
    if (split) {
        if (splitbysnp == TRUE) {
            FunciSNP:::theme_white()
            ggplot(dat, aes(x = R.squared)) + geom_histogram(binwidth = 0.05) + 
                geom_vline(xintercept = 0.5, linetype = 2) + 
                scale_x_continuous("Correlated SNPs R squared to Tag SNP (0-1)") + 
                scale_y_continuous("Total # of Correlated SNPs associated with tagSNP") + 
                opts(legend.position = "none", axis.text.y = theme_text(), 
                  axis.text.x = theme_text(angle = 90), title = "Distribution of correlated SNPs for each tagSNP\n", 
                  "at Rsquared values") + facet_wrap(chromosome ~ 
                tag.snp.id)
        placeholder<<<<<}
        else {
            tt <- count(df = dat, vars = "R.squared")
            tt <- na.omit(tt)
            ht <- range(tt[, "freq"])[2] * 1.2
            hh <- dat[, c("corr.snp.id", "R.squared")]
            hh <- na.omit(hh)
            hh.c <- count(round(hh$R.squared, digits = 1))
            dimnames(hh.c)[[1]] <- hh.c[, 1]
            k <- c(hh.c["0", 2], hh.c["0.1", 2], hh.c["0.2", 
                2], hh.c["0.3", 2], hh.c["0.4", 2], hh.c["0.5", 
                2], hh.c["0.6", 2], hh.c["0.7", 2], hh.c["0.8", 
                2], hh.c["0.9", 2], hh.c["1", 2])
            k[is.na(k)] <- 0
            plot(tt, xlim = c(0, 1), ylim = c(0, ht), pch = "*", 
                main = paste("Distribution of Correlated SNPs by Rsq values\n", 
                  "Total # of corr. SNPs: ", dim(dat)[1], "\n(with an Rsq value: ", 
                  sum(hh.c$freq), "; unique corr. SNPs: ", length(unique(hh$corr.snp.id)), 
                  ")", sep = ""), xlab = "R square values (0-1)", 
                ylab = "Number of correlated SNPs")
            abline(v = 0.1, lty = 2, col = "red")
            abline(v = 0.2, lty = 2, col = "red")
            abline(v = 0.3, lty = 2, col = "red")
            abline(v = 0.4, lty = 2, col = "red")
            abline(v = 0.5, col = "black")
            abline(v = 0.6, lty = 2, col = "green")
            abline(v = 0.7, lty = 2, col = "green")
            abline(v = 0.8, lty = 2, col = "green")
            abline(v = 0.9, lty = 2, col = "green")
            abline(h = ht * 0.9, col = "black", lty = 2)
            text(0.05, ht * 0.95, as.character(k[1] + k[2]))
            text(0.15, ht * 0.95, as.character(k[3]))
            text(0.25, ht * 0.95, as.character(k[4]))
            text(0.35, ht * 0.95, as.character(k[5]))
            text(0.45, ht * 0.95, as.character(k[6]))
            text(0.55, ht * 0.95, as.character(k[7]))
            text(0.65, ht * 0.95, as.character(k[8]))
            text(0.75, ht * 0.95, as.character(k[9]))
            text(0.85, ht * 0.95, as.character(k[10]))
            text(0.95, ht * 0.95, as.character(k[11]))
        placeholder<<<<<}
    placeholder<<<<<}
    if (tagSummary) {
        try(dir.create(path = paste("FunciSNP.", package.version("FunciSNP"), 
            "/plots", sep = ""), showWarnings = FALSE, recursive = TRUE), 
            silent = TRUE)
        require("ggplot2")
        theme_white()
        all.s <- try(subset(dat, R.squared >= rsq), silent = TRUE)
        all.ss <- try(subset(dat, R.squared < rsq), silent = TRUE)
        try(all.s$r.2 <- c("Yes"), silent = TRUE)
        try(all.ss$r.2 <- c("No"), silent = TRUE)
        if (nrow(all.s) > 0 && nrow(all.ss) > 0) {
            all <- try(rbind(all.s, all.ss), silent = TRUE)
        placeholder<<<<<}
        else {
            if (nrow(all.s) > 0 && nrow(all.ss) <= 0) {
                all <- all.s
            placeholder<<<<<}
            if (nrow(all.s) <= 0 && nrow(all.ss) > 0) {
                all <- all.ss
            placeholder<<<<<}
            if (nrow(all.s) <= 1 && nrow(all.ss) <= 0) {
                return()
            placeholder<<<<<}
        placeholder<<<<<}
        for (i in 1:length(summary(as.factor(all[, "bio.feature"])))) {
            bio <- names(summary(as.factor(all[, "bio.feature"])))
            tmp <- subset(all, bio.feature == bio[i])
            ggplot(tmp, aes(x = R.squared, fill = factor(r.2))) + 
                geom_histogram(binwidth = 0.05) + geom_vline(xintercept = rsq, 
                linetype = 2) + scale_x_continuous("Rsquare Values (0-1)", 
                limits = c(0, 1)) + scale_y_continuous("Total # of Surrogate SNPs associated with riskSNP") + 
                scale_fill_manual(values = c(Yes = "Red", No = "Black")) + 
                opts(legend.position = "none", axis.text.y = theme_text(), 
                  axis.text.x = theme_text(angle = 90), title = paste("riskSNP\nOverlapping: ", 
                    bio[i], sep = "")) + facet_wrap(~tag.snp.id)
            ggsave(file = paste("FunciSNP.", package.version("FunciSNP"), 
                "/plots/", bio[i], "_R2summary_riskSNP.pdf", 
                sep = ""))
            ggplot(tmp, aes(x = R.squared, y = distance.from.tag, 
                colour = r.2, size = factor(r.2))) + geom_point() + 
                geom_vline(xintercept = rsq, linetype = 2) + 
                scale_x_continuous("Rsquare Values (0-1)", limits = c(0, 
                  1)) + scale_y_continuous("Distance to Surrogate SNPs associated with riskSNP (bp)", 
                formatter = "comma") + scale_colour_manual(values = c(Yes = "Red", 
                No = "Black")) + scale_size_manual(values = c(Yes = 2, 
                No = 1)) + opts(legend.position = "none", axis.text.y = theme_text(), 
                axis.text.x = theme_text(angle = 90), title = paste("Distance between riskSNP\n", 
                  "and Surrogate SNP\nOverlapping: ", bio[i], 
                  sep = "")) + facet_wrap(~tag.snp.id)
            ggsave(file = paste("FunciSNP.", package.version("FunciSNP"), 
                "/plots/", bio[i], "_R2vsDist_riskSNP.pdf", sep = ""))
            cat("Finished plotting ", i, "/", length(bio), "\n")
        placeholder<<<<<}
        message("\n\nSee ", paste("FunciSNP.", package.version("FunciSNP"), 
            "/plots/", sep = ""), " folder in ", getwd(), " for all plots.\n\n")
    placeholder<<<<<}
    if (heatmap) {
        require("gplots")
        require("matlab")
        try(dir.create(path = paste("FunciSNP.", package.version("FunciSNP"), 
            "/plots", sep = ""), showWarnings = FALSE, recursive = TRUE), 
            silent = TRUE)
        all.s <- (table(subset(dat, R.squared >= rsq)[, "bio.feature"], 
            subset(dat, R.squared >= rsq)[, "tag.snp.id"]))
        all.s <- as.matrix(all.s)
        if (save) 
            png(filename = paste("FunciSNP.", package.version("FunciSNP"), 
                "/plots/FunciSNP_heatmap.png", sep = ""), bg = "white", 
                width = 3000, height = 3000)
        heatmap.2(all.s, na.rm = TRUE, scale = "none", col = jet.colors(max(all.s, 
            na.rm = T)), key = T, symkey = FALSE, density.info = "none", 
            trace = "none", Rowv = F, Colv = T, cexRow = 1, cexCol = 1, 
            keysize = 0.5, dendrogram = c("none"), main = paste("Rsquare values >= ", 
                rsq))
        if (save) 
            dev.off()
        if (save) 
            message("\nSee ", paste("FunciSNP.", package.version("FunciSNP"), 
                "/plots/", sep = ""), "folder in ", getwd(), 
                " for heatmap.\n\n")
    placeholder<<<<<}
    if (genomicSum) {
        if (rsq == 0) {
            dat.m <- melt(dat[, c(23:28)], measure.vars = c("Promoter", 
                "utr5", "Exon", "Intron", "utr3", "Intergenic"))
            t <- subset(dat.m, value == "NO")
            t$value <- "2.NO"
            tt <- subset(dat.m, value != "NO")
            tt$value <- "1.YES"
            dat.m <- rbind(t, tt)
            require(ggplot2)
            FunciSNP:::theme_white()
            ggplot(dat.m, aes(variable, fill = factor(value))) + 
                geom_bar() + opts(axis.text.y = theme_text(), 
                axis.text.x = theme_text(angle = 90), title = "Correlated SNPs distribution across genomic features") + 
                scale_fill_manual(values = c(`1.YES` = "Red", 
                  `2.NO` = "Black"), "Overlap") + scale_y_continuous("Total count of correlated SNPs")
        placeholder<<<<<}
        else {
            dat$r2 <- paste("R squared < ", rsq, sep = "")
            t <- subset(dat, R.squared >= rsq)
            t$r2 <- paste("R squared >= ", rsq, sep = "")
            dat <- rbind(t, subset(dat, R.squared < rsq))
            dat.m <- melt(dat[, c(23:29)], measure.vars = c("Promoter", 
                "utr5", "Exon", "Intron", "utr3", "Intergenic"))
            t <- subset(dat.m, value == "NO")
            t$value <- "2.NO"
            tt <- subset(dat.m, value != "NO")
            tt$value <- "1.YES"
            dat.m <- rbind(t, tt)
            require(ggplot2)
            FunciSNP:::theme_white()
            ggplot(dat.m, aes(variable, fill = factor(value))) + 
                geom_bar(position = "fill") + opts(axis.text.y = theme_text(), 
                axis.text.x = theme_text(angle = 90), title = paste("Distribution of corr. SNPs across genomic features\n", 
                  " at Rsquared cutoff of", rsq, sep = " ")) + 
                scale_fill_manual(values = c(`1.YES` = "Red", 
                  `2.NO` = "Black"), "Overlap") + scale_y_continuous("Percent of Total correlated SNPs at Rsquared cutoff") + 
                facet_wrap(~r2)
        placeholder<<<<<}
    placeholder<<<<<}
  placeholder<<<<<}
placeholder<<<<<}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 placeholder<<<<<}
\keyword{ ~kwd2 placeholder<<<<<}% __ONLY ONE__ keyword per line
