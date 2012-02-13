### R code from vignette source 'FunciSNP_vignette.Rnw'

###################################################
### code chunk number 1: packages
###################################################
#When package is offically posted in Bioconductor, uncomment next 2 lines.
#source("http://bioconductor.org/biocLite.R")
#biocLite("FunciSNP");
## Following two packages and options() are not required to run 'FunciSNP' but 
#will enhance the analysis experience.
#library(setwidth); ## Automatically set the value of options("width") when the 
#terminal emulator is resized
#library(colorout); ## colorize R output on terminal emulators
options(width=80);
##FunciSNP library and other related libraries needed.
library("org.Hs.eg.db");
library("gplots");
library("gtools");
library("ggplot2");
library("matlab");

library(FunciSNP);
package.version("FunciSNP");


###################################################
### code chunk number 2: FunciSNP_vignette.Rnw:96-102
###################################################
## Full path to the example GWAS SNP regions file for Glioblastoma 
#  (collected from SNPedia on Jan 2012)
glioma.snp <- file.path(system.file('extdata', package='FunciSNP'), 
dir(system.file('extdata',package='FunciSNP'), pattern='.snp$'));
gsnp <- read.delim(file=glioma.snp,sep=" ",header=FALSE);
gsnp;


###################################################
### code chunk number 3: FunciSNP_vignette.Rnw:108-117
###################################################
#glioma.snp;
## Full path to the example biological features BED files 
#  derived from the ENCODE project for Glioblastoma U-87 cell lines.
glioma.bio <- system.file('extdata',package='FunciSNP');
list.files(glioma.bio, pattern='.bed$');
nrsf.filename <- list.files(glioma.bio, pattern='.bed$')[2];
Nrsf <- read.delim(file=paste(glioma.bio, nrsf.filename,sep="/"), sep="\t",
        header=FALSE);
head(Nrsf);


###################################################
### code chunk number 4: FunciSNP_vignette.Rnw:120-132
###################################################
#glioma.bio;
## FunciSNP analysis, extracts correlated SNPs from the 
#  1000 genomes db ("ncbi" or "ebi") and finds overlaps between 
#  correlated SNP and biological features and then 
#  calculates LD (Rsquare, Dprime, distance, p-value).
## Depending on number of CPUs and internet connection, this step may take 
# some time. Please consider using a unix machine to access multiple cores.
# glioma <- FuncySNP(snp.regions.file=glioma.snp, 
#           bio.features.loc = glioma.bio, 
#           bio.features.TSS=FALSE);
# glioma;
# summary(glioma);


###################################################
### code chunk number 5: FunciSNP_vignette.Rnw:137-141
###################################################
data(glioma);
glioma;
summary(glioma);
class(glioma);


###################################################
### code chunk number 6: FunciSNP_vignette.Rnw:150-160
###################################################
glioma.anno <- FunciSNPAnnotateSummary(glioma);
class(glioma.anno);
gl.anno <- glioma.anno;
## remove rownames for this example section.
rownames(gl.anno) <- c(1:length(rownames(gl.anno)))
dim(gl.anno);
head(gl.anno); ## 
names(gl.anno);
summary(gl.anno[,c(1:18,20:28)]);
rm(gl.anno);


###################################################
### code chunk number 7: FunciSNP_vignette.Rnw:167-168
###################################################
FunciSNPtable(glioma.anno, rsq=0.5);


###################################################
### code chunk number 8: FunciSNP_vignette.Rnw:171-172
###################################################
FunciSNPtable(glioma.anno, rsq=0.5, geneSum=TRUE);


###################################################
### code chunk number 9: FunciSNP_vignette.Rnw:180-181
###################################################
FunciSNPsummaryOverlaps(glioma.anno)


###################################################
### code chunk number 10: FunciSNP_vignette.Rnw:185-186
###################################################
FunciSNPsummaryOverlaps(glioma.anno, rsq=0.5)


###################################################
### code chunk number 11: FunciSNP_vignette.Rnw:195-201
###################################################
rs6010620 <- FunciSNPidsFromSummary(glioma.anno, tagsnpid="rs6010620", 
        num.features=2, rsq=0.5)
summary(rs6010620);
dim(rs6010620);
class(rs6010620);
## See FunciSNPbed to visualize this data in a genome browser.


###################################################
### code chunk number 12: FunciSNP_vignette.Rnw:215-218
###################################################
pdf("glioma_dist.pdf")
FunciSNPplot(glioma.anno)
dev.off()


###################################################
### code chunk number 13: FunciSNP_vignette.Rnw:232-234
###################################################
FunciSNPplot(glioma.anno, splitbysnp=TRUE)
ggsave("glioma_dist_bysnp.pdf")


###################################################
### code chunk number 14: FunciSNP_vignette.Rnw:248-251
###################################################
pdf("glioma_genomic_sum_rcut.pdf")
FunciSNPplot(glioma.anno, rsq=0.5, genomicSum=TRUE, save=FALSE)
dev.off()


###################################################
### code chunk number 15: FunciSNP_vignette.Rnw:274-276
###################################################
## Following will output a series of plots for each biofeature at rsq=0.5
FunciSNPplot(glioma.anno, tagSummary=TRUE, rsq=0.5)


###################################################
### code chunk number 16: FunciSNP_vignette.Rnw:314-317
###################################################
## will output to current working directory.
FunciSNPbed(glioma.anno, rsq=0.5);
# FunciSNPbed(rs6010620, rsq=0.5);


###################################################
### code chunk number 17: sessionInfo
###################################################
sessionInfo()


