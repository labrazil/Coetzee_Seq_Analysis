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
### code chunk number 2: FunciSNP_vignette.Rnw:79-105
###################################################
## Full path to the example GWAS SNP regions file for Glioblastoma 
#  (collected from SNPedia on Jan 2012)
glioma.snp <- file.path(system.file('extdata', package='FunciSNP'), 
dir(system.file('extdata',package='FunciSNP'), pattern='.snp$'));
#glioma.snp;
## Full path to the example biological features BED files 
#  derived from the ENCODE project for Glioblastoma U-87 cell lines.
glioma.bio <- system.file('extdata',package='FunciSNP');
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

## The glioma example set above is called directly by  
data(glioma);
glioma;
summary(glioma);
class(glioma);


###################################################
### code chunk number 3: FunciSNP_vignette.Rnw:114-124
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
### code chunk number 4: FunciSNP_vignette.Rnw:132-134
###################################################
FunciSNPtable(glioma.anno, rsq=0.5);
FunciSNPtable(glioma.anno, rsq=0.5, geneSum=TRUE);


###################################################
### code chunk number 5: FunciSNP_vignette.Rnw:141-143
###################################################
FunciSNPsummaryOverlaps(glioma.anno)
FunciSNPsummaryOverlaps(glioma.anno, rsq=0.5)


###################################################
### code chunk number 6: FunciSNP_vignette.Rnw:152-158
###################################################
rs6010620 <- FunciSNPidsFromSummary(glioma.anno, tagsnpid="rs6010620", 
        num.features=2, rsq=0.5)
summary(rs6010620);
dim(rs6010620);
class(rs6010620);
## See FunciSNPbed to visualize this data in a genome browser.


###################################################
### code chunk number 7: FunciSNP_vignette.Rnw:172-175
###################################################
pdf("glioma_dist.pdf")
FunciSNPplot(glioma.anno)
dev.off()


###################################################
### code chunk number 8: FunciSNP_vignette.Rnw:188-190
###################################################
FunciSNPplot(glioma.anno, splitbysnp=TRUE)
ggsave("glioma_dist_bysnp.pdf")


###################################################
### code chunk number 9: FunciSNP_vignette.Rnw:204-207
###################################################
pdf("glioma_genomic_sum_rcut.pdf")
FunciSNPplot(glioma.anno, rsq=0.5, genomicSum=TRUE, save=FALSE)
dev.off()


###################################################
### code chunk number 10: FunciSNP_vignette.Rnw:226-228
###################################################
## Following will output a series of plots for each biofeature at rsq=0.5
FunciSNPplot(glioma.anno, tagSummary=TRUE, rsq=0.5)


###################################################
### code chunk number 11: FunciSNP_vignette.Rnw:260-263
###################################################
## will output to current working directory.
FunciSNPbed(glioma.anno, rsq=0.5);
# FunciSNPbed(rs6010620, rsq=0.5);


###################################################
### code chunk number 12: sessionInfo
###################################################
sessionInfo()


