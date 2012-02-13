pkgname <- "FunciSNP"
source(file.path(R.home("share"), "R", "examples-header.R"))
options(warn = 1)
library('FunciSNP')

assign(".oldSearch", search(), pos = 'CheckExEnv')
cleanEx()
nameEx("CorrGeno-class")
### * CorrGeno-class

flush(stderr()); flush(stdout())

### Name: CorrGeno-class
### Title: Class '"CorrGeno"'
### Aliases: CorrGeno-class
### Keywords: classes

### ** Examples

showClass("CorrGeno")



cleanEx()
nameEx("CorrelatedSNPs-class")
### * CorrelatedSNPs-class

flush(stderr()); flush(stdout())

### Name: CorrelatedSNPs-class
### Title: Class '"CorrelatedSNPs"'
### Aliases: CorrelatedSNPs-class AFR.D.prime<-,CorrelatedSNPs-method
###   AFR.D.prime,CorrelatedSNPs-method AFR.p.value<-,CorrelatedSNPs-method
###   AFR.p.value,CorrelatedSNPs-method
###   AFR.R.squared<-,CorrelatedSNPs-method
###   AFR.R.squared,CorrelatedSNPs-method
###   ALL.D.prime<-,CorrelatedSNPs-method ALL.D.prime,CorrelatedSNPs-method
###   ALL.p.value<-,CorrelatedSNPs-method ALL.p.value,CorrelatedSNPs-method
###   ALL.R.squared<-,CorrelatedSNPs-method
###   ALL.R.squared,CorrelatedSNPs-method
###   alt.allele<-,CorrelatedSNPs-method alt.allele,CorrelatedSNPs-method
###   AMR.D.prime<-,CorrelatedSNPs-method AMR.D.prime,CorrelatedSNPs-method
###   AMR.p.value<-,CorrelatedSNPs-method AMR.p.value,CorrelatedSNPs-method
###   AMR.R.squared<-,CorrelatedSNPs-method
###   AMR.R.squared,CorrelatedSNPs-method
###   ASN.D.prime<-,CorrelatedSNPs-method ASN.D.prime,CorrelatedSNPs-method
###   ASN.p.value<-,CorrelatedSNPs-method ASN.p.value,CorrelatedSNPs-method
###   ASN.R.squared<-,CorrelatedSNPs-method
###   ASN.R.squared,CorrelatedSNPs-method chr<-,CorrelatedSNPs-method
###   chr,CorrelatedSNPs-method EUR.D.prime<-,CorrelatedSNPs-method
###   EUR.D.prime,CorrelatedSNPs-method EUR.p.value<-,CorrelatedSNPs-method
###   EUR.p.value,CorrelatedSNPs-method
###   EUR.R.squared<-,CorrelatedSNPs-method
###   EUR.R.squared,CorrelatedSNPs-method
###   overlapping.features<-,CorrelatedSNPs-method
###   overlapping.features,CorrelatedSNPs-method
###   pop.genotype<-,CorrelatedSNPs-method
###   pop.genotype,CorrelatedSNPs-method position<-,CorrelatedSNPs-method
###   position,CorrelatedSNPs-method ref.allele<-,CorrelatedSNPs-method
###   ref.allele,CorrelatedSNPs-method snpid<-,CorrelatedSNPs-method
###   snpid,CorrelatedSNPs-method
### Keywords: classes

### ** Examples

showClass("CorrelatedSNPs")



cleanEx()
nameEx("FunciSNP-package")
### * FunciSNP-package

flush(stderr()); flush(stdout())

### Name: FunciSNP-package
### Title: Functional Identification of SNPs with Phenotype by Coincidence
###   with Chromatin Biofeatures
### Aliases: FunciSNP-package
### Keywords: package

### ** Examples

##
## Glioblastoma analysis using FunciSNP
##
## Full path to the example regions file for Glioblastoma 
#  (collected from SNPedia)
glioma.snp <- file.path(system.file('extdata',
  package='FunciSNP'),
  dir(system.file('extdata',package='FunciSNP'), 
  pattern='.snp$'));
 
## Full path to the example biological features BED files 
#  derived from the ENCODE project for Glioblastoma U-87 
#  cell lines.
glioma.bio <- system.file('extdata',package='FunciSNP');

## FunciSNP analysis, extracts correlated SNPs from the 
#  1000 genomes db ("ncbi") and finds overlaps between 
#  correlated SNP and biological features and then 
#  calculates LD (Rsquare, Dprime, distance, p-value).
# Do not run. Can take more than 5 min depending on internet connection and number of CPUs.
#glioma <- FunciSNP(snp.regions.file=glioma.snp, 
#  bio.features.loc = glioma.bio, bio.features.TSS=FALSE);

##
data(glioma);
class(glioma);
glioma;
summary(glioma);





cleanEx()
nameEx("FunciSNPAnnotateSummary")
### * FunciSNPAnnotateSummary

flush(stderr()); flush(stdout())

### Name: FunciSNPAnnotateSummary
### Title: Genomic Annotation of Func-y-SNPs.
### Aliases: FunciSNPAnnotateSummary
### Keywords: annotation

### ** Examples

data(glioma);
gl <- FunciSNPAnnotateSummary(glioma);
dim(gl)
head(gl)
names(gl)



cleanEx()
nameEx("FunciSNPbed")
### * FunciSNPbed

flush(stderr()); flush(stdout())

### Name: FunciSNPbed
### Title: Creates a BED file to view Func-y-SNPs in your favorite genome
###   browser
### Aliases: FunciSNPbed
### Keywords: ~kwd1 ~kwd2

### ** Examples

##
data(glioma);
glioma.anno <- FunciSNPAnnotateSummary(glioma);
FunciSNPbed(glioma.anno, rsq=0.9);
####
#Bed file "FunciSNP_results_rsq.0.9.bed" created successfully.
#(See folder: "/home/houtan/Downloads/")
#Total corSNP (RED):  15 
#Total tagSNP (BLK):  1 

#To view results, submit bed file as a
#  custom track in UCSC Genome Browser (genome.ucsc.edu), 

#Now have fun with your new Func-y SNPs!!
####






cleanEx()
nameEx("FunciSNPidsFromSummary")
### * FunciSNPidsFromSummary

flush(stderr()); flush(stdout())

### Name: FunciSNPidsFromSummary
### Title: coming soon.
### Aliases: FunciSNPidsFromSummary
### Keywords: ~kwd1 ~kwd2

### ** Examples

## coming soon



cleanEx()
nameEx("FunciSNPplot")
### * FunciSNPplot

flush(stderr()); flush(stdout())

### Name: FunciSNPplot
### Title: FunciSNPplot to visualize Func-y-SNP summary.
### Aliases: FunciSNPplot
### Keywords: ~kwd1 ~kwd2

### ** Examples

data(glioma)
gl <- FunciSNPAnnotateSummary(glioma)
FunciSNPplot(gl)
FunciSNPplot(gl, rsq=0, genomicSum=TRUE, save=FALSE)
FunciSNPplot(gl, rsq=0.5, genomicSum=TRUE, save=FALSE)
# DO NOT RUN
#FunciSNPplot(gl, tagSummary=TRUE, rsq=0.5)
#



cleanEx()
nameEx("FunciSNPsummaryOverlaps")
### * FunciSNPsummaryOverlaps

flush(stderr()); flush(stdout())

### Name: FunciSNPsummaryOverlaps
### Title: coming
### Aliases: FunciSNPsummaryOverlaps
### Keywords: ~kwd1 ~kwd2

### ** Examples

##coming soon.



cleanEx()
nameEx("FunciSNPtable")
### * FunciSNPtable

flush(stderr()); flush(stdout())

### Name: FunciSNPtable
### Title: Will output a summary report from FunciSNP at specified Rsquare
###   cut-offs.
### Aliases: FunciSNPtable
### Keywords: ~kwd1 ~kwd2

### ** Examples

data(glioma);
gl <- FunciSNPAnnotateSummary(glioma);
FunciSNPtable(gl, rsq=0.5);
FunciSNPtable(gl, rsq=0.5, geneSum=TRUE);



cleanEx()
nameEx("FuncySNP")
### * FuncySNP

flush(stderr()); flush(stdout())

### Name: FuncySNP
### Title: Functional Identification of SNPs with Phenotype by Coincidence
###   with Chromatin Biofeatures
### Aliases: FuncySNP glioma lincRNA.hg19
### Keywords: #kwd1 #kwd2

### ** Examples

##
## Glioblastoma analysis using FunciSNP
##
## Full path to the example regions file for Glioblastoma 
#  (collected from SNPedia)
glioma.snp <- file.path(system.file('extdata',
  package='FunciSNP'),
  dir(system.file('extdata',package='FunciSNP'), 
  pattern='.snp$'));
 
## Full path to the example biological features BED files 
#  derived from the ENCODE project for Glioblastoma U-87 
#  cell lines.
glioma.bio <- system.file('extdata',package='FunciSNP');

## FunciSNP analysis, extracts correlated SNPs from the 
#  1000 genomes db ("ncbi") and finds overlaps between 
#  correlated SNP and biological features and then 
#  calculates LD (Rsquare, Dprime, distance, p-value).
# Do not run. Can take more than 5 min depending on internet connection and number of CPUs.
#glioma <- FuncySNP(snp.regions.file=glioma.snp, 
#  bio.features.loc = glioma.bio, bio.features.TSS=FALSE);

##
data(glioma);
class(glioma);
glioma;
summary(glioma);





cleanEx()
nameEx("TSList-class")
### * TSList-class

flush(stderr()); flush(stdout())

### Name: TSList-class
### Title: Class '"TSList"'
### Aliases: TSList-class show,TSList-method summary,TSList-method
### Keywords: classes

### ** Examples

showClass("TSList")



cleanEx()
nameEx("TagSNP-class")
### * TagSNP-class

flush(stderr()); flush(stdout())

### Name: TagSNP-class
### Title: Class '"TagSNP"'
### Aliases: TagSNP-class AFR.overlapping.snps.geno,TagSNP-method
###   ALL.overlapping.snps.geno,TagSNP-method alt.allele<-,TagSNP-method
###   alt.allele,TagSNP-method AMR.overlapping.snps.geno,TagSNP-method
###   ASN.overlapping.snps.geno,TagSNP-method chr<-,TagSNP-method
###   chr,TagSNP-method correlated.snps<-,TagSNP-method
###   correlated.snps,TagSNP-method D.prime.corrsnps<-,TagSNP-method
###   D.prime.corrsnps,TagSNP-method
###   EUR.overlapping.snps.geno,TagSNP-method genotype<-,TagSNP-method
###   genotype,TagSNP-method overlapping.features<-,TagSNP-method
###   overlapping.features,TagSNP-method population<-,TagSNP-method
###   population,TagSNP-method position<-,TagSNP-method
###   position,TagSNP-method ref.allele<-,TagSNP-method
###   ref.allele,TagSNP-method R.squared.corrsnps<-,TagSNP-method
###   R.squared.corrsnps,TagSNP-method show,TagSNP-method
###   snpid<-,TagSNP-method snpid,TagSNP-method
### Keywords: classes

### ** Examples

showClass("TagSNP")



### * <FOOTER>
###
cat("Time elapsed: ", proc.time() - get("ptime", pos = 'CheckExEnv'),"\n")
grDevices::dev.off()
###
### Local variables: ***
### mode: outline-minor ***
### outline-regexp: "\\(> \\)?### [*]+" ***
### End: ***
quit('no')
