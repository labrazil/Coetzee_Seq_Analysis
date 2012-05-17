\name{FunciSNP-package}
\alias{FunciSNP-package}
\title{
Functional Identification of SNPs with Phenotype by Coincidence with Chromatin Biofeatures	
}
\description{
The package includes functions to identify and annotate putative functional SNPs using information derived from GWAS, 1000 genomes database, and sequences around peaks.
}
\details{
    \tabular{ll}{
Package: \tab FunciSNP\cr
             Type: \tab Package\cr
             Version: \tab 0.1.7\cr
             Date: \tab 2012-12-12\cr
             License: \tab GPL-3\cr
             LazyLoad: \tab yes\cr
}
}
\references{
	SG. Coetzee, SK. Rhie, BP. Berman, GA. Coetzee and H. Noushmehr, \emph{FunciSNP: An R/Bioconductor Tool Integrating Functional Non-coding Datasets with Genetic Association Studies to Identify Candidate Regulatory SNPs}, Nucleic Acids Research, In press, 2012 (doi:10.1093/nar/GKS542).
}
\author{
Simon Coetzee and Houtan Noushmehr, PhD
    

Maintainer:
    Simon G. Coetzee <scoetzee@gmail.com>
}

\seealso{
	\code{\link{FunciSNPplot}}, \code{\link{FunciSNPAnnotateSummary}}, \code{\link{FunciSNPtable}}, \code{\link{FunciSNPbed}}
}
\examples{
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


}
\keyword{ package }
