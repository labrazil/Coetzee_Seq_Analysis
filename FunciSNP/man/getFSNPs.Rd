\name{getFSNPs}
\alias{getFSNPs}
\alias{glioma}
\alias{lincRNA.hg19}
\title{
Functional Identification of SNPs with Phenotype by Coincidence with Chromatin Biofeatures	
}
\description{
Given a set of known tag-SNPs associated with a particular phenotype (e.g. disease, trait), and a set of available biological features (e.g. peaks derived from ChIP-seq experiments for phenotype), returns correlated SNPs (from the 1000 genomes db) which are in linkage disequilibrium (LD) to a known disease associated tag-SNP and overlaps chromatin biological features.  These identified correlated SNPs are characterized as putative functional SNPs for a particular trait. 
}
\usage{
getFSNPs(snp.regions.file, bio.features.loc = NULL,
                     built.in.biofeatures = TRUE,
                     par.threads=detectCores()/2,
                     verbose = par.threads < 2, method.p = "BH",
                     reduce.by.features = TRUE, search.window = 200000)

}
\arguments{
	\item{snp.regions.file}{
		path: Location of the regions file: Regions file is tab-deliminated and contains three elements per row. First element defines the genomic location of the tagSNP, 'chr:position' (e.g. 5:5420030). Second element contains the tagSNP name, 'rsID' (e.g. rs6010620). Third element defines the 'POPULATION' (ASN, EUR, AFR, ALL) where the tagSNP was identified (e.g. ASN, EUR, AFR, ALL).

		SNP Region file is imported and each row element (tagSNP element) is parsed for tagSNP name (rsXXXX), population (ASN, EUR, AFR, or ALL), and genomic location.  Genomic location is used to define the window size (see 'search.window' argument). See example file here: file.path(system.file('data',package='FunciSNP'), dir(system.file('data',package='FunciSNP'), pattern='.snp$')); 
	}
	\item{bio.features.loc}{
		path: Location of the biological features folder: Each biological feature for a particular genomic phenotype should be separated as individual BED files (tab deliminated file with chr, start and end). See UCSC for more information about BED formats \url{http://genome.ucsc.edu/FAQ/FAQformat.html#format1}. See example below. Default set to NULL.
	}
	\item{built.in.biofeatures}{
		logical: To include promoter regions, Encode DNaseI and CTCF sites as an additional biofeature in the analysis. Promoters defined as -1000 to +100 bp of a known TSS. File extracted on Feb. 9, 2012 from UCSC genome table browser. Default set to TRUE.
	}
	\item{par.threads}{
		an integer: Number of CPU cores to use for FunciSNP analysis. Default set at detectCores()/2. If par.threads > 1, then by default "verbose" = FALSE.	
	}
	\item{verbose}{
		logical: If set to TRUE, then regardless of par.threads value, all verbose message will output to terminal. If set to FALSE, no verbose message will output to terminal, except for warnings(). Default setting depends on number of 'par.threads' value.	
	}
	\item{method.p}{
		method: p-value correction (or adjustment) method (see ?p.adjust). Default set at "BH" (Benjamini & Hochberg (1995)).
	
	}
	\item{reduce.by.features}{
		logical: If set to TRUE, then only correlated SNPs overlapping biological features will be filtered and used to calculate Rsquared, Dprime, distance and p-value.  In addition, only these correlated SNPs will be used to generate plots and summary analysis.  If set to FALSE, all correlated SNPs regardless of overlap with biological features will have an associated Rsquare, Dprime, distance and p-value associated with the tag-SNP as defined by the 'search.window'.	
	}
	\item{search.window}{
		an integer: genomic window size used to extract all available correlated SNPs from the 1000 genomes db. The window size is centered around the tagSNP position as defined in the regions.file.
	}
}
\details{
This is the main funtion of FunciSNP. It will identify correlated SNPs which are in linkage disequilibrium (LD) to a known disease associated tagSNP. It will also determine if the correlated SNP in LD to the tagSNP overlaps a genomic biological feature. Correlated SNPs are directly imported from the current public release of the 1000 genomes database. 1000 genomes ftp servers available for the 1000 genomes public data: 1) National Center for Biotechnology Information (NCBI) \url{ftp://ftp-trace.ncbi.nih.gov/1000genomes/}; 2) European Bioinformatics Institute (EBI) \url{ftp://ftp.1000genomes.ebi.ac.uk/vol1/}.

Correlated SNPs in LD to a tagSNP and overlapping genomic biological features are known as putative functional SNPs (also defined as 'Func-y-SNP' elsewhere in the package.).  

%Using the location of a known tag-Single Nucleotide Polymorphisms (SNP) associated with a particular phenotype (e.g. disease, trait), FunciSNP will extract all available SNP from the 1000 Genomes Project within a selected window size (in bp, centered on each tag SNP). In addition, using the location of known biological features (e.g. statistically significant peaks identified by ChIP-seq, Chromatin ImmunoPrecipitation followed by genome-wide sequencing), correlated SNPs overlapping each peak is identified and further analyzed. The correlated SNPs overlapping a biological peak are then used to calculate the Rsquared, Dprime, p-value (as determined by fisher's exact test) and distance to the associated tag-SNP.
}
\value{
\item{TSList}{
	FunciSNP object. 
}
}
\references{
	SG. Coetzee, SK. Rhie, BP. Berman, GA. Coetzee and H. Noushmehr, \emph{FunciSNP: An R/Bioconductor Tool Integrating Functional Non-coding Datasets with Genetic Association Studies to Identify Candidate Regulatory SNPs}, Nucleic Acids Research, In press, 2012 (doi:10.1093/nar/GKS542).
}
\author{
	Simon G. Coetzee (maintainer: scoetzee@gmail.com); Houtan Noushmehr, PhD (houtan@usp.br)
}
\note{
    NA
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
#glioma <- getFSNPs(snp.regions.file=glioma.snp, 
#  bio.features.loc = glioma.bio);

##
data(glioma);
class(glioma);
glioma;
summary(glioma);


}
\keyword{ Functional }
\keyword{ SNP }% __ONLY ONE__ keyword per line
