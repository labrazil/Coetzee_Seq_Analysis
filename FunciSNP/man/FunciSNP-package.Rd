\name{FunciSNP}
\alias{FunciSNP}
\title{
Functional Identification of SNPs with Phenotype by Coincidence with Chromatin Biofeatures	
}
\description{
Given a set of known tag-SNPs associated with a particular phenotype (e.g. disease, trait), and a set of available biological features (e.g. peaks derived from ChIP-seq experiments for phenotype), returns correlated SNPs (from the 1000 genomes db) which are in linkage disequilibrium (LD) to a known disease associated tag-SNP and overlaps chromatin biological features.  These identified correlated SNPs are characterized as putative functional SNPs for a particular trait. 
}
\usage{
FunciSNP(snp.regions.file, bio.features.loc = NULL, primary.server = "ncbi", 
par.threads = detectCores()/2, verbose = par.threads < 2, 
method.p = "BH", reduce.by.features = TRUE, 
search.window = 2e+05)

	detectCores()
	#function implemented from library(parallel) to identify number of 
	#CPUs running on current R session.
}
\arguments{
	\item{snp.regions.file}{
		Location of the regions file: Regions file is tab-deliminated and contains three elements per row. First element defines the genomic location of the tagSNP, 'chr:position' (e.g. 5:5420030). Second element contains the tagSNP name, 'rsID' (e.g. rs6010620). Third element defines the 'POPULATION' (ASN, EUR, AFR, ALL) where the tagSNP was identified (e.g. ASN, EUR, AFR, ALL).

		SNP Region file is imported and each row element (tagSNP element) is parsed for tagSNP name (rsXXXX), population (ASN, EUR, AFR, or ALL), and genomic location.  Genomic location is used to define the window size (see 'search.window' argument). See example file here: file.path(system.file('data',package='FunciSNP'), dir(system.file('data',package='FunciSNP'), pattern='.snp$')); 
	}
	\item{bio.features.loc}{
		Location of the biological features folder: Each biological features for a particular genomic phenotype should be stored as separate BED files (tab deliminated file with chr, start and end). See UCSC for more information about BED formats \url{http://genome.ucsc.edu/FAQ/FAQformat.html#format1}. See example below.
	}
	\item{primary.server}{
		1000 genomes ftp server: Two 'File Transfer Protocol' (FTP) mirrors are available for the 1000 genomes public data: 1) National Center for Biotechnology Information (NCBI); 2) European Bioinformatics Institute (EBI).  To use the NCBI server (\url{ftp://ftp-trace.ncbi.nih.gov/1000genomes/}), use "ncbi". To use the EBI server (\url{ftp://ftp.1000genomes.ebi.ac.uk/vol1/}), use "ebi". Default set at "ncbi".
	}
	\item{par.threads}{
		an integer: Number of CPU cores to use for FunciSNP analysis. Default set at detectCores()/2. If par.threads > 1, then by default "verbose" = FALSE.	}
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
Identification of correlated SNPs (from the 1000 genomes db) which are in Linkage Disequilibrium (LD) to a known disease associated tag-SNP and which the identified correlated SNP overlaps chromatin biological features. Correlated SNP in LD to a tag-SNP and overlapping chromatin biological features are known as putative functional SNP.  
Using the location of a known tag-Single Nucleotide Polymorphisms (SNP) associated with a particular phenotype (e.g. disease, trait), FunciSNP will extract all available SNP from the 1000 Genomes Project within a selected window size (in bp, centered on each tag SNP). In addition, using the location of known biological features (e.g. statistically significant peaks identified by ChIP-seq, Chromatin ImmunoPrecipitation followed by genome-wide sequencing), correlated SNPs overlapping each peak is identified and further analyzed. The correlated SNPs overlapping a biological peak are then used to calculate the Rsquared, Dprime, p-value (as determined by fisher's exact test) and distance to the associated tag-SNP.
}
\value{
##%%  #Describe the value returned
##%%  If it is a LIST, use
##%%  \item{comp1 }{Description of 'comp1'}
##%%  \item{comp2 }{Description of 'comp2'}
##%% ...
}
\references{
	Coetzee SG et al. submitted for review. 2012
}
\author{
	Simon Coetzee, Houtan Noushmehr
}
\note{
}


\seealso{
	\code{\link{FunciSNPplot}}
}
\examples{
##
## Glioblastoma analysis using FunciSNP
##
## Full path to the example regions file for Glioblastoma 
#  (collected from SNPedia)
glioma.snp <- file.path(system.file('data',
  package='FunciSNP'),
  dir(system.file('data',package='FunciSNP'), 
  pattern='.snp$'));
 
## Full path to the example biological features BED files 
#  derived from the ENCODE project for Glioblastoma U-87 
#  cell lines.
glioma.bio <- system.file('extdata',package='FunciSNP');

## FunciSNP analysis, extracts correlated SNPs from the 
#  1000 genomes db ("ncbi") and finds overlaps between 
#  correlated SNP and biological features and then 
#  calculates LD (Rsquare, Dprime, distance, p-value).
glioma <- FunciSNP(snp.regions.file=glioma.snp, 
  bio.features.loc = glioma.bio, bio.features.TSS=FALSE);

##
data(glioma);
data(glioma_anno);

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ #kwd1 }
\keyword{ #kwd2 }% __ONLY ONE__ keyword per line
