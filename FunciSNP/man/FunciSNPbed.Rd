\name{FunciSNPbed}
\alias{FunciSNPbed}

\title{
Creates a BED file to view Func-y-SNPs in your favorite genome browser
}
\description{
FunciSNPbed will output a BED file to a specified folder. The BED file is in standard UCSC Genome Browser format (\url{http://genome.ucsc.edu/FAQ/FAQformat}). Each tagSNP is colored black and each Func-y-SNP is colored red.

}
\usage{
FunciSNPbed(dat, rsq, path = getwd(), filename = NULL)
}

\arguments{
  \item{dat}{
             FunciSNP data.frame: \code{dat} is a data.frame object from FunciSNPAnnotateSummary. Need to run \code{\link{FunciSNPAnnotateSummary}} first.
}
  \item{rsq}{
            an interger (0-1): \code{rsq} is the Rsquared cutoff used to subset.
}
  \item{path}{
           a character: \code{path} is the path to the folder where to save the BED file. Default to getwd() or current working directory.
}
  \item{filename}{
           a character: \code{filename} is the name of the BED file. If NULL, filename is 'FunciSNP_results_rsq.RSQ value.bed'
}
}
\details{
FunciSNPbed outputs a unique BED file which can be used to view in any genomic browser compatible with BED formats. To learn more about BED formats, see UCSC Genome Browser FAQ (\url{http://genome.ucsc.edu/FAQ/FAQformat}). Each tagSNP which is in LD to a corresponding Func-y-SNP overlapping at least one biofeature is colored black, while the Func-y-SNP is colored red. The initial position is provided by the first tagSNP and the first linked Func-y-SNP. We recommend using UCSC genome browser to view your BED files. This is useful so you can view all public and private tracks in relation to FunciSNP results.
}
\value{
BED file is outputed as a tab-deliminated file in the specified 'path' folder. See example below.
}
\references{
	Coetzee SG et al. submitted for review. 2012
}
\author{
	Simon Coetzee, Houtan Noushmehr
}
\note{
NA
}



\seealso{
	\code{\link{getFSNPs}}, \code{\link{FunciSNPplot}}, \code{\link{FunciSNPAnnotateSummary}}, \code{\link{FunciSNPtable}}, \code{\link{FunciSNPbed}}
}
\examples{
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



}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
