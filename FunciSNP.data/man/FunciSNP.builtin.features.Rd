\name{FunciSNP.builtin.features}
\alias{FunciSNP.data}
\alias{FunciSNP.builtin.features}
\alias{ctcf_only}
\alias{encode_dnase1_only}
\alias{encode_dnase1_with_ctcf}
\alias{encode_faire}
\alias{known_gene_promoters}
\docType{data}
\title{Builtin biological features used by FunciSNP package}
\description{
    These data were taken from the UCSC table browser and Encode website. Data
retrieved on March of 2012. They include CTCF sites do not overlap a known
Encode DNAseI site, Encode DNAseI site containing a CTCF site, Encode DNAseI not
containing a CTCF site, ENCODE FAIRE clusters, promoters of known genes defined
by -1000 and +100bp of a known transcription start site (TSS). These data were
created to be used in the vignette for the \pkg{FunciSNP} pacakge.

}
\usage{

ctcf_only
encode_dnase1_only
encode_dnase1_with_ctcf
encode_faire
known_gene_promoters

}
\details{
	Each data set was extracted from the ENCODE and/or UCSC table browser on
March of 2012. We retained only autosomal and sex linked chromosomes, removing
unknown and random chromosomes. See reference for detail information on how
these bed files were generated. Objects are stored as RangedData.
}
\seealso{
	\code{\link[FunciSNP.data]{FunciSNP.builtin.features}}
	\code{\link[FunciSNP.data]{lincRNA.hg19.rda}}
	\code{\link[FunciSNP.data]{refseqgenes.rda}}
}
\references{
SG. Coetzee, SK. Rhie, BP. Berman, GA. Coetzee and H. Noushmehr, FunciSNP: An
R/Bioconductor Tool Integrating Functional Non-coding Datasets with Genetic
Association Studies to Identify Candidate Regulatory SNPs., Nucleic Acids
Research, In press, 2012 (doi:10.1093/nar/gks542).
}
\examples{
dim(ctcf_only);
dim(encode_dnase1_only);
dim(encode_dnase1_with_ctcf);
dim(encode_faire);
dim(known_gene_promoters);
class(ctcf_only);
}
\keyword{datasets}
