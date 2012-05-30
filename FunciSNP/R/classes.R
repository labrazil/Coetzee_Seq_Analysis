setClass("CorrGeno",
         representation(SnpMatrix="SnpMatrix",
                        populations="list"
                        )
         )

setClass("CorrelatedSNPs",
         representation(chromosome="character",
                        position="integer",
                        snpid="character",
                        ref.allele="character",
                        alt.allele="character",
                        overlapping.features="GRanges",
                        genotype="CorrGeno"
                        )
         )

setClass("TagSNP",
         representation(chromosome="character",
                        position="integer",
                        snpid="character",
                        population="character",
                        ref.allele="character",
                        alt.allele="character",
#                        overlapping.features="GRanges",
                        genotype="SnpMatrix",
                        yafsnp.rsq="dgCMatrix",
                        yafsnp.dprime="dgCMatrix",
                        yafsnp.pvalue="list",
                        correlated.snps="CorrelatedSNPs"
                        )
         )
setClass("TSList",
         contains = "List",
         representation(snp.data="list",
                        summary.data="data.frame"
                        ),
         prototype(elementType = "TagSNP"
                   ),
         )
### TSList Constructor
TSList <- function(...) {
  list.data <- list(...)
  if(length(list.data) == 0L) {
    snp.data <- new("TagSNP")
  } else {
    if(length(list.data) == 1L && is.list(list.data[[1L]]))
      list.data <- list.data[[1L]]
    if(!all(sapply(list.data, is, "TagSNP")))
      stop("all elements in '...' must be TagSNP objects")
    list.data <- do.call("c", list.data)
  }
  summary.data <- AnnotateSummary(list.data)
  x <- new("TSList", snp.data=list.data, summary.data=summary.data)
  return(x)
}


setGeneric("chr", function(x) standardGeneric("chr"))
setMethod("chr", "TagSNP", function(x) x@chromosome)
setGeneric("chr<-", function(x, value) standardGeneric("chr<-"))
setReplaceMethod("chr", "TagSNP", function(x, value) {
                 x@chromosome <- value
                 x
         })
setGeneric("position", function(x) standardGeneric("position"))
setMethod("position", "TagSNP", function(x) x@position)
setGeneric("position<-", function(x, value) standardGeneric("position<-"))
setReplaceMethod("position", "TagSNP", function(x, value) {
                 x@position <- value
                 x
         })
setGeneric("snpid", function(x) standardGeneric("snpid"))
setMethod("snpid", "TagSNP", function(x) x@snpid)
setGeneric("snpid<-", function(x, value) standardGeneric("snpid<-"))
setReplaceMethod("snpid", "TagSNP", function(x, value) {
                 x@snpid <- value
                 x
         })
setGeneric("population", function(x) standardGeneric("population"))
setMethod("population", "TagSNP", function(x) x@population)
setGeneric("population<-", function(x, value) standardGeneric("population<-"))
setReplaceMethod("population", "TagSNP", function(x, value) {
                 x@population <- value
                 x
         })
setGeneric("ref.allele", function(x) standardGeneric("ref.allele"))
setMethod("ref.allele", "TagSNP", function(x) x@ref.allele)
setGeneric("ref.allele<-", function(x, value) standardGeneric("ref.allele<-"))
setReplaceMethod("ref.allele", "TagSNP", function(x, value) {
                 x@ref.allele <- value
                 x
         })
setGeneric("alt.allele", function(x) standardGeneric("alt.allele"))
setMethod("alt.allele", "TagSNP", function(x) x@alt.allele)
setGeneric("alt.allele<-", function(x, value) standardGeneric("alt.allele<-"))
setReplaceMethod("alt.allele", "TagSNP", function(x, value) {
                 x@alt.allele <- value
                 x
         })
#setGeneric("overlapping.features", function(x) standardGeneric("overlapping.features"))
#setMethod("overlapping.features", "TagSNP", function(x) x@overlapping.features)
#setGeneric("overlapping.features<-", function(x, value) standardGeneric("overlapping.features<-"))
#setReplaceMethod("overlapping.features", "TagSNP", function(x, value) {
#                 x@overlapping.features <- value
#                 x
#         })
setGeneric("genotype", function(x) standardGeneric("genotype"))
setMethod("genotype", "TagSNP", function(x) x@genotype)
setGeneric("genotype<-", function(x, value) standardGeneric("genotype<-"))
setReplaceMethod("genotype", "TagSNP", function(x, value) {
                 x@genotype <- value
                 x
         })
setGeneric("yafsnp.rsq", function(x) standardGeneric("yafsnp.rsq"))
setMethod("yafsnp.rsq", "TagSNP", function(x) x@yafsnp.rsq)
setGeneric("yafsnp.rsq<-", function(x, value) standardGeneric("yafsnp.rsq<-"))
setReplaceMethod("yafsnp.rsq", "TagSNP", function(x, value) {
                 x@yafsnp.rsq <- value
                 x
         })
setGeneric("yafsnp.dprime", function(x) standardGeneric("yafsnp.dprime"))
setMethod("yafsnp.dprime", "TagSNP", function(x) x@yafsnp.dprime)
setGeneric("yafsnp.dprime<-", function(x, value) standardGeneric("yafsnp.dprime<-"))
setReplaceMethod("yafsnp.dprime", "TagSNP", function(x, value) {
                 x@yafsnp.dprime <- value
                 x
         })
setGeneric("yafsnp.pvalue", function(x) standardGeneric("yafsnp.pvalue"))
setMethod("yafsnp.pvalue", "TagSNP", function(x) x@yafsnp.pvalue)
setGeneric("yafsnp.pvalue<-", function(x, value) standardGeneric("yafsnp.pvalue<-"))
setReplaceMethod("yafsnp.pvalue", "TagSNP", function(x, value) {
                 x@yafsnp.pvalue <- value
                 x
         })
setGeneric("correlated.snps", function(x) standardGeneric("correlated.snps"))
setMethod("correlated.snps", "TagSNP", function(x) x@correlated.snps)
setGeneric("correlated.snps<-", function(x, value) standardGeneric("correlated.snps<-"))
setReplaceMethod("correlated.snps", "TagSNP", function(x, value) {
                 x@correlated.snps <- value
                 x
         })



##### CorrelatedSNPs methods and generics
setMethod("chr", "CorrelatedSNPs", function(x) x@chromosome)
setReplaceMethod("chr", "CorrelatedSNPs", function(x, value) {
                 x@chromosome <- value
                 x
         })
setMethod("position", "CorrelatedSNPs", function(x) x@position)
setReplaceMethod("position", "CorrelatedSNPs", function(x, value) {
                 x@position <- value
                 x
         })
setMethod("snpid", "CorrelatedSNPs", function(x) x@snpid)
setReplaceMethod("snpid", "CorrelatedSNPs", function(x, value) {
                 x@snpid <- value
                 x
         })
setMethod("ref.allele", "CorrelatedSNPs", function(x) x@ref.allele)
setReplaceMethod("ref.allele", "CorrelatedSNPs", function(x, value) {
                 x@ref.allele <- value
                 x
         })
setMethod("alt.allele", "CorrelatedSNPs", function(x) x@alt.allele)
setReplaceMethod("alt.allele", "CorrelatedSNPs", function(x, value) {
                 x@alt.allele <- value
                 x
         })
setGeneric("overlapping.features", function(x) standardGeneric("overlapping.features"))
setGeneric("overlapping.features<-", function(x, value) standardGeneric("overlapping.features<-"))
setMethod("overlapping.features", "CorrelatedSNPs", function(x) x@overlapping.features)
setReplaceMethod("overlapping.features", "CorrelatedSNPs", function(x, value) {
                 x@overlapping.features <- value
                 x
         })
setGeneric("pop.genotype", function(x, population) standardGeneric("pop.genotype"))
setGeneric("pop.genotype<-", function(x, value) standardGeneric("pop.genotype<-"))
setMethod("pop.genotype",
          "CorrelatedSNPs",
          function(x, population="ALL"){
            snps <- x@genotype@SnpMatrix
            populations <- x@genotype@populations
            x <- snps[row.names(snps) %in% populations[[population]]]
            return(x)
          }
          )
setReplaceMethod("pop.genotype", "CorrelatedSNPs", function(x, value) {
                 x@genotype@SnpMatrix <- value
                 x
          })

setGeneric("ALL.overlapping.snps.geno", function(object) standardGeneric("ALL.overlapping.snps.geno"))
setMethod("ALL.overlapping.snps.geno",
          "TagSNP",
          function(object){
            overlapping.snps <- unique(names(overlapping.features(correlated.snps(object))))
            y <- pop.genotype(correlated.snps(object), "ALL")[, c(overlapping.snps, snpid(object))]
            return(y)
          })
setGeneric("AMR.overlapping.snps.geno", function(object) standardGeneric("AMR.overlapping.snps.geno"))
setMethod("AMR.overlapping.snps.geno",
          "TagSNP",
          function(object){
            overlapping.snps <- unique(names(overlapping.features(correlated.snps(object))))
            y <- pop.genotype(correlated.snps(object), "AMR")[, c(overlapping.snps, snpid(object))]
            return(y)
          })
setGeneric("ASN.overlapping.snps.geno", function(object) standardGeneric("ASN.overlapping.snps.geno"))
setMethod("ASN.overlapping.snps.geno",
          "TagSNP",
          function(object){
            overlapping.snps <- unique(names(overlapping.features(correlated.snps(object))))
            y <- pop.genotype(correlated.snps(object), "ASN")[, c(overlapping.snps, snpid(object))]
            return(y)
          })
setGeneric("AFR.overlapping.snps.geno", function(object) standardGeneric("AFR.overlapping.snps.geno"))
setMethod("AFR.overlapping.snps.geno",
          "TagSNP",
          function(object){
            overlapping.snps <- unique(names(overlapping.features(correlated.snps(object))))
            y <- pop.genotype(correlated.snps(object), "AFR")[, c(overlapping.snps, snpid(object))]
            return(y)
          })
setGeneric("EUR.overlapping.snps.geno", function(object) standardGeneric("EUR.overlapping.snps.geno"))
setMethod("EUR.overlapping.snps.geno",
          "TagSNP",
          function(object){
            overlapping.snps <- unique(names(overlapping.features(correlated.snps(object))))
            y <- pop.genotype(correlated.snps(object), "EUR")[, c(overlapping.snps, snpid(object))]
            return(y)
          })
setMethod("show",
          signature=signature(object="TagSNP"),
          function(object){
            cat("class:", class(object), "\n")
            cat("Tag SNP: ", snpid(object), "\n")
            if(dim(elementMetadata(overlapping.features(correlated.snps(object))))[1] > 0) {
              cat(length(unique(names(overlapping.features(correlated.snps(object))))), 
                  "Surrounding SNPs are overlapping: \n")
              cat(unique(elementMetadata(overlapping.features(correlated.snps(object)))[, "feature"]), "\n")
            } else {
              cat("No Surrounding SNPs are overlapping any biofeatures\n")
            }
          })
setMethod("summary",
          signature=signature(object="TSList"),
          function(object){
            cat("TagSNP List with ", length(object@snp.data), " Tag SNPs and \n",
                length(unique(object@summary.data$corr.snp.id)), "nearby, ",
                "potentially correlated SNPs, that overlap at least one biofeature \n")
              cat("Number of potentially correlated SNPs", "\noverlapping at least x",
                  "biofeatures, per Tag SNP at a specified R squared\n")
            r.squared.values <- c(0.10, 0.50, 0.90)
            bio.features.counting <- list()
            for(i in r.squared.values) {
              summ.overlaps <- FunciSNPsummaryOverlaps(object@summary.data, rsq=i)
              bio.features.counting[[ as.character(paste("R squared: ", 
                                                         i, 
                                                         " in ", 
                                                         nrow(summ.overlaps) - 1, 
                                                         " Tag SNPs with a total of ", 
                                                         sep="")) ]] <- summ.overlaps
            }
            print(bio.features.counting, quote=F)
          })
setMethod("show",
          signature=signature(object="TSList"),
          function(object){
            cat("TagSNP List with ", length(object@snp.data), " Tag SNPs and \n",
                length(unique(object@summary.data$corr.snp.id)), "nearby, ",
                "potentially correlated SNPs, that overlap at least one biofeature \n")
            r.squared.values <- c(0.10, 0.50, 0.90)
            funcitable <- list()
            for(i in r.squared.values) {
              summ.overlaps <- FunciSNPtable(object@summary.data, rsq=i)
              funcitable[[ as.character(paste("R squared: ", 
                                                         i, 
                                                         sep="")) ]] <- summ.overlaps
            }
            print(funcitable, quote=F)
          })
#setMethod("summary", "TagSNP", function(object) {
#          data.frame(chromosome=chromosome(object), position=position(object),
#                     snpid=snpid(object), population=population(object),
#                     ref.allele=ref.allele(object), alt.allele=alt.allele(object),
#                     overlapping.features=lapply(overlapping.features(object), function(x) length(x)),
#                     genotype=genotype(object), R.squared.CorrelatedSNPs=R.squared.CorrelatedSNPs(object),
#                     D.prime.CorrelatedSNPs=D.prime.CorrelatedSNPs(object))
#         }
#)
