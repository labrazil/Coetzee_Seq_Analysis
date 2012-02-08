setClass("CorrGeno",
         representation(SnpMatrix="SnpMatrix",
                        populations="list"
                        )
         )

setClass("CorrelatedSNPs",
         representation(chromosome="integer",
                        position="integer",
                        snpid="character",
                        ref.allele="character",
                        alt.allele="character",
                        overlapping.features="GRanges",
                        genotype="CorrGeno",
                        ALL.R.squared="matrix",
                        AFR.R.squared="matrix",
                        AMR.R.squared="matrix",
                        ASN.R.squared="matrix",
                        EUR.R.squared="matrix",
                        ALL.D.prime="matrix",
                        AFR.D.prime="matrix",
                        AMR.D.prime="matrix",
                        ASN.D.prime="matrix",
                        EUR.D.prime="matrix",
                        ALL.p.value="list",
                        AFR.p.value="list",
                        AMR.p.value="list",
                        ASN.p.value="list",
                        EUR.p.value="list"
                        )
         )

setClass("TagSNP",
         representation(chromosome="integer",
                        position="integer",
                        snpid="character",
                        population="character",
                        ref.allele="character",
                        alt.allele="character",
                        overlapping.features="GRanges",
                        genotype="SnpMatrix",
                        R.squared.corrsnps="dgCMatrix",
                        D.prime.corrsnps="dgCMatrix",
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
  summary.data <- FunciSNPAnnotateSummary(list.data)
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
setGeneric("overlapping.features", function(x) standardGeneric("overlapping.features"))
setMethod("overlapping.features", "TagSNP", function(x) x@overlapping.features)
setGeneric("overlapping.features<-", function(x, value) standardGeneric("overlapping.features<-"))
setReplaceMethod("overlapping.features", "TagSNP", function(x, value) {
                 x@overlapping.features <- value
                 x
         })
setGeneric("genotype", function(x) standardGeneric("genotype"))
setMethod("genotype", "TagSNP", function(x) x@genotype)
setGeneric("genotype<-", function(x, value) standardGeneric("genotype<-"))
setReplaceMethod("genotype", "TagSNP", function(x, value) {
                 x@genotype <- value
                 x
         })
setGeneric("R.squared.corrsnps", function(x) standardGeneric("R.squared.corrsnps"))
setMethod("R.squared.corrsnps", "TagSNP", function(x) x@R.squared.corrsnps)
setGeneric("R.squared.corrsnps<-", function(x, value) standardGeneric("R.squared.corrsnps<-"))
setReplaceMethod("R.squared.corrsnps", "TagSNP", function(x, value) {
                 x@R.squared.corrsnps <- value
                 x
         })
setGeneric("D.prime.corrsnps", function(x) standardGeneric("D.prime.corrsnps"))
setMethod("D.prime.corrsnps", "TagSNP", function(x) x@D.prime.corrsnps)
setGeneric("D.prime.corrsnps<-", function(x, value) standardGeneric("D.prime.corrsnps<-"))
setReplaceMethod("D.prime.corrsnps", "TagSNP", function(x, value) {
                 x@D.prime.corrsnps <- value
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
setMethod("overlapping.features", "CorrelatedSNPs", function(x) x@overlapping.features)
setReplaceMethod("overlapping.features", "CorrelatedSNPs", function(x, value) {
                 x@overlapping.features <- value
                 x
         })
setGeneric("genotype", function(object, populations) standardGeneric("genotype"))
setGeneric("genotype<-", function(object, value) standardGeneric("genotype<-"))
setMethod("genotype",
          "CorrelatedSNPs",
          function(object, populations="ALL"){
            x <- object@genotype@SnpMatrix[row.names(object@genotype@SnpMatrix)
                                           %in% object@genotype@populations[[population]]]
            return(x)
          }
          )
setReplaceMethod("genotype", "CorrelatedSNPs", function(object, value) {
                 object@genotype@SnpMatrix <- value
                 object
          })

setGeneric("ALL.R.squared", function(x) standardGeneric("ALL.R.squared"))
setGeneric("ALL.R.squared<-", function(x, value) standardGeneric("ALL.R.squared<-"))
setMethod("ALL.R.squared", "CorrelatedSNPs", function(x) x@ALL.R.squared)
setReplaceMethod("ALL.R.squared", "CorrelatedSNPs", function(x, value) {
                 x@ALL.R.squared <- value
                 x
         })
setGeneric("AFR.R.squared", function(x) standardGeneric("AFR.R.squared"))
setGeneric("AFR.R.squared<-", function(x, value) standardGeneric("AFR.R.squared<-"))
setMethod("AFR.R.squared", "CorrelatedSNPs", function(x) x@AFR.R.squared)
setReplaceMethod("AFR.R.squared", "CorrelatedSNPs", function(x, value) {
                 x@AFR.R.squared <- value
                 x
         })
setGeneric("AMR.R.squared", function(x) standardGeneric("AMR.R.squared"))
setGeneric("AMR.R.squared<-", function(x, value) standardGeneric("AMR.R.squared<-"))
setMethod("AMR.R.squared", "CorrelatedSNPs", function(x) x@AMR.R.squared)
setReplaceMethod("AMR.R.squared", "CorrelatedSNPs", function(x, value) {
                 x@AMR.R.squared <- value
                 x
         })
setGeneric("ASN.R.squared", function(x) standardGeneric("ASN.R.squared"))
setGeneric("ASN.R.squared<-", function(x, value) standardGeneric("ASN.R.squared<-"))
setMethod("ASN.R.squared", "CorrelatedSNPs", function(x) x@ASN.R.squared)
setReplaceMethod("ASN.R.squared", "CorrelatedSNPs", function(x, value) {
                 x@ASN.R.squared <- value
                 x
         })
setGeneric("EUR.R.squared", function(x) standardGeneric("EUR.R.squared"))
setGeneric("EUR.R.squared<-", function(x, value) standardGeneric("EUR.R.squared<-"))
setMethod("EUR.R.squared", "CorrelatedSNPs", function(x) x@EUR.R.squared)
setReplaceMethod("EUR.R.squared", "CorrelatedSNPs", function(x, value) {
                 x@EUR.R.squared <- value
                 x
         })
setGeneric("ALL.D.prime", function(x) standardGeneric("ALL.D.prime"))
setGeneric("ALL.D.prime<-", function(x, value) standardGeneric("ALL.D.prime<-"))
setMethod("ALL.D.prime", "CorrelatedSNPs", function(x) x@ALL.D.prime)
setReplaceMethod("ALL.D.prime", "CorrelatedSNPs", function(x, value) {
                 x@ALL.D.prime <- value
                 x
         })
setGeneric("AFR.D.prime", function(x) standardGeneric("AFR.D.prime"))
setGeneric("AFR.D.prime<-", function(x, value) standardGeneric("AFR.D.prime<-"))
setMethod("AFR.D.prime", "CorrelatedSNPs", function(x) x@AFR.D.prime)
setReplaceMethod("AFR.D.prime", "CorrelatedSNPs", function(x, value) {
                 x@AFR.D.prime <- value
                 x
         })
setGeneric("AMR.D.prime", function(x) standardGeneric("AMR.D.prime"))
setGeneric("AMR.D.prime<-", function(x, value) standardGeneric("AMR.D.prime<-"))
setMethod("AMR.D.prime", "CorrelatedSNPs", function(x) x@AMR.D.prime)
setReplaceMethod("AMR.D.prime", "CorrelatedSNPs", function(x, value) {
                 x@AMR.D.prime <- value
                 x
         })
setGeneric("ASN.D.prime", function(x) standardGeneric("ASN.D.prime"))
setGeneric("ASN.D.prime<-", function(x, value) standardGeneric("ASN.D.prime<-"))
setMethod("ASN.D.prime", "CorrelatedSNPs", function(x) x@ASN.D.prime)
setReplaceMethod("ASN.D.prime", "CorrelatedSNPs", function(x, value) {
                 x@ASN.D.prime <- value
                 x
         })
setGeneric("EUR.D.prime", function(x) standardGeneric("EUR.D.prime"))
setGeneric("EUR.D.prime<-", function(x, value) standardGeneric("EUR.D.prime<-"))
setMethod("EUR.D.prime", "CorrelatedSNPs", function(x) x@EUR.D.prime)
setReplaceMethod("EUR.D.prime", "CorrelatedSNPs", function(x, value) {
                 x@EUR.D.prime <- value
                 x
         })
setGeneric("ALL.p.value", function(x) standardGeneric("ALL.p.value"))
setGeneric("ALL.p.value<-", function(x, value) standardGeneric("ALL.p.value<-"))
setMethod("ALL.p.value", "CorrelatedSNPs", function(x) x@ALL.p.value)
setReplaceMethod("ALL.p.value", "CorrelatedSNPs", function(x, value) {
                 x@ALL.p.value <- value
                 x
         })
setGeneric("AFR.p.value", function(x) standardGeneric("AFR.p.value"))
setGeneric("AFR.p.value<-", function(x, value) standardGeneric("AFR.p.value<-"))
setMethod("AFR.p.value", "CorrelatedSNPs", function(x) x@AFR.p.value)
setReplaceMethod("AFR.p.value", "CorrelatedSNPs", function(x, value) {
                 x@AFR.p.value <- value
                 x
         })
setGeneric("AMR.p.value", function(x) standardGeneric("AMR.p.value"))
setGeneric("AMR.p.value<-", function(x, value) standardGeneric("AMR.p.value<-"))
setMethod("AMR.p.value", "CorrelatedSNPs", function(x) x@AMR.p.value)
setReplaceMethod("AMR.p.value", "CorrelatedSNPs", function(x, value) {
                 x@AMR.p.value <- value
                 x
         })
setGeneric("ASN.p.value", function(x) standardGeneric("ASN.p.value"))
setGeneric("ASN.p.value<-", function(x, value) standardGeneric("ASN.p.value<-"))
setMethod("ASN.p.value", "CorrelatedSNPs", function(x) x@ASN.p.value)
setReplaceMethod("ASN.p.value", "CorrelatedSNPs", function(x, value) {
                 x@ASN.p.value <- value
                 x
         })
setGeneric("EUR.p.value", function(x) standardGeneric("EUR.p.value"))
setGeneric("EUR.p.value<-", function(x, value) standardGeneric("EUR.p.value<-"))
setMethod("EUR.p.value", "CorrelatedSNPs", function(x) x@EUR.p.value)
setReplaceMethod("EUR.p.value", "CorrelatedSNPs", function(x, value) {
                 x@EUR.p.value <- value
                 x
         })
setGeneric("ALL.overlapping.snps.geno", function(object) standardGeneric("ALL.overlapping.snps.geno"))
setMethod("ALL.overlapping.snps.geno",
          "TagSNP",
          function(object){
            overlapping.snps <- unique(names(overlapping.features(correlated.snps(object))))
            y <- genotype(correlated.snps(x), "ALL")[, c(overlapping.snps, snpid(object))]
            return(y)
          })
setGeneric("AMR.overlapping.snps.geno", function(object) standardGeneric("AMR.overlapping.snps.geno"))
setMethod("AMR.overlapping.snps.geno",
          "TagSNP",
          function(object){
            overlapping.snps <- unique(names(overlapping.features(correlated.snps(object))))
            y <- genotype(correlated.snps(x), "AMR")[, c(overlapping.snps, snpid(object))]
            return(y)
          })
setGeneric("ASN.overlapping.snps.geno", function(object) standardGeneric("ASN.overlapping.snps.geno"))
setMethod("ASN.overlapping.snps.geno",
          "TagSNP",
          function(object){
            overlapping.snps <- unique(names(overlapping.features(correlated.snps(object))))
            y <- genotype(correlated.snps(x), "ASN")[, c(overlapping.snps, snpid(object))]
            return(y)
          })
setGeneric("AFR.overlapping.snps.geno", function(object) standardGeneric("AFR.overlapping.snps.geno"))
setMethod("AFR.overlapping.snps.geno",
          "TagSNP",
          function(object){
            overlapping.snps <- unique(names(overlapping.features(correlated.snps(object))))
            y <- genotype(correlated.snps(x), "AFR")[, c(overlapping.snps, snpid(object))]
            return(y)
          })
setGeneric("EUR.overlapping.snps.geno", function(object) standardGeneric("EUR.overlapping.snps.geno"))
setMethod("EUR.overlapping.snps.geno",
          "TagSNP",
          function(object){
            overlapping.snps <- unique(names(overlapping.features(correlated.snps(object))))
            y <- genotype(correlated.snps(x), "EUR")[, c(overlapping.snps, snpid(object))]
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
setMethod("show",
          signature=signature(object="TSList"),
          function(object){
            cat("TagSNP List with ", length(object), " Tag SNPs and \n"
                length(unique(object@summary.data$corr.snp.id)), "nearby, ",
                "potentially correlated SNPs, that overlap at least one biofeature \n")
              cat("Number of potentially correlated SNPs overlapping at least x",
                  " number of biofeatures at an R squared of\n")
            for(i in c(0.20, 0.50, 0.60, 0.70, 0.80, 0.90)) {
              summ.overlaps <- FunciSNPsummaryOverlaps(object@summary.data, rsq=i)
              cat("at least ", i)
              cat(summ.overlaps)
            }
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
