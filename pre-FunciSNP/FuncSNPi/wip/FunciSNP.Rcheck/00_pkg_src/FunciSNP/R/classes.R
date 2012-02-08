setClass("corrsnps",
         representation(chromosome="integer",
                        position="integer",
                        snpid="character",
                        ref.allele="character",
                        alt.allele="character",
                        overlapping.features="GRanges",
                        ALL.geno="SnpMatrix",
                        AFR.geno="SnpMatrix",
                        AMR.geno="SnpMatrix",
                        ASN.geno="SnpMatrix",
                        EUR.geno="SnpMatrix",
                        ALL.geno.overlapping.snps="SnpMatrix",
                        AFR.geno.overlapping.snps="SnpMatrix",
                        AMR.geno.overlapping.snps="SnpMatrix",
                        ASN.geno.overlapping.snps="SnpMatrix",
                        EUR.geno.overlapping.snps="SnpMatrix",
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

setClass("tagsnp",
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
                        correlated.snps="corrsnps"
                        )
         )

setGeneric("chromosome", function(x) standardGeneric("chromosome"))
setMethod("chromosome", "tagsnp", function(x) x@chromosome)
setGeneric("chromosome<-", function(x, value) standardGeneric("chromosome<-"))
setReplaceMethod("chromosome", "tagsnp", function(x, value) {
                 x@chromosome <- value
                 x
         })
setGeneric("position", function(x) standardGeneric("position"))
setMethod("position", "tagsnp", function(x) x@position)
setGeneric("position<-", function(x, value) standardGeneric("position<-"))
setReplaceMethod("position", "tagsnp", function(x, value) {
                 x@position <- value
                 x
         })
setGeneric("snpid", function(x) standardGeneric("snpid"))
setMethod("snpid", "tagsnp", function(x) x@snpid)
setGeneric("snpid<-", function(x, value) standardGeneric("snpid<-"))
setReplaceMethod("snpid", "tagsnp", function(x, value) {
                 x@snpid <- value
                 x
         })
setGeneric("population", function(x) standardGeneric("population"))
setMethod("population", "tagsnp", function(x) x@population)
setGeneric("population<-", function(x, value) standardGeneric("population<-"))
setReplaceMethod("population", "tagsnp", function(x, value) {
                 x@population <- value
                 x
         })
setGeneric("ref.allele", function(x) standardGeneric("ref.allele"))
setMethod("ref.allele", "tagsnp", function(x) x@ref.allele)
setGeneric("ref.allele<-", function(x, value) standardGeneric("ref.allele<-"))
setReplaceMethod("ref.allele", "tagsnp", function(x, value) {
                 x@ref.allele <- value
                 x
         })
setGeneric("alt.allele", function(x) standardGeneric("alt.allele"))
setMethod("alt.allele", "tagsnp", function(x) x@alt.allele)
setGeneric("alt.allele<-", function(x, value) standardGeneric("alt.allele<-"))
setReplaceMethod("alt.allele", "tagsnp", function(x, value) {
                 x@alt.allele <- value
                 x
         })
setGeneric("overlapping.features", function(x) standardGeneric("overlapping.features"))
setMethod("overlapping.features", "tagsnp", function(x) x@overlapping.features)
setGeneric("overlapping.features<-", function(x, value) standardGeneric("overlapping.features<-"))
setReplaceMethod("overlapping.features", "tagsnp", function(x, value) {
                 x@overlapping.features <- value
                 x
         })
setGeneric("genotype", function(x) standardGeneric("genotype"))
setMethod("genotype", "tagsnp", function(x) x@genotype)
setGeneric("genotype<-", function(x, value) standardGeneric("genotype<-"))
setReplaceMethod("genotype", "tagsnp", function(x, value) {
                 x@genotype <- value
                 x
         })
setGeneric("R.squared.corrsnps", function(x) standardGeneric("R.squared.corrsnps"))
setMethod("R.squared.corrsnps", "tagsnp", function(x) x@R.squared.corrsnps)
setGeneric("R.squared.corrsnps<-", function(x, value) standardGeneric("R.squared.corrsnps<-"))
setReplaceMethod("R.squared.corrsnps", "tagsnp", function(x, value) {
                 x@R.squared.corrsnps <- value
                 x
         })
setGeneric("D.prime.corrsnps", function(x) standardGeneric("D.prime.corrsnps"))
setMethod("D.prime.corrsnps", "tagsnp", function(x) x@D.prime.corrsnps)
setGeneric("D.prime.corrsnps<-", function(x, value) standardGeneric("D.prime.corrsnps<-"))
setReplaceMethod("D.prime.corrsnps", "tagsnp", function(x, value) {
                 x@D.prime.corrsnps <- value
                 x
         })
setGeneric("correlated.snps", function(x) standardGeneric("correlated.snps"))
setMethod("correlated.snps", "tagsnp", function(x) x@correlated.snps)
setGeneric("correlated.snps<-", function(x, value) standardGeneric("correlated.snps<-"))
setReplaceMethod("correlated.snps", "tagsnp", function(x, value) {
                 x@correlated.snps <- value
                 x
         })



##### corrsnps methods and generics
setMethod("chromosome", "corrsnps", function(x) x@chromosome)
setReplaceMethod("chromosome", "corrsnps", function(x, value) {
                 x@chromosome <- value
                 x
         })
setMethod("position", "corrsnps", function(x) x@position)
setReplaceMethod("position", "corrsnps", function(x, value) {
                 x@position <- value
                 x
         })
setMethod("snpid", "corrsnps", function(x) x@snpid)
setReplaceMethod("snpid", "corrsnps", function(x, value) {
                 x@snpid <- value
                 x
         })
setMethod("ref.allele", "corrsnps", function(x) x@ref.allele)
setReplaceMethod("ref.allele", "corrsnps", function(x, value) {
                 x@ref.allele <- value
                 x
         })
setMethod("alt.allele", "corrsnps", function(x) x@alt.allele)
setReplaceMethod("alt.allele", "corrsnps", function(x, value) {
                 x@alt.allele <- value
                 x
         })
setMethod("overlapping.features", "corrsnps", function(x) x@overlapping.features)
setReplaceMethod("overlapping.features", "corrsnps", function(x, value) {
                 x@overlapping.features <- value
                 x
         })
setGeneric("ALL.geno", function(x) standardGeneric("ALL.geno"))
setGeneric("ALL.geno<-", function(x, value) standardGeneric("ALL.geno<-"))
setMethod("ALL.geno", "corrsnps", function(x) x@ALL.geno)
setReplaceMethod("ALL.geno", "corrsnps", function(x, value) {
                 x@ALL.geno <- value
                 x
         })
setGeneric("AFR.geno", function(x) standardGeneric("AFR.geno"))
setGeneric("AFR.geno<-", function(x, value) standardGeneric("AFR.geno<-"))
setMethod("AFR.geno", "corrsnps", function(x) x@AFR.geno)
setReplaceMethod("AFR.geno", "corrsnps", function(x, value) {
                 x@AFR.geno <- value
                 x
         })
setGeneric("AMR.geno", function(x) standardGeneric("AMR.geno"))
setGeneric("AMR.geno<-", function(x, value) standardGeneric("AMR.geno<-"))
setMethod("AMR.geno", "corrsnps", function(x) x@AMR.geno)
setReplaceMethod("AMR.geno", "corrsnps", function(x, value) {
                 x@AMR.geno <- value
                 x
         })
setGeneric("ASN.geno", function(x) standardGeneric("ASN.geno"))
setGeneric("ASN.geno<-", function(x, value) standardGeneric("ASN.geno<-"))
setMethod("ASN.geno", "corrsnps", function(x) x@ASN.geno)
setReplaceMethod("ASN.geno", "corrsnps", function(x, value) {
                 x@ASN.geno <- value
                 x
         })
setGeneric("EUR.geno", function(x) standardGeneric("EUR.geno"))
setGeneric("EUR.geno<-", function(x, value) standardGeneric("EUR.geno<-"))
setMethod("EUR.geno", "corrsnps", function(x) x@EUR.geno)
setReplaceMethod("EUR.geno", "corrsnps", function(x, value) {
                 x@EUR.geno <- value
                 x
         })

setGeneric("ALL.R.squared", function(x) standardGeneric("ALL.R.squared"))
setGeneric("ALL.R.squared<-", function(x, value) standardGeneric("ALL.R.squared<-"))
setMethod("ALL.R.squared", "corrsnps", function(x) x@ALL.R.squared)
setReplaceMethod("ALL.R.squared", "corrsnps", function(x, value) {
                 x@ALL.R.squared <- value
                 x
         })
setGeneric("AFR.R.squared", function(x) standardGeneric("AFR.R.squared"))
setGeneric("AFR.R.squared<-", function(x, value) standardGeneric("AFR.R.squared<-"))
setMethod("AFR.R.squared", "corrsnps", function(x) x@AFR.R.squared)
setReplaceMethod("AFR.R.squared", "corrsnps", function(x, value) {
                 x@AFR.R.squared <- value
                 x
         })
setGeneric("AMR.R.squared", function(x) standardGeneric("AMR.R.squared"))
setGeneric("AMR.R.squared<-", function(x, value) standardGeneric("AMR.R.squared<-"))
setMethod("AMR.R.squared", "corrsnps", function(x) x@AMR.R.squared)
setReplaceMethod("AMR.R.squared", "corrsnps", function(x, value) {
                 x@AMR.R.squared <- value
                 x
         })
setGeneric("ASN.R.squared", function(x) standardGeneric("ASN.R.squared"))
setGeneric("ASN.R.squared<-", function(x, value) standardGeneric("ASN.R.squared<-"))
setMethod("ASN.R.squared", "corrsnps", function(x) x@ASN.R.squared)
setReplaceMethod("ASN.R.squared", "corrsnps", function(x, value) {
                 x@ASN.R.squared <- value
                 x
         })
setGeneric("EUR.R.squared", function(x) standardGeneric("EUR.R.squared"))
setGeneric("EUR.R.squared<-", function(x, value) standardGeneric("EUR.R.squared<-"))
setMethod("EUR.R.squared", "corrsnps", function(x) x@EUR.R.squared)
setReplaceMethod("EUR.R.squared", "corrsnps", function(x, value) {
                 x@EUR.R.squared <- value
                 x
         })
setGeneric("ALL.D.prime", function(x) standardGeneric("ALL.D.prime"))
setGeneric("ALL.D.prime<-", function(x, value) standardGeneric("ALL.D.prime<-"))
setMethod("ALL.D.prime", "corrsnps", function(x) x@ALL.D.prime)
setReplaceMethod("ALL.D.prime", "corrsnps", function(x, value) {
                 x@ALL.D.prime <- value
                 x
         })
setGeneric("AFR.D.prime", function(x) standardGeneric("AFR.D.prime"))
setGeneric("AFR.D.prime<-", function(x, value) standardGeneric("AFR.D.prime<-"))
setMethod("AFR.D.prime", "corrsnps", function(x) x@AFR.D.prime)
setReplaceMethod("AFR.D.prime", "corrsnps", function(x, value) {
                 x@AFR.D.prime <- value
                 x
         })
setGeneric("AMR.D.prime", function(x) standardGeneric("AMR.D.prime"))
setGeneric("AMR.D.prime<-", function(x, value) standardGeneric("AMR.D.prime<-"))
setMethod("AMR.D.prime", "corrsnps", function(x) x@AMR.D.prime)
setReplaceMethod("AMR.D.prime", "corrsnps", function(x, value) {
                 x@AMR.D.prime <- value
                 x
         })
setGeneric("ASN.D.prime", function(x) standardGeneric("ASN.D.prime"))
setGeneric("ASN.D.prime<-", function(x, value) standardGeneric("ASN.D.prime<-"))
setMethod("ASN.D.prime", "corrsnps", function(x) x@ASN.D.prime)
setReplaceMethod("ASN.D.prime", "corrsnps", function(x, value) {
                 x@ASN.D.prime <- value
                 x
         })
setGeneric("EUR.D.prime", function(x) standardGeneric("EUR.D.prime"))
setGeneric("EUR.D.prime<-", function(x, value) standardGeneric("EUR.D.prime<-"))
setMethod("EUR.D.prime", "corrsnps", function(x) x@EUR.D.prime)
setReplaceMethod("EUR.D.prime", "corrsnps", function(x, value) {
                 x@EUR.D.prime <- value
                 x
         })
setGeneric("ALL.p.value", function(x) standardGeneric("ALL.p.value"))
setGeneric("ALL.p.value<-", function(x, value) standardGeneric("ALL.p.value<-"))
setMethod("ALL.p.value", "corrsnps", function(x) x@ALL.p.value)
setReplaceMethod("ALL.p.value", "corrsnps", function(x, value) {
                 x@ALL.p.value <- value
                 x
         })
setGeneric("AFR.p.value", function(x) standardGeneric("AFR.p.value"))
setGeneric("AFR.p.value<-", function(x, value) standardGeneric("AFR.p.value<-"))
setMethod("AFR.p.value", "corrsnps", function(x) x@AFR.p.value)
setReplaceMethod("AFR.p.value", "corrsnps", function(x, value) {
                 x@AFR.p.value <- value
                 x
         })
setGeneric("AMR.p.value", function(x) standardGeneric("AMR.p.value"))
setGeneric("AMR.p.value<-", function(x, value) standardGeneric("AMR.p.value<-"))
setMethod("AMR.p.value", "corrsnps", function(x) x@AMR.p.value)
setReplaceMethod("AMR.p.value", "corrsnps", function(x, value) {
                 x@AMR.p.value <- value
                 x
         })
setGeneric("ASN.p.value", function(x) standardGeneric("ASN.p.value"))
setGeneric("ASN.p.value<-", function(x, value) standardGeneric("ASN.p.value<-"))
setMethod("ASN.p.value", "corrsnps", function(x) x@ASN.p.value)
setReplaceMethod("ASN.p.value", "corrsnps", function(x, value) {
                 x@ASN.p.value <- value
                 x
         })
setGeneric("EUR.p.value", function(x) standardGeneric("EUR.p.value"))
setGeneric("EUR.p.value<-", function(x, value) standardGeneric("EUR.p.value<-"))
setMethod("EUR.p.value", "corrsnps", function(x) x@EUR.p.value)
setReplaceMethod("EUR.p.value", "corrsnps", function(x, value) {
                 x@EUR.p.value <- value
                 x
         })
setGeneric("ALL.geno.overlapping.snps", function(x) standardGeneric("ALL.geno.overlapping.snps"))
setGeneric("ALL.geno.overlapping.snps<-", function(x, value) standardGeneric("ALL.geno.overlapping.snps<-"))
setMethod("ALL.geno.overlapping.snps", "corrsnps", function(x) x@ALL.geno.overlapping.snps)
setReplaceMethod("ALL.geno.overlapping.snps", "corrsnps", function(x, value) {
                 x@ALL.geno.overlapping.snps <- value
                 x
         })
setGeneric("AFR.geno.overlapping.snps", function(x) standardGeneric("AFR.geno.overlapping.snps"))
setGeneric("AFR.geno.overlapping.snps<-", function(x, value) standardGeneric("AFR.geno.overlapping.snps<-"))
setMethod("AFR.geno.overlapping.snps", "corrsnps", function(x) x@AFR.geno.overlapping.snps)
setReplaceMethod("AFR.geno.overlapping.snps", "corrsnps", function(x, value) {
                 x@AFR.geno.overlapping.snps <- value
                 x
         })
setGeneric("AMR.geno.overlapping.snps", function(x) standardGeneric("AMR.geno.overlapping.snps"))
setGeneric("AMR.geno.overlapping.snps<-", function(x, value) standardGeneric("AMR.geno.overlapping.snps<-"))
setMethod("AMR.geno.overlapping.snps", "corrsnps", function(x) x@AMR.geno.overlapping.snps)
setReplaceMethod("AMR.geno.overlapping.snps", "corrsnps", function(x, value) {
                 x@AMR.geno.overlapping.snps <- value
                 x
         })
setGeneric("ASN.geno.overlapping.snps", function(x) standardGeneric("ASN.geno.overlapping.snps"))
setGeneric("ASN.geno.overlapping.snps<-", function(x, value) standardGeneric("ASN.geno.overlapping.snps<-"))
setMethod("ASN.geno.overlapping.snps", "corrsnps", function(x) x@ASN.geno.overlapping.snps)
setReplaceMethod("ASN.geno.overlapping.snps", "corrsnps", function(x, value) {
                 x@ASN.geno.overlapping.snps <- value
                 x
         })
setGeneric("EUR.geno.overlapping.snps", function(x) standardGeneric("EUR.geno.overlapping.snps"))
setGeneric("EUR.geno.overlapping.snps<-", function(x, value) standardGeneric("EUR.geno.overlapping.snps<-"))
setMethod("EUR.geno.overlapping.snps", "corrsnps", function(x) x@EUR.geno.overlapping.snps)
setReplaceMethod("EUR.geno.overlapping.snps", "corrsnps", function(x, value) {
                 x@EUR.geno.overlapping.snps <- value
                 x
         })

setMethod("summary", "tagsnp", function(object) {
          data.frame(chromosome=chromosome(object), position=position(object),
                     snpid=snpid(object), population=population(object),
                     ref.allele=ref.allele(object), alt.allele=alt.allele(object),
                     overlapping.features=lapply(overlapping.features(object), function(x) length(x)),
                     genotype=genotype(object), R.squared.corrsnps=R.squared.corrsnps(object),
                     D.prime.corrsnps=D.prime.corrsnps(object))
         }
)
