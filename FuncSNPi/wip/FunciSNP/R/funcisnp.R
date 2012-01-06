library("Rsamtools")
library("rtracklayer")
library("GGtools")
library("parallel")
source("/home/scoetzee/scripts/Coetzee_Seq_Analysis/FuncSNPi/wip/FunciSNP/R/classes.R")

ReadRegionsFile <- function(regions.file) {
  # Reads a tab seperated regions file in the form
  # chr:start-end snp_name    ethnicity   chr#
  #
  # returns the variables for snp.range, snp.name, and snp.ethno
  snp.regions <- read.table(regions.file)
  snp.region.split <- unlist(strsplit(as.vector(snp.regions[,1]), ":"))

  snp.chromosome <- grep("^\\d{0,1}[0-9X-Y]$", snp.region.split, value = TRUE,
                         perl = TRUE)

  m <- regexpr("^\\d{1,}-", snp.region.split)
  snp.region.start <- substr(snp.region.split,
                             m,
                             m + attr(m, "match.length") - 1)
  snp.region.start <- as.vector(na.omit(as.numeric(gsub("-$",
                                                        "", 
                                                        snp.region.start))))

  m <- regexpr("-\\d{1,}$", snp.region.split)
  snp.region.end <- substr(snp.region.split,
                           m,
                           m + attr(m, "match.length") - 1)
  snp.region.end <- as.vector(na.omit(as.numeric(gsub("^-",
                                                      "",
                                                      snp.region.end))))

  snp.name <- as.character(snp.regions[, 2])
  snp.ethno <- as.character(snp.regions[, 3])
  snp.regions <- data.frame(snp.chromosome, 
                            snp.region.start, 
                            snp.region.end, 
                            snp.name, 
                            snp.ethno, 
                            stringsAsFactors = FALSE)
  return(snp.regions)
}

ServerCheck <- function(primary.server, verbose=TRUE) {
  ncbi <- "ftp://ftp-trace.ncbi.nih.gov/1000genomes/"
  ebi <- "ftp://ftp.1000genomes.ebi.ac.uk/vol1/"
  test.file <- "ftp/release/20110521/phase1_integrated_calls.20101123.ALL.panel"
  if(primary.server == "ebi"){
    primary.server <- ebi
    secondary.server <- ncbi
    primary.server.name <- "ebi"
    secondary.server.name <- "ncbi"
  } else {
    primary.server <- ncbi
    secondary.server <- ebi
    primary.server.name <- "ncbi"
    secondary.server.name <- "ebi"
  }
  if(verbose) {
    message("trying ", primary.server.name, " as 1000 genomes server\n")
  }
  server.up <- url(paste(primary.server, test.file, sep=""), open = "rt")
  server.error <- inherits(server.up, "try-error")
  if(server.error) {
    warning(primary.server.name, " failed \ntrying ",  secondary.server.name,
            " as 1000 genomes server")
    rm(server.error)
    server.up <- url(paste(secondary.server, test.file, sep=""), open = "rt")
    server.error <- inherits(server.up, "try-error")
    if(server.error) {
      warning(secondary.server.name, " failed")
      close(server.up)
      stop("Neither EBI nor NCBI mirrors were found")
    } else {
      if(verbose) {
        message("OK using ", secondary.server.name, " : ", secondary.server)
      }
      close(server.up)
      return(secondary.server)
    }
  } else {
    if(verbose) {
      message("OK using ", primary.server.name, " : ", primary.server)
    }
    close(server.up)
    return(primary.server)
  }
}

CreateTagSNP <- function(tag.snp.name) {
  tag.id <- grep("^rs",
                 unlist(strsplit(tag.snp.name, ":")),
                 value=T,
                 perl = TRUE)
  tag.population <- grep("[A-Za-z]$", 
                         unlist(strsplit(tag.snp.name, ":")), 
                         value = TRUE, 
                         perl = TRUE)
  assign(tag.snp.name,
         new("tagsnp", 
             snpid=as.character(tag.id), 
             population=toupper(as.character(tag.population))))
}
CreateCorrelatedSNPs<- function(tag.snp.name, snp.list, primary.server,
                                snp.region, bio.features.file, populations,
                                verbose, method.p, reduce.by.features) {
#  tag.snp.name <<- tag.snp.name
#  snp.list <<- snp.list
#  primary.server <<- primary.server
#  snp.region <<- snp.region
#  bio.features.file <<- bio.features.file
#  populations <<- populations
#  stop("stopped already")
  tag.snp.complete <- PullInVariants(tag.snp.name, snp.list, primary.server,
                                     snp.region, populations, verbose)
  overlapping.features(correlated.snps(tag.snp.complete)) <-
    unlist(GRangesList(lapply(as.character(bio.features.file),
                              FilterByFeatures,
                              tag.snp.name=tag.snp.name,
                              tag.snp.complete=tag.snp.complete,
                              verbose=verbose)))
  if(reduce.by.features) {
    overlapping.snps <- unique(names(overlapping.features(correlated.snps(tag.snp.complete))))
    ALL.geno.overlapping.snps(correlated.snps(tag.snp.complete)) <- ALL.geno(correlated.snps(tag.snp.complete))[, c(overlapping.snps, snpid(tag.snp.complete))]
    AFR.geno.overlapping.snps(correlated.snps(tag.snp.complete)) <- AFR.geno(correlated.snps(tag.snp.complete))[, c(overlapping.snps, snpid(tag.snp.complete))]
    AMR.geno.overlapping.snps(correlated.snps(tag.snp.complete)) <- AMR.geno(correlated.snps(tag.snp.complete))[, c(overlapping.snps, snpid(tag.snp.complete))]
    ASN.geno.overlapping.snps(correlated.snps(tag.snp.complete)) <- ASN.geno(correlated.snps(tag.snp.complete))[, c(overlapping.snps, snpid(tag.snp.complete))]
    EUR.geno.overlapping.snps(correlated.snps(tag.snp.complete)) <- EUR.geno(correlated.snps(tag.snp.complete))[, c(overlapping.snps, snpid(tag.snp.complete))]
    tag.snp.complete <<- tag.snp.complete
    tag.snp.complete <- FilteredLDTesting(tag.snp.complete, verbose)
    if(verbose) message("Calculating p-value for ", tag.snp.name)
    ALL.p.value(correlated.snps(tag.snp.complete)) <-
      ChiSquaredPvalue(ALL.geno.overlapping.snps(correlated.snps(tag.snp.complete)),
                       snpid(tag.snp.complete), method.p)
    AFR.p.value(correlated.snps(tag.snp.complete)) <-
      ChiSquaredPvalue(AFR.geno.overlapping.snps(correlated.snps(tag.snp.complete)),
                       snpid(tag.snp.complete), method.p)
    AMR.p.value(correlated.snps(tag.snp.complete)) <-
      ChiSquaredPvalue(AMR.geno.overlapping.snps(correlated.snps(tag.snp.complete)),
                       snpid(tag.snp.complete), method.p)
    ASN.p.value(correlated.snps(tag.snp.complete)) <-
      ChiSquaredPvalue(ASN.geno.overlapping.snps(correlated.snps(tag.snp.complete)),
                       snpid(tag.snp.complete), method.p)
    EUR.p.value(correlated.snps(tag.snp.complete)) <-
      ChiSquaredPvalue(EUR.geno.overlapping.snps(correlated.snps(tag.snp.complete)),
                       snpid(tag.snp.complete), method.p)
  } else {
    tag.snp.complete <- LDTesting(tag.snp.complete, verbose)
    if(verbose) message("Calculating p-value for ", tag.snp.name)
    ALL.p.value(correlated.snps(tag.snp.complete)) <-
      ChiSquaredPvalue(ALL.geno(correlated.snps(tag.snp.complete)),
                       snpid(tag.snp.complete), method.p)
    AFR.p.value(correlated.snps(tag.snp.complete)) <-
      ChiSquaredPvalue(AFR.geno(correlated.snps(tag.snp.complete)),
                       snpid(tag.snp.complete), method.p)
    AMR.p.value(correlated.snps(tag.snp.complete)) <-
      ChiSquaredPvalue(AMR.geno(correlated.snps(tag.snp.complete)),
                       snpid(tag.snp.complete), method.p)
    ASN.p.value(correlated.snps(tag.snp.complete)) <-
      ChiSquaredPvalue(ASN.geno(correlated.snps(tag.snp.complete)),
                       snpid(tag.snp.complete), method.p)
    EUR.p.value(correlated.snps(tag.snp.complete)) <-
      ChiSquaredPvalue(EUR.geno(correlated.snps(tag.snp.complete)),
                       snpid(tag.snp.complete), method.p)
  }
  message("Tag SNP ", snpid(tag.snp.complete), " has completed")
  return(tag.snp.complete)
  gc(verbose = FALSE)
}

CreatePopulations <- function(primary.server="ncbi") {
  onek.genome.server <- ServerCheck(primary.server, verbose = FALSE)
  manifest <- read.delim(paste(onek.genome.server,
                               "ftp/release/20110521/phase1_integrated_calls.20101123.ALL.panel", sep=""),
                         sep="\t", header = FALSE)
  for(i in c("AFR", "ASN", "EUR", "AMR", "ALL")) {
    ifelse(i != "ALL",
           assign(i, as.character(subset(manifest, V3==i)[,1])),
           assign(i, as.character(manifest[, 1])))
  }
  populations <- list(AFR=AFR, AMR=AMR, ASN=ASN, EUR=EUR, ALL=ALL)
  return(populations)
}

FunciSNP <- function(snp.regions.file, bio.features.loc = NULL,
                     primary.server="ncbi", par.threads=detectCores()/2,
                     verbose = par.threads < 2, method.p = "BH", reduce.by.features = TRUE ) {
  message("
          ####################################
          ##                                ##
          ##      Welcome to FunciSNP       ##
          ##                                ##
          ####################################
          ::args used::", "\n",
          "          verbose:                          ",
          verbose, "\n",
          "          cores in use:                     ",
          par.threads, "\n",
          "          snp.regions.file:                 ",
          as.character(snp.regions.file), "\n",
          "          p-value adjustment by:            ",
          method.p)
  if(identical(bio.features.loc, NULL)) {
    message("          Bio Features:                     no biofeatures selected")
  } else {
    message("          Bio Features:                     ",
            gsub(".bed$", ", ", list.files(bio.features.loc,
                                           pattern="*.bed$",
                                           full.names=FALSE)))
  }
  if(primary.server == "ebi" ||  primary.server == "ncbi") {
    message("          Primary Server:                   ", primary.server)
  } else {
    stop("please select \"ebi\" or \"ncbi\" as your primary server")
  }
  snp.region <- ReadRegionsFile(snp.regions.file)
  message("          Number of TagSNPs Interrogated:   ", nrow(snp.region))
  options(mc.cores=par.threads)
  populations <- CreatePopulations("ncbi")
  if(identical(bio.features.loc, NULL)) {
    bio.features.file <- NULL
  } else {
    bio.features.file <- list.files(bio.features.loc, pattern="*.bed$",
                                    full.names = TRUE)
  }
  tag.snp.names <- paste(snp.region$snp.name, ":", snp.region$snp.ethno, sep="")

  snp.list <- lapply(tag.snp.names, CreateTagSNP)
  names(snp.list) <- tag.snp.names
  snp.list <- mclapply(tag.snp.names, CreateCorrelatedSNPs,
                       snp.list=snp.list,
                       primary.server=primary.server,
                       snp.region=snp.region,
                       bio.features.file=bio.features.file,
                       populations=populations,
                       method.p=method.p,
                       verbose=verbose,
                       reduce.by.features=reduce.by.features) 
  names(snp.list) <- tag.snp.names
  return(snp.list)
}

PullInVariants <- function(tag.snp.name, snp.list, primary.server, snp.region,
                           populations, verbose = TRUE) {
#  tag.snp.name <<- tag.snp.name
#  snp.list <<- snp.list
#  primary.server <<- primary.server
#  snp.region <<- snp.region
#  populations <<- populations
#  verbose <<- verbose
#  stop("no error, just stopping")
  tag.id <- snpid(snp.list[[tag.snp.name]])
  window.size <- prettyNum(as.numeric(snp.region$snp.region.end[snp.region$snp.name == tag.id]) - as.numeric(snp.region$snp.region.start[snp.region$snp.name == tag.id]),
                           big.mark=",", 
                           scientific = FALSE)
  message("\nloading ", tag.id, " window size: ", window.size, " bp")
  onek.genome.server <- ServerCheck(primary.server, verbose = FALSE)
  variants.reference <- 
    paste(onek.genome.server, 
          "/ftp/release/20101123/interim_phase1_release/ALL.chr", 
          snp.region$snp.chromosome[snp.region$snp.name == tag.id][1], 
          ".phase1.projectConsensus.genotypes.vcf.gz", sep="")

  param <- GRanges(seqnames=snp.region$snp.chromosome[snp.region$snp.name == tag.id][1],
                   IRanges(snp.region$snp.region.start[snp.region$snp.name == tag.id][1],
                           snp.region$snp.region.end[snp.region$snp.name == tag.id][1]))

  kgeno <- TabixFile(variants.reference)
  tabix.header <- 
    strsplit(headerTabix(variants.reference)$header, split="\t")[[6]]
  if(verbose) message("scanning tabix file for ", tag.id)
  tabix.file <- 
    strsplit(scanTabix(kgeno, param=param)[[1]], split="\t")
  if(verbose) message("creating snpMatrix for ", tag.id)
  snps.geno <- vcf2sm(kgeno, gr=param, nmetacol=9L)
  snps.support <- 
    t(as.data.frame(tabix.file, stringsAsFactors = FALSE)[1:5, ])
  colnames(snps.support) <- tabix.header[1:5]
  row.names(snps.support) <- NULL
  tag.support <- snps.support[(snps.support[, 3] %in% tag.id)]
  snps.support <- snps.support[!(snps.support[,3] %in% tag.id), ]


  corr.snps <- new("corrsnps",
                   chromosome=as.integer(snps.support[, 1]),
                   position=as.integer(snps.support[, 2]),
                   snpid=ifelse(snps.support[, 3] == ".", 
                                paste("chr", snps.support[, 1], ":",
                                      snps.support[, 2], sep=""),
                                as.character(snps.support[, 3])),
                   ref.allele=as.character(snps.support[, 4]),
                   alt.allele=as.character(snps.support[, 5]),
                   ALL.geno=snps.geno,
                   AFR.geno=snps.geno[row.names(snps.geno) %in% populations$AFR],
                   AMR.geno=snps.geno[row.names(snps.geno) %in% populations$AMR],
                   ASN.geno=snps.geno[row.names(snps.geno) %in% populations$ASN],
                   EUR.geno=snps.geno[row.names(snps.geno) %in% populations$EUR]
                   )

  correlated.snps(snp.list[[tag.snp.name]]) <- corr.snps
  position(snp.list[[tag.snp.name]]) <- as.integer(tag.support[2])
  ref.allele(snp.list[[tag.snp.name]]) <- as.character(tag.support[4])
  alt.allele(snp.list[[tag.snp.name]]) <- as.character(tag.support[5])
  tag.ethno <- substr(tag.snp.name, nchar(tag.snp.name) - 2, nchar(tag.snp.name))
  tag.population <- getElement(populations, tag.ethno)
  ##questionable line break
  genotype(snp.list[[tag.snp.name]]) <- snps.geno[, tag.id][row.names(snps.geno)
                                                            %in% tag.population]

  return(snp.list[[tag.snp.name]])
}

FilterByFeatures <- function(bio.features.file = NULL, tag.snp.name,
                             tag.snp.complete, verbose = TRUE) {
#  bio.features.file <<- bio.features.file
#  tag.snp.name <<- tag.snp.name
#  tag.snp.complete <<- tag.snp.complete
  if(verbose) message("Filtering ", snpid(tag.snp.complete), " against ",
                      c(substr(grep(".bed$",
                                    unlist(strsplit(bio.features.file, "/")), 
                                    value=TRUE), 
                               1, nchar(grep(".bed$",
                                             unlist(strsplit(bio.features.file,
                                                             "/")),
                                             value=TRUE))-4)))
  close.snp.ranges <- 
    GRanges(seqnames=paste(
                           "chr",
                           as.character(chromosome(correlated.snps(tag.snp.complete))), 
                           sep=""), 
            ranges=(IRanges(
                            start=as.integer(position(correlated.snps(tag.snp.complete))), 
                            width=1)))
  names(close.snp.ranges) <- snpid(correlated.snps(tag.snp.complete))

  if(identical(bio.features.file, NULL)) {
    snps.included <- close.snp.ranges
  } else {
    bio.features.file.interval <- import(bio.features.file, asRangedData = FALSE)
    elementMetadata(bio.features.file.interval) <- NULL
    elementMetadata(bio.features.file.interval)[, "feature"] <-
      c(substr(grep(".bed$",
                    unlist(strsplit(bio.features.file, "/")), value=TRUE),
               1,
               nchar(grep(".bed$",
                          unlist(strsplit(bio.features.file, "/")),
                          value=TRUE))-4))
    snps.included <-
      !is.na(match(close.snp.ranges, bio.features.file.interval) > 0)
    ranges.included <-
      !is.na(match(bio.features.file.interval, close.snp.ranges) > 0)
    peaks.with.overlapping.snps <-
      subset(bio.features.file.interval, ranges.included)
    snps.overlapping.peaks <-
      subset(close.snp.ranges, snps.included)
    overlaps <-
      findOverlaps(peaks.with.overlapping.snps,
                   snps.overlapping.peaks,
                   select="all")
    snps.included <-
      lapply(queryHits(overlaps), function(x) peaks.with.overlapping.snps[x])

    ##write method to call a specific snp that occurs under two peaks 
    ##use the following technique test00[names(test00)=="rs123456"]
  }

  if(length(snps.included) >= 1) {
    snps.included <- unlist(GRangesList(snps.included))
    names(snps.included) <- names(snps.overlapping.peaks)
    if(identical(overlapping.features(correlated.snps(tag.snp.complete)),
                 GRanges())) {
      overlapping.features(correlated.snps(tag.snp.complete)) <-
        snps.included
    } else {
      overlapping.features(correlated.snps(tag.snp.complete)) <-
        c(overlapping.features(correlated.snps(tag.snp.complete)),
          snps.included)
    }
    return(overlapping.features(correlated.snps(tag.snp.complete)))
  } else {
    feature.name <-
      c(substr(grep(".bed$", unlist(strsplit(bio.features.file, "/")),
                    value=TRUE),
               1,
               nchar(grep(".bed$", unlist(strsplit(bio.features.file, "/")),
                          value=TRUE))-4))
    ## need to put the following in log files
    warning("There is no overlap for: \n",
            "\tTag SNP: \t\t", snpid(tag.snp.complete), "\n", 
            "\tbiofeature: \t\t", feature.name, "\n", 
            "\tpopulation: \t\t", population(tag.snp.complete));
    return(overlapping.features(correlated.snps(tag.snp.complete)))
  }
  message(snpid(tag.snp.complete), " has ", sum(unique(names(overlapping.features(correlated.snps(tag.snp.complete))))), " nearby SNPs overlapping with feature ", feature.name)
}

LDTesting <- function(tag.snp.complete, verbose = TRUE) {
  if(verbose) message("Calculating R^2 and D' for ", snpid(tag.snp.complete))
  snp.name.chosen <- snpid(tag.snp.complete)
  #  snp.name.chosen <- grep("^rs", unlist(strsplit(snp.list.name, ":")), value=T)

  ALL.R.squared(correlated.snps(tag.snp.complete)) <-
    ld(ALL.geno(correlated.snps(tag.snp.complete))[, snp.name.chosen],
       ##questionable line break
       ALL.geno(correlated.snps(tag.snp.complete))[, !(colnames(ALL.geno(correlated.snps(tag.snp.complete))))
                                                   %in% snp.name.chosen], stats="R.squared", depth=100000)
  AFR.R.squared(correlated.snps(tag.snp.complete)) <-
    ld(AFR.geno(correlated.snps(tag.snp.complete))[, snp.name.chosen],
       AFR.geno(correlated.snps(tag.snp.complete))[, !(colnames(AFR.geno(correlated.snps(tag.snp.complete))))
                                                   %in% snp.name.chosen], stats="R.squared", depth=100000)
  AMR.R.squared(correlated.snps(tag.snp.complete)) <-
    ld(AMR.geno(correlated.snps(tag.snp.complete))[, snp.name.chosen],
       AMR.geno(correlated.snps(tag.snp.complete))[, !(colnames(AMR.geno(correlated.snps(tag.snp.complete))))
                                                   %in% snp.name.chosen], stats="R.squared", depth=100000)
  ASN.R.squared(correlated.snps(tag.snp.complete)) <-
    ld(ASN.geno(correlated.snps(tag.snp.complete))[, snp.name.chosen],
       ASN.geno(correlated.snps(tag.snp.complete))[, !(colnames(ASN.geno(correlated.snps(tag.snp.complete))))
                                                   %in% snp.name.chosen], stats="R.squared", depth=100000)
  EUR.R.squared(correlated.snps(tag.snp.complete)) <-
    ld(EUR.geno(correlated.snps(tag.snp.complete))[, snp.name.chosen],
       EUR.geno(correlated.snps(tag.snp.complete))[, !(colnames(EUR.geno(correlated.snps(tag.snp.complete))))
                                                   %in% snp.name.chosen], stats="R.squared", depth=100000)


  ALL.D.prime(correlated.snps(tag.snp.complete)) <-
    ld(ALL.geno(correlated.snps(tag.snp.complete))[, snp.name.chosen],
       ALL.geno(correlated.snps(tag.snp.complete))[, !(colnames(ALL.geno(correlated.snps(tag.snp.complete))))
                                                   %in% snp.name.chosen], stats="D.prime", depth=100000)
  AFR.D.prime(correlated.snps(tag.snp.complete)) <-
    ld(AFR.geno(correlated.snps(tag.snp.complete))[, snp.name.chosen],
       AFR.geno(correlated.snps(tag.snp.complete))[, !(colnames(AFR.geno(correlated.snps(tag.snp.complete))))
                                                   %in% snp.name.chosen], stats="D.prime", depth=100000)
  AMR.D.prime(correlated.snps(tag.snp.complete)) <-
    ld(AMR.geno(correlated.snps(tag.snp.complete))[, snp.name.chosen],
       AMR.geno(correlated.snps(tag.snp.complete))[, !(colnames(AMR.geno(correlated.snps(tag.snp.complete))))
                                                   %in% snp.name.chosen], stats="D.prime", depth=100000)
  ASN.D.prime(correlated.snps(tag.snp.complete)) <-
    ld(ASN.geno(correlated.snps(tag.snp.complete))[, snp.name.chosen],
       ASN.geno(correlated.snps(tag.snp.complete))[, !(colnames(ASN.geno(correlated.snps(tag.snp.complete))))
                                                   %in% snp.name.chosen], stats="D.prime", depth=100000)
  EUR.D.prime(correlated.snps(tag.snp.complete)) <-
    ld(EUR.geno(correlated.snps(tag.snp.complete))[, snp.name.chosen],
       EUR.geno(correlated.snps(tag.snp.complete))[, !(colnames(EUR.geno(correlated.snps(tag.snp.complete))))
                                                   %in% snp.name.chosen], stats="D.prime", depth=100000)



  R.squared.corrsnps(tag.snp.complete) <-
    ld(eval(parse(text=(paste(population(tag.snp.complete),
                              ".geno(correlated.snps(tag.snp.complete))",
                              sep="")))), stats="R.squared", depth=2000)
  D.prime.corrsnps(tag.snp.complete) <-
    ld(eval(parse(text=(paste(population(tag.snp.complete),
                              ".geno(correlated.snps(tag.snp.complete))",
                              sep="")))), stats="D.prime", depth=2000)

  return(tag.snp.complete)
}

FilteredLDTesting <- function(tag.snp.complete, verbose = TRUE) {
  if(ncol(eval(parse(text=(paste(population(tag.snp.complete),
                              ".geno.overlapping.snps(correlated.snps(tag.snp.complete))",
                              sep=""))))) > 1) {
  if(verbose) message("Calculating R^2 and D' for ", snpid(tag.snp.complete))
  snp.name.chosen <- snpid(tag.snp.complete)
  #  snp.name.chosen <- grep("^rs", unlist(strsplit(snp.list.name, ":")), value=T)

  R.squared.corrsnps(tag.snp.complete) <-
    ld(eval(parse(text=(paste(population(tag.snp.complete),
                              ".geno.overlapping.snps(correlated.snps(tag.snp.complete))",
                              sep="")))), stats="R.squared", depth=2000)
  D.prime.corrsnps(tag.snp.complete) <-
    ld(eval(parse(text=(paste(population(tag.snp.complete),
                              ".geno.overlapping.snps(correlated.snps(tag.snp.complete))",
                              sep="")))), stats="D.prime", depth=2000)
  return(tag.snp.complete)
  } else {
    return(tag.snp.complete)
  }
}

ChiSquaredPvalue <- function(tag.snp.complete, tag.snp.id, method.p) {
#  tag.snp.complete <<- tag.snp.complete
#  tag.snp.id <<- tag.snp.id
#  stop("stopped")
  snp.list <- lapply(colnames(tag.snp.complete),
                     function(x) as(tag.snp.complete[, x], "character"))
  names(snp.list) <- colnames(tag.snp.complete)
  genotype.table.snps <- lapply(snp.list,
                                function(x) table(unlist(snp.list[tag.snp.id]),
                                                  unlist(x),
                                                  dnn = (c("tag", "corr"))))

  pre.genotype.table.snps <<- genotype.table.snps
  genotype.table.snps <-
    lapply(pre.genotype.table.snps, function(x) {
           #   message("\n", x)
           x <- as.matrix(x)
           if((dim(x) > 3)[1]) {
             x <- x[1:3, 1:ncol(x)]
           }
           if((dim(x) > 3)[2]) {
             x <- x[1:nrow(x), 1:3]
           }
           xnames <- dimnames(x)
           #   message(x)
           #   message(xnames)

           #message(xnames)
           ##rows
           #message(1)
           if(dim(x)[1] == 1 && length(grep("A/A",dimnames(x)[[1]])) >= 1) {
             x <- rbind("A/A"=x[1, ], "A/B"=0, "B/B"=0)
             xnames$tag <- c("A/A", "A/B", "B/B")
             dimnames(x) <- xnames
           }
           #message(1)
           if(dim(x)[1] == 1 && length(grep("A/B", dimnames(x)[[1]])) >= 1) {
             x <- rbind("A/A"=0, "A/B"=x[1, ], "B/B"=0)
             xnames$tag <- c("A/A", "A/B", "B/B")
             dimnames(x) <- xnames
           }
           #message(1)
           if(dim(x)[1] == 1 && length(grep("B/B", dimnames(x)[[1]])) >= 1) {
             x <- rbind("A/A"=0, "A/B"=0, "B/B"=x[1, ])
             xnames$tag <- c("A/A", "A/B", "B/B")
             dimnames(x) <- xnames
           }
           #message(1)

           if(dim(x)[1] == 2 && length(grep("A/A", dimnames(x)[[1]])) < 1) {
             x <- rbind("A/A"=0, as.matrix(x[1:2, ]))
             xnames$tag <- c("A/A", "A/B", "B/B")
             dimnames(x) <- xnames
           }
           #message(1)
           if(dim(x)[1] == 2 && length(grep("A/B", dimnames(x)[[1]])) < 1) {
             x <- rbind("A/A"=as.matrix(x[1, ]), "A/B"=0, "B/B"=as.matrix(x[2, ]))
             xnames$tag <- c("A/A", "A/B", "B/B")
             dimnames(x) <- xnames
           }
           #message(1)
           if(dim(x)[1] == 2 && length(grep("B/B", dimnames(x)[[1]])) < 1) {
             x <- rbind(as.matrix(x[1:2, ]), "B/B"=0)
             xnames$tag <- c("A/A", "A/B", "B/B")
             dimnames(x) <- xnames
           }
           #message(1)
           ##columns                                
           if(dim(x)[2] == 1 && length(grep("A/A", dimnames(x)[[2]])) >= 1) {
             x <- cbind("A/A"=x[, 1], "A/B"=0, "B/B"=0)
             xnames$corr <- c("A/A", "A/B", "B/B")
             dimnames(x) <- xnames
           }
           #message(1)
           if(dim(x)[2] == 1 && length(grep("A/B", dimnames(x)[[2]])) >= 1) {
             x <- cbind("A/A"=0, "A/B"=x[, 1], "B/B"=0)
             xnames$corr <- c("A/A", "A/B", "B/B")
             dimnames(x) <- xnames
           }
           #message(1)
           if(dim(x)[2] == 1 && length(grep("B/B", dimnames(x)[[2]])) >= 1) {
             x <- cbind("A/A"=0, "A/B"=0, "B/B"=x[, 1])
             xnames$corr <- c("A/A", "A/B", "B/B")
             dimnames(x) <- xnames
           }

           #message(1)
           if(dim(x)[2] == 2 && length(grep("A/A", dimnames(x)[[2]])) < 1) {
             x <- cbind("A/A"=0, x[, 1:2])
             xnames$corr <- c("A/A", "A/B", "B/B")
             dimnames(x) <- xnames
           }
           #message(1)
           if(dim(x)[2] == 2 && length(grep("A/B", dimnames(x)[[2]])) < 1) {
             x <- cbind("A/A"=x[, 1], "A/B"=0, "B/B"=x[, 2])
             xnames$corr <- c("A/A", "A/B", "B/B")
             dimnames(x) <- xnames
           }
           #message(1)
           if(dim(x)[2] == 2 && length(grep("B/B", dimnames(x)[[2]])) < 1) {
             x <- cbind(x[, 1:2], "B/B"=0)
             xnames$corr <- c("A/A", "A/B", "B/B")
             dimnames(x) <- xnames
           }
           #message(1)
           return(x)
                                }
  )
  OR <- ld(tag.snp.complete[, tag.snp.id],
           tag.snp.complete[, !colnames(tag.snp.complete) %in% tag.snp.id],
           stats="OR", depth=10000)
  p.e <- lapply(colnames(tag.snp.complete), function(x, tag.snp.id) {
                if(!(identical(x, tag.snp.id))) {
                  genotype.table.snps[[x]][2, 2] * OR[, x]/(1 + OR[, x])
  }}, tag.snp.id)
  q.e <- lapply(colnames(tag.snp.complete), function(x, tag.snp.id) {
                if(!(identical(x, tag.snp.id))) {
                  genotype.table.snps[[x]][2, 2] * 1/(1 + OR[, x])
                  }}, tag.snp.id)
  names(p.e) <- names(snp.list)
  names(q.e) <- names(snp.list)
  q.e <- unlist(q.e)
  p.e <- unlist(p.e)
#  q.e <<- q.e
#  p.e <<- p.e
#  snp.list <<- snp.list
#  genotype.table.snps <<- genotype.table.snps
  haplotype.table.snps <-
    lapply(names(snp.list), function(x, genotype.table.snps, p.e, q.e) {
           A <- sum(2 * genotype.table.snps[[x]][1, 1],
                    genotype.table.snps[[x]][1, 2],
                    genotype.table.snps[[x]][2, 1],
                    p.e[x], na.rm = TRUE)
           B <- sum(2 * genotype.table.snps[[x]][1, 3],
                    genotype.table.snps[[x]][1, 2],
                    genotype.table.snps[[x]][2, 3],
                    q.e[x], na.rm = TRUE)
           C <- sum(2 * genotype.table.snps[[x]][3, 1],
                    genotype.table.snps[[x]][3, 2],
                    genotype.table.snps[[x]][2, 1],
                    q.e[x], na.rm = TRUE)
           D <- sum(2 * genotype.table.snps[[x]][3, 3],
                    genotype.table.snps[[x]][3, 2],
                    genotype.table.snps[[x]][2, 3],
                    p.e[x], na.rm = TRUE)
           return(matrix(c(A, B, C, D), ncol = 2, byrow = TRUE))
                  },
                  genotype.table.snps, p.e, q.e)
  names(haplotype.table.snps) <- names(snp.list)
  raw.p <- lapply(haplotype.table.snps, function(x) {
                  if(sum(x, na.rm = TRUE) > 0) {
                    return(suppressWarnings(fisher.test(x)$p.value)) ## done b/c many tables have cells with < 5 values
                  } else {
                    return(NA)
                  }})
  adj.p <- p.adjust(raw.p, method=method.p)
  names(raw.p) <- names(snp.list)
  names(adj.p) <- names(snp.list)
  return(list(raw.p=raw.p, adj.p=adj.p))
}

