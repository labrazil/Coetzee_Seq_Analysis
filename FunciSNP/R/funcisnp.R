ReadRegionsFile <- function(regions.file, search.window=200000) {
  # Reads a tab seperated regions file in the form
  # chr:loc snp_name    ethnicity
  # 8:130685457 rs4295627 EUR
  # returns the variables for snp.range, snp.name, and snp.ethno
  snp.regions <- read.table(regions.file)
  snp.region.split <- unlist(strsplit(as.vector(snp.regions[,1]), ":"))

  snp.chromosome <- as.character(sapply(strsplit(as.vector(snp.regions[,1]), ":"), function(x) x[1]))
  snp.loc <- as.numeric(sapply(strsplit(as.vector(snp.regions[,1]), ":"), function(x) x[2]))


  snp.region.start <- round(snp.loc - search.window/2)
  snp.region.end <- round(snp.loc + search.window/2)

  snp.name <- as.character(snp.regions[, 2])
  snp.ethno <- as.character(snp.regions[, 3])
  snp.regions <- data.frame(snp.chromosome,
                            snp.loc,
                            snp.region.start,
                            snp.region.end,
                            snp.name,
                            snp.ethno,
                            stringsAsFactors = FALSE)
  snp.regions$snp.ethno <- toupper(snp.regions$snp.ethno)
  tag.snp.all <- subset(snp.regions, snp.ethno == "ALL")
  snp.regions <- subset(snp.regions, snp.ethno != "ALL")
  if(dim(tag.snp.all)[1] > 0){
    for(i in 1:dim(tag.snp.all)[1]){
      x <- tag.snp.all[i, ]
      tag.snp.afr <- x
      tag.snp.afr["snp.ethno"] <- "AFR"
      tag.snp.amr <- x
      tag.snp.amr["snp.ethno"] <- "AMR"
      tag.snp.asn <- x
      tag.snp.asn["snp.ethno"] <- "ASN"
      tag.snp.eur <- x
      tag.snp.eur["snp.ethno"] <- "EUR"
      tag.snp.each <- rbind(x,
                            tag.snp.afr,
                            tag.snp.amr,
                            tag.snp.asn,
                            tag.snp.eur)
      snp.regions <- rbind(snp.regions, tag.snp.each)
    }
  }
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
  server.up <- try(url(paste(primary.server, test.file, sep=""), open = "rt"), 
                   silent=T)
  server.error <- inherits(server.up, "try-error")
  if(server.error) {
    if(verbose) warning(primary.server.name, " failed \ntrying ",  
                        secondary.server.name,
                        " as 1000 genomes server")
    rm(server.error)
    server.up <- try(url(paste(secondary.server, test.file, sep=""), 
                         open = "rt"), silent=T)
    server.error <- inherits(server.up, "try-error")
    if(server.error) {
      if(verbose) warning(secondary.server.name, " failed")
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
         new("TagSNP",
             snpid=as.character(tag.id),
             population=toupper(as.character(tag.population))))
}
CreateCorrelatedSNPs<- function(tag.snp.name, snp.list, primary.server,
                                snp.region, bio.features.file, populations,
                                verbose, method.p, reduce.by.features, 
                                window.size, par.threads=1) {
  message()
  timestamp()
  tag.snp.complete <- try(PullInVariants(tag.snp.name, snp.list, primary.server,
                                         snp.region, populations, verbose, 
                                         window.size, par.threads), silent=TRUE)
  tag.snp.error <- inherits(tag.snp.complete, "try-error")
  if(tag.snp.error){
    if(identical(length(grep("not in 1000 genomes data",tag.snp.complete[[1]])),
                             as.integer(0))){
      while(tag.snp.error) {
        tag.snp.complete <- try(PullInVariants(tag.snp.name, snp.list,
                                               primary.server,
                                               snp.region, populations, verbose,
                                               window.size, par.threads), 
                                               silent=TRUE)
        tag.snp.error <- inherits(tag.snp.complete, "try-error")

      }
    } else {
      if(par.threads > 1) {
        message("\n #### The Tag SNP ", tag.snp.name, 
                " seems to be unavailable from the current ",
           "1000 genomes data \n", 
           " please check this Tag SNP on the 1000 genomes browser: \n",
           "http://browser.1000genomes.org/")
        stop("\n #### The Tag SNP ", tag.snp.name, 
             " seems to be unavailable from the current ",
             "1000 genomes data \n", 
             " please check this Tag SNP on the 1000 genomes browser: \n",
             "http://browser.1000genomes.org/")
      } else {
        stop("\n #### The Tag SNP ", tag.snp.name, 
             " seems to be unavailable from the current ",
             "1000 genomes data \n", 
             " please check this Tag SNP on the 1000 genomes browser: \n",
             "http://browser.1000genomes.org/")
      }
    }
  }
  if(!(tag.snp.error)){
    overlapping.features(correlated.snps(tag.snp.complete)) <-
      IRanges::unlist(GRangesList(lapply(as.character(bio.features.file),
                                         FilterByFeatures,
                                         tag.snp.name=tag.snp.name,
                                         tag.snp.complete=tag.snp.complete,
                                         verbose=verbose)))
    if(reduce.by.features) {
      overlapping.snps <- unique(names(overlapping.features(correlated.snps(tag.snp.complete))))
      ALL.overlapping.snps.geno(correlated.snps(tag.snp.complete)) <- ALL.geno(correlated.snps(tag.snp.complete))[, c(overlapping.snps, snpid(tag.snp.complete))]
      AFR.overlapping.snps.geno(correlated.snps(tag.snp.complete)) <- AFR.geno(correlated.snps(tag.snp.complete))[, c(overlapping.snps, snpid(tag.snp.complete))]
      AMR.overlapping.snps.geno(correlated.snps(tag.snp.complete)) <- AMR.geno(correlated.snps(tag.snp.complete))[, c(overlapping.snps, snpid(tag.snp.complete))]
      ASN.overlapping.snps.geno(correlated.snps(tag.snp.complete)) <- ASN.geno(correlated.snps(tag.snp.complete))[, c(overlapping.snps, snpid(tag.snp.complete))]
      EUR.overlapping.snps.geno(correlated.snps(tag.snp.complete)) <- EUR.geno(correlated.snps(tag.snp.complete))[, c(overlapping.snps, snpid(tag.snp.complete))]
      tag.snp.complete <- FilteredLDTesting(tag.snp.complete, verbose)
      if(verbose) message("Calculating p-value for ", tag.snp.name)
      ALL.p.value(correlated.snps(tag.snp.complete)) <-
        ChiSquaredPvalue(ALL.overlapping.snps.geno(correlated.snps(tag.snp.complete)),
                         snpid(tag.snp.complete), method.p)
      AFR.p.value(correlated.snps(tag.snp.complete)) <-
        ChiSquaredPvalue(AFR.overlapping.snps.geno(correlated.snps(tag.snp.complete)),
                         snpid(tag.snp.complete), method.p)
      AMR.p.value(correlated.snps(tag.snp.complete)) <-
        ChiSquaredPvalue(AMR.overlapping.snps.geno(correlated.snps(tag.snp.complete)),
                         snpid(tag.snp.complete), method.p)
      ASN.p.value(correlated.snps(tag.snp.complete)) <-
        ChiSquaredPvalue(ASN.overlapping.snps.geno(correlated.snps(tag.snp.complete)),
                         snpid(tag.snp.complete), method.p)
      EUR.p.value(correlated.snps(tag.snp.complete)) <-
        ChiSquaredPvalue(EUR.overlapping.snps.geno(correlated.snps(tag.snp.complete)),
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
    message("\nTag SNP ", snpid(tag.snp.complete), " has completed")
    return(tag.snp.complete)
  } else {
    tag.id <- snpid(snp.list[[tag.snp.name]])
    message("Tag SNP ", tag.id, " has an error; skipping ahead")
    return(NULL)
  }
  gc(verbose = FALSE)
}

CreatePopulations <- function(primary.server="ncbi") {
  onek.genome.server <- ServerCheck(primary.server, verbose = FALSE)
  manifest <- read.delim(paste(onek.genome.server,
                               "ftp/release/20110521/phase1_integrated_calls.20101123.ALL.panel", sep=""),
                         sep="\t", header = FALSE)
  for(i in c("AFR", "ASN", "EUR", "AMR", "ALL")) {
    ifelse(i != "ALL",
           assign(i, as.character(subset(manifest, manifest[, 3]==i)[,1])),
           assign(i, as.character(manifest[, 1])))
  }
  populations <- list(AFR=AFR, AMR=AMR, ASN=ASN, EUR=EUR, ALL=ALL)
  return(populations)
}

FunciSNP <- function(snp.regions.file, bio.features.loc = NULL,
                     par.threads=detectCores()/2,
                     verbose = par.threads < 2, method.p = "BH", 
                     reduce.by.features = TRUE, search.window = 200000 ) {
  message("\n",
	  "| | _  |  _  _ __  _    _|_ _ \n",
	  "|^|(/_ | (_ (_)|||(/_    |_(_)\n",
"\n",
	  " __             __    _ \n",
	  "|_    __  _  o (_ |\\||_)\n",
	  "|  |_|| |(_  | __)| ||  \n",
	  "\nVersion: ", package.version("FunciSNP"),"\n",
          "System: ",
          if(Sys.info()[["sysname"]] != "Windows") {
            Sys.info()[["sysname"]]
          } else {
            paste(Sys.info()[["sysname"]], " :-( ; parallel code not in effect, reverting to 1 core", sep="")
          }, "\n",
          "::args used::", "\n",
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
            length(list.files(bio.features.loc,
                              pattern="*.bed$",
                              full.names=FALSE)), ": ",
            gsub(".bed$", ", ", list.files(bio.features.loc,
                                           pattern="*.bed$",
                                           full.names=FALSE)))
  }
  primary.server <- sample(c("ncbi", "ebi"), size=1)
  snp.region <- ReadRegionsFile(snp.regions.file, search.window)
  message("          Number of TagSNPs Interrogated:   ", nrow(snp.region), 
          " representing ", length(unique(snp.region$snp.name)), " unique tag SNPs")
  options(mc.cores=par.threads)
  populations <- CreatePopulations(primary.server)
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
                       reduce.by.features=reduce.by.features,
                       window.size=search.window,
                       par.threads=par.threads)
  names(snp.list) <- tag.snp.names
  return(snp.list)
}


PullInVariants <- function(tag.snp.name, snp.list, primary.server, snp.region,
                           populations, verbose = TRUE, window.size, par.threads=1) {
  tag.id <- snpid(snp.list[[tag.snp.name]])
  window.size <- prettyNum(window.size, big.mark=",", scientific = FALSE)
  if(verbose) message("loading ", tag.id, " window size: ", window.size, " bp")
  primary.server <- sample(c("ncbi", "ebi"), size=1, prob=c(1,2))
  onek.genome.server <- ServerCheck(primary.server, verbose = FALSE)
  variants.reference <-
    paste(onek.genome.server,
          "/ftp/release/20101123/interim_phase1_release/ALL.chr",
          snp.region$snp.chromosome[snp.region$snp.name == tag.id][1],
          ".phase1.projectConsensus.genotypes.vcf.gz", sep="")

  param <- GRanges(seqnames=snp.region$snp.chromosome[snp.region$snp.name == tag.id][1],
                   IRanges(snp.region$snp.region.start[snp.region$snp.name == tag.id][1],
                           snp.region$snp.region.end[snp.region$snp.name == tag.id][1]))
  primes <- c(2, 3, 5, 7, 11, 13, 17, 19, 23, 29,
              31, 37, 41, 43, 47, 53, 59, 61, 67, 71,
              73, 79, 83, 89, 97, 101, 103, 107, 109, 113, 
              127, 131, 137, 139, 149, 151, 157, 163, 167, 173,
              179, 181, 191, 193, 197, 199, 211, 223, 227, 229,
              233, 239, 241, 251, 257, 263, 269, 271, 277, 281,
              283, 293, 307, 311, 313, 317, 331, 337, 347, 349,
              353, 359, 367, 373, 379, 383, 389, 397, 401, 409,
              419, 421, 431, 433, 439, 443, 449, 457, 461, 463,
              467, 479, 487, 491, 499, 503, 509, 521, 523, 541,
              547, 557, 563, 569, 571, 577, 587, 593, 599, 601,
              607, 613, 617, 619, 631, 641, 643, 647, 653, 659,
              661, 673, 677, 683, 691, 701, 709, 719, 727, 733,
              739, 743, 751, 757, 761, 769, 773, 787, 797, 809,
              811, 821, 823, 827, 829, 839, 853, 857, 859, 863,
              877, 881, 883, 887, 907, 911, 919, 929, 937, 941,
              947, 953, 967, 971, 977, 983, 991, 997, 1009, 1013)

  primes <- primes[4:(4+par.threads)]
  wait.time <- sample(primes, size=1)

  kgeno <- try(TabixFile(variants.reference), silent = TRUE)
  while(inherits(kgeno, "try-error")){
    Sys.sleep(wait.time)
    kgeno <- try(TabixFile(variants.reference), silent = TRUE)
  }

  tabix.header <-
    try(strsplit(headerTabix(variants.reference)$header, split="\t")[[6]], silent = TRUE)
  while(inherits(tabix.header, "try-error")){
    Sys.sleep(wait.time)
  tabix.header <-
    try(strsplit(headerTabix(variants.reference)$header, split="\t")[[6]], silent = TRUE)
  }

  if(verbose) message("scanning tabix file for ", tag.id)
  
  tabix.file <-
    try(strsplit(scanTabix(kgeno, param=param)[[1]], split="\t"), silent = TRUE)
  while(inherits(tabix.file, "try-error")) {
    message("Delay Connecting to ", primary.server, ", waiting ", wait.time,
            " seconds to try SNP ", tag.id, " again")
    Sys.sleep(wait.time)
    tabix.file <-
      try(strsplit(scanTabix(kgeno, param=param)[[1]], split="\t"), silent = TRUE)
    wait.time <- sample(primes, size=1)
    if(primary.server == "ncbi") {
      primary.server <- "ebi"
    } else {
      primary.server <- "ncbi"
    }
    onek.genome.server <- ServerCheck(primary.server, verbose = FALSE)
    variants.reference <-
      paste(onek.genome.server,
            "/ftp/release/20101123/interim_phase1_release/ALL.chr",
            snp.region$snp.chromosome[snp.region$snp.name == tag.id][1],
            ".phase1.projectConsensus.genotypes.vcf.gz", sep="")
    kgeno <- TabixFile(variants.reference)
  }
  snps.support <- lapply(tabix.file, function(x) x[1:5])
  snps.support <-
    t(as.data.frame(snps.support, stringsAsFactors = FALSE))
  
  timer <- 0
  while(match(tag.id, snps.support[, 3], nomatch=0) == 0 && timer < 5) {
    tabix.file <-
      try(strsplit(scanTabix(kgeno, param=param)[[1]], split="\t"), silent = TRUE)
    while(inherits(tabix.file, "try-error")) {
      message("Delay Connecting to ", primary.server, ", waiting ", wait.time,
              " seconds to try SNP ", tag.id, " again")
      Sys.sleep(wait.time)
      tabix.file <-
        try(strsplit(scanTabix(kgeno, param=param)[[1]], split="\t"), silent = TRUE)
      wait.time <- sample(primes, size=1)
      if(primary.server == "ncbi") {
        primary.server <- "ebi"
      } else {
        primary.server <- "ncbi"
      }
      onek.genome.server <- ServerCheck(primary.server, verbose = FALSE)
      variants.reference <-
        paste(onek.genome.server,
              "/ftp/release/20101123/interim_phase1_release/ALL.chr",
              snp.region$snp.chromosome[snp.region$snp.name == tag.id][1],
              ".phase1.projectConsensus.genotypes.vcf.gz", sep="")
      kgeno <- TabixFile(variants.reference)
    }
    snps.support <- lapply(tabix.file, function(x) x[1:5])
    snps.support <-
      t(as.data.frame(snps.support, stringsAsFactors = FALSE))
    timer <- timer + 1
  }


  colnames(snps.support) <- tabix.header[1:5]
  row.names(snps.support) <- NULL

  if(verbose) message("creating snpMatrix for ", tag.id)
  snps.geno <- try(vcf2sm(kgeno, gr=param, nmetacol=9L), silent = TRUE)
  while(inherits(snps.geno, "try-error")) {
    message("Delay Connecting to ", primary.server, ", waiting ", wait.time,
            " seconds to try SNP ", tag.id, " again")
    Sys.sleep(wait.time)
    snps.geno <- try(vcf2sm(kgeno, gr=param, nmetacol=9L), silent = TRUE)
    wait.time <- sample(primes, size=1)
    onek.genome.server <- ServerCheck(primary.server, verbose = FALSE)
    variants.reference <-
      paste(onek.genome.server,
            "/ftp/release/20101123/interim_phase1_release/ALL.chr",
            snp.region$snp.chromosome[snp.region$snp.name == tag.id][1],
            ".phase1.projectConsensus.genotypes.vcf.gz", sep="")
    kgeno <- TabixFile(variants.reference)
  }

  if((match(tag.id, snps.support[, 3], nomatch=0) == 0) || (dim(snps.support)[1] != dim(snps.geno)[2])) {
    #message(tag.id)
    #message(dim(snps.support))
    stop("Tag SNP id ", tag.id, " not in 1000 genomes data, look for alternative ID")
  }
  tag.support <- snps.support[(snps.support[, 3] %in% tag.id)]
  snps.support <- snps.support[!(snps.support[, 3] %in% tag.id), ]

  corr.snps <- new("CorrelatedSNPs",
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
  chr(snp.list[[tag.snp.name]]) <- as.integer(tag.support[1])
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
                           as.character(chr(correlated.snps(tag.snp.complete))),
                           sep=""),
            ranges=(IRanges(
                            start=as.integer(position(correlated.snps(tag.snp.complete))),
                            width=1)),
            snp.id=snpid(correlated.snps(tag.snp.complete)))
  names(close.snp.ranges) <- elementMetadata(close.snp.ranges)[, "snp.id"]
  if(identical(bio.features.file, NULL)) {
    snps.included <- close.snp.ranges
  } else {
    bio.features.file.interval <- import(bio.features.file, asRangedData = FALSE)
    bio.features.file.interval <- IRanges::sort(bio.features.file.interval)
    elementMetadata(bio.features.file.interval) <- NULL
    elementMetadata(bio.features.file.interval)[, "feature"] <-
      c(substr(grep(".bed$",
                    unlist(strsplit(bio.features.file, "/")), value=TRUE),
               1,
               nchar(grep(".bed$",
                          unlist(strsplit(bio.features.file, "/")),
                          value=TRUE))-4))
#    snps.included <-
#      !is.na(match(close.snp.ranges, bio.features.file.interval) > 0)
#    message("there")
#    ranges.included <-
#      !is.na(match(bio.features.file.interval, close.snp.ranges) > 0)
#    message("done")
#    stop("here")
#    peaks.with.overlapping.snps <-
#      subset(bio.features.file.interval, ranges.included)
#    snps.overlapping.peaks <-
#      subset(close.snp.ranges, snps.included)
    overlaps <-
      findOverlaps(bio.features.file.interval,
                   close.snp.ranges,
                   select="all")
    snps.included <-
      lapply(queryHits(overlaps), function(x) bio.features.file.interval[x])

    ##write method to call a specific snp that occurs under two peaks
    ##use the following technique test00[names(test00)=="rs123456"]
  }

  if(length(snps.included) >= 1) {
    snps.included <- IRanges::unlist(GRangesList(snps.included))
    names(snps.included) <- names(close.snp.ranges)[subjectHits(overlaps)]
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

  corr.snp.depth <- (dim(ALL.geno)[2]) - 1
  ALL.R.squared(correlated.snps(tag.snp.complete)) <-
    ld(ALL.geno(correlated.snps(tag.snp.complete))[, snp.name.chosen],
       ##questionable line break
       ALL.geno(correlated.snps(tag.snp.complete))[, !(colnames(ALL.geno(correlated.snps(tag.snp.complete))))
                                                   %in% snp.name.chosen], stats="R.squared", depth=corr.snp.depth)
  AFR.R.squared(correlated.snps(tag.snp.complete)) <-
    ld(AFR.geno(correlated.snps(tag.snp.complete))[, snp.name.chosen],
       AFR.geno(correlated.snps(tag.snp.complete))[, !(colnames(AFR.geno(correlated.snps(tag.snp.complete))))
                                                   %in% snp.name.chosen], stats="R.squared", depth=corr.snp.depth)
  AMR.R.squared(correlated.snps(tag.snp.complete)) <-
    ld(AMR.geno(correlated.snps(tag.snp.complete))[, snp.name.chosen],
       AMR.geno(correlated.snps(tag.snp.complete))[, !(colnames(AMR.geno(correlated.snps(tag.snp.complete))))
                                                   %in% snp.name.chosen], stats="R.squared", depth=corr.snp.depth)
  ASN.R.squared(correlated.snps(tag.snp.complete)) <-
    ld(ASN.geno(correlated.snps(tag.snp.complete))[, snp.name.chosen],
       ASN.geno(correlated.snps(tag.snp.complete))[, !(colnames(ASN.geno(correlated.snps(tag.snp.complete))))
                                                   %in% snp.name.chosen], stats="R.squared", depth=corr.snp.depth)
  EUR.R.squared(correlated.snps(tag.snp.complete)) <-
    ld(EUR.geno(correlated.snps(tag.snp.complete))[, snp.name.chosen],
       EUR.geno(correlated.snps(tag.snp.complete))[, !(colnames(EUR.geno(correlated.snps(tag.snp.complete))))
                                                   %in% snp.name.chosen], stats="R.squared", depth=corr.snp.depth)


  ALL.D.prime(correlated.snps(tag.snp.complete)) <-
    ld(ALL.geno(correlated.snps(tag.snp.complete))[, snp.name.chosen],
       ALL.geno(correlated.snps(tag.snp.complete))[, !(colnames(ALL.geno(correlated.snps(tag.snp.complete))))
                                                   %in% snp.name.chosen], stats="D.prime", depth=corr.snp.depth)
  AFR.D.prime(correlated.snps(tag.snp.complete)) <-
    ld(AFR.geno(correlated.snps(tag.snp.complete))[, snp.name.chosen],
       AFR.geno(correlated.snps(tag.snp.complete))[, !(colnames(AFR.geno(correlated.snps(tag.snp.complete))))
                                                   %in% snp.name.chosen], stats="D.prime", depth=corr.snp.depth)
  AMR.D.prime(correlated.snps(tag.snp.complete)) <-
    ld(AMR.geno(correlated.snps(tag.snp.complete))[, snp.name.chosen],
       AMR.geno(correlated.snps(tag.snp.complete))[, !(colnames(AMR.geno(correlated.snps(tag.snp.complete))))
                                                   %in% snp.name.chosen], stats="D.prime", depth=corr.snp.depth)
  ASN.D.prime(correlated.snps(tag.snp.complete)) <-
    ld(ASN.geno(correlated.snps(tag.snp.complete))[, snp.name.chosen],
       ASN.geno(correlated.snps(tag.snp.complete))[, !(colnames(ASN.geno(correlated.snps(tag.snp.complete))))
                                                   %in% snp.name.chosen], stats="D.prime", depth=corr.snp.depth)
  EUR.D.prime(correlated.snps(tag.snp.complete)) <-
    ld(EUR.geno(correlated.snps(tag.snp.complete))[, snp.name.chosen],
       EUR.geno(correlated.snps(tag.snp.complete))[, !(colnames(EUR.geno(correlated.snps(tag.snp.complete))))
                                                   %in% snp.name.chosen], stats="D.prime", depth=corr.snp.depth)



  R.squared.corrsnps(tag.snp.complete) <-
    ld(eval(parse(text=(paste(population(tag.snp.complete),
                              ".geno(correlated.snps(tag.snp.complete))",
                              sep="")))), stats="R.squared", depth=corr.snp.depth)
  D.prime.corrsnps(tag.snp.complete) <-
    ld(eval(parse(text=(paste(population(tag.snp.complete),
                              ".geno(correlated.snps(tag.snp.complete))",
                              sep="")))), stats="D.prime", depth=corr.snp.depth)

  return(tag.snp.complete)
}

FilteredLDTesting <- function(tag.snp.complete, verbose = TRUE) {
  if(ncol(eval(parse(text=(paste(population(tag.snp.complete),
                              ".overlapping.snps.geno(correlated.snps(tag.snp.complete))",
                              sep=""))))) > 1) {
  if(verbose) message("Calculating R^2 and D' for ", snpid(tag.snp.complete))
  snp.name.chosen <- snpid(tag.snp.complete)
  #  snp.name.chosen <- grep("^rs", unlist(strsplit(snp.list.name, ":")), value=T)
  corr.snp.depth <- (dim(eval(parse(text=(paste(population(tag.snp.complete),
                              ".overlapping.snps.geno(correlated.snps(tag.snp.complete))",
                              sep="")))))[2]) - 1
  R.squared.corrsnps(tag.snp.complete) <-
    ld(eval(parse(text=(paste(population(tag.snp.complete),
                              ".overlapping.snps.geno(correlated.snps(tag.snp.complete))",
                              sep="")))), stats="R.squared", depth=corr.snp.depth)
  D.prime.corrsnps(tag.snp.complete) <-
    ld(eval(parse(text=(paste(population(tag.snp.complete),
                              ".overlapping.snps.geno(correlated.snps(tag.snp.complete))",
                              sep="")))), stats="D.prime", depth=corr.snp.depth)
  return(tag.snp.complete)
  } else {
    return(tag.snp.complete)
  }
}

ChiSquaredPvalue <- function(tag.snp.complete, tag.snp.id, method.p) {
  snp.list <- lapply(colnames(tag.snp.complete),
                     function(x) as(tag.snp.complete[, x], "character"))
  names(snp.list) <- colnames(tag.snp.complete)
  genotype.table.snps <- lapply(snp.list,
                                function(x) table(unlist(snp.list[tag.snp.id]),
                                                  unlist(x),
                                                  dnn = (c("tag", "corr"))))

  pre.genotype.table.snps <- genotype.table.snps
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
  raw.p <- unlist(raw.p)
  adj.p <- p.adjust(raw.p, method=method.p)
  names(raw.p) <- names(snp.list)
  names(adj.p) <- names(snp.list)
  return(list(raw.p=raw.p, adj.p=adj.p))
}

SNPSummary <- function(snp.list) {
  if(length(overlapping.features(correlated.snps(snp.list))) > 0) {
    id.matrix.overlap <- as.data.frame(as.matrix(names(overlapping.features(correlated.snps(snp.list)))))
    id.matrix.overlap$pos <- dimnames(id.matrix.overlap)[[1]]
    id.matrix.complete <- as.data.frame(as.matrix(snpid(correlated.snps(snp.list))))
    id.matrix.complete$pos <- dimnames(id.matrix.complete)[[1]]
    id.matrix <- merge(id.matrix.overlap, id.matrix.complete, by.x = "V1", by.y = "V1", all.x = TRUE)
    dimnames(id.matrix)[[2]] <- c("corr.snp.id", "overlap.pos", "complete.pos")
    id.matrix$overlap.pos <- as.numeric(id.matrix$overlap.pos)
    id.matrix$complete.pos <- as.numeric(id.matrix$complete.pos)
    id.matrix <- id.matrix[order(id.matrix$overlap.pos), ]

    elementMetadata(overlapping.features(correlated.snps(snp.list)))[, "corr.snp.id"] <- snpid(correlated.snps(snp.list))[id.matrix$complete.pos]
    elementMetadata(overlapping.features(correlated.snps(snp.list)))[, "corr.snp.position"] <- position(correlated.snps(snp.list))[id.matrix$complete.pos]
    elementMetadata(overlapping.features(correlated.snps(snp.list)))[, "tag.snp.id"] <- snpid(snp.list)
    elementMetadata(overlapping.features(correlated.snps(snp.list)))[, "tag.snp.position"] <- position(snp.list)
    elementMetadata(overlapping.features(correlated.snps(snp.list)))[, "D.prime"] <- D.prime.corrsnps(snp.list)[, snpid(snp.list)][snpid(correlated.snps(snp.list))[id.matrix$complete.pos]]
    elementMetadata(overlapping.features(correlated.snps(snp.list)))[, "R.squared"] <- R.squared.corrsnps(snp.list)[, snpid(snp.list)][snpid(correlated.snps(snp.list))[id.matrix$complete.pos]]
    elementMetadata(overlapping.features(correlated.snps(snp.list)))[, "p.value"] <- 
      eval(parse(text=(paste(population(snp.list),
                             ".p.value(correlated.snps(snp.list))[[2]][snpid(correlated.snps(snp.list))[id.matrix$complete.pos]]",
                             sep=""))))
    elementMetadata(overlapping.features(correlated.snps(snp.list)))[, "distance.from.tag"] <- position(correlated.snps(snp.list))[id.matrix$complete.pos] - position(snp.list)
    elementMetadata(overlapping.features(correlated.snps(snp.list)))[, "population.count"] <- length(genotype(snp.list))
    elementMetadata(overlapping.features(correlated.snps(snp.list)))[, "population"] <- population(snp.list)
    return(overlapping.features(correlated.snps(snp.list)))
  } else {
    return(NULL)
  }
}

FunciSNPAnnotateSummary <- function(snp.list, verbose=TRUE) {
  if(identical(snp.list, NULL)) {
    return(NULL)
  } else{
    
    summary.snp.list <- lapply(snp.list, SNPSummary)
    summary.snp.list <- summary.snp.list[!sapply(summary.snp.list, is.null)]
    summary.snp.list <- IRanges::unlist(GRangesList(summary.snp.list))
    
    names(summary.snp.list) <-
      paste(names(summary.snp.list),
            elementMetadata(summary.snp.list)[, "feature"],
            sep=".")
    ### Taking only unique entries from summary.snp.list
    ### There can be duplicate row names when bio features have overlapping peaks
    ### especially when there are replicates

    summary.snp.list <- summary.snp.list[ which(!(duplicated(names(summary.snp.list)))), ]

    summary.snp.list <- as.data.frame(summary.snp.list)
    summary.snp.list$width <- NULL
    summary.snp.list$strand <- NULL
    
    colnames(summary.snp.list) <- c("chromosome",
                                    "bio.feature.start",
                                    "bio.feature.end",
                                    "bio.feature",
                                    "corr.snp.id",
                                    "corr.snp.position",
                                    "tag.snp.id",
                                    "tag.snp.position",
                                    "D.prime",
                                    "R.squared",
                                    "p.value",
                                    "distance.from.tag",
                                    "population.count",
                                    "population")
    
    gr.corr.snp.loc <- GRanges(seqnames=summary.snp.list$chromosome,
                               ranges=IRanges(start=summary.snp.list$corr.snp.position,
                                              width=1),
                               snpid=rownames(summary.snp.list)
                               )
    x <- strsplit(as.character(summary.snp.list$chromosome), split="chr")
    
    summary.snp.list$chromosome <- sapply(x, "[", 2)
    
    rd.corr.snp.loc <- RangedData(space=summary.snp.list$chromosome,
                                  ranges=IRanges(start=summary.snp.list$corr.snp.position,
                                                 width=1),
                                  snpid=rownames(summary.snp.list)
                                  )
    rownames(rd.corr.snp.loc) <- rd.corr.snp.loc$snpid
    rd.corr.snp.loc$snpid <- NULL


    data(TSS.human.GRCh37)
    data(lincRNA.hg19)
    
    ##nearest linc RNAs
    cat("Adding lincRNA")
    nearest.RNA <-
      annotatePeakInBatch(myPeakList = rd.corr.snp.loc,
                          AnnotationData = lincRNA,
                          output="nearestStart")
    summary.snp.list$nearest.lincRNA.ID <- NA
    summary.snp.list[nearest.RNA$peak, ]$nearest.lincRNA.ID <-
      nearest.RNA$feature
    summary.snp.list$nearest.lincRNA.ID <-
      as.factor(summary.snp.list$nearest.lincRNA.ID)

    summary.snp.list$nearest.lincRNA.distancetoFeature <- NA
    summary.snp.list[nearest.RNA$peak, ]$nearest.lincRNA.distancetoFeature <-
      nearest.RNA$distancetoFeature

    summary.snp.list$nearest.lincRNA.coverage <- NA
    summary.snp.list[nearest.RNA$peak, ]$nearest.lincRNA.coverage <-
      nearest.RNA$insideFeature
    summary.snp.list$nearest.lincRNA.coverage <-
      as.factor(summary.snp.list$nearest.lincRNA.coverage)
    cat(" ... done\n")
    ##nearest TSS (conanical gene)
    cat("Adding gene annotations\n\n")
    nearest.TSS <- annotatePeakInBatch(myPeakList = rd.corr.snp.loc, 
                                       AnnotationData = TSS.human.GRCh37, 
                                       output="nearestStart")
    nearest.TSS <- addGeneIDs(nearest.TSS, 
                              "org.Hs.eg.db", 
                              IDs2Add = c("symbol", "refseq"), 
                              silence = TRUE)
    nearest.TSS <- as(nearest.TSS, "GRanges")
    nearest.TSS <- nearest.TSS[order(elementMetadata(nearest.TSS)[, "peak"]), ]
    summary.snp.list <- summary.snp.list[order(row.names(summary.snp.list)), ]

    summary.snp.list$nearest.TSS.GeneSymbol <- elementMetadata(nearest.TSS)[, "symbol"]
    summary.snp.list$nearest.TSS.GeneSymbol <- as.factor(summary.snp.list$nearest.TSS.GeneSymbol)

    summary.snp.list$nearest.TSS.refseq <- NA
    summary.snp.list$nearest.TSS.refseq <- elementMetadata(nearest.TSS)[, "refseq"]
    summary.snp.list$nearest.TSS.refseq <- as.factor(summary.snp.list$nearest.TSS.refseq)

    summary.snp.list$nearest.TSS.ensembl <- NA
    summary.snp.list$nearest.TSS.ensembl <- elementMetadata(nearest.TSS)[, "feature"]
    summary.snp.list$nearest.TSS.ensembl <- as.factor(summary.snp.list$nearest.TSS.ensembl)

    summary.snp.list$nearest.TSS.coverage <- NA
    summary.snp.list$nearest.TSS.coverage <- elementMetadata(nearest.TSS)[, "insideFeature"]
    summary.snp.list$nearest.TSS.coverage <- as.factor(summary.snp.list$nearest.TSS.coverage)

    summary.snp.list$nearest.TSS.distancetoFeature <- NA
    summary.snp.list$nearest.TSS.distancetoFeature <- elementMetadata(nearest.TSS)[, "distancetoFeature"]
    #cat(" ... done\n")
    ## overlap genomic features (intergenic, utr5, utr3, intron, exon
    gr.corr.snp.loc <- gr.corr.snp.loc[order(elementMetadata(gr.corr.snp.loc)[,"snpid"]),]
    cat("\nAdding genomic annotations")
    gf.overlaps <- locateVariants(gr.corr.snp.loc, TxDb.Hsapiens.UCSC.hg19.knownGene)
    #cat(" ... done")
    genomic.feature <-as.character(gf.overlaps$Location)
    queryRow <-(gf.overlaps$queryHits)
    ddd <-(cbind(queryRow, genomic.feature)) ## used for identifying intergenic

    ## create set columns with null values ('NO')
    summary.snp.list$Promoter <- "NO"
    summary.snp.list$utr5 <- "NO"
    summary.snp.list$Exon <- "NO"
    summary.snp.list$Intron <- "NO"
    summary.snp.list$utr3 <- "NO"
    summary.snp.list$Intergenic <- "NO"

    ## promoter defined
    promoter.state <- subset(summary.snp.list, (nearest.TSS.distancetoFeature < 100) & (nearest.TSS.distancetoFeature > -1000))
    if(dim(promoter.state)[1] > 0) summary.snp.list[rownames(promoter.state), ]$Promoter <- "YES"
    summary.snp.list$Promoter <- as.factor(summary.snp.list$Promoter)
    ## utr5 defined
    utr5.rows <- as.numeric(subset(ddd, genomic.feature=="5'UTR", select="queryRow"))
    if(isTRUE(length(unique(utr5.rows)) > 0)) summary.snp.list[utr5.rows,"utr5"] <- "YES"; summary.snp.list$utr5 <- as.factor(summary.snp.list$utr5)
    ## exon defined
    exon.rows <- as.numeric(subset(ddd, genomic.feature=="coding")[,1])
    if(isTRUE(length(unique(exon.rows)) > 0)) summary.snp.list[exon.rows,"Exon"] <- "YES"; summary.snp.list$Exon <- as.factor(summary.snp.list$Exon)
    ## intron defined
    intron.rows <- as.numeric(subset(ddd, genomic.feature=="intron")[,1])
    if(isTRUE(length(unique(intron.rows)) > 0)) summary.snp.list[unique(intron.rows),"Intron"] <- "YES"; summary.snp.list$Intron <- as.factor(summary.snp.list$Intron)
    ## utr3 defined
    utr3.rows <- as.numeric(subset(ddd, genomic.feature=="3'UTR")[,1])
    if(isTRUE(length(unique(utr3.rows)) > 0)) summary.snp.list[utr3.rows,"utr3"] <- "YES"; summary.snp.list$utr3 <- as.factor(summary.snp.list$utr3)


    ## intergenic defined
    intergenic.rows <- as.numeric(subset(ddd, genomic.feature=="intergenic")[,1])
    if(isTRUE(length(unique(intergenic.rows)) > 0)) summary.snp.list[intergenic.rows,"Intergenic"] <- "YES"; summary.snp.list$Intergenic <- as.factor(summary.snp.list$Intergenic)
    promoter.intergenic.rows <- dimnames(subset(summary.snp.list, Intergenic=="YES" & Promoter=="YES"))[[1]]
    if(isTRUE(length(promoter.intergenic.rows) > 0)) summary.snp.list[promoter.intergenic.rows,"Intergenic"] <- "NO";
    cat(" ... done\n\nNow do the FunciStuff!\n");
    #rm("lincRNA"); ## remove object after annotation
    #rm("TSS.human.GRCh37"); ## remove object after annotation
    return(summary.snp.list)

  }
}

FunciSNPsummaryOverlaps <- function(dat, rsq=0) {
  dat <- subset(dat, R.squared >= rsq)
  tag.snps.with.overlaps <- unique(as.character(dat$tag.snp.id))

  tag.snp.features <- lapply(tag.snps.with.overlaps, function(x) {
                               overlap.counts <- count(as.character(dat[dat$tag.snp.id == x, ]$corr.snp.id))
                               overlap.counts$x <- as.character(overlap.counts$x)
                               max.freq <- max(overlap.counts$freq)
                               range.freq <- c(1:(max.freq+1))

                               z <- t(ldply(range.freq, function(x, overlap.counts) 
                                            { 
                                              y <- dim(overlap.counts[which(overlap.counts$freq == x), ])[[1]]
                                              z <- data.frame(SNPs=y)
                                              return(z)
                                            }, overlap.counts=overlap.counts))

                               colnames(z) <- range.freq
                               return(z)
                             })
  names(tag.snp.features) <- tag.snps.with.overlaps
  tag.snp.features <- lapply(tag.snp.features, function(x) {
                               y <- cumsum(x)
                               total <- rep(y[length(y)],length(y))
                               x + total - y
                             })
  columnnames <- c(1:max(unlist(lapply(tag.snp.features, 
                   function(x) dim(x)[[2]]))))
  overlap.counts <- do.call("rbind", lapply(tag.snp.features, "[", columnnames))
  overlap.counts[is.na(overlap.counts)] <- 0
  colnames(overlap.counts) <- paste("bio.",as.character(columnnames),sep="")
  return(overlap.counts)
}


FunciSNPidsFromSummary <- function(dat, tagsnpid, num.features, rsq=0) {
  dat <- subset(dat, R.squared >= rsq)
  overlap.counts <- 
        count(
        as.character(
        dat[dat$tag.snp.id == tagsnpid, ]$corr.snp.id))
  corr.overlapping.xfeatures <- 
                  subset(overlap.counts, freq >= num.features, select=x)
  corr.overlapping.xfeatures <- 
                  as.character(as.list(corr.overlapping.xfeatures)$x)
  summary.corr.snps <- 
                adply(corr.overlapping.xfeatures, 1, function(x, dat) {
                            dat[which(dat$corr.snp.id == x), ]
                             }, dat)
  summary.corr.snps$X1 <- NULL
  summary.corr.snps <- subset(summary.corr.snps, tag.snp.id == tagsnpid)
  rownames(summary.corr.snps) <- 
             paste(summary.corr.snps$tag.snp.id, ":", 
             summary.corr.snps$population, ".", 
             summary.corr.snps$corr.snp.id, ".", 
             summary.corr.snps$bio.feature, sep="")
  return(summary.corr.snps)
}




FunciSNPtable <- function(dat, rsq, geneSum = FALSE) {
  if(!geneSum){
    total.tagSNPs <-length(unique(dat[,"tag.snp.id"]))  
    total.1kSNPs  <-length(unique(dat[,"corr.snp.id"]))
    total.feature  <-length(unique(subset(dat, R.squared>rsq)[,"bio.feature"]))   

    total.tagSNPs.cutoff <- 
                     length(unique(subset(dat, R.squared>rsq)[,"tag.snp.id"]))
    total.1kSNPs.cutoff  <- 
                     length(unique(subset(dat, R.squared>rsq)[,"corr.snp.id"]))
    total.feature.cutoff  <- 
                     length(unique(subset(dat, R.squared>rsq)[,"bio.feature"]))

    total.summary.snp.list <- matrix(c(total.tagSNPs,total.1kSNPs,total.feature, 
                                       total.tagSNPs.cutoff,total.1kSNPs.cutoff,
                                       total.feature.cutoff), 
                                       nrow = 3, ncol=2, byrow=FALSE,
                                     dimnames = list(c("tagSNPs", 
                                                       "1kSNPs",
                                                       "bio.features"),
                                                     c("Total", 
                                                       paste("R.squared.cuff.",
                                                       rsq,sep=""))))
    total.summary.snp.list <- as.data.frame(total.summary.snp.list)
    total.summary.snp.list$Percent <- 
            round((total.summary.snp.list[,2]/total.summary.snp.list[,1])*100,2)
    return(total.summary.snp.list);
  } else {
    ###new function try out
    dat.s <- subset(dat, R.squared > rsq)
    y <- colSums(table(dat.s$corr.snp.id, dat.s$nearest.TSS.GeneSymbol))
    y[y==0] <- NA
    y <- na.omit(y)
    y <- as.matrix(names(y))
    y <- as.data.frame(y)
    dimnames(y)[[2]] <- "Gene_Names"
    return(y)
  }
}


FunciSNPplot <- function (dat, rsq = 0, split = FALSE, splitbysnp = FALSE, 
                          tagSummary = FALSE, heatmap = FALSE, 
                          genomicSum = FALSE) 
{
  if(sum(c(split,tagSummary,heatmap,genomicSum)) == 0){
     split = TRUE;
  }
require(ggplot2)
  if(split){
    if(splitbysnp == TRUE){
      FunciSNP:::theme_white()
      ggplot(dat, aes(x = R.squared)) + geom_histogram(binwidth = 0.05) + 
      geom_vline(xintercept = 0.5, linetype = 2) + 
      scale_x_continuous("Correlated SNPs R squared to Tag SNP (0-1)") + 
      scale_y_continuous("Total # of Correlated SNPs associated with tagSNP") + 
      opts(legend.position = "none", axis.text.y = theme_text(), 
           axis.text.x = theme_text(angle = 90), 
           title = "Distribution of correlated SNPs for each tagSNP\n",
                   "at Rsquared values") + 
      facet_wrap(chromosome ~ tag.snp.id)
    }else{
      tt <- count(df = dat, vars = "R.squared")
      tt <- na.omit(tt)
      ht <- range(tt[, "freq"])[2]*1.2
      hh <- dat[,c("corr.snp.id","R.squared")]
      hh <- na.omit(hh)
      hh.c <- count(round(hh$R.squared,digits = 1))
      dimnames(hh.c)[[1]] <- hh.c[,1]
      k <- c(hh.c["0",2],hh.c["0.1",2],hh.c["0.2",2],hh.c["0.3",2],
           hh.c["0.4",2],hh.c["0.5",2],hh.c["0.6",2],hh.c["0.7",2],
           hh.c["0.8",2],hh.c["0.9",2],hh.c["1",2])
      k[is.na(k)] <- 0;

      plot(tt, 
           xlim = c(0, 1), 
           ylim = c(0, ht), 
           pch = "*", 
           main = paste("Distribution of Correlated SNPs by Rsq values\n",
			"Total # of corr. SNPs: ",dim(dat)[1],
			"\n(with an Rsq value: ", sum(hh.c$freq),
			"; unique corr. SNPs: ", 
                        length(unique(hh$corr.snp.id)),")", 
                        sep = ""),
           xlab = "R square values (0-1)", 
           ylab = "Number of correlated SNPs")
      abline(v = 0.1, lty = 2, col = "red")
      abline(v = 0.2, lty = 2, col = "red")
      abline(v = 0.3, lty = 2, col = "red")
      abline(v = 0.4, lty = 2, col = "red")
      abline(v = 0.5, col = "black")
      abline(v = 0.6, lty = 2, col = "green")
      abline(v = 0.7, lty = 2, col = "green")
      abline(v = 0.8, lty = 2, col = "green")
      abline(v = 0.9, lty = 2, col = "green")
      abline(h = ht*.90, col = "black", lty = 2)
      text(0.05, ht*.95, as.character(k[1]+k[2]))
      text(0.15, ht*.95, as.character(k[3]))
      text(0.25, ht*.95, as.character(k[4]))
      text(0.35, ht*.95, as.character(k[5]))
      text(0.45, ht*.95, as.character(k[6]))
      text(0.55, ht*.95, as.character(k[7]))
      text(0.65, ht*.95, as.character(k[8]))
      text(0.75, ht*.95, as.character(k[9]))
      text(0.85, ht*.95, as.character(k[10]))
      text(0.95, ht*.95, as.character(k[11]))
    }
  }
  if(tagSummary){
    ## directory used to store plots
    try(dir.create(path=paste("FunciSNP.",package.version("FunciSNP"),"/plots",
    sep=""), showWarnings = FALSE, recursive=TRUE), silent=TRUE) 
    require("ggplot2")

    ### ggplot2 plots#####

    theme_white()

    all.s <- try(subset(dat, R.squared >= rsq), silent = TRUE)
    all.ss <- try(subset(dat, R.squared < rsq), silent = TRUE)
    try(all.s$r.2 <- c("Yes"), silent = TRUE)
    try(all.ss$r.2 <- c("No"), silent = TRUE)
    if(nrow(all.s) > 0 && nrow(all.ss) > 0) {
      all <- try(rbind(all.s, all.ss), silent = TRUE)
    } else {
      if(nrow(all.s) > 0 && nrow(all.ss) <= 0) {
        all <- all.s
      }
      if(nrow(all.s) <= 0 && nrow(all.ss) > 0) {
        all <- all.ss
      }
      if(nrow(all.s) <= 1 && nrow(all.ss) <= 0) {
        return()
      }
    }
    for( i in 1:length(summary(as.factor(all[,"bio.feature"]))) ){
      bio <- names(summary(as.factor(all[,"bio.feature"])))
      tmp <- subset(all, bio.feature==bio[i])

      ## plot r.2 values
      ggplot(tmp, aes(x=R.squared, fill=factor(r.2))) + 
      geom_histogram(binwidth=0.05) + 
      geom_vline(xintercept = rsq, linetype=2) +
      scale_x_continuous("Rsquare Values (0-1)", limits=c(0,1)) + 
      scale_y_continuous("Total # of Surrogate SNPs associated with riskSNP") + 
      scale_fill_manual(values = c("Yes" = "Red", "No" = "Black")) +
      opts(legend.position = "none", axis.text.y = theme_text(), 
           axis.text.x = theme_text(angle=90), 
           title = paste("riskSNP\nOverlapping: ", bio[i], sep="")) + 
      facet_wrap(~ tag.snp.id)

      ggsave(file=paste("FunciSNP.",package.version("FunciSNP"),"/plots/",
             bio[i],"_R2summary_riskSNP.pdf",sep=""))

      ## plot r.2 vs. distance values
      ggplot(tmp, aes(x=R.squared, y=distance.from.tag, colour=r.2, 
             size=factor(r.2))) + 
      geom_point() + 
      geom_vline(xintercept = rsq, linetype=2) +
      #geom_abline(intercept = 0, slope = 1) +
      scale_x_continuous("Rsquare Values (0-1)", limits=c(0,1)) + 
      scale_y_continuous(
           "Distance to Surrogate SNPs associated with riskSNP (bp)",
           formatter="comma") + 
      scale_colour_manual(values = c("Yes" = "Red", "No" = "Black")) +
      scale_size_manual(values = c("Yes" = 2, "No" = 1)) +
      opts(legend.position = "none", axis.text.y = theme_text(), 
           axis.text.x = theme_text(angle=90), 
           title = paste("Distance between riskSNP\n",
           "and Surrogate SNP\nOverlapping: ", bio[i], sep="")) + 
      facet_wrap(~ tag.snp.id)
      ggsave(file=paste("FunciSNP.",package.version("FunciSNP"),"/plots/",
             bio[i],"_R2vsDist_riskSNP.pdf",sep=""))
      cat("Finished plotting ", i, "/",length(bio), "\n")
    }
    message("\n\nSee ",
            paste("FunciSNP.",package.version("FunciSNP"),"/plots/",sep=""),
            " folder in ", getwd(), " for all plots.\n\n")
  }
  if(heatmap){
    require("gplots")
    require('matlab')
    ## directory used to store plots
    try(dir.create(path=paste("FunciSNP.",package.version("FunciSNP"),"/plots",
        sep=""), showWarnings = FALSE, recursive=TRUE), silent=TRUE)
    all.s<-(table( subset(dat,R.squared>=rsq)[,"bio.feature"], 
            subset(dat,R.squared>=rsq)[,"tag.snp.id"] ))
    all.s <- as.matrix(all.s)
    png(filename=paste("FunciSNP.",package.version("FunciSNP"),
        "/plots/FunciSNP_heatmap.png",sep=""), bg = "white")
    heatmap.2(
              all.s,
              na.rm=TRUE,
              scale="none",
              #RowSideColor=probe.cc,
              #ColSideColors=cc.col,
              col=jet.colors(max(all.s,na.rm=T)/2),
              #col=redgreen(75),
              key=T,
              symkey=FALSE,
              density.info="none",
              trace="none",
              Rowv=F,
              Colv=T,
              cexRow=1,
              cexCol=1,
              keysize=1,
              dendrogram=c("col"),
              main = paste("Rsquare values >= ",rsq)
              #labCol=NULL
              )
    dev.off()
    message("\nSee ",paste("FunciSNP.",package.version("FunciSNP"),
        "/plots/",sep=""), "folder in ", getwd()," for heatmap.\n\n")

  }
  if(genomicSum){
    if(rsq==0){
      dat.m <- melt(dat[,c(23:28)], 
               measure.vars=c("Promoter", 
                              "utr5", 
                              "Exon",
                              "Intron",
                              "utr3",
                              "Intergenic"))

      t <- subset(dat.m, value=="NO")
      t$value <- "2.NO"
      tt <- subset(dat.m, value!="NO")
      tt$value <- "1.YES"
      dat.m <- rbind(t,tt)


      require(ggplot2)
      FunciSNP:::theme_white()
      ggplot(dat.m, aes(variable, fill=factor(value))) + 
      geom_bar() +
      opts(axis.text.y = theme_text(), 
           axis.text.x = theme_text(angle = 90), 
           title = "Correlated SNPs distribution across genomic features") +
scale_fill_manual(values = c("1.YES" = "Red", "2.NO" = "Black"), "Overlap") +
#scale_x_continuous("Genomic Features") +
scale_y_continuous("Total count of correlated SNPs")

    } else {

      dat$r2 <- paste("R squared < ", rsq,sep="")
      t <- subset(dat, R.squared>=rsq)
      t$r2 <- paste("R squared >= ", rsq,sep="")
      dat <- rbind(t, subset(dat, R.squared<rsq)) 
      dat.m <- melt(dat[,c(23:29)], 
                    measure.vars=c("Promoter", 
                                   "utr5", 
                                   "Exon",
                                   "Intron",
                                   "utr3",
                                   "Intergenic"))

      t <- subset(dat.m, value=="NO")
      t$value <- "2.NO"
      tt <- subset(dat.m, value!="NO")
      tt$value <- "1.YES"
      dat.m <- rbind(t,tt)

      require(ggplot2)
      FunciSNP:::theme_white()
      ggplot(dat.m, aes(variable, fill=factor(value))) + 
      geom_bar(position="fill") +
      opts(axis.text.y = theme_text(), 
           axis.text.x = theme_text(angle = 90), 
           title = paste("Distribution of corr. SNPs across genomic features\n",
                         " at Rsquared cutoff of", rsq, sep=" ")) +
scale_fill_manual(values = c("1.YES" = "Red", "2.NO" = "Black"), "Overlap") +
#scale_x_continuous("Genomic Features") +
scale_y_continuous("Percent of Total correlated SNPs at Rsquared cutoff") +
facet_wrap(~ r2)
    }

  }
}

### generic functions used above ####
vlookup <- function(val, df, col){
  df[df[1] == val, col][1]
}


theme_white <- function() {
  require(ggplot2)
  theme_update (
   plot.background = theme_blank(),
   panel.background=theme_rect(colour="black", size=1),
   axis.text.x= theme_text(colour="black",vjust= 1, size=12),
   axis.text.y= theme_text(colour="black",hjust=1, size=12),
   axis.title.x =theme_text(colour="black",face="bold", size=12),
   axis.title.y =theme_text(colour="black",face="bold", angle = 90, size=12)
  )
}

yapply <- function(X,FUN, ...) { 
  index <- seq(length.out=length(X)) 
  namesX <- names(X) 
  if(is.null(namesX)) {
    namesX <- rep(NA,length(X))
  }
  FUN <- match.fun(FUN) 
  fnames <- names(formals(FUN)) 
  if(!"INDEX" %in% fnames) { 
    formals(FUN) <- append(formals(FUN), alist(INDEX=))
  } 
  if(!"NAMES" %in% fnames) { 
    formals(FUN) <- append(formals(FUN), alist(NAMES=))
  } 
  mapply(FUN,X,INDEX=index, NAMES=namesX,MoreArgs=list(...)) 
}

