## FunciSNP Code
## Author: Simon G. Coetzee; Houtan Noushmehr, PhD
## scoetzee@gmail.com; houtan@usp.br
## 310.570.2362
## All rights reversed.

ReadRegionsFile <- function(regions.file, search.window=200000) {
  # Reads a tab seperated regions file in the form
  # chr:loc snp_name    ethnicity
  # 8:130685457 rs4295627 EUR
  # returns the variables for snp.range, snp.name, and snp.ethno
  snp.regions <- read.table(regions.file)
  snp.region.split <- unlist(strsplit(as.vector(snp.regions[,1]), ":"))

  snp.chromosome <- as.character(sapply(strsplit(as.vector(snp.regions[,1]),
                                                 ":"), function(x) x[1]));
  snp.loc <- as.numeric(sapply(strsplit(as.vector(snp.regions[,1]), ":"),
                               function(x) x[2]));


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
    apollo <- "ftp://asclepius.hsc.usc.edu/1000genomes/"
    test.file <- "ftp/release/20110521/phase1_integrated_calls.20101123.ALL.panel"
    if(primary.server == "ebi"){
        primary.server <- ebi
        secondary.server <- ncbi
        primary.server.name <- "ebi"
        secondary.server.name <- "ncbi"
    } else if(primary.server == "ncbi"){
        primary.server <- ncbi
        secondary.server <- ebi
        primary.server.name <- "ncbi"
        secondary.server.name <- "apollo"
    } else if(primary.server == "apollo"){
        primary.server <- apollo
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
  if(verbose) {message(); timestamp()}
  tag.snp.complete <- try(PullInVariants(tag.snp.name, snp.list, primary.server,
                                         snp.region, populations, verbose,
                                         window.size, par.threads), silent=TRUE)
  tag.snp.error <- inherits(tag.snp.complete, "try-error")
  if(tag.snp.error){
    message(tag.snp.complete)
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
      overlapping.snps <-
          unique(names(overlapping.features(correlated.snps(tag.snp.complete))))
      tag.snp.complete <- FilteredLDTesting(tag.snp.complete, verbose)
      if(verbose) message("Calculating p-value for ", tag.snp.name)
      yafsnp.pvalue(tag.snp.complete) <- ChiSquaredPvalue(tag.snp.complete, snpid(tag.snp.complete), method.p)
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
    AFR <- NULL
    ASN <- NULL
    EUR <- NULL
    AMR <- NULL
    ALL <- NULL
    onek.genome.server <- ServerCheck(primary.server, verbose = FALSE)
    manifest <- read.delim(paste(onek.genome.server,
                                 "ftp/release/20110521/phase1_integrated_calls.20101123.ALL.panel",
                                 sep=""), sep="\t", header = FALSE)
    for(i in c("AFR", "ASN", "EUR", "AMR", "ALL")) {
        ifelse(i != "ALL",
               assign(i, as.character(subset(manifest, manifest[, 3]==i)[,1])),
               assign(i, as.character(manifest[, 1])))
    }
    populations <- list(AFR=AFR, AMR=AMR, ASN=ASN, EUR=EUR, ALL=ALL)
    return(populations)
}

getFSNPs <- function(snp.regions.file, bio.features.loc = NULL,
                     built.in.biofeatures = TRUE,
                     par.threads=detectCores()/2,
                     verbose = par.threads < 2, method.p = "BH",
                     search.window = 200000, primary.server="ebi" ) {
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
            paste(Sys.info()[["sysname"]], " :-( ; parallel code not in effect,
                  reverting to 1 core", sep="")
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
    message("          Bio Features:                     no biofeatures
            selected")
  } else {
    message("          Bio Features:                     ",
            if(built.in.biofeatures) {
              (length(list.files(bio.features.loc,
                                 pattern="*.bed$",
                                 full.names=FALSE))) + 5
            } else {
              length(list.files(bio.features.loc,
                                pattern="*.bed$",
                                full.names=FALSE))
            }, ": ",
            gsub(".bed$", ", ", list.files(bio.features.loc,
                                           pattern="*.bed$",
                                           full.names=FALSE)),

            if(built.in.biofeatures) {
            builtins <- gsub(".bed$", ", ", list.files(system.file('extdata',
                                                                   package='FunciSNP.data'),
                                              pattern="*.bed$",
                                              full.names=FALSE))
            builtins[length(builtins)] <- sub(", $", "", builtins[length(builtins)])
            builtins
            }
            )
  }
  snp.region <- ReadRegionsFile(snp.regions.file, search.window)
  message("          Number of TagSNPs Interrogated:   ", nrow(snp.region),
          " representing ", length(unique(snp.region$snp.name)), " unique tagSNPs")
          options(mc.cores=par.threads)
          populations <- CreatePopulations(primary.server)
          if(identical(bio.features.loc, NULL)) {
            bio.features.file <- NULL
            if(built.in.biofeatures) bio.features.file <-
              list.files(system.file('extdata',package='FunciSNP.data'),
                         pattern="known.bed$", full.names = TRUE)
          } else {
            bio.features.file <- list.files(bio.features.loc, pattern="*.bed$",
                                            full.names = TRUE)
            if(built.in.biofeatures) {
              bio.features.file.known <-
                list.files(system.file('extdata',package='FunciSNP.data'),
                           pattern="known.bed$", full.names = TRUE)
              bio.features.file <- c(bio.features.file, bio.features.file.known)
            }
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
                               reduce.by.features=TRUE,
                               window.size=search.window,
                               par.threads=par.threads)
          names(snp.list) <- tag.snp.names
          snp.list <- TSList(snp.list)
          return(snp.list)
}


PullInVariants <- function(tag.snp.name, snp.list, primary.server, snp.region,
                           populations, verbose = TRUE, window.size,
                           par.threads=1) {

  tag.id <- snpid(snp.list[[tag.snp.name]])
  window.size <- prettyNum(window.size, big.mark=",", scientific = FALSE)
  if(verbose) message("loading ", tag.id, " window size: ", window.size, " bp")
  onek.genome.server <- ServerCheck(primary.server, verbose = FALSE)
  variants.reference <-
    paste(onek.genome.server,
          "/ftp/release/20110521/ALL.chr",
          snp.region$snp.chromosome[snp.region$snp.name == tag.id][1],
          ".phase1_release_v3.20101123.snps_indels_svs.genotypes.vcf.gz", sep="")

  param <- GRanges(seqnames=snp.region$snp.chromosome[snp.region$snp.name ==
                   tag.id][1],
  IRanges(snp.region$snp.region.start[snp.region$snp.name ==
          tag.id][1],
  snp.region$snp.region.end[snp.region$snp.name ==
                            tag.id][1]))
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
    try(strsplit(headerTabix(variants.reference)$header, split="\t"),
        silent = TRUE)
  tabix.header <- tabix.header[[length(tabix.header)]]
  while(inherits(tabix.header, "try-error")){
    Sys.sleep(wait.time)
    tabix.header <-
      try(strsplit(headerTabix(variants.reference)$header, split="\t"),
          silent = TRUE)
    tabix.header <- tabix.header[[length(tabix.header)]]
  }

  if(verbose) message("scanning tabix file for ", tag.id)

  chunk <- try(scanTabix(kgeno, param=param), silent = TRUE)
  while(inherits(chunk, "try-error")) {
    message("Delay Connecting to ", primary.server, ", waiting ", wait.time,
            " seconds to try SNP ", tag.id, " again")
    Sys.sleep(wait.time)
    chunk <- try(scanTabix(kgeno, param=param), silent = TRUE)
    wait.time <- sample(primes, size=1)
    if(primary.server == "ncbi") {
        primary.server <- "ebi"
    } else if(primary.server == "ebi"){
        primary.server <- "apollo"
    } else if(primary.server == "apollo"){
        primary.server <- "ncbi"
    }
    onek.genome.server <- ServerCheck(primary.server, verbose = FALSE)
    variants.reference <-
        paste(onek.genome.server,
              "/ftp/release/20110521/ALL.chr",
              snp.region$snp.chromosome[snp.region$snp.name == tag.id][1],
              ".phase1_release_v3.20101123.snps_indels_svs.genotypes.vcf.gz", sep="")
    kgeno <- TabixFile(variants.reference)
  }


  ##### modified code from GGtools vcf2sm and parsVCFrec to work with remote tabix files
  if (Rsamtools:::isOpen(kgeno)) Rsamtools:::close.TabixFile(kgeno)
  Rsamtools:::open.TabixFile(kgeno)
  #  tabix.header <-
  #    try(strsplit(headerTabix(kgeno)$header, split="\t"),
  #        silent = TRUE)
#  tabix.header <- tabix.header[[length(tabix.header)]]
#  while(inherits(tabix.header, "try-error")){
#    Sys.sleep(wait.time)
#    tabix.header <-
#      try(strsplit(headerTabix(kgeno)$header, split="\t"),
#          silent = TRUE)
#    tabix.header <- tabix.header[[length(tabix.header)]]
#  }
  sampids <- tabix.header[10:length(tabix.header)]
  out <- list()
  for (i in 1:length(chunk)) {
    if (length(chunk[[i]]) == 0) next
    out[[i]] = lapply(chunk[[i]], function(rec) {
                      vec <- strsplit(rec, "\t")[[1]]
                      if((nchar(vec[4]) == 1) && (nchar(vec[5]) == 1)) {
                        meta <- vec[1:9]
                        calls <- vec[-c(1:9)]
                        nalt <- strsplit(calls, "")
                        nums <- lapply(nalt, "[", c(1,3))  # extract the call components
                        hasmiss <- which(sapply(nums, function(x) any(x == ".")))
                        nalt <- sapply(nums, function(x) 2-sum(x=="0"))  # this is correct only for diallelic locus; note in doc
                        if (length(hasmiss)>0) nalt[hasmiss] <- -1
                        nalt <- nalt+1
                        if (meta[3] == "." ) meta[3] <- paste("chr", meta[1], ":", meta[2], sep="")
                        x <- list(chr=meta[1], id=meta[3], loc=meta[2], ref=meta[4], alt=meta[5],
                                calls=as.raw(nalt))
                        return(x)
                      } else {
                        return(NULL)
                      }
          })
  }
  out <- unlist(out, recursive=FALSE)
  if (length(out) == 0) return(NULL)
  out.filter <- lapply(out, is.null)
  out.filter <- unlist(out.filter)
  out.filter <- !out.filter
  out <- out[out.filter]
  rsid <- sapply(out, "[[", "id")
  nsnp <- length(out)
  mat <- matrix(as.raw(0), nrow=length(sampids), ncol=nsnp)
  for (i in 1:nsnp) mat[,i] = out[[i]]$calls
  rownames(mat) <- sampids
  colnames(mat) <- rsid
  Rsamtools:::close.TabixFile(kgeno)
  snps.geno <- new("SnpMatrix", mat)
  snps.support <- NULL
  snps.support <- data.frame(CHROM=sapply(out, "[[", "chr"), POS=sapply(out, "[[", "loc"), ID=sapply(out, "[[", "id"), REF=sapply(out, "[[", "ref"), ALT=sapply(out, "[[", "alt"), stringsAsFactors = FALSE)
  ## End Code Adpated from GGtools
  timer <- 0
  while(match(tag.id, snps.support[, 3], nomatch=0) == 0 && timer < 5) {
    chunk <- try(scanTabix(kgeno, param=param), silent = TRUE)
    while(inherits(chunk, "try-error")) {
      message("Delay Connecting to ", primary.server, ", waiting ", wait.time,
              " seconds to try SNP ", tag.id, " again")
      Sys.sleep(wait.time)
      chunk <- try(scanTabix(kgeno, param=param), silent = TRUE)
      wait.time <- sample(primes, size=1)
      if(primary.server == "ncbi") {
        primary.server <- "ebi"
      } else {
        primary.server <- "ncbi"
      }
      onek.genome.server <- ServerCheck(primary.server, verbose = FALSE)
      variants.reference <-
        paste(onek.genome.server,
              "/ftp/release/20110521/ALL.chr",
              snp.region$snp.chromosome[snp.region$snp.name == tag.id][1],
              ".phase1_release_v3.20101123.snps_indels_svs.genotypes.vcf.gz", sep="")
      kgeno <- TabixFile(variants.reference)
    }
    ##### modified code from GGtools vcf2sm and parsVCFeec to work with remote tabix files
    if (Rsamtools:::isOpen(kgeno)) Rsamtools:::close.TabixFile(kgeno)
    Rsamtools:::open.TabixFile(kgeno)
    tabix.header <-
      try(strsplit(headerTabix(kgeno)$header, split="\t"),
          silent = TRUE)
    tabix.header <- tabix.header[[length(tabix.header)]]
    while(inherits(tabix.header, "try-error")){
      Sys.sleep(wait.time)
      tabix.header <-
        try(strsplit(headerTabix(kgeno)$header, split="\t"),
            silent = TRUE)
      tabix.header <- tabix.header[[length(tabix.header)]]
    }
    sampids <- tabix.header[10:length(tabix.header)]
    rsid <- sapply(out, "[[", "id")
    nsnp <- length(out)
    mat <- matrix(as.raw(0), nrow=length(sampids), ncol=nsnp)
    out <- list()
    for (i in 1:length(chunk)) {
      if (length(chunk[[i]]) == 0) next
      out[[i]] = lapply(chunk[[i]], function(rec) {
                        vec <- strsplit(rec, "\t")[[1]]
                        if((nchar(vec[4]) == 1) && (nchar(vec[5]) == 1)) {
                          meta <- vec[1:9]
                          calls <- vec[-c(1:9)]
                          nalt <- strsplit(calls, "")
                          nums <- lapply(nalt, "[", c(1,3))  # extract the call components
                          hasmiss <- which(sapply(nums, function(x) any(x == ".")))
                          nalt <- sapply(nums, function(x) 2-sum(x=="0"))  # this is correct only for diallelic locus; note in doc
                          if (length(hasmiss)>0) nalt[hasmiss] <- -1
                          nalt <- nalt+1
                          chr <- meta[1]
                          id <- meta[3]
                          loc <- meta[2]
                          if (id == "." ) id <- paste("chr", chr, ":", loc, sep="")
                          x <- list(chr=chr, id=id, loc=loc, ref=meta[4], alt=meta[5], depth=meta[8],
                                    calls=as.raw(nalt), support=meta[1:5])
                          return(x)
                        } else {
                          return(NULL)
                        }
            })
    }
    out <- unlist(out, recursive=FALSE)
    if (length(out) == 0) return(NULL)
    out.filter <- lapply(out, is.null)
    out.filter <- unlist(out.filter)
    out.filter <- !out.filter
    out <- out[out.filter]
    rsid <- sapply(out, "[[", "id")
    nsnp <- length(out)
    mat <- matrix(as.raw(0), nrow=length(sampids), ncol=nsnp)
    for (i in 1:nsnp) mat[,i] = out[[i]]$calls
    rownames(mat) <- sampids
    colnames(mat) <- rsid
    Rsamtools:::close.TabixFile(kgeno)
    snps.geno <- new("SnpMatrix", mat)
    snps.support <- sapply(out, "[[", "support")
    ## End Code Adpated from GGtools
    snps.support <-
      t(as.data.frame(snps.support, stringsAsFactors = FALSE))
    timer <- timer + 1
  }

  colnames(snps.support) <- tabix.header[1:5]
  row.names(snps.support) <- NULL
  if((match(tag.id, snps.support[, 3], nomatch=0) == 0) || (dim(snps.support)[1]
                                                            != dim(snps.geno)[2])) {
    #message(tag.id)
    #message(dim(snps.support))
    stop("TagSNP id ", tag.id, " not in 1000 genomes data, look for alternative
         ID")
  }
  tag.support <- snps.support[(snps.support[, 3] %in% tag.id), ]
  snps.support <- snps.support[!(snps.support[, 3] %in% tag.id), ]

  tag.ethno <- substr(tag.snp.name, nchar(tag.snp.name) - 2, nchar(tag.snp.name))
  tag.population <- getElement(populations, tag.ethno)

  corr.snps <- new("CorrelatedSNPs",
                   chromosome=as.character(snps.support[, 1]),
                   position=as.integer(snps.support[, 2]),
                   snpid=ifelse(snps.support[, 3] == ".",
                                paste("chr", snps.support[, 1], ":",
                                      snps.support[, 2], sep=""),
                                as.character(snps.support[, 3])),
                   ref.allele=as.character(snps.support[, 4]),
                   alt.allele=as.character(snps.support[, 5]),
                   genotype=new("CorrGeno", 
                                SnpMatrix=snps.geno,
                                populations=populations
                                )
                   )

  correlated.snps(snp.list[[tag.snp.name]]) <- corr.snps
  chr(snp.list[[tag.snp.name]]) <- as.character(tag.support[1])
  position(snp.list[[tag.snp.name]]) <- as.integer(tag.support[2])
  ref.allele(snp.list[[tag.snp.name]]) <- as.character(tag.support[4])
  alt.allele(snp.list[[tag.snp.name]]) <- as.character(tag.support[5])
  ##questionable line break
  genotype(snp.list[[tag.snp.name]]) <- snps.geno[, tag.id][row.names(snps.geno)
                                                            %in% tag.population]
  return(snp.list[[tag.snp.name]])
}

FilterByFeatures <- function(bio.features.file = NULL, tag.snp.name,
                             tag.snp.complete, verbose = TRUE) {

  ending.in.bed <- c(substr(grep(".bed$",
                                 unlist(strsplit(bio.features.file, "/")), value=TRUE),
                            1,
                            nchar(grep(".bed$",
                                       unlist(strsplit(bio.features.file, "/")),
                                       value=TRUE))-4))

  if(verbose) message("Filtering ", snpid(tag.snp.complete), " against ", ending.in.bed[[length(ending.in.bed)]])
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
    bio.features.file.interval <- import.bed(bio.features.file, asRangedData =
                                         FALSE)
    bio.features.file.interval <- IRanges::sort(bio.features.file.interval)
    elementMetadata(bio.features.file.interval) <- NULL
    elementMetadata(bio.features.file.interval)[, "feature"] <- ending.in.bed[[length(ending.in.bed)]]
    overlaps <- findOverlaps(bio.features.file.interval,
                             close.snp.ranges,
                             select="all")
    snps.included <- lapply(queryHits(overlaps), function(x) bio.features.file.interval[x])
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
            "\ttagSNP: \t\t", snpid(tag.snp.complete), "\n",
            "\tbiofeature: \t\t", feature.name, "\n",
            "\tpopulation: \t\t", population(tag.snp.complete));
    return(overlapping.features(correlated.snps(tag.snp.complete)))
  }
  message(snpid(tag.snp.complete), " has ",
          sum(unique(names(overlapping.features(correlated.snps(tag.snp.complete))))),
          " nearby SNPs overlapping with feature ", feature.name)
}

FilteredLDTesting <- function(tag.snp.complete, verbose = TRUE) {
  if(ncol(eval(parse(text=(paste(population(tag.snp.complete),
                                 ".overlapping.snps.geno(tag.snp.complete)",
                                 sep=""))))) > 1) {
    if(verbose) message("Calculating R\u00B2 and D' for ", snpid(tag.snp.complete))
    snp.name.chosen <- snpid(tag.snp.complete)
    corr.snp.depth <- (dim(eval(parse(text=(paste(population(tag.snp.complete),
                                                  ".overlapping.snps.geno(tag.snp.complete)",
                                                  sep="")))))[2]) - 1
    yafsnp.rsq(tag.snp.complete) <-
      ld(eval(parse(text=(paste(population(tag.snp.complete),
                                ".overlapping.snps.geno(tag.snp.complete)",
                                sep="")))),
         stats="R.squared",
         depth=corr.snp.depth)
    yafsnp.dprime(tag.snp.complete) <-
      ld(eval(parse(text=(paste(population(tag.snp.complete),
                                ".overlapping.snps.geno(tag.snp.complete)",
                                sep="")))), stats="D.prime", depth=corr.snp.depth)
    return(tag.snp.complete)
  } else {
    return(tag.snp.complete)
  }
}

ChiSquaredPvalue <- function(tag.snp.complete, tag.snp.id, method.p) {
tag.snp.complete <- eval(parse(text=(paste(population(tag.snp.complete),
                                ".overlapping.snps.geno(tag.snp.complete)",
                                sep=""))))
corr.snp.depth <- (dim(tag.snp.complete)[2]) - 1
  snp.list <- lapply(colnames(tag.snp.complete),
                     function(x) as(tag.snp.complete[, x], "character"))
  names(snp.list) <- colnames(tag.snp.complete)

  genotype.table.snps <- lapply(colnames(tag.snp.complete), function(x) {
         freq <- col.summary(tag.snp.complete[, x])[, c("P.AA", "P.AB", "P.BB")]
         tag.freq <- col.summary(tag.snp.complete[, tag.snp.id])[, c("P.AA", "P.AB", "P.BB")]
         genotype.table <- table(unlist(snp.list[tag.snp.id]),
                                 unlist(snp.list[x]),
                                 dnn = (c("tag", "corr")))
         if(identical(min(tag.freq), 0)) {
           fix.geno.table <- array(dim=c(3,3))
           colnames(fix.geno.table) <- c("A/A", "A/B", "B/B")
           rownames(fix.geno.table) <- c("A/A", "A/B", "B/B")
           if(identical((tag.freq)$P.AA, 0)) {
             fix.geno.table[1, ] <- 0
           } else {
             fix.geno.table[1, 1:ncol(genotype.table)] <- genotype.table["A/A", ]
           }
           if(identical((tag.freq)$P.AB, 0)) {
             fix.geno.table[2, ] <- 0
           } else {
             fix.geno.table[2, 1:ncol(genotype.table)] <- genotype.table["A/B", ]
           }
           if(identical((tag.freq)$P.BB, 0)) {
             fix.geno.table[3, ] <- 0
           } else {
             fix.geno.table[3, 1:ncol(genotype.table)] <- genotype.table["B/B", ]
           }
           genotype.table <- fix.geno.table
         }
         if(identical(min(freq), 0)) {
           fix.geno.table <- array(dim=c(3,3))
           colnames(fix.geno.table) <- c("A/A", "A/B", "B/B")
           rownames(fix.geno.table) <- c("A/A", "A/B", "B/B")
           if(identical((freq)$P.AA, 0)) {
             fix.geno.table[, 1] <- 0
           } else {
             fix.geno.table[1:nrow(genotype.table), 1] <- genotype.table[, "A/A"]
           }
           if(identical((freq)$P.AB, 0)) {
             fix.geno.table[, 2] <- 0
           } else {
             fix.geno.table[1:nrow(genotype.table), 2] <- genotype.table[, "A/B"]
           }
           if(identical((freq)$P.BB, 0)) {
             fix.geno.table[, 3] <- 0
           } else {
             fix.geno.table[1:nrow(genotype.table), 3] <- genotype.table[, "B/B"]
           }
           genotype.table <- fix.geno.table
         }
         return(genotype.table)
                                })
  names(genotype.table.snps) <- colnames(tag.snp.complete)
  OR <- ld(tag.snp.complete[, tag.snp.id],
           tag.snp.complete[, !colnames(tag.snp.complete) %in% tag.snp.id],
           stats="OR", depth=corr.snp.depth)
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
                    x <- round(x)
                    ## done b/c many tables have cells with < 5 values
                    return(suppressWarnings(fisher.test(x)$p.value))
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
    id.matrix.overlap <-
      as.data.frame(as.matrix(names(overlapping.features(correlated.snps(snp.list)))))
    id.matrix.overlap$pos <- dimnames(id.matrix.overlap)[[1]]
    id.matrix.complete <-
      as.data.frame(as.matrix(snpid(correlated.snps(snp.list))))
    id.matrix.complete$pos <- dimnames(id.matrix.complete)[[1]]
    id.matrix <- merge(id.matrix.overlap, id.matrix.complete, by.x = "V1", by.y
                       = "V1", all.x = TRUE)
    dimnames(id.matrix)[[2]] <- c("corr.snp.id", "overlap.pos", "complete.pos")
    id.matrix$overlap.pos <- as.numeric(id.matrix$overlap.pos)
    id.matrix$complete.pos <- as.numeric(id.matrix$complete.pos)
    id.matrix <- id.matrix[order(id.matrix$overlap.pos), ]

    elementMetadata(overlapping.features(correlated.snps(snp.list)))[,"corr.snp.id"] <- snpid(correlated.snps(snp.list))[id.matrix$complete.pos]
    elementMetadata(overlapping.features(correlated.snps(snp.list)))[,"corr.snp.position"] <- position(correlated.snps(snp.list))[id.matrix$complete.pos]
    elementMetadata(overlapping.features(correlated.snps(snp.list)))[,"tag.snp.id"] <- snpid(snp.list)
    elementMetadata(overlapping.features(correlated.snps(snp.list)))[,"tag.snp.position"] <- position(snp.list)
    elementMetadata(overlapping.features(correlated.snps(snp.list)))[,"D.prime"] <- yafsnp.dprime(snp.list)[, snpid(snp.list)][snpid(correlated.snps(snp.list))[id.matrix$complete.pos]]
    elementMetadata(overlapping.features(correlated.snps(snp.list)))[,"R.squared"] <- yafsnp.rsq(snp.list)[, snpid(snp.list)][snpid(correlated.snps(snp.list))[id.matrix$complete.pos]]
    elementMetadata(overlapping.features(correlated.snps(snp.list)))[,"p.value"] <- yafsnp.pvalue(snp.list)[[2]][snpid(correlated.snps(snp.list))[id.matrix$complete.pos]]
    elementMetadata(overlapping.features(correlated.snps(snp.list)))[,"distance.from.tag"] <- position(correlated.snps(snp.list))[id.matrix$complete.pos] - position(snp.list)
    elementMetadata(overlapping.features(correlated.snps(snp.list)))[,"population.count"] <- length(genotype(snp.list))
    elementMetadata(overlapping.features(correlated.snps(snp.list)))[,"population"] <- population(snp.list)
    return(overlapping.features(correlated.snps(snp.list)))
  } else {
    return(NULL)
  }
}



FunciSNPAnnotateSummary <- function(snp.list){
  return(snp.list@summary.data)
}

AnnotateSummary <- function(snp.list, verbose=TRUE) {
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
                ### There can be duplicate row names when bio features have overlapping
                ### peaks especially when there are replicates

                summary.snp.list <- summary.snp.list[
                                                     which(!(duplicated(names(summary.snp.list)))), ]

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

                data(TSS.human.GRCh37, package='ChIPpeakAnno')
                #    data(refseqgenes, package='FunciSNP')
                #    data(lincRNA.hg19, package='FunciSNP')

                ##nearest linc RNAs
                cat("Putative Functional SNPs identified!!\nAnnotation will begin\n~~\n")
                cat("Adding lincRNA")
                #    lincRNA <- system.file('extdata/annotation/lincRNA.hg19.rda',package='FunciSNP')
                #    load(lincRNA)
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
                cat("Adding gene annotations")
                #    refseqgenes <- system.file('extdata/annotation/refseqgenes.rda', package='FunciSNP')
                #    load(refseqgenes)
               txdb <<- TxDb.Hsapiens.UCSC.hg19.knownGene 

                nearest.TSS <- annotatePeakInBatch(myPeakList = rd.corr.snp.loc,
                                                   AnnotationData = refseqgenes,
                                                   output="nearestStart")
                nearest.TSS <- as(nearest.TSS, "GRanges")
                nearest.TSS <- nearest.TSS[order(elementMetadata(nearest.TSS)[, "peak"]), ]
                summary.snp.list <- summary.snp.list[order(row.names(summary.snp.list)), ]

                summary.snp.list$nearest.TSS.refseq <- NA
                summary.snp.list$nearest.TSS.refseq <- elementMetadata(nearest.TSS)[, "feature"]
                summary.snp.list$nearest.TSS.refseq <-
                        as.factor(summary.snp.list$nearest.TSS.refseq)

                summary.snp.list$nearest.TSS.GeneSymbol <- NA
                summary.snp.list$nearest.TSS.GeneSymbol <- refseqgenes$genesymbol[match(summary.snp.list$nearest.TSS.refseq, refseqgenes$names)]
                summary.snp.list$nearest.TSS.GeneSymbol <-
                        as.factor(summary.snp.list$nearest.TSS.GeneSymbol)

                summary.snp.list$nearest.TSS.ensembl <- NA
                summary.snp.list$nearest.TSS.ensembl <- refseqgenes$ensembl[match(summary.snp.list$nearest.TSS.refseq, refseqgenes$names)]
                summary.snp.list$nearest.TSS.ensembl <-
                        as.factor(summary.snp.list$nearest.TSS.ensembl)

                summary.snp.list$nearest.TSS.coverage <- NA
                summary.snp.list$nearest.TSS.coverage <- elementMetadata(nearest.TSS)[, "insideFeature"]
                summary.snp.list$nearest.TSS.coverage <-
                        as.factor(summary.snp.list$nearest.TSS.coverage)
                summary.snp.list$nearest.TSS.distancetoFeature <- NA
                summary.snp.list$nearest.TSS.distancetoFeature <-
                        elementMetadata(nearest.TSS)[, "distancetoFeature"]
                cat(" ... done\n")
                ## overlap genomic features (intergenic, utr5, utr3, intron, exon
                chroms <- append(unlist(lapply(as.character(1:22), function(x) {paste("chr", x, sep="")})), 'chrX')
                cat("\nAdding genomic annotations")
                myenv <- new.env()
                gr.corr.snp.loc <<- gr.corr.snp.loc
                gf.overlaps.utr5 <<- locateVariants(query=gr.corr.snp.loc[order(elementMetadata(gr.corr.snp.loc)[,"snpid"]),],
                                                    subject=txdb,
                                                    region=FiveUTRVariants(),
                                                    cache=myenv)
                gf.overlaps.utr3 <<- locateVariants(query=gr.corr.snp.loc[order(elementMetadata(gr.corr.snp.loc)[,"snpid"]),],
                                                    subject=txdb,
                                                    region=ThreeUTRVariants(),
                                                    cache=myenv)
                gf.overlaps <<- locateVariants(query=gr.corr.snp.loc[order(elementMetadata(gr.corr.snp.loc)[,"snpid"]),],
                                               subject=txdb,
                                               region=AllVariants(),
                                               cache=myenv)
                #cat(" ... done")
                #genomic.feature <- as.character(gf.overlaps$Location)
                #queryRow <-(gf.overlaps$queryHits)
                genomic.feature <- as.character(elementMetadata(gf.overlaps)[, "LOCATION" ])
                queryRow <- elementMetadata(gf.overlaps)[, "QUERYID"]
                ddd <-(cbind(queryRow, genomic.feature)) ## used for identifying intergenic
                ## create set columns with null values ('NO')
                summary.snp.list$Promoter <- "NO"
                summary.snp.list$utr5 <- "NO"
                summary.snp.list$Exon <- "NO"
                summary.snp.list$Intron <- "NO"
                summary.snp.list$utr3 <- "NO"
                summary.snp.list$Intergenic <- "NO"

                ## promoter defined
                #    promoter.state <- subset(summary.snp.list, (nearest.TSS.distancetoFeature <
                #                                                100) & (nearest.TSS.distancetoFeature > -1000))
                #    if(dim(promoter.state)[1] > 0) summary.snp.list[rownames(promoter.state),
                #
                # Promoter 2000 bp down 200 bp up
                promoter.rows <- as.numeric(subset(ddd, genomic.feature=="promoter")[,1])
                if(isTRUE(length(unique(promoter.rows)) > 0)){
                        summary.snp.list[promoter.rows,"Promoter"] <- "YES"; 
                        summary.snp.list$Promoter <- as.factor(summary.snp.list$Promoter)
                }
                ## utr5 defined
                utr5.rows <- as.numeric(subset(ddd, genomic.feature=="fiveUTR")[,1])
                if(isTRUE(length(unique(utr5.rows)) > 0)){
                        summary.snp.list[utr5.rows,"utr5"] <- "YES"; 
                        summary.snp.list$utr5 <- as.factor(summary.snp.list$utr5)
                }
                ## exon defined
                exon.rows <- as.numeric(subset(ddd, genomic.feature=="coding")[,1])
                if(isTRUE(length(unique(exon.rows)) > 0)){
                        summary.snp.list[exon.rows,"Exon"] <- "YES";
                        summary.snp.list$Exon <- as.factor(summary.snp.list$Exon)
                }
                ## intron defined
                intron.rows <- as.numeric(subset(ddd, genomic.feature=="intron")[,1])
                if(isTRUE(length(unique(intron.rows)) > 0)){
                        summary.snp.list[intron.rows,"Intron"] <- "YES";
                        summary.snp.list$Intron <- as.factor(summary.snp.list$Intron)
                }
                ## utr3 defined
                utr3.rows <- as.numeric(subset(ddd, genomic.feature=="threeUTR")[,1])
                if(isTRUE(length(unique(utr3.rows)) > 0)){
                        summary.snp.list[utr3.rows,"utr3"] <- "YES";
                        summary.snp.list$utr3 <- as.factor(summary.snp.list$utr3)
                }

                ## intergenic defined
                intergenic.rows <- as.numeric(subset(ddd, genomic.feature=="intergenic")[,1])
                if(isTRUE(length(unique(intergenic.rows)) > 0)){
                        summary.snp.list[intergenic.rows,"Intergenic"] <- "YES";
#                        summary.snp.list[which(summary.snp.list$Promoter == "YES"), "Intergenic"] <- "NO"
                        summary.snp.list$Intergenic <- as.factor(summary.snp.list$Intergenic)
                }
#                promoter.intergenic.rows <- dimnames(subset(summary.snp.list,
#                                                            Intergenic=="YES" & Promoter=="YES"))[[1]]
#                if(isTRUE(length(promoter.intergenic.rows) > 0)){
#                    summary.snp.list[promoter.intergenic.rows,"Intergenic"] <- "NO";
#                }
                cat(" ... done\n\nNow do the Funci Dance!\n");
                return(summary.snp.list)

        }
}

bedColors <- function(dat, rsq=0, filename, filepath) {
    jet.colors <- colorRampPalette(c("#00007F", "blue", "#007FFF", "cyan",
                                     "#7FFF7F", "yellow", "#FF7F00", "red", "#7F0000"))
    dat <- dat[which(dat$R.squared >= rsq), ]
    tag.snps.with.overlaps <- unique(as.character(dat$tag.snp.id))
    corr.snp.counts <- lapply(tag.snps.with.overlaps, function(x) {
                              count(as.character(dat[dat$tag.snp.id == x,]$corr.snp.id))})
    corr.snp.counts <- do.call("rbind", corr.snp.counts)
    max.freq <- max(corr.snp.counts$freq)
  #  corr.snp.counts$score <- ((corr.snp.counts$freq) / (max.freq)) * 1000
  colors <- t(col2rgb(jet.colors(max.freq)))
  corr.snp.counts$color <- unlist(lapply(corr.snp.counts$freq, function(x) {
                                         do.call(paste, c(as.list(colors[x, ]), sep=","))}))
  png(filename=paste(filepath, filename, sep="/"), width=800, height=200)
  z <- seq(-1, 1, length = 200)
  n <- max.freq
  image(matrix(z, ncol = 1), col = jet.colors(n), 
        xaxt = "n", yaxt = "n", main = "Overlap count Key")
  box()
  par(usr = c(0, n, 0, n))
  axis(1, at = c(0:n))
  dev.off()
  return(corr.snp.counts)
}


FunciSNPsummaryOverlaps <- function(dat, rsq=0) {
  dat <- dat[which(dat$R.squared >= rsq), ]
  if(dim(dat)[[1]] == 0) {
    x <- paste("At rsq > ", rsq, " no potentially correlated SNPs remain", sep="")
    return(x)
  } else {
  tag.snps.with.overlaps <- unique(as.character(dat$tag.snp.id))
  tag.snp.features <- lapply(tag.snps.with.overlaps, function(x) {
                             overlap.counts <- count(as.character(dat[dat$tag.snp.id ==
                                                                  x, ]$corr.snp.id))
                      overlap.counts$x <- as.character(overlap.counts$x)
                      max.freq <- max(overlap.counts$freq)
                      range.freq <- c(1:(max.freq))

                      z <- t(ldply(range.freq, function(x, overlap.counts)
                                   {
                                     y <-
                                       dim(overlap.counts[which(overlap.counts$freq
                                                                == x), ])[[1]]
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
  overlap.counts <- rbind(overlap.counts, colSums(overlap.counts))
  rownames(overlap.counts)[length(rownames(overlap.counts))] <- 
    "TOTAL # 1kgSNPs"
  return(overlap.counts)
  }
}


FunciSNPidsFromSummary <- function(dat, tagsnpid=NULL, num.features, rsq=0) {
  dat <- dat[which(dat$R.squared>=rsq), ]
  if(identical(tagsnpid, NULL)) {
    dat.sum <- FunciSNPsummaryOverlaps(dat=dat, rsq=rsq)
    tag.snps <- dat.sum[1:nrow(dat.sum)-1, num.features]
    tag.snps <- names(tag.snps[which(tag.snps > 0)])
    tagsnpid <- tag.snps
  }
  summary.corr.snps.list <- NULL
  for(i in tagsnpid) {
    overlap.counts <-
      count(
            as.character(
                         dat[which(dat$tag.snp.id == i), ]$corr.snp.id))
    corr.overlapping.xfeatures <- overlap.counts[which(overlap.counts$freq >=
                                                       num.features), ]
#  subset(overlap.counts, freq >= num.features, select=x)
  corr.overlapping.xfeatures <-
    as.character(as.list(corr.overlapping.xfeatures)$x)
  summary.corr.snps <-
    adply(corr.overlapping.xfeatures, 1, function(x, dat) {
          dat[which(dat$corr.snp.id == x), ]
                                                       }, dat)
  summary.corr.snps$X1 <- NULL
  summary.corr.snps <- summary.corr.snps[which(summary.corr.snps$tag.snp.id==i),
                                         ]
  summary.corr.snps.list <- rbind(summary.corr.snps.list, summary.corr.snps)
  }
  rownames(summary.corr.snps.list) <-
    paste(summary.corr.snps.list$tag.snp.id, ":",
          summary.corr.snps.list$population, ".",
          summary.corr.snps.list$corr.snp.id, ".",
          summary.corr.snps.list$bio.feature, sep="")
  return(summary.corr.snps.list)
}




FunciSNPtable <- function(dat, rsq, geneSum = FALSE) {
  if(!geneSum){
    total.tagSNPs <- length(unique(dat[,"tag.snp.id"]))
    total.1kSNPs <- length(unique(dat[,"corr.snp.id"]))
    total.feature <- length(unique(dat[,"bio.feature"]))

    total.tagSNPs.cutoff <-
      length(unique(dat[which(dat$R.squared>rsq),
                    "tag.snp.id"]))
    total.1kSNPs.cutoff  <-
      length(unique(dat[which(dat$R.squared>rsq),
                    "corr.snp.id"]))
    total.feature.cutoff  <-
      length(unique(dat[which(dat$R.squared>rsq),
                    "bio.feature"]))

    total.summary.snp.list <- matrix(c(total.tagSNPs,total.1kSNPs,total.feature,
                                       total.tagSNPs.cutoff,total.1kSNPs.cutoff,
                                       total.feature.cutoff),
                                     nrow = 3, ncol=2, byrow=FALSE,
                                     dimnames = list(c("tagSNPs",
                                                       "1K SNPs",
                                                       "Biofeatures"),
                                                     c("Total",
                                                       paste("R.sq>=",
                                                             rsq,sep=""))))
    total.summary.snp.list <- as.data.frame(total.summary.snp.list)
    total.summary.snp.list$Percent <-
      round((total.summary.snp.list[,2]/total.summary.snp.list[,1])*100,2)
    return(total.summary.snp.list);
  } else {
    ###new function try out
    dat.s <- dat[which(dat$R.squared > rsq), ]
    y <- colSums(table(dat.s$corr.snp.id, dat.s$nearest.TSS.GeneSymbol))
    y[y==0] <- NA
    y <- na.omit(y)
    y <- as.matrix(names(y))
    y <- as.data.frame(y)
    dimnames(y)[[2]] <- "Gene_Names"
    return(y)
  }
}

FunciSNPbed <- function(dat, rsq, path=getwd(), filename=NULL) {
  if(identical(filename,NULL)){
    filename <- paste("FunciSNP_results_rsq.",rsq,".bed",sep="")
  }else{
    filename <- filename
  }
  key.filename <- paste(filename, "colorkey", "png", sep=".")
  key.path <- path
  ###new function try out
  #d.s <- subset(dat, R.squared > rsq)
  d.s <- dat[which(dat$R.squared>rsq), ]
  d.cor <- d.s[ which(!(duplicated(d.s[,"corr.snp.id"]))), ]
  d.tag <- d.s[ which(!(duplicated(d.s[,"tag.snp.id"]))), ]
  d.tag <- (d.tag[,c(1,8,8,7,10,8,8)])
  d.cor <- (d.cor[,c(1,6,6,5,10,6,6)])
  d.cor$chromosome <- paste("chr",d.cor$chromosome,sep="")
  d.tag$chromosome <- paste("chr",d.tag$chromosome,sep="")
  d.cor[,2] <- d.cor[,2]-1
  d.tag[,2] <- d.tag[,2]-1
  d.cor[,6] <- d.cor[,6]-1
  d.tag[,6] <- d.tag[,6]-1
  d.cor$strand <- "+"
  d.tag$strand <- "+"
  d.cor$color <- NA
  d.tag$color <- "0,0,0"
  d.cor <- d.cor[,c(1:5,8,6:7,9)]
  d.tag <- d.tag[,c(1:5,8,6:7,9)]
  dimnames(d.cor)[[2]] <- c("chr", "snp.pos.s", "snp.pos.e", "snp.id",
                            "rsquare", "strand", "snp.pos.s", "snp.pos.e", "color")
  dimnames(d.tag)[[2]] <- dimnames(d.cor)[[2]]

  bed.colors <- bedColors(d.s, filename = key.filename, filepath = key.path)
  d.cor$color <- bed.colors$color[ match(d.cor$snp.id, bed.colors[ ,1])]

  d.cor$rsquare <- round(d.cor$rsquare, digits=4);
  d.tag$rsquare <- round(d.tag$rsquare, digits=4);
  d.cor$snp.id <- paste(d.cor$snp.id, "--", d.cor$rsquare, sep="")  
  y <- rbind(d.tag, d.cor); 
  con <- file(paste(path,filename,sep="/"),open="wt")
  writeLines(paste("browser position chr",d.s[1,1],":",
                   if(d.s[1,6]<d.s[1,8]){    
                     paste(d.s[1,6]-500,"-",d.s[1,8]+500,sep="")
                   }else{
                     paste(d.s[1,8]-500,"-",d.s[1,6]+500,sep="")
                   }
                   ,
                   "\ntrack name=\"FunciSNP_results\" description=\"Funci{SNP} Results : Rsquare cut-off at ", rsq, " (ver. ", 
                   package.version("FunciSNP"),")\" visibility=3 itemRgb=\"On\"", sep=""), con)
  write.table(y, row.names=F, col.names=F, sep="\t", file=con, quote=F)
  close(con)
  message("####\nBed file \"", filename, 
          "\" created successfully.\n(See folder: \"", path,"\")")
  cat("Total corSNP (RED): ", dim(d.cor)[1],"\nTotal tagSNP (BLK): ",
      dim(d.tag)[1],"\n")
  message("\nTo view results, submit bed file as a\ncustom track in ", 
          "UCSC Genome Browser (genome.ucsc.edu),", 
          " \n\nNow have fun with your new YAFSNPs",
          if(dim(d.cor)[1]>1){"s"}else{""},"!!\n####");

}


FunciSNPplot <- function (dat, rsq = 0, split = FALSE, splitbysnp = FALSE, 
                          tagSummary = FALSE, heatmap = FALSE, heatmap.key = FALSE,
                          genomicSum = FALSE, save = FALSE, pathplot=getwd(),
                          text.size=10, save.width=7, save.height=7) 
{
  save.width <- save.width * 25.4
  save.height <- save.height * 25.4
#  setThemeWhite(size = text.size)
#  fsnptheme <- theme_set(theme_bw(base_size = text.size))
  fsnptheme <- theme_bw(base_size = text.size) + theme(axis.text.x = element_text(angle = 90, size = text.size * 0.8, hjust=1))
  theme_set(fsnptheme)
#  require(scales)
  if(sum(c(split,tagSummary,heatmap,genomicSum)) == 0){
    split = TRUE;
  }
  if(save){
    try(dir.create(path=paste(pathplot, "/FunciSNP", "/plots",
                              sep=""), showWarnings = FALSE, recursive=TRUE), silent=TRUE) 
  }
  if(split){
    if(splitbysnp){
      p <- ggplot(dat, aes(x = R.squared)) +
      geom_histogram(binwidth = 0.05) + 
      geom_vline(xintercept = 0.5, linetype = 2) + 
      ggtitle("Distribution of 1kgSNPs for each tagSNP\nat R\u00B2 values") +
      scale_x_continuous(
                         "1kgSNPs R\u00B2 to tagSNP (0-1)") + 
         scale_y_continuous(
                            "Total # of 1kgSNPs associated with tagSNP") + 
         theme(legend.position = "none") + 
         facet_wrap(chromosome ~ tag.snp.id)
         if(save){
           ggplot2::ggsave(filename=paste(pathplot, "/FunciSNP",
                             "/plots/Distribution_for_each_tagSNP.pdf",
                             sep=""),
                  plot=p,
             dpi = 600,
             width = save.width, 
             height = save.height,
             units = "mm")
                  
         }else{
           return(p)
         }
    } else {
      tt <- count(df = dat, vars = "R.squared")
      tt <- na.omit(tt)
      ht <- range(tt[, "freq"])[2]*1.2
      hh <- dat[,c("corr.snp.id","R.squared")]
      hh <- na.omit(hh)
      hh.c <- count(round(hh$R.squared,digits = 1))
      dimnames(hh.c)[[1]] <- hh.c[,1]
      k <- c(hh.c["0",2], hh.c["0.1",2], hh.c["0.2",2],
             hh.c["0.3",2], hh.c["0.4",2], hh.c["0.5",2],
             hh.c["0.6",2], hh.c["0.7",2], hh.c["0.8",2],
             hh.c["0.9",2], hh.c["1",2])
      k[is.na(k)] <- 0;
      if(save){
        pdf(file=paste(pathplot, "/FunciSNP",
                       "/plots/Distribution_for_all_tagSNP.pdf",
                       sep=""), width=10, height=10)
      }
      plot(tt, 
           xlim = c(0, 1), 
           ylim = c(0, ht), 
           pch = "*", 
           main = paste(
                        "Distribution of 1kgSNPs by R\u00B2 values\n",
                        "Total # of 1kgSNPs: ",dim(dat)[1],
                        "\n(with an Rsq value: ", sum(hh.c$freq),
                        "; unique 1kgSNPs: ", 
                        length(unique(hh$corr.snp.id)),")", 
                        sep = ""),
           xlab = "R\u00B2 values (0-1)", 
           ylab = "Number of 1kgSNPs")
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
      if(save){
        dev.off() 
      }	
    }
  }
  if(tagSummary){
    try(dir.create(path=paste(pathplot, "/FunciSNP", "/plots",
                              sep=""), showWarnings = FALSE, recursive=TRUE), silent=TRUE)

    ### ggplot2 plots#####


    all.s <- try(dat[which(dat$R.squared >= rsq), ], silent = TRUE)
    all.ss <- try(dat[which(dat$R.squared < rsq), ], silent = TRUE)
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
      tmp <- all[which(all$bio.feature == bio[i]), ]
      ## plot r.2 values
      p <- ggplot(tmp, aes(x=R.squared, fill=factor(r.2))) + 
      geom_histogram(binwidth=0.05) + 
      geom_vline(xintercept = rsq, linetype=2) +
      scale_x_continuous("R\u00B2 Values (0-1)", limits=c(0,1)) + 
      scale_y_continuous("Total # of 1kgSNPs associated with riskSNP") + 
      scale_fill_manual(values = c("Yes" = "Red", "No" = "Black")) +
      ggtitle(paste("Distribution of 1kgSNPs R\u00B2",
                         "\ndivided by tagSNP & Overlapping biofeature:\n ", 
                         bio[i], sep="")) + 
      theme(legend.position = "none") +
         facet_wrap(chromosome ~ tag.snp.id)

         ggplot2::ggsave(filename=paste(pathplot, "/FunciSNP", "/plots/",
                           bio[i],"_R2summary_riskSNP.pdf",sep=""),
                plot=p,
             dpi = 600,
             width = save.width, 
             height = save.height,
             units = "mm")


         ## plot r.2 vs. distance values
         p <- ggplot(tmp, aes(x=R.squared, y=distance.from.tag, colour=r.2, 
                         size=factor(r.2))) + 
         geom_point() + 
         geom_vline(xintercept = rsq, linetype=2) +
         #geom_abline(intercept = 0, slope = 1) +
         scale_x_continuous("R\u00B2 Values (0-1)", limits=c(0,1)) + 
         scale_y_continuous("Distance to 1kgSNPs associated with tagSNP (bp)", labels = comma_format()) + 
         scale_colour_manual(values = 
                             c("Yes" = "Red", "No" = "Black")) +
         scale_size_manual(values = c("Yes" = 2, "No" = 1)) +
         ggtitle(paste("Distance between tagSNP ",
                            "and 1kgSNP\nOverlapping biofeature: ", 
                            bio[i], sep="")) +
         theme(legend.position = "none") +
         facet_wrap(chromosome ~ tag.snp.id)
         ggplot2::ggsave(filename=paste(pathplot, "/FunciSNP", "/plots/",
                           bio[i],"_R2vsDist_riskSNP.pdf",sep=""),
                plot=p,
             dpi = 600,
             width = save.width, 
             height = save.height,
             units = "mm")
         cat("Finished plotting ", i, "/",length(bio), "\n")
    }
    message("\n\nSee ",
            paste("FunciSNP","/plots/",sep=""),
            " folder in ", pathplot, " for all plots.\n\n")
  }
  if(heatmap){


    all.s <- table( dat[which(dat$R.squared>=rsq),"bio.feature"], 
                  dat[which(dat$R.squared>=rsq) ,"tag.snp.id"] )
#    rownames(all.s) <- paste(rownames(all.s), "\n(n=", rowSums(all.s), ")", sep="")
#    colnames(all.s) <- paste(colnames(all.s), "\n(n=", colSums(all.s), ")", sep="")
    x <- as.matrix(all.s)
#    dd.col <- as.dendrogram(hclust(dist(x)))
#    col.ord <- order.dendrogram(dd.col)
    dd.row <- as.dendrogram(hclust(dist(t(x))))
    row.ord <- order.dendrogram(dd.row)
    xx <- all.s[1:(nrow(all.s)), row.ord]
    xx_names <- attr(xx, "dimnames")
    df <- as.data.frame(xx)
    colnames(df) <- xx_names[[2]]
    df$sig <- xx_names[[1]]
    df$sig <- with(df, factor(sig, levels=sig, ordered=T))
    mdf <- reshape::melt(df, id.vars="sig")
    mdf$value <- as.numeric(mdf$value)
    all.s <- mdf
    if(isTRUE(heatmap.key)) {
      plot.here <- ggplot(all.s, aes(variable, sig, label=value)) +
      geom_tile(aes(fill=value), color="gray60") +
      scale_fill_gradient(low="white",
                          high="palevioletred4",
                          guide=guide_colorbar(direction = "horizontal", barheight = .5, "SNP count"),
                          "# of potentially correlated SNPs") +
      geom_text(size = text.size * 0.2) +
      labs(x = "", y = "") +
      ggtitle(paste("tagSNP vs Biofeature\n1kgSNP with R\u00B2 >= ", rsq, sep="")) +
      theme(axis.ticks = element_blank(),
#           axis.text.x = theme_text(angle=90, hjust = 1),
#           axis.text.y = theme_text(hjust=1),
           panel.background = element_rect(fill="white", colour="white"),
           legend.position = c(0,1),
           legend.justification=c(0,0))
      #    heatmap.2(
      #              all.s,
      #              na.rm=TRUE,
      #              scale="none",
      #              col=rev(terrain.colors(max(all.s, na.rm=TRUE))),
      #              key=TRUE,
      #              symkey=FALSE,
      #              density.info="none",
      #              trace="none",
      #              xlab="tagSNP",
      #              ylab="Biofeature",
      #              cellnote=hm.notes,
      #              colsep=c(1:(ncol(all.s)-1)),
      #              rowsep=c(1:(nrow(all.s)-1)),
      #              sepwidth=c(0.01, 0.01),
      #              sepcolor="black",
      #              notecol="black",
      #              notecex=0.75, 
      #              Rowv=FALSE,
      #              Colv=TRUE,
      #              cexRow=1,
      #              cexCol=1,
      #              keysize=0.75,
      #              dendrogram=c("none"),
      #              main = paste(
      #                           "tagSNP vs Biofeature\n1kgSNP with", 
      #                           "R\u00B2 >=", 
      #                           " ", rsq, sep=""))
    } else {
      plot.here <- ggplot(all.s, aes(variable, sig, label=value)) +
      geom_tile(aes(fill=value), color="gray60") +
      scale_fill_gradient(low="white",
                          high="palevioletred4",
                          guide=guide_colorbar(direction = "horizontal", barheight = .5, ""),
                          "# of potentially correlated SNPs") +
      labs(x = "", y = "") +
      ggtitle(paste("tagSNP vs Biofeature\n1kgSNP with R\u00B2 >= ", rsq, sep="")) +
      theme(axis.ticks = element_blank(),
#           axis.text.x = theme_text(angle=90, hjust = 1),
#           axis.text.y = theme_text(hjust=1),
           panel.background = element_rect(fill="white", colour="white"),
           legend.position = c(0,1),
           legend.justification=c(0,0))
      #    heatmap.2(
      #              all.s,
      #              na.rm=TRUE,
      #              scale="none",
      #              col=rev(terrain.colors(max(all.s, na.rm=TRUE))),
      #              key=TRUE,
      #              symkey=FALSE,
      #              density.info="none",
#              trace="none",
#              xlab="tagSNP",
#              ylab="Biofeature",
#              colsep=c(1:(ncol(all.s)-1)),
#              rowsep=c(1:(nrow(all.s)-1)),
#              sepwidth=c(0.01, 0.01),
#              sepcolor="black",
#              Rowv=FALSE,
#              Colv=TRUE,
#              cexRow=1,
#              cexCol=1,
#              keysize=0.75,
#              dendrogram=c("none"),
#              main = paste(
#                           "tagSNP vs Biofeature\n1kgSNP with ", 
#                           "R\u00B2 >=", 
#                           " ", rsq, sep="")
#              )
    }
    ### reverse matrix/dataframe x <- x[nrow(x):1, ]
    if(save) {
      ggplot2::ggsave(filename=paste(pathplot, "/FunciSNP",
                            "/plots/FunciSNP_heatmap.eps", sep=""),
             plot=plot.here, bg = "white",
             dpi = 600,
             width = save.width, 
             height = save.height,
             units = "mm")
      #	message("\nSee ",paste("FunciSNP.",package.version("FunciSNP"),
      #	 "/plots/",sep=""), "folder in ", pathplot," for heatmap.\n\n")
    } else {
      show(plot.here)
    }

  }
  if(genomicSum){
    if(rsq==0){
      dat.m <- reshape::melt(dat[,c(23:28)], 
                    measure.vars=c("Promoter", 
                                   "utr5", 
                                   "Exon",
                                   "Intron",
                                   "utr3",
                                   "Intergenic"))

      t <- dat.m[which(dat.m$value=="NO"), ]
      t$value <- "2.NO"
      tt <- dat.m[which(dat.m$value!="NO"), ]
      tt$value <- "1.YES"
      dat.m <- rbind(t,tt)


      qd <- ggplot(dat.m, aes(variable, fill=factor(value))) + 
      geom_bar() +
      ggtitle("1kgSNPs distribution across Genomic Features") +
      theme(axis.text.x = element_text(angle = 90, size = text.size*.8, hjust = 1)) +
      guides(fill = guide_legend(keywidth = .5, keyheight = 1)) +
      scale_x_discrete("") +
      scale_y_continuous("Total count of 1kgSNPs") +
      scale_fill_manual(values = c("1.YES" = "Red", "2.NO" = "Black"),
                        "Overlap")
      if(save){
        ggplot2::ggsave(filename=paste(pathplot, "/FunciSNP",
                          "/plots/Genomic_Summary_All.pdf", sep=""),
               plot=qd,
               dpi = 600,
               width = save.width, 
               height = save.height,
               units = "mm")
      }else{
        return(qd)
      }

    } else {

      dat$r2 <- paste("All 1kgSNPs", sep="")
      t <- dat[which(dat$R.squared >= rsq), ]
      t$r2 <- paste("R\u00B2 >= ", rsq, sep="")
      dat <- rbind(t, dat) 
      dat.m <- reshape::melt(dat[,c(23:29)], 
                    measure.vars=c("Promoter", 
                                   "utr5", 
                                   "Exon",
                                   "Intron",
                                   "utr3",
                                   "Intergenic"))

      t <- dat.m[which(dat.m$value=="NO"), ]
      t$value <- "2.NO"
      tt <- dat.m[which(dat.m$value!="NO"), ]
      tt$value <- "1.YES"
      dat.m <- rbind(t,tt)

      plot.title = paste("Distribution of 1kgSNP SNPs across Genomic Features\n",
                         " at R\u00B2 cut-off of", rsq, sep=" ")
      qp<-ggplot(dat.m, aes(variable, fill=factor(value))) + 
      geom_bar(position="fill") +
      ggtitle(plot.title) +
      theme(axis.text.x = element_text(angle = 90, size = text.size*.8, hjust = 1)) +
      guides(fill = guide_legend(keywidth = .5, keyheight = 1)) +
      scale_x_discrete("") +
      scale_fill_manual(values = c("1.YES" = "Red", "2.NO" = "Black"),"Overlap") +
      scale_y_continuous("Percent of Total 1kgSNPs at R\u00B2 cut-off") +
      facet_wrap(~ r2)
         if(save){
           ggplot2::ggsave(filename=paste(pathplot, "/FunciSNP",
                             "/plots/Genomic_Summary_by_rsq.", rsq, ".pdf",
                             sep=""),
                  plot=qp,
             dpi = 600,
             width = save.width, 
             height = save.height,
             units = "mm")
         }else{
           return(qp)
         }
    }

  }
}

### generic functions used above ####
vlookup <- function(val, df, col){
  df[df[1] == val, col][1]
}

setThemeWhite <- function(size = 10) {
theme.white <- theme_set(theme_grey(base_size = size))
theme.white <- theme_update(plot.background = element_blank(),
                            panel.background = element_rect(colour="black", size=.5))
theme_set(theme.white)
}

#theme_white <- function(size=10) {
#  require(ggplot2)
#  theme_grey() <- theme_update (
#                plot.background = theme_blank(),
#                panel.background=theme_rect(colour="black", size=1)
#              #  axis.text.x=theme_text(colour="black",vjust=1, angle=90),
#              #  axis.text.y=theme_text(colour="black",hjust=1),
#              #  axis.title.x=theme_text(colour="black",face="bold"),
#              #  axis.title.y=theme_text(colour="black",face="bold", angle = 90)
#                )
#  theme_grey(base_size = size)
#}

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

## FunciSNP Code                                                                 
## Author: Simon G. Coetzee; Houtan Noushmehr, PhD                               
## scoetzee@gmail.com; houtan@usp.br                                         
## 310.570.2362                                                                  
## All rights reversed.

