#!/usr/bin/env R

library("Rsamtools")
library("genetics")
library("rtracklayer")
## TODO:::
## verbose options for everything
##
ReadRegionsFile <- function(regions.file) {
    # Reads a tab seperated regions file in the form
    # chr:start-end snp_name    ethnicity   chr#
    #
    # returns the variables for snp.range, snp.name, and snp.ethno
    snp.regions <- read.table(regions.file)
    snp.region.split <- unlist(strsplit(as.vector(snp.regions[,1]), ":"))
    
    snp.chromosome <- grep("^\\d{0,1}[0-9X-Y]$", snp.region.split, value = TRUE, perl = TRUE)

	m <- regexpr("^\\d{1,}-", snp.region.split)
	snp.region.start <- substr(snp.region.split, m, m + attr(m, "match.length") - 1)
	snp.region.start <- as.vector(na.omit(as.numeric(gsub("-$", "", snp.region.start))))

	m <- regexpr("-\\d{1,}$", snp.region.split)
	snp.region.end <- substr(snp.region.split, m, m + attr(m, "match.length") - 1)
	snp.region.end <- as.vector(na.omit(as.numeric(gsub("^-", "", snp.region.end))))
	
	snp.name <- snp.regions[, 2]
	snp.ethno <- snp.regions[, 3]
    snp.regions <<- data.frame(snp.chromosome, 
                              snp.region.start, 
                              snp.region.end, 
                              snp.name, 
                              snp.ethno, 
                              stringsAsFactors = FALSE)
}

LoopOverRiskSNPs <- function(snp.names = snp.regions$snp.name, 
                             ethno = c("AFR", "ASN", "EUR", "AMR", "ALL"), 
                             bio.features.loc, preferred.1000genomes.server) {
    
    ncbi <- "ftp://ftp-trace.ncbi.nih.gov/1000genomes/"
    ebi <- "ftp://ftp.1000genomes.ebi.ac.uk/vol1/"
    manifest.file <- "ftp/release/20110521/phase1_integrated_calls.20101123.ALL.panel"

    server.up <- url(paste(ncbi, manifest.file, sep=""))
    if(exists("server.up")) {
        try(open(server.up), silent = TRUE)
        print(1)
        if(isOpen(server.up)) {
            print(2)
            manifest <<- read.delim(paste(ncbi, manifest.file, sep=""), sep="\t", header = FALSE)
        } else {
            print(3)
            try(rm(server.up), silent = TRUE)
            server.up <- url(paste(ebi, manifest.file, sep=""))
            if(exists("server.up")) {
                print(4)
                try(open(server.up), silent = TRUE)
                if(isOpen(server.up)) {
                    print(5)
                    manifest <<- read.delim(paste(ebi, manifest.file, sep=""), sep="\t", header = FALSE)
                } else {
                    print(7)
                    stop("Neither EBI nor NCBI mirrors for the 1000 genomes project were found")
                }
            } else {
                print(6)
                stop("Neither EBI nor NCBI mirrors for the 1000 genomes project were found")
            }
        }
    }

    bio.features.file <<- list.files(bio.features.loc, pattern="*.bed$", full.names = TRUE)
    bio.features <<- list.files(bio.features.loc, pattern="*.bed$", full.names = FALSE)
    ifelse(isOpen(server.up),  
           manifest <<- read.delim(paste(ncbi, manifest.file, sep=""), sep="\t", header = FALSE),
           ifelse(isOpen(url(paste(ebi, manifest.file), "rt")),
                  manifest <<- read.delim(paste(ebi, manifest.file, sep=""), sep="\t", header = FALSE),
                  stop("Neither EBI nor NCBI mirrors for the 1000 genomes project were found")))
    
    for (j in snp.names) {
        for (i in ethno) {
            PullInVariants(i, j)
#            print(paste("completed PullInVariants", i, j, sep=" "))
            bio.features.count <- 0
            for (h in bio.features.file) {
                cat(paste("began searching for overlap with new biofeature: ", h, sep=""))
                bio.features.count <<- bio.features.count + 1
                FilterByFeatures(h)
                print("began LDTesting")
                LDTesting(snp.geno)
                try(save(snp.ld.frame, file="risk.snp.Rda"), silent = TRUE)
                try(write.table(snp.ld.frame, file="risk.snp.txt", sep="\t", quote = FALSE, row.names = FALSE), silent = TRUE)
            }
        }
    }
write.table(snp.ld.frame, file="risk.snp.txt", sep="\t", quote = FALSE, row.names = FALSE)
}

    
PullInVariants <- function(ethno.chosen, snp.names.chosen) {
    ethno.chosen <<- ethno.chosen
    ifelse (ethno.chosen != "ALL",
            ethno.sample.set <- as.character(subset(manifest, V3==ethno.chosen)[,1]),
            ethno.sample.set <- as.character(manifest[,1]))
#    print(ethno.sample.set)
    snp.names.chosen <<- snp.names.chosen
    print(snp.names.chosen)
    intermediate.vcf <- paste(snp.names.chosen, ".vcf", sep="")
    variants.file <- paste(intermediate.vcf, ".gz", sep="") 
    if (file.exists(variants.file) == FALSE) {
        ## the tabix file is being created for virtually every single loop
        ## this is most likely b/c of poor naming of intermediate and variants
        ## file.  Find a solution
        system(paste("tabix -hf ", "ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20101123/interim_phase1_release/ALL.chr", snp.regions$snp.chromosome[snp.regions$snp.name == snp.names.chosen][1], ".phase1.projectConsensus.genotypes.vcf.gz", " ", snp.regions$snp.chromosome[snp.regions$snp.name == snp.names.chosen][1], ":", snp.regions$snp.region.start[snp.regions$snp.name == snp.names.chosen][1], "-", snp.regions$snp.region.end[snp.regions$snp.name == snp.names.chosen][1], " > ", intermediate.vcf, " && ", "bgzip -f ", intermediate.vcf, " && ", "tabix -hf ", intermediate.vcf, ".gz" , sep=""))
        subset.variants.file <<- unpackVcf(scanVcf(variants.file))[[1]]
    }
    
    cat("check 1 \n")
    
    ifelse(ethno.chosen != "ALL",
           genotype.data <- subset.variants.file$GENO$GT[, c(ethno.sample.set)],
           genotype.data <- subset.variants.file$GENO$GT)
    
    variants.data <- data.frame(subset.variants.file$CHROM, 
                                subset.variants.file$POS, 
                                subset.variants.file$ID, 
                                subset.variants.file$REF, 
                                subset.variants.file$ALT,
                                genotype.data,
                                stringsAsFactors = FALSE)
#    summary(variants.data[,1:10]) 
    core.names <-  c("CHROM", "POS", "ID", "REF", "ALT")
    colnames(variants.data)[1:5] <- core.names
    variants.data.length <- length(variants.data)
    
    variants.data[, 6:variants.data.length] <- 
        ifelse(variants.data[, 6:variants.data.length] == "0|0", 
               paste(variants.data$REF, "/", variants.data$REF, sep=""),
               ifelse(variants.data[, 6:variants.data.length] == "0|1",
                      paste(variants.data$REF, "/", variants.data$ALT, sep=""),
                      ifelse(variants.data[, 6:variants.data.length] == "1|0",
                             paste(variants.data$ALT, "/", variants.data$REF, sep=""),
                             ifelse(variants.data[, 6:variants.data.length] == "1|1",
                                    paste(variants.data$ALT, "/", variants.data$ALT
                                          , sep=""),
                                    NA))))
    variants.data$ID <- 
        ifelse(variants.data$ID == ".", variants.data$POS, variants.data$ID)
    variants.data <<- variants.data
#    summary(variants.data[,1:10])
}

FilterByFeatures <- function(features.file) {
    close.snp.ranges <- 
        GRanges(seqnames=paste(
                               "chr",
                               as.character(variants.data$CHROM), 
                               sep=""), 
                ranges=(IRanges(
                                start=as.integer(variants.data$POS), 
                                width=1)))
#    print(close.snp.ranges)
    names(close.snp.ranges) <- variants.data$ID
#    print(close.snp.ranges)
    features.file.interval <<- import(features.file, asRangedData = FALSE)
    snps.included <- !is.na(match(close.snp.ranges, features.file.interval) > 0)
    if (sum(snps.included) >= 1) {
        snps.included <- subset(variants.data, snps.included)

        snp.geno <- data.frame(t(snps.included[, c(6:length(snps.included))]))
        colnames(snp.geno) <- snps.included$ID
        snp.geno <<- snp.geno[, colSums(is.na(snp.geno))<nrow(snp.geno)]
        if (snp.names.chosen %in% colnames(snp.geno)) {
            print(paste(snp.names.chosen, ", the risk snp, already overlaps with the feature", features.file, sep=""))
        } else {
            temp <- data.frame(t(subset(variants.data, ID==snp.names.chosen))[6:length(variants.data), ])
            colnames(temp) <- snp.names.chosen
#            print(temp)
            snp.geno <- cbind(snp.geno, temp)
            rm(temp)
        }
        for (i in 1:ncol(snp.geno)) {
            snp.geno[, i] <- genotype(snp.geno[, i])
        }
        snp.geno <<- snp.geno
    } else {
       try(rm(snp.geno), silent = TRUE)
       cat(paste("There is no overlap for: \n",
                    "\tRisk SNP: \t\t", snp.names.chosen, "\n", 
                    "\tbiofeature: \t\t", features.file, "\n", 
                    "\tracial/ethnic group: \t", ethno.chosen, "\n", sep=""))
    }
    close.snp.ranges <<- close.snp.ranges
}

LDTesting <- function(snps) {
    try(snp.ld <- LD(snps), silent = TRUE)
    if (exists("snp.ld")) {
        d <- na.omit(as.data.frame(snp.ld$"D"[, c(snp.names.chosen)]))
        colnames(d) <- "d"
        D.prime <- na.omit(as.data.frame(snp.ld$"D'"[, c(snp.names.chosen)]))
        colnames(D.prime) <- "D.prime"
        r <- na.omit(as.data.frame(snp.ld$"r"[, c(snp.names.chosen)]))
        colnames(r) <- "r"
        R.squared <- na.omit(as.data.frame(snp.ld$"R^2"[, c(snp.names.chosen)]))
        colnames(R.squared) <- "R.squared"
        n <- na.omit(as.data.frame(snp.ld$"n"[, c(snp.names.chosen)]))
        colnames(n) <- "n"
        chi.squared <- na.omit(as.data.frame(snp.ld$"X^2"[, c(snp.names.chosen)]))
        colnames(chi.squared) <- "chi.squared"
        p.val <- na.omit(as.data.frame(snp.ld$"P-value"[, c(snp.names.chosen)]))
        colnames(p.val) <- "p.val"
        snp.names.cor <- dimnames(p.val)[[1]]
        if (exists("snp.ld.frame") == FALSE) {
            snp.ld.frame <- data.frame(snp.names.chosen, snp.names.cor, d, D.prime, r, R.squared, n, chi.squared, p.val)
            snp.ld.frame$feature <- unlist(strsplit(bio.features[bio.features.count], "\\.bed")) 
            snp.ld.frame$loc.feature <- paste(seqnames(features.file.interval[c(as.matrix(findOverlaps(features.file.interval, close.snp.ranges[snp.names.cor]))[,1])]), ":", start(features.file.interval[c(as.matrix(findOverlaps(features.file.interval, close.snp.ranges[snp.names.cor]))[,1])]), "-", end(features.file.interval[c(as.matrix(findOverlaps(features.file.interval, close.snp.ranges[snp.names.cor]))[,1])]), sep="")
            ## loc.feature has "blank" areas where the only thing visible is :-, meaning that for some reason,
            ## features.file.interval[c(as.matrix(findOverlaps(features.file.interval, close.snp.ranges[snp.names.cor]))[,1])]
            ## is not being found sometimes
            snp.ld.frame$eth.risk <- snp.regions$snp.ethno[snp.regions$snp.name == snp.names.chosen][1]
            snp.ld.frame$eth.chosen <- ethno.chosen
            snp.ld.frame$chr <- snp.regions$snp.chromosome[snp.regions$snp.name == snp.names.chosen][1]
            snp.ld.frame$loc.risk <- subset(variants.data, ID==snp.names.chosen)$POS
            snp.ld.frame$loc.cor <- subset(variants.data, variants.data$ID %in% snp.names.cor)$POS
            snp.ld.frame$dist <- abs(snp.ld.frame$loc.risk - snp.ld.frame$loc.cor)
        } else {
            print(1)
            tmp.snp.ld.frame <- data.frame(snp.names.chosen, snp.names.cor, d, D.prime, r, R.squared, n, chi.squared, p.val)
            print(2)
            tmp.snp.ld.frame$feature <- unlist(strsplit(bio.features[bio.features.count], "\\.bed")) 
            print(3)
            tmp.snp.ld.frame$loc.feature <- paste(seqnames(features.file.interval[c(as.matrix(findOverlaps(features.file.interval, close.snp.ranges[snp.names.cor]))[,1])]), ":", start(features.file.interval[c(as.matrix(findOverlaps(features.file.interval, close.snp.ranges[snp.names.cor]))[,1])]), "-", end(features.file.interval[c(as.matrix(findOverlaps(features.file.interval, close.snp.ranges[snp.names.cor]))[,1])]), sep="")
            ## loc.feature has "blank" areas where the only thing visible is :-, meaning that for some reason,
            ## features.file.interval[c(as.matrix(findOverlaps(features.file.interval, close.snp.ranges[snp.names.cor]))[,1])]
            ## is not being found sometimes
            print(4)
            tmp.snp.ld.frame$eth.risk <- snp.regions$snp.ethno[snp.regions$snp.name == snp.names.chosen][1]
            print(5)
            tmp.snp.ld.frame$eth.chosen <- ethno.chosen
            print(6)
            tmp.snp.ld.frame$chr <- snp.regions$snp.chromosome[snp.regions$snp.name == snp.names.chosen][1]
            print(7)
            tmp.snp.ld.frame$loc.risk <- subset(variants.data, ID==snp.names.chosen)$POS
            print(8)
            tmp.snp.ld.frame$loc.cor <- subset(variants.data, variants.data$ID %in% snp.names.cor)$POS
            print(9)
            tmp.snp.ld.frame$dist <- abs(tmp.snp.ld.frame$loc.risk - tmp.snp.ld.frame$loc.cor)
            print(10)
            snp.ld.frame <- rbind(snp.ld.frame, tmp.snp.ld.frame)
        }
        rownames(snp.ld.frame) <- NULL
        snp.ld.frame <<- snp.ld.frame
    } else {
        print(paste("Fewer than two of the snps grouped with ", snp.names.chosen, ", that overlap with ", unlist(strsplit(bio.features[bio.features.count], "\\.bed")), " for the ", ethno.chosen, " racial/ethnic group, have more than one allele", sep=""))
    }
}
