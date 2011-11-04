#!/usr/bin/env R

library("Rsamtools")
library("genetics")
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
    
    snp.chromosome <- grep("^\\d{1,2}$", snp.region.split, value=T, perl=T)

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
                              stringsAsFactors=F)
}

ChooseFeaturesFiles <- function(folder, verbose = TRUE) {
    # Reads the directory of the Features files, to check for bed files
    # and arranges them into a list
    bio.features.file <<- list.files(folder, pattern="*.bed$", full.names = TRUE)
    bio.features <<- list.files(folder, pattern="*.bed$", full.names = FALSE)
    if(verbose)
        cat("You have chosen", length(bio.features.file), "features\n", 
            print(bio.features), "\n", sep=" ")
}

LoopOverRiskSNPs <- function(snp.names = snp.regions$snp.name, 
                             ethno = c("AFR", "ASN", "EUR", "AMR", "ALL"), 
                             bio.features.file = bio.features.file) {
    manifest <<- read.delim("ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20110521/phase1_integrated_calls.20101123.ALL.panel", sep="\t",header=F)
    snp.names <<- snp.names
    ethno <<- ethno
    for (j in snp.names) {
        for (i in ethno) {
            ifelse (i != "ALL",
                    ethno.sample.set <<- as.character(subset(manifest, V3==i)[,1]),
                    ethno.sample.set <<- as.character(manifest[,1]))
            PullInVariants(i, j)
            print(paste("completed PullInVariants", i, j, sep=" "))
            for (h in bio.features.file) {
                print("began searching for biofeatures")
                bio.features.count <<- h
                FilterByFeatures(h)
                LDTesting(snp.geno)
                save(snp.ld.frame, file="risk.snp.Rda")
            }
        }
    }
write.table(snp.ld.frame, file="risk.snp.txt", sep="\t", quote=F, row.names=F)
}

    
PullInVariants <- function(ethno.chosen, snp.names.chosen) {
    ethno.chosen <<- ethno.chosen
    snp.names.chosen <<- snp.names.chosen
    intermediate.vcf <- paste(snp.names.chosen, "_", ethno.chosen, ".vcf", sep="")
    print(paste( "ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20101123/interim_phase1_release/ALL.chr", snp.regions$snp.chromosome[snp.regions$snp.name == snp.names.chosen][1], ".phase1.projectConsensus.genotypes.vcf.gz", sep=""))
    system(paste("tabix -hf ", "ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20101123/interim_phase1_release/ALL.chr", snp.regions$snp.chromosome[snp.regions$snp.name == snp.names.chosen][1], ".phase1.projectConsensus.genotypes.vcf.gz", " ", snp.regions$snp.chromosome[snp.regions$snp.name == snp.names.chosen][1], ":", snp.regions$snp.region.start[snp.regions$snp.name == snp.names.chosen][1], "-", snp.regions$snp.region.end[snp.regions$snp.name == snp.names.chosen][1], " > ", intermediate.vcf, " && ", "bgzip -f ", intermediate.vcf, " && ", "tabix -hf ", intermediate.vcf, ".gz" , sep=""))
    
    variants.file <- paste(intermediate.vcf, ".gz", sep="") 
    
    subset.variants.file <- unpackVcf(scanVcf(variants.file))[[1]]
    cat("check 1 \n") 
    variants.data <- data.frame(subset.variants.file$CHROM, 
                                subset.variants.file$POS, 
                                subset.variants.file$ID, 
                                subset.variants.file$REF, 
                                subset.variants.file$ALT,
                                subset.variants.file$GENO$GT[, c(ethno.sample.set)],
                                stringsAsFactors = FALSE)
    cat("check 2 \n")
    core.names <-  c("CHROM", "POS", "ID", "REF", "ALT")
    
    colnames(variants.data)[1:5] <- core.names
    
    variants.data.length <<- length(variants.data)

    cat("check 3 \n")
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
    cat("check 4 \n")
    variants.data$ID <- 
        ifelse(variants.data$ID == ".", variants.data$POS, variants.data$ID)
    variants.data <<- variants.data
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
    names(close.snp.ranges) <- variants.data$ID
    features.file.interval <<- import(features.file, asRangedData = FALSE)
    snps.included <- !is.na(match(close.snp.ranges, features.file.interval) > 0)
    if (sum(snps.included) >= 1) {
        snps.included <- subset(variants.data, snps.included)

        snp.geno <- data.frame(t(snps.included[, c(6:variants.data.length)]))
        colnames(snp.geno) <- snps.included$ID
        snp.geno <<- snp.geno[, colSums(is.na(snp.geno))<nrow(snp.geno)] #removes NA cols
        if (snp.names.chosen %in% colnames(snp.geno)) {
            print(paste(snp.names.chosen, ", the risk snp, already overlaps with the feature", features.file, sep=""))
        }
        else {
            temp <- data.frame(t(subset(variants.data, ID==snp.names.chosen))[6:variants.data.length, ])
            colnames(temp) <- snp.names.chosen
            snp.geno <- cbind(snp.geno, temp)
            rm(temp)
        }
        for (i in 1:ncol(snp.geno)) {
            snp.geno[, i] <- genotype(snp.geno[, i])
        }
        snp.geno <<- snp.geno
    }
    else {
        print(paste("There is no overlap for: \n",
                    "SNP: ", snp.names.chosen, "\n", 
                    "biofeature: ", features.file, "\n", 
                    "racial/ethnic group: ", ethno.chosen, sep=""))
    }
}

LDTesting <- function(snps) {
    snp.ld <- snps
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
    if (exists(as.character(substitute(snp.ld.frame))) == FALSE) {
        snp.ld.frame <- data.frame(snp.names.chosen, snp.names.cor, d, D.prime, r, R.squared, n, chi.squared, p.val)
        snp.ld.frame$feature <- unlist(strsplit(bio.features[bio.features.count], "\\.bed")) 
        snp.ld.frame$loc.feature <- paste(seqnames(subsetByOverlaps(features.file.interval, test[snp.names.corr])), ":", start(subsetByOverlaps(features.file.interval, test[snp.names.corr])), "-", end(subsetByOverlaps(features.file.interval, test[snp.names.corr])), sep="")
        snp.ld.frame$eth.risk <- snp.regions$snp.ethno[snp.regions$snp.name == snp.names.chosen][1]
        snp.ld.frame$eth.chosen <- ethno.chosen
        snp.ld.frame$chr <- snp.regions$snp.chromosome[snp.regions$snp.name == snp.names.chosen][1]
        snp.ld.frame$loc.risk <- subset(variants.data, ID==snp.names.chosen)$POS
        snp.ld.frame$loc.cor <- subset(variants.data, variants.data$ID %in% snp.names.cor)$POS
        snp.ld.frame$dist <- abs(snp.ld.frame$loc.risk - snp.ld.frame$loc.cor)
    }
    else {
        tmp.snp.ld.frame <- data.frame(snp.names.chosen, snp.names.cor, d, D.prime, r, R.squared, n, chi.squared, p.val)
        tmp.snp.ld.frame$feature <- unlist(strsplit(bio.features[bio.features.count], "\\.bed")) 
        tmp.snp.ld.frame$loc.feature <- paste(seqnames(subsetByOverlaps(features.file.interval, test[snp.names.corr])), ":", start(subsetByOverlaps(features.file.interval, test[snp.names.corr])), "-", end(subsetByOverlaps(features.file.interval, test[snp.names.corr])), sep="")
        tmp.snp.ld.frame$eth.risk <- snp.regions$snp.ethno[snp.regions$snp.name == snp.names.chosen][1]
        tmp.snp.ld.frame$eth.chosen <- ethno.chosen
        tmp.snp.ld.frame$chr <- snp.regions$snp.chromosome[snp.regions$snp.name == snp.names.chosen][1]
        tmp.snp.ld.frame$loc.risk <- subset(variants.data, ID==snp.names.chosen)$POS
        tmp.snp.ld.frame$loc.cor <- subset(variants.data, variants.data$ID %in% snp.names.cor)$POS
        tmp.snp.ld.frame$dist <- abs(snp.ld.frame$loc.risk - snp.ld.frame$loc.cor)
        snp.ld.frame <- rbind(snp.ld.frame, tmp.snp.ld.frame)
    }
    rownames(snp.ld.frame) <- NULL
    snp.ld.frame <<- snp.ld.frame
}
