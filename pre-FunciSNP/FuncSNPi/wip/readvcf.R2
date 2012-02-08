#!/usr/bin/env R

library("VariantAnnotation")
library("rtracklayer")
library("LDheatmap")
library("genetics")
x <- "/home/scoetzee/scripts/Coetzee_Seq_Analysis/FuncSNPi/wip/PCa_riskSNP_withEthnicity_2MB.txt"

ReadRegionsFile <- function(x) {
    # Reads a tab seperated regions file in the form
    # chr:start-end snp_name    ethnicity   chr#
    #
    # returns the variables for snp.range, snp.name, and snp.ethno
    snp.regions <- read.table(x)
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
    
    snp.regions <- data.frame(snp.chromosome, 
                              snp.region.start, 
                              snp.region.end, 
                              snp.name, 
                              snp.ethno, 
                              stringsAsFactors=F)
}

ChooseFeaturesFiles <- function(x, verbose = TRUE) {
    # Reads the directory of the Features files, to check for bed files
    # and arranges them into a list
    features <- list.files("x", pattern="*.bed$")
    if(verbose)
        cat("You have chosen", length(features), "features\n", 
            print(features), "\n", sep="")
}

ExtractGenomicWindow <- function(
features.file <- "/home/scoetzee/data/AR_hg19.bed"

input.ethno <- "AFR"

variants.file <- paste("/media/bigboy/shared_data/public/SNP/1000_genomes_",
                       input.ethno, "_20100804.genotypes.vcf.gz", sep="")
#region.size <- Iranges(start=start, end=end)

##gr is the genomic range where you want to look for SNPs, seqnames="chrom"
SNP.range <- GRanges(seqnames=as.character(c(8)), 
                     ranges=IRanges(start=128400000, 
                                    end=128900000))

subset.variants.file <- unpackVcf(scanVcf(variants.file, param=SNP.range), 
                                  variants.file)
subset.variants.file <- subset.variants.file[[1]]
variants.data <- cbind(
                       as.data.frame(
                                     cbind(subset.variants.file$CHROM, 
                                           subset.variants.file$POS, 
                                           subset.variants.file$ID, 
                                           subset.variants.file$REF, 
                                           subset.variants.file$ALT, 
                                           subset.variants.file$INFO$AF)), 
                       as.data.frame(subset.variants.file$GENO$GT))
variants.data.length <- length(variants.data)
variants.data[, 7:variants.data.length] <- 
    sapply(variants.data[, 7:variants.data.length], as.character)
colnames(variants.data)[1:6] <- c("CHROM", "POS", "ID", "REF", "ALT", "AF")
variants.data[, 7:variants.data.length] <- 
    ifelse(variants.data[, 7:variants.data.length] == "0|0", 
           paste(variants.data$REF, "/", variants.data$REF, sep=""),
           ifelse(variants.data[, 7:variants.data.length] == "0|1",:w
                  paste(variants.data$REF, "/", variants.data$ALT, sep=""),
                  ifelse(variants.data[, 7:variants.data.length] == "1|0",
                         paste(variants.data$ALT, "/", variants.data$REF, sep=""),
                         ifelse(variants.data[, 7:variants.data.length] == "1|1",
                                paste(variants.data$ALT, "/", variants.data$ALT
                                      , sep=""),
                                NA))))

variants.data$ID <- 
    ifelse(variants.data$ID == ".", 
           paste("POS:", variants.data$CHROM, ":", variants.data$POS, sep=""),
           variants.data$ID)
snp.ranges <- 
    GRanges(seqnames=paste(
                           "chr",
                           as.character(variants.data$CHROM), 
                           sep=""), 
            ranges=(IRanges(
                            start=as.integer(variants.data$POS), 
                            width=1)))
names(snp.ranges) <- variants.data$ID
bedfile <- import(features.file, asRangedData = FALSE)
snps.included <- snp.ranges %in% bedfile
snps.included <- subset(variants.data, snps.included)
snp.geno <- data.frame(t(snps.included[, c(7:variants.data.length)]))
colnames(snp.geno) <- snps.included$ID
snp.geno <- snp.geno[, colSums(is.na(snp.geno))<nrow(snp.test)] #removes NA cols

#overlapping.snps <- subsetByOverlaps(snp.ranges, bedfile) # maxgap is an option
                                                          # that allows for an
                                                          # interval in the bed
                                                          # file to be made
                                                          # larger by by maxgap
                                                          # on either side of
                                                          # the interval


snp <- t(variants.data[, c(3,7:variants.data.length)]) ## creates the data format for LDsnp (mostly; there is a problem with unnamed snps)
