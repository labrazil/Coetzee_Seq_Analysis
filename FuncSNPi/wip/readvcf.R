#!/usr/bin/env R

#### FunciSNP - the Functional Integration of SNPs with Phenotype
####           by Coincidence with Chromatin Biofeatures
#### Copyright (C) 2011 Simon G. Coetzee and Houtan Noushmehr
####
#### This program is free software: you can redistribute it and/or modify
#### it under the terms of the GNU General Public License as published by
#### the Free Software Foundation, either version 3 of the License, or
#### (at your option) any later version.
####
#### This program is distributed in the hope that it will be useful,
#### but WITHOUT ANY WARRANTY; without even the implied warranty of
#### MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#### GNU General Public License for more details.
####
#### You should have received a copy of the GNU General Public License
#### along with this program.  If not, see <http://www.gnu.org/licenses/>.

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
    snp.regions <- data.frame(snp.chromosome, 
                              snp.region.start, 
                              snp.region.end, 
                              snp.name, 
                              snp.ethno, 
                              stringsAsFactors = FALSE)
    snp.regions
}


ServerCheck <- function(primary.server) {
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
    cat("trying", primary.server.name, "as 1000 genomes server\n")
    server.up <- try(url(paste(primary.server, test.file, sep=""), open = "rt"), silent = TRUE)
    server.error <- inherits(server.up, "try-error")
    if(server.error) {
        cat(primary.server.name, "failed \ntrying",  secondary.server.name, "as 1000 genomes server\n")
        rm(server.error)
        server.up <- try(url(paste(secondary.server, test.file, sep=""), open = "rt"), silent = TRUE)
        server.error <- inherits(server.up, "try-error")
        if(server.error) {
            cat(secondary.server.name, "failed \n")
            close(server.up)
            stop("Neither EBI nor NCBI mirrors for the 1000 genomes project were found")
        } else {
            cat("OK using", secondary.server.name, ":", secondary.server, "\n")
            close(server.up)
            secondary.server
        }
    } else {
        cat("OK using", primary.server.name, ":", primary.server, "\n")
        close(server.up)
        primary.server
    }
}







FunciSNP <- function(ethno = c("AFR", "ASN", "EUR", "AMR", "ALL"), 
                     bio.features.loc, snp.regions.file, R.squared.cutoff = 0, primary.server, make.plots = TRUE) {
    cat("
        ####################################
        ##                                ##
        ##      Welcome to FunciSNP       ##
        ##                                ##
        ####################################
        args used:\n
        snp.regions.file:   ", as.character(snp.regions.file), "\n",
        "       ethno:            ", ethno, "\n",
        "       bio.features.loc: ", list.files(bio.features.loc, pattern="*.bed$", full.names = FALSE), "\n",
        "       Make Plots?:      ", make.plots, "\n")
        if(primary.server == "ebi" ||  primary.server == "ncbi") {
            cat("you have selected", primary.server, "as your primary server \n")
        } else {
            stop("please select ebi or ncbi as your primary server")
        }
        if(R.squared.cutoff < 1 || R.squared.cutoff > 0) {
            cat("you have selected" , R.squared.cutoff, "as your prefered R squared cut off used for summary reports and plots \n")
        } else {
            stop("please select an R squared cutoff between 0-1")
        }
    cat("\n ### WARNING:: existing vcf.gz files will be removed from", getwd(), "/funcisnp/",
        "\n ### and new ones will be placed there, press
 ### <Esc> or <Ctrl-C> to cancel if you do not want this to happen after the
 ### 10 second countdown:", "\n")
        timer <- 10
        pb <- txtProgressBar(min = 0, max = timer, style = 2)
        for(i in 1:timer) {
            Sys.sleep(1.0)
            setTxtProgressBar(pb, i)
        }
        close(pb)

        ### Not sure if we want to make this silent or error out if permission is denied. I use plot folder to store graphs.
        
        try(dir.create(path="funcisnp_results", showWarnings = FALSE), silent=TRUE) ## directory used to store plots
        try(dir.create(path="funcisnp_results/plots", showWarnings = FALSE), silent=TRUE) ## directory used to store plots
        try(dir.create(path="funcisnp_results/tmp", showWarnings = FALSE), silent=TRUE) ## directory used to store temp files
        try(dir.create(path="funcisnp_results/log", showWarnings = FALSE), silent=TRUE) ## directory used to store log files
        try(dir.create(path="funcisnp_results/tables", showWarnings = FALSE),silent=TRUE) ## directory used to store summary tables

        snp.region <- ReadRegionsFile(snp.regions.file)
        snp.names <- snp.region$snp.name
        system("find funcisnp_results/ -iname \\*vcf.gz\\* | xargs rm", ignore.stdout = TRUE, ignore.stderr = TRUE)

        onek.genome.server <- ServerCheck(primary.server)
        manifest.read <- read.delim(paste(onek.genome.server, "ftp/release/20110521/phase1_integrated_calls.20101123.ALL.panel", sep=""), sep="\t", header = FALSE)
        bio.features.file <<- list.files(bio.features.loc, pattern="*.bed$", full.names = TRUE)
        bio.features <<- list.files(bio.features.loc, pattern="*.bed$", full.names = FALSE)

        snp.ld.frame <- NULL
        snp.names.count <- 0

        for(j in as.character(snp.names)) {
            snp.names.count <- snp.names.count + 1
            cat(paste("\n#####\n","Pulling In Variants File for SNP: ", j, "\n", sep=""))
            for(i in as.character(ethno)) {
                variants.data <- PullInVariants(i, j, snp.region, primary.server, manifest.read)
                cat(paste("subsetting Variants File for racial/ethnic group: ", i, "\n", sep="")) 
                bio.features.count <- 0
                for(h in bio.features.file) {
                    bio.features.count <- bio.features.count + 1
                    cat(paste("time: ", format(Sys.time(), "%H:%M:%S"), "\tRisk SNP: ", j, ", ", i, " (", snp.names.count, "/", length(snp.names), ")\n>>>>> (", bio.features.count, "/", length(bio.features), ")", "\t", bio.features[bio.features.count], "\n", sep=""))
                    snp.geno <- FilterByFeatures(h, j, i, variants.data)
                    snp.ld.frame <- LDTesting(snp.geno, j, i, snp.ld.frame, bio.features.count, snp.region, variants.data)
                    try(save(snp.ld.frame, file="funcisnp_results/snp_table.Rda"), silent = TRUE)
                    try(write.table(snp.ld.frame, file="funcisnp_results/snp_table.txt", sep="\t", quote = FALSE, row.names = FALSE), silent = TRUE)
                }
            }
            file.remove(paste("funcisnp_results/", j, ".vcf.gz", sep=""))
            file.remove(paste("funcisnp_results/", j, ".vcf.gz.tbi", sep=""))
            summary(snp.ld.frame)
        }
        write.table(snp.ld.frame, file="funcisnp_results/snp_table.txt", sep="\t", quote = FALSE, row.names = FALSE)
        if(make.plots == TRUE) {
                snp.ld.frame <<- snp.ld.frame
                FunciSNPSummary(R.squared.cutoff, snp.ld.frame)
                FunciSNPPlot(R.squared.cutoff, snp.ld.frame)
                FunciSNPHeatmap(R.squared.cutoff, snp.ld.frame)
                cat("##########################################\nFunciSNP is complete!\nSee 'funcisnp_results/tables/' and 'funcisnp_results/plots/' folder for overall results\nSee 'funcisnp_results/snp_table' for complete results.\n##########################################\n")
        } else {
            FunciSNPSummary(R.2=R.squared.cutoff, dat = snp.ld.frame)
            cat("##########################################\nFunciSNP is complete!\nSee 'funcisnp_results/tables/' folder for overall results\nSee 'funcisnp_results/snp_table' for complete results.\n##########################################")
        }
}


PullInVariants <- function(ethno.chosen, snp.names.chosen, snp.regions, primary.server, manifest) {
    ifelse(ethno.chosen != "ALL",
           ethno.sample.set <- as.character(subset(manifest, V3==ethno.chosen)[,1]),
           ethno.sample.set <- as.character(manifest[,1]))
    intermediate.vcf <- paste("funcisnp_results/", snp.names.chosen, ".vcf", sep="")
    variants.file <- paste(intermediate.vcf, ".gz", sep="") 
    if(file.exists(variants.file) == FALSE) {
        onek.genome.server <- ServerCheck(primary.server)
        system(paste("tabix -hf ", onek.genome.server, "/ftp/release/20101123/interim_phase1_release/ALL.chr", snp.regions$snp.chromosome[snp.regions$snp.name == snp.names.chosen][1], ".phase1.projectConsensus.genotypes.vcf.gz", " ", snp.regions$snp.chromosome[snp.regions$snp.name == snp.names.chosen][1], ":", snp.regions$snp.region.start[snp.regions$snp.name == snp.names.chosen][1], "-", snp.regions$snp.region.end[snp.regions$snp.name == snp.names.chosen][1], " > ", intermediate.vcf, " && ", "bgzip -f ", intermediate.vcf, " && ", "tabix -hf ", intermediate.vcf, ".gz" , sep=""))
        subset.variants.file <<- unpackVcf(scanVcf(variants.file))[[1]]
    }
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
    variants.data
}

FilterByFeatures <- function(features.file, snp.names.chosen, ethno.chosen, variants.data) {
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
    if(sum(snps.included) >= 1) {
        snps.included <- subset(variants.data, snps.included)
        snp.geno <- data.frame(t(snps.included[, c(6:length(snps.included))]))
        snp.geno <- as.data.frame(snp.geno[, colSums(is.na(snp.geno))<nrow(snp.geno)])
        colnames(snp.geno) <- snps.included$ID
        if(snp.names.chosen %in% colnames(snp.geno)) {
            cat(snp.names.chosen, ", the risk snp, already overlaps with the feature", features.file, "\n");
        } else {
            temp <- data.frame(t(subset(variants.data, ID==snp.names.chosen))[6:length(variants.data), ])
            colnames(temp) <- snp.names.chosen
            snp.geno <- cbind(snp.geno, temp)
            rm(temp)
        }
        for(i in 1:ncol(snp.geno)) {
            snp.geno[, i] <- genotype(snp.geno[, i])
        }
        close.snp.ranges <<- close.snp.ranges
        snp.geno
    } else {
        ## need to put the following in log files
        #        cat(paste("There is no overlap for: \n",
        #                  "\tRisk SNP: \t\t", snp.names.chosen, "\n", 
        #                  "\tbiofeature: \t\t", features.file, "\n", 
        #                  "\tracial/ethnic group: \t", ethno.chosen, "\n", sep=""));
        NULL
        close.snp.ranges <<- close.snp.ranges
    }
}

LDTesting <- function(snps, snp.names.chosen, ethno.chosen, snp.ld.frame, bf.count, snp.regions, var.data) {
    try(snp.ld <- LD(snps), silent = TRUE)
    if(exists("snp.ld")) {
        if(snp.names.chosen %in% colnames(snp.ld$"D")) {
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
            if(identical(snp.ld.frame, NULL)) {
                snp.ld.frame <- data.frame(snp.names.chosen, snp.names.cor, d, D.prime, r, R.squared, n, chi.squared, p.val)
                snp.ld.frame$feature <- unlist(strsplit(bio.features[bf.count], "\\.bed")) 
                snp.ld.frame$loc.feature <- paste(seqnames(features.file.interval[c(as.matrix(findOverlaps(features.file.interval, close.snp.ranges[snp.names.cor]))[,1])]), ":", start(features.file.interval[c(as.matrix(findOverlaps(features.file.interval, close.snp.ranges[snp.names.cor]))[,1])]), "-", end(features.file.interval[c(as.matrix(findOverlaps(features.file.interval, close.snp.ranges[snp.names.cor]))[,1])]), sep="")
                snp.ld.frame$eth.risk <- snp.regions$snp.ethno[snp.regions$snp.name == snp.names.chosen][1]
                snp.ld.frame$eth.chosen <- ethno.chosen
                snp.ld.frame$chr <- snp.regions$snp.chromosome[snp.regions$snp.name == snp.names.chosen][1]
                snp.ld.frame$loc.risk <- subset(var.data, ID==snp.names.chosen)$POS
                snp.ld.frame$loc.cor <- subset(var.data, var.data$ID %in% snp.names.cor)$POS
                snp.ld.frame$dist <- abs(snp.ld.frame$loc.risk - snp.ld.frame$loc.cor)
            } else {
                tmp.snp.ld.frame <- data.frame(snp.names.chosen, snp.names.cor, d, D.prime, r, R.squared, n, chi.squared, p.val)
                tmp.snp.ld.frame$feature <- unlist(strsplit(bio.features[bf.count], "\\.bed")) 
                tmp.snp.ld.frame$loc.feature <- paste(seqnames(features.file.interval[c(as.matrix(findOverlaps(features.file.interval, close.snp.ranges[snp.names.cor]))[,1])]), ":", start(features.file.interval[c(as.matrix(findOverlaps(features.file.interval, close.snp.ranges[snp.names.cor]))[,1])]), "-", end(features.file.interval[c(as.matrix(findOverlaps(features.file.interval, close.snp.ranges[snp.names.cor]))[,1])]), sep="")
                tmp.snp.ld.frame$eth.risk <- snp.regions$snp.ethno[snp.regions$snp.name == snp.names.chosen][1]
                tmp.snp.ld.frame$eth.chosen <- ethno.chosen
                tmp.snp.ld.frame$chr <- snp.regions$snp.chromosome[snp.regions$snp.name == snp.names.chosen][1]
                tmp.snp.ld.frame$loc.risk <- subset(var.data, ID==snp.names.chosen)$POS
                tmp.snp.ld.frame$loc.cor <- subset(var.data, var.data$ID %in% snp.names.cor)$POS
                tmp.snp.ld.frame$dist <- abs(tmp.snp.ld.frame$loc.risk - tmp.snp.ld.frame$loc.cor)
                snp.ld.frame <- rbind(snp.ld.frame, tmp.snp.ld.frame)
            }
            rownames(snp.ld.frame) <- NULL
            snp.ld.frame
        } else {
            ## need to put the following in log files
            #cat("risk snp", snp.names.chosen, "is not variant within the", ethno.chosen, "group \n")
            snp.ld.frame
        }
    } else {
        ## need to put the following in log files
        #cat(paste("\n", "Fewer than two of the snps grouped with ", snp.names.chosen, ", that overlap with ", unlist(strsplit(bio.features[bf.count], "\\.bed")), " for the ", ethno.chosen, " racial/ethnic group, have more than one allele", "\n", sep=""))
        snp.ld.frame
    }
}



FunciSNPSummary <- function(R.2, dat) {

    total.tagSNPs <-length(unique(dat[,"snp.names.chosen"]))  
    total.1kSNPs  <-length(unique(dat[,"snp.names.cor"]))
    total.feature  <-length(unique(subset(dat, R.squared>R.2)[,"feature"]))   

    total.tagSNPs.cutoff <-length(unique(subset(dat, R.squared>R.2)[,"snp.names.chosen"]))
    total.1kSNPs.cutoff  <-length(unique(subset(dat, R.squared>R.2)[,"snp.names.cor"]))
    total.feature.cutoff  <-length(unique(subset(dat, R.squared>R.2)[,"feature"]))

    total.dat <- matrix(c(total.tagSNPs,total.1kSNPs,total.feature, total.tagSNPs.cutoff,total.1kSNPs.cutoff,total.feature.cutoff), nrow = 3, ncol=2, byrow=FALSE,
                        dimnames = list(c("tagSNPs", "1kSNPs","bio.features"),
                                        c("Total", paste("R.squared.cuff.",R.2,sep=""))))
    total.dat <- as.data.frame(total.dat)
    total.dat$Percent <- round((total.dat[,2]/total.dat[,1])*100,2)
    write.table(total.dat, file="funcisnp_results/tables/summary.txt", sep="\t", row.names=T, col.names=T, quote=F)
    return(total.dat);
}

FunciSNPPlot <- function(R.2, dat) {
    require("ggplot2")
    cat("\n\nPlotting will begin ")
    cat("using R square cut off of", R.2, "\n")

    ### ggplot2 plots#####
    theme_white <- function() {

        theme_update (
                      plot.background = theme_blank(),
                      panel.background=theme_rect(colour="black", size=1),
                      axis.text.x= theme_text(colour="black",vjust= 1, size=12),
                      axis.text.y= theme_text(colour="black",hjust=1, size=12),
                      axis.title.x =theme_text(colour="black",face="bold", size=12),
                      axis.title.y =theme_text(colour="black",face="bold", angle = 90, size=12)
                      )
    }
    theme_white()

    all.s <- try(subset(dat, R.squared >= R.2), silent = TRUE)
    all.ss <- try(subset(dat, R.squared < R.2), silent = TRUE)
    try(all.s$r2 <- c("Yes"), silent = TRUE)
    try(all.ss$r2 <- c("No"), silent = TRUE)
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
    for( i in 1:length(summary(as.factor(all[,"feature"]))) ){

        bio <- names(summary(as.factor(all[,"feature"])))
        ## plot r.2 values
        tmp <- subset(all, feature==bio[i])
        p.all <- ggplot(tmp, aes(x=R.squared, fill=factor(r2)))
        p.all + 
        geom_histogram() + 
        geom_vline(xintercept = R.2, linetype=2) +
        scale_x_continuous("Rsquare Values (0-1)", limits=c(0,1)) + 
        scale_y_continuous("Total # of Surrogate SNPs associated with riskSNP") + 
        scale_fill_manual(values = c("Yes" = "Red", "No" = "Black")) +
        opts(legend.position = "none", axis.text.y = theme_text(), axis.text.x = theme_text(angle=90), title = paste("riskSNP\nOverlapping: ", bio[i], sep="")) + 
        facet_wrap(~ snp.names.chosen)
        ggsave(file=paste("funcisnp_results/plots/",bio[i],"_R2summary_riskSNP.pdf",sep=""))

        ## plot r.2 vs. distance values
        p.all.d <- ggplot(tmp, aes(x=R.squared, y=dist, colour=r2, size=factor(r2)))
        p.all.d + 
        geom_point() + 
        geom_vline(xintercept = R.2, linetype=2) +
        #geom_abline(intercept = 0, slope = 1) +
        scale_x_continuous("Rsquare Values (0-1)", limits=c(0,1)) + 
        scale_y_continuous("Distance to Surrogate SNPs associated with riskSNP (bp)", formatter="comma") + 
        scale_colour_manual(values = c("Yes" = "Red", "No" = "Black")) +
        scale_size_manual(values = c("Yes" = 2, "No" = 1)) +
        opts(legend.position = "none", axis.text.y = theme_text(), axis.text.x = theme_text(angle=90), title = paste("Distance between riskSNP\nand Surrogate SNP\nOverlapping: ", bio[i], sep="")) + 
        facet_wrap(~ snp.names.chosen)
        ggsave(file=paste("funcisnp_results/plots/",bio[i],"_R2vsDist_riskSNP.pdf",sep=""))
        cat("Finished plotting ", i, "/",length(bio), "\n")
    }
    cat("\n\nMaking plots is finished, see 'funcisnp_results/plots/' folder in ", getwd(), "\n\n")
}
FunciSNPHeatmap <- function(R.2, dat) {
require("gplots")
require('matlab')
    cat("Heatmap generation will begin")
    cat(" using R square cut off of ", R.2,"\n")
    all.s<-(table( subset(dat,R.squared>=R.2)[,"feature"], subset(dat,R.squared>=R.2)[,"snp.names.chosen"] ))
    all.s <- as.matrix(all.s)
    png(filename="funcisnp_results/plots/heatmap.png", bg = "white")
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
              main = paste("Rsquare values >= ",R.2)
              #labRow=c("ACH_OV","AR","AR_nonOV","AR_OV", "Sahu_AR","Sahu_AR_siFoxA1", "Sahu_FoxA1", "Sahu_FoxA1_siFoxA1", "AR_DHT", "AR_DHT_siFoxA1", "FoxA1", "FoxA1_DHT", "K4me1","K4me1_DHT", "K4me1_siFoxA1","K4me1_DH!_siFoxA1","K27ac_DHT","K27ac_DHT_siFoxA1", "Dnase_DHT", "Dnase","Dnase_broad1","Dnase_broad2","Dnase_narrow1","Dnase_narrow2"),
              #labCol=NULL
              )
    dev.off()
    cat("\n\nHeatmap plot is complete, see 'funcisnp_results/plots/' folder in ", getwd(),"\n\n")
}


