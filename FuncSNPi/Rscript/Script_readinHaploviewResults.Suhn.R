## R script to read in all Haploview data into a datamatrix for downstream analysis.
#starting in this directory "/home/houtana/Documents/SNP_Analysis/PCa_analysis/Ver4/results/hg19/2MB"
## change the following to match your own path. then copy from here to there
setwd("/home/houtana/Documents/SNP_Analysis/PCa_analysis/Ver5/results/hg19/2MB")



## copy from here ###
## copy from here ###
## copy from here ###
dat <- vector("list", length(list.files()));
data <- vector("list", length(list.files()));
for(i in 1:length(list.files())){ ## for each tagSNP
	tmp <- list.files(paste(list.files()[i],"/haploview_results/",sep=""), pattern=".LD$")
	#names(dat)[i] <- paste(unlist(strsplit(tmp[i], ".", fixed = TRUE))[4],unlist(strsplit(tmp[i], ".", fixed = TRUE))[2],sep="_")
	#names(data)[i] <- paste(unlist(strsplit(tmp[i], ".", fixed = TRUE))[4],unlist(strsplit(tmp[i], ".", fixed = TRUE))[2],sep="_")	
	dat[[i]] <- vector("list", length(tmp));	
		for(j in 1:length(tmp)){ ## for each data within a tagSNP
		

		
			tmp.name <- unlist(strsplit(tmp[j], ".", fixed = TRUE))
		    print(tmp.name);
			#names(dat[[i]])[j] <- tmp[j]
			dat[[i]][[j]] <- read.delim(file=paste(list.files()[i],"/haploview_results/",tmp[j],sep=""),sep="\t", check.names=T);
            print(dat[[i]][[j]]);
				if(dim(dat[[i]][[j]])[1] == 0){
					dat[[i]][[j]][1,] <- NA
				}
        }
}
			
			dat[[i]][[j]]$tagSNP <- tmp.name[2]
			dat[[i]][[j]]$bioFeature <- tmp.name[5]
			dat[[i]][[j]]$location <- tmp.name[4]
			dat[[i]][[j]]$ethnicity <- tmp.name[3]
			dat[[i]][[j]]$genome1000 <- tmp.name[1]
			
		        #cat("Time: ", date(), " - Finished loading data: for ", paste(tmp.name[1],".",tmp.name[2],".",tmp.name[3],".",tmp.name[4],".",tmp.name[5],".",tmp.name[6], sep=""), "\n");
		}## for each biological tagSNP
		
		data[[i]] <- rbind(dat[[i]][[1]], dat[[i]][[2]]);
		
	    for(h in 3:length(tmp)){
		data[[i]] <- rbind(data[[i]], dat[[i]][[h]]);
	    }
    cat("Finished riskSNP: ", list.files()[i], " - Time: ", date(), "\n");
} ## for each tagSNP

riskSNP <- rbind(data[[1]], data[[2]]);
		
    for(k in 3:length(list.files())){
	riskSNP <- rbind(riskSNP, data[[k]]);
	cat(k,", ");
   }
   
   riskSNP$TagTRUE1 <- c(as.character(riskSNP[,1])==as.character(riskSNP[,10]))
   riskSNP$TagTRUE2 <- c(as.character(riskSNP[,2])==as.character(riskSNP[,10]))
## copy above ##   
## copy above ##   
## copy above ##   
## copy above ##   
## copy above ##
   save(riskSNP, file="/home/houtana/Documents/SNP_Analysis/PCa_analysis/Ver5/results/riskSNP_data.rda")
   
      setwd("/home/houtana/Documents/SNP_Analysis/PCa_analysis/Ver5/results/")
   load(file="/home/houtana/Documents/SNP_Analysis/PCa_analysis/Ver5/results/riskSNP_data.rda")
   ## run the following to find out which riskSNP has a correlation values more than 0.5
(subset(riskSNP, TagTRUE2==TRUE & r.2>0.5))
(subset(riskSNP, TagTRUE1==TRUE & r.2>0.5))

val<- read.delim(file="/home/houtana/Documents/SNP_Analysis/PCa_analysis/data/Enhancer_Clone_Assay_JohnLai_2009.txt", sep="\t")
riskSNP$ID <- c("byHoutan")
val$ID <- c("byJohn")

dat <- merge(riskSNP, val, by.x="tagSNP", by.y = "riskID",all = T)

###############
###############plots for Genomic Annotations (TSS, Exons, introns, etc.)
###############

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

c.r <- ggplot(riskSNP, aes(y=value, x=factor(genomic_type), fill=variable));
c.r + geom_bar(position=dodge) + scale_fill_manual(values=c("Peaks" = "RED", "Background" = "darkgrey")) + opts(axis.text.x = theme_text(angle=90), title = "SNP Characterization") + scale_y_continuous("Raw Counts") + facet_grid(type ~ Class.s, scales="free")


ggsave(file="results/genomicAnnotations_Raw_ver2.pdf")

