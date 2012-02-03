## reads all LD files in haploview_results for each biological feature/riskSNP.
## need to be in directory where chrN.rsXXXX.XX:XXXX-XXXX is located.

dat <- vector("list", length(list.files()));
data <- vector("list", length(list.files()));
for(i in 1:length(list.files())){ ## for each tagSNP
	tmp <- list.files(paste(list.files()[i],"/haploview_results/",sep=""), pattern="*.LD$")
	#names(dat)[i] <- paste(unlist(strsplit(tmp[i], ".", fixed = TRUE))[4],unlist(strsplit(tmp[i], ".", fixed = TRUE))[2],sep="_")
	#names(data)[i] <- paste(unlist(strsplit(tmp[i], ".", fixed = TRUE))[4],unlist(strsplit(tmp[i], ".", fixed = TRUE))[2],sep="_")	
	dat[[i]] <- vector("list", length(tmp));	
		for(j in 1:length(tmp)){ ## for each data within a tagSNP
		

		
			tmp.name <- unlist(strsplit(tmp[j], ".", fixed = TRUE))
		
			#names(dat[[i]])[j] <- tmp[j]
			dat[[i]][[j]] <- read.delim(file=paste(list.files()[i],"/haploview_results/",tmp.name[1],".",tmp.name[2],".",tmp.name[3],".",tmp.name[4],".",tmp.name[5],".LD", sep=""), sep="\t", check.names=T);
				if(dim(dat[[i]][[j]])[1] == 0){
					dat[[i]][[j]][1,] <- NA
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

save(riskSNP, file="riskSNP_data.rda")
