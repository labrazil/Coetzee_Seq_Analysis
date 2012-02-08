h3k4me2 <- import("~/Desktop/H3K4me2_sample.bed", asRangedData=F)
faire <- import("~/Desktop/faire.bed", asRangedData=F)
elementMetadata(h3k4me2)[, "name"] <- as.numeric(values(h3k4me2)$name)
elementMetadata(faire)[, "name"] <- as.numeric(values(faire)$name)
overlap <- subsetByOverlaps(faire, h3k4me2, maxgap=300)
strong <- overlap[values(overlap)$name > 6]
weak <- overlap[values(overlap)$name <= 8]
