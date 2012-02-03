library(ggplot2)
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

##

load(file="riskSNP.rda")

all <- rbind(
	subset(riskSNP,TagTRUE1==TRUE),
	subset(riskSNP,TagTRUE2==TRUE)
)

## plot r.2 values
p.all <- ggplot(all, aes(x=r.2))
p.all + 
	geom_histogram() + 
	geom_vline(xintercept = 0.5, linetype=2) +
	scale_x_continuous("Rsquare Values (0-1)") + 
	scale_y_continuous("Total # of Surrogate SNPs associated with riskSNP") + 
	opts(legend.position = "none", axis.text.y = theme_text(), axis.text.x = theme_text(angle=90), title = "riskSNP associated with Prostate Cancer diseases\nOverlapping") + 
	facet_wrap(~ tagSNP)
ggsave(file="plots/Prostate_Cancer_R2summary_riskSNP_ver1.pdf")

## plot r.2 vs. distance values
p.all.d <- ggplot(all, aes(x=r.2, y=Dist))
p.all.d + 
	geom_point() + 
	geom_vline(xintercept = 0.5, linetype=2) +
	#geom_abline(intercept = 0, slope = 1) +
	scale_x_continuous("Rsquare Values (0-1)") + 
	scale_y_continuous("Distance to Surrogate SNPs associated with riskSNP (bp)") + 
	opts(legend.position = "none", axis.text.y = theme_text(), axis.text.x = theme_text(angle=90), title = "Distance between riskSNP associated with Prostate Cancer diseases\nand Surrogate SNP\nOverlapping") + 
	facet_wrap(~ tagSNP)
ggsave(file="plots/Prostate_Cancer_MCF7_allpeaks_R2vsDist_riskSNP_ver1.pdf")

