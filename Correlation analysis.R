library(cocor)
library(openxlsx)
library(Hmisc)
library(ggplot2)
library(psych)
library(reshape2)
library(readr)
library(corrplot)
library(scales)

# Get correlation and p values
rel = read.csv("Grape_phylum_RA_difference.csv", row.names = 1)
bac = read.csv("Grape_bileacid_change.csv", row.names = 1)
selcor= cor(rel, bac, method="spearman", use = "pairwise.complete.obs")
co = reshape2::melt(selcor)
write.csv(co, "corrcoef3.csv")
cortest = corr.test(rel, bac, method = "spearman", ci=FALSE, adjust="none")
familyp = reshape2::melt(cortest$p)
write.csv(familyp, "corrcoef3_p.csv")

# Generate heatmap AFTER COMBINING the correlation and p values
cordata = read.csv("corrcoef_spearman.csv")

ggplot(cordata, aes(Var2, Var1))+geom_tile(aes(fill = value)) +
geom_text(aes(fill = value, label = round(value, 2)), size = 2) + 
scale_fill_gradient2(low = muted("darkred"), mid = "white", high = muted("midnightblue"), midpoint = 0) + 
theme(panel.grid.major.x=element_blank(),panel.grid.minor.x=element_blank(),panel.grid.major.y=element_blank(), panel.grid.minor.y=element_blank(),panel.background=element_rect(fill="white"), axis.text.x = element_text(angle=90, hjust = 1,vjust=1,size = 9),plot.title = element_text(size=16,face="bold"),axis.text.y = element_text(size = 9)) +
ggtitle("Correlation Plot") +
theme(legend.title=element_text(size=12)) +
scale_x_discrete(name="")+scale_y_discrete(name="") +
labs(fill="Corr. Coef.", size = 1)  

cordata$Var2 = as.character(cordata$Var2)
cordata$Var2 = factor(cordata$Var2, levels=unique(cordata$Var2))
cordata$stars = cut(cordata$p, breaks=c(-Inf, 0.001, 0.01, 0.05, Inf), label=c("***", "**", "*", ""))
p = ggplot(aes(Var2, Var1, fill=value), data=cordata)
SHFig = p + geom_tile() + scale_fill_gradient2(low="#D7191C", mid="white", high="#2C7BB6")+
geom_text(aes(label=stars), color="black", size=5) + 
labs(y=NULL, x=NULL, fill="Corr. Coef.") + theme_bw() +
theme(axis.text.x=element_text(angle = 90, hjust = 1, vjust=1,size = 9),plot.title = element_text(size=20),axis.text.y = element_text(size = 9)) +
ggtitle("Spearman Correlation")

ggsave("Grape_bileacid_RA_spearman.png")
