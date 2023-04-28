library(phyloseq)
library(ggplot2)
library(openxlsx)
library(dplyr)
library(vegan)
library(data.table)
library(microbiome)

#Map file
# Use encoding = UTF-8 in argument.
map = import_qiime_sample_data(mapfilename = 'Grape_mapping2.txt')
map[["Patient"]] = factor(map[["Patient"]])

#ASV counts file
ASV_count_only = read.csv("Grape_counts.csv", header=T, row.names=1)
ASV_count_only = as.matrix(ASV_count_only)
ASV_count_only = otu_table(ASV_count_only, taxa_are_rows = TRUE)

#Rank file
tax = read.csv("OTU Rank.csv", header=T, row.names=1)
# FOR PATHWAY ANALYSIS: don't need tax file
tax = as.matrix(tax)
tax_final = tax_table(tax)

#Tree file
tree = read_tree("GrapeStudy_Aug2021_tree.tree")

#Combine and filter
ASVphylo = merge_phyloseq(tax_final,ASV_count_only,map, tree)
ASVphyloF10 = filter_taxa(ASVphylo, function(x) sum(x > 3) > (0.1*length(x)), TRUE)
ASVphyloF10.rarefied = rarefy_even_depth(ASVphyloF10, rngseed=1, sample.size=0.9*min(sample_sums(ASVphyloF10)), replace=F)
rarecurve(t(otu_table(ASVphylo)),step=50,cex=0.5)

#Plot abundance
plot_bar(ASVphyloF10.rarefied, fill="phylum") + facet_wrap(~Group, scales = "free_x", nrow = 1)

#Plot alpha
alpha_meas = c("Chao1", "shannon")
alphalabels = c("Baseline", "Intervention")
p = plot_richness(ASVphyloF10.rarefied, "Group", measures=alpha_meas)
p + geom_boxplot(data=p$ASVphyloF10.rarefied, aes(x=Group, y=value, color=NULL), alpha=0.1) +
theme(axis.text.x = element_text(angle = 90, hjust = 1, size=11,color="Black", face="bold"))+
scale_x_discrete(name=" ", labels= alphalabels)+
theme(panel.border = element_blank(), panel.grid.minor=element_blank(), panel.background = element_blank(), panel.grid.major = element_blank(), axis.line = element_line(color = 'black'))
alpha_div<-estimate_richness(ASVphyloF10, measures=c("Shannon", "Chao1"))
write.csv (as.data.frame(alpha_div),file="Grape_alpha_diversity.csv")

#Plot beta, repeatable for different methods
uwunifrac_dist = phyloseq::distance(ASVphyloF10, method='unifrac', weighted=F)
ordination_uw_uf = ordinate(ASVphyloF10,method='PCoA', distance = uwunifrac_dist)
plot_ordination(ASVphyloF10,ordination_bray, color='Group') + geom_point(size=4.5) + 
scale_color_manual(name='Group', labels = alphalabels, values=c('Blue', 'Orange')) + ggtitle("PCoA: Bray-Curtis") +
theme_bw() + theme(axis.text.x=element_text(angle=0, hjust=1,size=9,color='darkred')) +
theme(axis.text.y=element_text(angle=0, hjust=1,size=9,color='Black'))+ 
theme(axis.title=element_text(size=10,color='Black',face='bold')) +
theme(legend.title=element_text(color='white',face='bold',size=0)) +
theme(legend.text=element_text(color='Black', size=10,face='bold'))+
theme(legend.position=c(0.8,0.9))

#Get beta statistic
betastat = as(sample_data(ASVphyloF10), "data.frame")
adonis(uwunifrac_dist~Group, data=betastat)

#DEseq2, repeatable for multiple tax ranks
library(DESeq2)
alpha = 0.05
gm_mean = function(x, na.rm=TRUE){+ exp(sum(log(x[x>10]),na.rm=na.rm)/length(x))}
genusASVF10 = tax_glom(ASVphyloF10, taxrank = 'genus')
phylumASVF10 = tax_glom(ASVphyloF10, taxrank = 'phylum')
diagadds_gen = phyloseq_to_deseq2(genusASVF10, ~Patient+Group) 
##If paired: control patient and ask for group by doing phyloseq_to_deseq2(genusASVF10, ~Patient+Group)
geomeans_gen = apply(counts(diagadds_gen),1,gm_mean)
diagadds_gen = estimateSizeFactors(diagadds_gen, geoMeans = geomeans_gen)
diagadds_gen = DESeq(diagadds_gen,fitType='local')
resTT_gen = results(diagadds_gen)
sigtab_gen = resTT_gen[which(resTT_gen$pvalue < alpha),]
sigtab_gen_bind = cbind(as(sigtab_gen,"data.frame"), as(tax_table(ASVphyloF10)[rownames(sigtab_gen),],"matrix"))
write.csv(sigtab_gen_bind, "sigtab_gentest0827.csv")

#Plot DEseq2, repeatable for multiple tax ranks
taxforplot = tax_table(ASVphyloF10)[rownames(sigtab_gen),]
#Genus = 6
sigtab_gen$stars = cut(sigtab_gen$pvalue, c(-Inf,0.001,0.01,0.05,Inf), label=c("***", "**", "*", ""))
genusforplot = taxforplot@.Data[,6]
x_gen = tapply(sigtab_gen$log2FoldChange, genusforplot, function(x) max(x))
sigtab_gen$Genus = factor(as.character(genusforplot), levels=names(x_gen))
#Phylum = 2
phylumforplot = taxforplot@.Data[,2]
x_phy = tapply(sigtab_gen$log2FoldChange, phylumforplot, function(x) max(x))
sigtab_gen$Phylum = factor(as.character(phylumforplot), levels=names(x_phy))
sigtab_gen_df = as(sigtab_gen, "data.frame")

ggplot(sigtab_gen_df, aes(x=Genus, y=log2FoldChange, color=Phylum)) + geom_point(size=4) +
geom_text(aes(label=stars), color='black',size=4) +
theme(legend.text=element_text(color='Black', size=8)) +
theme(legend.title=element_text(color='black',face='bold',size=14)) +
theme(axis.text.y = element_text(angle=0,hjust=1,size=8,face='bold',color='Black'))+
theme(axis.title=element_text(size=14,color='Black', face='bold')) +
theme(axis.text.x=element_text(angle=-90,hjust=0,size=8,face='bold',vjust=0.5),panel.grid.major=element_blank(),panel.grid.minor=element_blank(),panel.border=element_rect(color='black', fill=NA,size=1),panel.background=element_blank())

#Extract and plot relative abundances
phylumtable = abundances(phylumASVF10, "compositional")
taxonomy = tax_table(phylumASVF10)
write.csv(phylumtable, "phylum_table.csv")
write.csv(taxonomy, "taxonomy.csv")
phylumASVF10_RA = transform_sample_counts(phylumASVF10, function(x) x/sum(x))
plot_bar(phylumASVF10_RA, fill="phylum") + ggtitle("Relative Abundance %") +
  scale_fill_brewer(palette="Paired") + theme_bw() + facet_wrap(~Group, scales="free_x", nrow=1)
