###Load Packages
library(phyloseq)
library(ggplot2)
library(scales)
library(grid)
library(vegan)
library(RColorBrewer)
library(DESeq2)
library(reshape2)
library(dplyr)
library(rio)
library(data.table)
library(VennDiagram)
library(microbiome)
library(goeveg)
library(venneuler)

###Set working directory
setwd("~/Desktop/Nina_Collab/")
list.files()

##DADA2
load(file="Nina_Collab_DADA2.RData")
samdf <- read.table("Nina_Collab_mappingfile.txt")
colnames(samdf)<- c("SampleID","Animal","Litter","BirthCage","TreatmentCage","TreatmentPlan")
mappingfile <- samdf[,-1]
rownames(mappingfile) <- samdf[,1]
physeq <- phyloseq(otu_table(seqtab.nochim, taxa_are_rows=FALSE), 
                   sample_data(mappingfile), 
                   tax_table(taxa.plus))
library("ape")
random_tree = rtree(ntaxa(physeq), rooted=TRUE, tip.label=taxa_names(physeq))
physeq1 = merge_phyloseq(physeq, mappingfile, random_tree)
save.image(file="Nina_Collab_physeq1.RData")

###Load Data
load(file="Nina_Collab_physeq1.RData")
theme_set(theme_bw())

###Remove Singletons
physeq1_trim<-prune_taxa(taxa_sums(physeq1) > 1, physeq1)
physeq1_trim

###Reads per sample
Analysis<-physeq1_trim
sdt = data.table(as(sample_data(Analysis), "data.frame"),
                 TotalReads = sample_sums(Analysis), keep.rownames = TRUE)
setnames(sdt, "rn", "SampleID")
pSeqDepth = ggplot(sdt, aes(TotalReads)) + geom_histogram() + ggtitle("Sequencing Depth")
pSeqDepth
sdt

Analysis = subset_samples(Analysis, Animal != "F4-182F")

##Combined Weighted Unifrac PCoA
Analysis<-physeq1_trim
GP1=transform_sample_counts(Analysis, function(x) 1E6 * x/sum(x))
orduW = ordinate(GP1, "PCoA","unifrac",weighted=TRUE)
pW = plot_ordination(GP1, orduW, color = "TreatmentPlan", title = "Weighted Unifrac") + geom_point(size=6) + scale_colour_manual(values=c("blue","red"))
pW + theme(plot.title = element_text(size=18))+ theme(legend.text=element_text(size=14)) + theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank()) #+ylim(-0.1,0.1) +xlim(-0.2,0.0)

###Permanova by TreatmentPlan
set.seed(42)
analysis_unifrac_weighted<-phyloseq::distance(Analysis,method = "unifrac", weighted=TRUE)
sampledf<- data.frame(sample_data(Analysis))
adonis(analysis_unifrac_weighted ~ TreatmentPlan, data = sampledf)

##Unweighted Unifrac PCoA
GP1=transform_sample_counts(Analysis, function(x) 1E6 * x/sum(x))
orduU = ordinate(GP1, "PCoA","unifrac",weighted=FALSE)
pU = plot_ordination(GP1, orduU, color = "TreatmentPlan", title = "Unweighted Unifrac") + geom_point(size=6) + scale_colour_manual(values=c("blue","red"))
pU + theme(plot.title = element_text(size=18))+ theme(legend.text=element_text(size=14)) + theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank())  #+ ylim(-0.2,0.1) +xlim(-0.30,0.30) 

###Permanova by TreatmentPlan
set.seed(42)
analysis_unifrac_unweighted<-phyloseq::distance(Analysis,method = "unifrac", weighted=FALSE)
sampledf<- data.frame(sample_data(Analysis))
adonis(analysis_unifrac_unweighted ~ TreatmentPlan, data = sampledf)


##Differential Abundance between oil and tam

###Change ASV IDs
Analysis<-physeq1_trim
taxa_names(Analysis) <- paste0("Seq", seq(ntaxa(Analysis)))

###Condense by taxonomic level
ps.phylum<-tax_glom(Analysis,taxrank="Phylum")
ps.class<-tax_glom(Analysis,taxrank="Class")
ps.order<-tax_glom(Analysis,taxrank="Order")
ps.family<-tax_glom(Analysis,taxrank="Family")

###Differential Abundance by Phylum
dds = phyloseq_to_deseq2(ps.phylum, ~ TreatmentPlan)
dds = DESeq(dds, test="Wald", fitType="local")
res_WD = results(dds, cooksCutoff = FALSE)
alpha = 0.05
dds
sigtab_WD = res_WD[which(res_WD$padj < alpha), ]
diff.abund = cbind(as(sigtab_WD, "data.frame"), as(tax_table(ps.phylum)[rownames(sigtab_WD), ], "matrix"))
diff.abund
write.csv(diff.abund,"diff.abund.phylum.csv")

###Plot counts for ASVs determined differentially abundant
data <- plotCounts(dds, "Seq4446", intgroup=c("TreatmentPlan"), returnData=TRUE)
plot<-ggplot(data, aes(x=TreatmentPlan, y=count)) +
  scale_y_log10() + geom_point(size=3) + geom_line() + theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank())
print(plot + ggtitle("Phylum Proteobacteria"))

data <- plotCounts(dds, "Seq6100", intgroup=c("TreatmentPlan"), returnData=TRUE)
plot<-ggplot(data, aes(x=TreatmentPlan, y=count)) +
  scale_y_log10() + geom_point(size=3) + geom_line() + theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank())
print(plot + ggtitle("Phylum Bacteroidetes"))

###Differential Abundance by Class
dds = phyloseq_to_deseq2(ps.class, ~ TreatmentPlan)
dds = DESeq(dds, test="Wald", fitType="local")
res_WD = results(dds, cooksCutoff = FALSE)
alpha = 0.05
dds
sigtab_WD = res_WD[which(res_WD$padj < alpha), ]
diff.abund = cbind(as(sigtab_WD, "data.frame"), as(tax_table(ps.class)[rownames(sigtab_WD), ], "matrix"))
diff.abund
write.csv(diff.abund,"diff.abund.class.csv")

###Plot counts for ASVs determined differentially abundant
data <- plotCounts(dds, "Seq4446", intgroup=c("TreatmentPlan"), returnData=TRUE)
plot<-ggplot(data, aes(x=TreatmentPlan, y=count)) +
  scale_y_log10() + geom_point(size=3) + geom_line() + theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank())
print(plot + ggtitle("Class Betaproteobacteria"))

data <- plotCounts(dds, "Seq6100", intgroup=c("TreatmentPlan"), returnData=TRUE)
plot<-ggplot(data, aes(x=TreatmentPlan, y=count)) +
  scale_y_log10() + geom_point(size=3) + geom_line() + theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank())
print(plot + ggtitle("Class Bacteroidia"))

###Differential Abundance by Order
dds = phyloseq_to_deseq2(ps.order, ~ TreatmentPlan)
dds = DESeq(dds, test="Wald", fitType="local")
res_WD = results(dds, cooksCutoff = FALSE)
alpha = 0.05
dds
sigtab_WD = res_WD[which(res_WD$padj < alpha), ]
diff.abund = cbind(as(sigtab_WD, "data.frame"), as(tax_table(ps.order)[rownames(sigtab_WD), ], "matrix"))
diff.abund
write.csv(diff.abund,"diff.abund.order.csv")

###Plot counts for ASVs determined differentially abundant
data <- plotCounts(dds, "Seq4446", intgroup=c("TreatmentPlan"), returnData=TRUE)
plot<-ggplot(data, aes(x=TreatmentPlan, y=count)) +
  scale_y_log10() + geom_point(size=3) + geom_line() + theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank())
print(plot + ggtitle("Order Burkholderiales"))

data <- plotCounts(dds, "Seq6100", intgroup=c("TreatmentPlan"), returnData=TRUE)
plot<-ggplot(data, aes(x=TreatmentPlan, y=count)) +
  scale_y_log10() + geom_point(size=3) + geom_line() + theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank())
print(plot + ggtitle("Order Bacteroidales"))

###Differential Abundance by Family
dds = phyloseq_to_deseq2(ps.family, ~ TreatmentPlan)
dds = DESeq(dds, test="Wald", fitType="local")
res_WD = results(dds, cooksCutoff = FALSE)
alpha = 0.05
dds
sigtab_WD = res_WD[which(res_WD$padj < alpha), ]
diff.abund = cbind(as(sigtab_WD, "data.frame"), as(tax_table(ps.family)[rownames(sigtab_WD), ], "matrix"))
diff.abund
write.csv(diff.abund,"diff.abund.family.csv")

###Plot counts for ASVs determined differentially abundant
data <- plotCounts(dds, "Seq765", intgroup=c("TreatmentPlan"), returnData=TRUE)
plot<-ggplot(data, aes(x=TreatmentPlan, y=count)) +
  scale_y_log10() + geom_point(size=3) + geom_line() + theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank())
print(plot + ggtitle("Family Bacteroidaceae"))

data <- plotCounts(dds, "Seq3111", intgroup=c("TreatmentPlan"), returnData=TRUE)
plot<-ggplot(data, aes(x=TreatmentPlan, y=count)) +
  scale_y_log10() + geom_point(size=3) + geom_line() + theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank())
print(plot + ggtitle("Family Prevotellaceae"))

data <- plotCounts(dds, "Seq4446", intgroup=c("TreatmentPlan"), returnData=TRUE)
plot<-ggplot(data, aes(x=TreatmentPlan, y=count)) +
  scale_y_log10() + geom_point(size=3) + geom_line() + theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank())
print(plot + ggtitle("Family Sutterellaceae"))














##Rarefaction Curves
set.seed(42)
psdata<-Analysis

calculate_rarefaction_curves <- function(psdata, measures, depths) {
  require('plyr') # ldply
  require('reshape2') # melt
  
  estimate_rarified_richness <- function(psdata, measures, depth) {
    if(max(sample_sums(psdata)) < depth) return()
    psdata <- prune_samples(sample_sums(psdata) >= depth, psdata)
    
    rarified_psdata <- rarefy_even_depth(psdata, depth, verbose = FALSE)
    
    alpha_diversity <- estimate_richness(rarified_psdata, measures = measures)
    
    # as.matrix forces the use of melt.array, which includes the Sample names (rownames)
    molten_alpha_diversity <- melt(as.matrix(alpha_diversity), varnames = c('Sample', 'Measure'), value.name = 'Alpha_diversity')
    
    molten_alpha_diversity
  }
  
  names(depths) <- depths # this enables automatic addition of the Depth to the output by ldply
  rarefaction_curve_data <- ldply(depths, estimate_rarified_richness, psdata = psdata, measures = measures, .id = 'Depth', .progress = ifelse(interactive(), 'text', 'none'))
  
  # convert Depth from factor to numeric
  rarefaction_curve_data$Depth <- as.numeric(levels(rarefaction_curve_data$Depth))[rarefaction_curve_data$Depth]
  
  rarefaction_curve_data
}

rarefaction_curve_data <- calculate_rarefaction_curves(psdata, c('Observed'), rep(c(100,1000,2000,4000,6000,8000,10000,12000,14000,16000,18000,20000,22000), each = 50))

summary(rarefaction_curve_data)
###Summarize Alpha Diversity
rarefaction_curve_data_summary <- ddply(rarefaction_curve_data, c('Depth', 'Sample', 'Measure'), summarise, Alpha_diversity_mean = mean(Alpha_diversity), Alpha_diversity_sd = sd(Alpha_diversity))
###Add Sample Data
rarefaction_curve_data_summary_verbose <- merge(rarefaction_curve_data_summary, data.frame(sample_data(psdata)), by.x = 'Sample', by.y = 'row.names')
###Plot
ggplot(
  data = rarefaction_curve_data_summary_verbose,
  mapping = aes(
    x = Depth,
    y = Alpha_diversity_mean,
    ymin = Alpha_diversity_mean - Alpha_diversity_sd,
    ymax = Alpha_diversity_mean + Alpha_diversity_sd,
    colour = TreatmentPlan,
    group = Sample
  )
) + geom_line(
) + geom_pointrange(
) + facet_wrap(
  facets = ~ Measure,
  scales = 'free_y'
) + geom_point(size=1) + scale_colour_manual(values=c("blue","red")) + ggtitle("Alpha Rarefaction") + theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank()) 

#Bar plots
###Change ASV IDs
Analysis<-physeq1_trim
taxa_names(Analysis) <- paste0("Seq", seq(ntaxa(Analysis)))

###Normalize dataset
subset=transform_sample_counts(Analysis,function(x)x/sum(x))

###Condense by taxonomic level
ps.class<-tax_glom(subset,taxrank="Class")

###Subset Data
oil.class <- ps.class %>%
  subset_samples(TreatmentPlan == "oil")
oil.class

tam.class <- ps.class %>%
  subset_samples(TreatmentPlan == "tam")
tam.class


###Percent Abundance
otu.table<-otu_table(oil.class)
write.csv(otu.table,"oil.class.otu.table.csv")
tax.table<-tax_table(oil.class)
write.csv(tax.table,"oil.class.tax.table.csv")

otu.table<-otu_table(tam.class)
write.csv(otu.table,"tam.class.otu.table.csv")
tax.table<-tax_table(tam.class)
write.csv(tax.table,"tam.class.tax.table.csv")
