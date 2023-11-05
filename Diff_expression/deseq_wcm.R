setwd("/projectnb/mullenl/alqassar/diff_expression_wcm/DEseq/")

library(tidyr)
library(plyr)
library(dplyr)
library(DESeq2)
library(ggplot2)
library(cowplot)
library(readr)
library(RColorBrewer)
library(gplots)
library(knitr)
library(plotly)
library(vegan)
library(kableExtra)
library(reshape2)
library(prettydoc)
library(VennDiagram)


#library(affycoretools)
#library(arrayQualityMetrics)
#library(genefilter)
#library(DESeq)

setwd("/projectnb/mullenl/alqassar/diff_expression_wcm/")
countData <- read.table("wcm_iso_counts.txt")
names(countData) = c("adult", "late_stage_larva", "early_stage_larva")

treat = c("adult","larvae","larvae")
experimental_design = data.frame(names(countData), treat)
names(experimental_design) = c("sample","treat")

#summing across the counts for genes and seeing how many total counts we have for each sample
totalCounts=colSums(countData)
totalCounts
barplot(totalCounts, ylab="raw counts", main = "WCM total counts")
# looks good

dds<-DESeqDataSetFromMatrix(countData=countData, colData=experimental_design, design=~treat) #can only test for the main effects of treatment
dds = DESeq(dds)
results = results(dds)
summary(results)

#Let's check to see if we set up our contrast correctly. We should have the treatment condition first and the control second in the log2 fold change (MLE) output. 
head(results)

# write out file of normalized counts
norm.counts = counts(dds, normalized = TRUE) # these are the counts DESeq uses
write.csv(norm.counts, "normalized_iso_counts_wcm.csv") #these are all counts, not considering treatment comparisons

#Now do an rlogged transformation, which is useful various unsupervised clustering analyses. Be sure the include blind = TRUE as it doesn't normalize the data with any priors from our experimental design. 
rlogged = rlogTransformation(dds, blind = TRUE)


res_wcm <- results(dds, contrast=c("treat","larvae","adult"))
head(res_wcm)
#how many FDR < 10%?
table(res_wcm$padj<0.01) # we actually get more DEGs when genotype is included in model.
# 0.1=17643
# 0.05=14838
# 0.01=10837
summary(res_wcm)
write.table(res_wcm, file="wcm_iso_dds_results.txt", quote=F, sep="\t")

# Heatmap of overall expression, sample by distance.
sampleDists <- as.matrix(dist(t(assay(rlogged))))
heatmap.2(as.matrix(sampleDists), key=F, trace="none",
          col=colorpanel(100, "black", "white"),
          margin=c(10, 10))

# Individual gene heatmaps
#Now plotting the z scores, as heatmap2 creates nice clean clusters by doing this. Upregulation indicated by warmer colors, downregulation indicated by cooler colors.

norm_counts = read.csv("normalized_iso_counts_wcm.csv")
hm = read.table("wcm_iso_dds_results.txt", header=TRUE) %>% 
  tibble::rownames_to_column("X") %>%
  filter(padj < 0.01) %>% # only want the most DEGs
  select(X) %>%
  merge(norm_counts, by.x = "X", by.y = "X")  # turn into a countdatafile
row.names(hm) = hm$X
hm$X = NULL

head(hm)

## Turning into z-score table
hm.z = data.matrix(hm)
hm.z = sweep(hm.z, 1L, rowMeans(hm.z), check.margin = FALSE)
hm.z.sx = apply(hm.z, 1L, sd)
hm.z = sweep(hm.z, 1L, hm.z.sx, "/", check.margin = FALSE)
hm.z = data.matrix(hm.z)

colour = colorRampPalette(rev(c("chocolate1","#FEE090","grey10", "cyan3","cyan")))(100)
# heatmap.2(hm.z, col = colour, Rowv = TRUE, Colv = TRUE, scale = "row", 
          dendrogram = "both",
          trace = "none", 
          margin = c(5,15))

# heatmap.2(hm.z, col = colour, Rowv = TRUE, Colv = TRUE, scale = "column", 
          dendrogram = "both",
          trace = "none", 
          margin = c(5,15))

heatmap.2(hm.z, col = colour, Rowv = TRUE, Colv = TRUE, scale = "column", 
          dendrogram = "both",
          trace = "none", 
          margin = c(6,20),
          cexCol=0.7)
dim(hm.z)
# Make GO input file
#sym
res = read.table("wcm_iso_dds_results.txt")
head(res)
res$isogroup=row.names(res)

go_input_wcm = res %>%
  mutate(mutated_p = -log(pvalue)) %>%
  mutate(mutated_p_updown = ifelse(log2FoldChange < 0, mutated_p*-1, mutated_p*1)) %>%
  select(isogroup, mutated_p_updown) %>%
  na.omit()

nrow(go_input_wcm)
head(go_input_wcm)
colnames(go_input_wcm) <- c("gene", "pval")
head(go_input_wcm)
write.csv(go_input_wcm, file="wcm_iso_GO.csv", quote=F, row.names=FALSE)

### Now GO pipeline 
## see GO_MWU_wcm.R





