# Adapted by Jasmine Alqassar from Mikhail Matz's GO_MWU.R https://github.com/z0on/GO_MWU

setwd("/projectnb/mullenl/alqassar/diff_expression_wcm/DEseq/")

# Edit these to match your data file names: 
goAnnotations="wcm_gene2go.tab" # two-column, tab-delimited, one line per gene, multiple GO terms separated by semicolon. If you have multiple lines per gene, use nrify_GOtable.pl prior to running this script.
goDatabase="go.obo" # download from http://www.geneontology.org/GO.downloads.ontology.shtml
input="wcm_iso_GO.csv" # two columns of comma-separated values: gene id, continuous measure of significance. To perform standard GO enrichment analysis based on Fisher's exact test, use binary measure (0 or 1, i.e., either sgnificant or not).
goDivision="MF" # either MF, or BP, or CC
source("gomwu.functions.R")

# ------------- Calculating stats
# It might take a few minutes for MF and BP. Do not rerun it if you just want to replot the data with different cutoffs, go straight to gomwuPlot. If you change any of the numeric values below, delete the files that were generated in previos runs first.

gomwuStats(input, goDatabase, goAnnotations, goDivision,
           perlPath="perl", # replace with full path to perl executable if it is not in your system's PATH already
           largest=0.1,  # a GO category will not be considered if it contains more than this fraction of the total number of genes
           smallest=5,   # a GO category should contain at least this many genes to be considered
           clusterCutHeight=0.25, # threshold for merging similar (gene-sharing) terms. See README for details.
           #	Alternative="g" # by default the MWU test is two-tailed; specify "g" or "l" of you want to test for "greater" or "less" instead. 
           #	Module=TRUE,Alternative="g" # un-remark this if you are analyzing a SIGNED WGCNA module (values: 0 for not in module genes, kME for in-module genes). In the call to gomwuPlot below, specify absValue=0.001 (count number of "good genes" that fall into the module)
           #	Module=TRUE # un-remark this if you are analyzing an UNSIGNED WGCNA module 
) 

# do not continue if the printout shows that no GO terms pass 10% FDR.
# RESULT: 448 GO terms at 10% FDR

quartz()
mf_wcm_results=gomwuPlot(input,goAnnotations,goDivision,
                                    absValue=-log(0.05,10),  # genes with the measure value exceeding this will be counted as "good genes". This setting is for signed log-pvalues. Specify absValue=0.001 if you are doing Fisher's exact test for standard GO enrichment or analyzing a WGCNA module (all non-zero genes = "good genes").
                                    #absValue=1, # un-remark this if you are using log2-fold changes
                                    level1=0.05, # FDR threshold for plotting. Specify level1=1 to plot all GO categories containing genes exceeding the absValue.
                                    level2=0.01, # FDR cutoff to print in regular (not italic) font.
                                    level3=0.001, # FDR cutoff to print in large bold font.
                                    txtsize=1,    # decrease to fit more on one page, or increase (after rescaling the plot so the tree fits the text) for better "word cloud" effect
                                    treeHeight=0.5, # height of the hierarchical clustering tree
                                    #	colors=c("dodgerblue2","firebrick1","skyblue2","lightcoral") # these are default colors, un-remar and change if needed
)
# text representation of results, with actual adjusted p-values
mf_hot_brown_host_results[[1]]
write.csv(mf_hot_brown_host_results, "mf_oculina_host_brown_heat_results.csv")

goAnnotations="/projectnb/mullenl/alqassar/diff_expression_wcm/DEseq/wcm_gene2go.tab" # two-column, tab-delimited, one line per gene, multiple GO terms separated by semicolon. If you have multiple lines per gene, use nrify_GOtable.pl prior to running this script.
goDatabase="/projectnb/mullenl/alqassar/diff_expression_wcm/DEseq/go.obo" # download from http://www.geneontology.org/GO.downloads.ontology.shtml
input="/projectnb/mullenl/alqassar/diff_expression_wcm/DEseq/wcm_iso_GO.csv" # two columns of comma-separated values: gene id, continuous measure of significance. To perform standard GO enrichment analysis based on Fisher's exact test, use binary measure (0 or 1, i.e., either sgnificant or not).
goDivision="MF" # either MF, or BP, or CC
source("/projectnb/mullenl/alqassar/diff_expression_wcm/DEseq/gomwu.functions.R")

goDivision="BP" 
gomwuStats(input, goDatabase, goAnnotations, goDivision,
           perlPath="/usr/bin/perl", # replace with full path to perl executable if it is not in your system's PATH already
           largest=0.1,  # a GO category will not be considered if it contains more than this fraction of the total number of genes
           smallest=5,   # a GO category should contain at least this many genes to be considered
           clusterCutHeight=0.25, # threshold for merging similar (gene-sharing) terms. See README for details.
           #	Alternative="g" # by default the MWU test is two-tailed; specify "g" or "l" of you want to test for "greater" or "less" instead. 
           #	Module=TRUE,Alternative="g" # un-remark this if you are analyzing a SIGNED WGCNA module (values: 0 for not in module genes, kME for in-module genes). In the call to gomwuPlot below, specify absValue=0.001 (count number of "good genes" that fall into the module)
           #	Module=TRUE # un-remark this if you are analyzing an UNSIGNED WGCNA module 
) 
# do not continue if the printout shows that no GO terms pass 10% FDR.
# RESULT:  1688 GO terms at 10% FDR

BP_results=gomwuPlot(input,goAnnotations,goDivision,
                     absValue=-log(0.05,10),  # genes with the measure value exceeding this will be counted as "good genes". This setting is for signed log-pvalues. Specify absValue=0.001 if you are doing Fisher's exact test for standard GO enrichment or analyzing a WGCNA module (all non-zero genes = "good genes").
                     #absValue=1, # un-remark this if you are using log2-fold changes
                     level1=0.001, # FDR threshold for plotting. Specify level1=1 to plot all GO categories containing genes exceeding the absValue.
                     level2=0.0001, # FDR cutoff to print in regular (not italic) font.
                     level3=0.00001, # FDR cutoff to print in large bold font.
                     txtsize=0.75,    # decrease to fit more on one page, or increase (after rescaling the plot so the tree fits the text) for better "word cloud" effect
                     treeHeight=0.5, # height of the hierarchical clustering tree
                     #	colors=c("dodgerblue2","firebrick1","skyblue2","lightcoral") # these are default colors, un-remar and change if needed
)

