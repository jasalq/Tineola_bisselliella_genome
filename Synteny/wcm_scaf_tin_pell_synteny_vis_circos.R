# Script written by Jasmine Alqassar 2024 for circular visulization of synteny alignment data from Satsuma 
# Inspiration taken from: 
# https://github.com/rstewa03/Pieris_macdunnoughii_genome/blob/8eaafd3a8e4153a6015b5e589d3f39b42b595164/SI/Pmac07_Synteny.R/Pmac07_Synteny.R
# https://github.com/annaorteu/Hypolimnas_genome_Wchr_evolution/blob/44b5d865eba05ae2ac8fbd6da54663e1100acc2f/scripts/Satsuma_comparisons_multiple_species_against_each_other_SuppFigure.R
# to make a karyotype file, use samtools -faidx [fasta] and take the first and second column from this file into an excel sheet and place a column before the second column with 0s, then add headers" Chr, Start, End
# load necessary packages  
library(tidyverse)
library(circlize)
library(RColorBrewer)

data_dir <- "/projectnb/mullenl/alqassar/synteny/circos/"
ref_karyotype <- "karyotype_tin_pell" # ref = your comparison genome 
qry_karyotype <- "karyotype_wcm_scaf" #qry = your genome 
synteny_data <- "satsuma_scaff_summary.chained.out"
setwd(data_dir)


Ref_pal <- read.delim(paste0(data_dir,ref_karyotype), sep = "\t", header = TRUE) 
colnames(Ref_pal) <- c("Chr1", "name", "Start", "End", "fill", "species", "size", "color1")


for(i in 1:nrow(Ref_pal)){
  Ref_pal$Chr[i] <-  strsplit(Ref_pal$Chr1[i], ' ') [[1]][3]}


Ref_pal = Ref_pal %>%
  select(Chr,Start,End) %>%
  relocate("Chr", .before = "Start")

Ref_pal$color <- NULL #creating an empty column to fill later


Ref_pal = Ref_pal %>%
  mutate(color = rep(brewer.pal(n=12, "Paired"), length.out = length(Chr))) #Paired is the name of a colorblindness friendly palette 
#to display colorblind friendly palettes to choose from: display.brewer.all(colorblindFriendly = TRUE) 

Qry_pal <- read.delim(paste0(data_dir,qry_karyotype), sep = "\t", header = TRUE) 

Qry_pal$color <- NULL #creating an empty column to fill later

Qry_pal = Qry_pal %>%
  mutate(color = "#FFFFFF")

# before visualizing in an ideogram you need to add different prefixes to the chromosomes
Ref_pal[ ,1] = paste0("Tinea_chr", Ref_pal[, 1])
Pal1 <- union(Ref_pal, Qry_pal) #combine these files 

synteny_alignment <- read.delim(paste0(data_dir, synteny_data), sep = "\t", header = FALSE)
colnames(synteny_alignment) <- c("ref_chrom_temp", "Start_ref", "End_ref", "qry_chrom_temp", 
                                 "Start_qry", "End_qry", "V6", "V7")

synteny_alignment = synteny_alignment %>%
  select(-V6, -V7) 

for(i in 1:nrow(synteny_alignment)){
  synteny_alignment$Ref[i] <-  strsplit(synteny_alignment$ref_chrom_temp[i], '_')[[1]][7]
  synteny_alignment$Qry[i] <-  strsplit(synteny_alignment$qry_chrom_temp[i], '_')[[1]][2]}

synteny_alignment = synteny_alignment %>%
  select(-ref_chrom_temp, -qry_chrom_temp) %>%
  relocate("Ref", .before = "Start_ref") %>%
  relocate("Qry", .before = "Start_qry") %>%
  filter(!grepl('mitochondrion', Ref)) %>%
  filter(!grepl('scaffold', Ref)) %>%
  drop_na() 

synteny_alignment[ ,1] = paste0("Tinea_chr", synteny_alignment[, 1])
synteny_alignment[ ,4] = paste0("Tineola_scaf", synteny_alignment[, 4])

# circos visualization
#svg(file="wcm_pseudochr_tinpell_circos_May15.svg") #making an empty svg file to be filled 
png(filename = "scaf_tin_pell_May20_highres2.png",
    width = 2000, height = 2000, units = "px")
circos.initializeWithIdeogram(Pal1,plotType = NULL)
circos.track(ylim = c(0, 1), panel.fun = function(x, y) {
  chr = CELL_META$sector.index
  xlim = CELL_META$xlim
  ylim = CELL_META$ylim
  circos.rect(xlim[1], 0, xlim[2], 0.5, col = Pal1$color[Pal1$Chr == chr])
  circos.text(mean(xlim), 1, chr, cex = 2, facing = "clockwise",
              niceFacing = TRUE,  adj = c(0.4, 0.5)) }
  , bg.border = NA)


#text(-0.9, -0.8, "Tinea pellionella\ngenome")
#text(0.9, 0.8, "Tineola bisselliella\ngenome")


for ( i in c(1:dim(synteny_alignment)[1])){
  reference <- synteny_alignment$Ref[i]
  query <- synteny_alignment$Qry[i]
  RS <- synteny_alignment$Start_ref[i]
  RE <- synteny_alignment$End_ref[i]
  QS <- synteny_alignment$Start_qry[i]
  QE <- synteny_alignment$End_qry[i]
  
  circos.link(reference, c(RS,RE),
              query, c(QS, QE), col = Pal1$color[Pal1$Chr == reference]
  )}


circos.clear()
dev.off()





# then on your computer convert SVG file 
# rsvg-convert -d 1000 -p 1000 file.svg > file.png 




#Notes from Joe:
#DotPlot(object = Hm.sops.sct, features = row3,group.by='alpha',dot.min=0.1,col.min=0)+ theme(axis.text.x = element_text(angle = 30, hjust=1))+scale_color_batlow(reverse=T)+scale_x_discrete(labels=row3names)
#dev.copy(pdf,'DotPlot.alpha.row3.pdf',width=20,height=4)
#dev.off()

